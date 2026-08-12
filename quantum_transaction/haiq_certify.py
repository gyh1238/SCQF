"""
Statevector certification of the zone-local law.
================================================
`haiq_zone.sample_zone` claims to reproduce the accepted branch of the
circuit exactly rather than to approximate it.  This module checks that
claim on a real Qiskit circuit built to the manuscript's specification:

  * one code register per local UE, prepared uniformly (q_0);
  * one controlled-R_Y per candidate with theta = 2 arccos sqrt(g), Eq.(theta);
  * an RB occupancy counter and an AP admission counter, each compared
    against its limit and recorded in a violation flag, then uncomputed;
  * acceptance = (all cost qubits |0>) AND (all violation flags |0>).

The conditional distribution of the code registers on that acceptance event
is compared with exp(lambda J_z) restricted to F_z.  Agreement to statevector
precision is what licenses using the classical sampler at region sizes no
present-day processor can run.

The comparison is a distributional identity, so it is independent of how the
counters are decomposed; the arithmetic here uses ripple incrementers and
exact threshold patterns rather than the QFT adders of Sec. V-C, which cost
differently but realize the same phi_z.
"""

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import Statevector

from haiq_instance import make_instance, utility_scale
from haiq_partition import partition_aps
from haiq_zone import build_zone, enumerate_zone, _per_ue_weights


def _increment(qc, ctrl, cnt):
    """Controlled +1 on a little-endian counter (cnt[0] = LSB)."""
    for j in reversed(range(len(cnt))):
        qc.mcx([ctrl] + list(cnt[:j]), cnt[j])


def _decrement(qc, ctrl, cnt):
    for j in range(len(cnt)):
        qc.mcx([ctrl] + list(cnt[:j]), cnt[j])


def _pattern_mcx(qc, reg, value, target):
    """Fire `target` exactly when `reg` holds `value` (little-endian)."""
    zeros = [q for j, q in enumerate(reg) if not (value >> j) & 1]
    for q in zeros:
        qc.x(q)
    qc.mcx(list(reg), target)
    for q in zeros:
        qc.x(q)


def build_zone_circuit(zone, lam):
    """Circuit realizing q_0, the utility rotations, and the strict flags."""
    n = zone.n_ue
    widths = zone.code_widths()
    assert all(w == 1 for w in widths), "certification uses 1-qubit codes"

    code = QuantumRegister(n, "code")
    cost = QuantumRegister(n, "cost")
    w_cnt = max(1, int(np.ceil(np.log2(n + 1))))
    cnt = QuantumRegister(w_cnt, "cnt")
    rb_flag = QuantumRegister(len(zone.rbs), "rbf")
    ap_flag = QuantumRegister(len(zone.aps), "apf")
    # a UE whose candidate count is not a power of two leaves unused
    # codewords; Sec. IV-A marks each of them with a validity flag so that
    # they receive zero accepted probability
    novel = [k for k in range(n) if len(zone.cand[k]) < 2 ** widths[k]]
    val_flag = QuantumRegister(max(1, len(novel)), "valf")
    qc = QuantumCircuit(code, cost, cnt, rb_flag, ap_flag, val_flag)

    # ---- q_0 : uniform over valid codewords -----------------------------
    qc.h(code)

    # ---- validity flags for unused codewords ----------------------------
    for fi, k in enumerate(novel):
        for v in range(len(zone.cand[k]), 2 ** widths[k]):
            _pattern_mcx(qc, [code[k]], v, val_flag[fi])

    # ---- utility rotations, Eq. (theta) ---------------------------------
    g = _per_ue_weights(zone, lam)
    for k in range(n):
        for c in range(len(zone.cand[k])):
            theta = 2 * np.arccos(np.sqrt(np.clip(g[k][c], 0.0, 1.0)))
            if theta == 0.0:
                continue
            if c == 0:
                qc.x(code[k])
            qc.cry(theta, code[k], cost[k])
            if c == 0:
                qc.x(code[k])
    qc.barrier()

    # ---- RB access limit -------------------------------------------------
    for fi, r in enumerate(zone.rbs):
        users = [(k, int(np.where(zone.cand[k] == r)[0][0]))
                 for k in range(n) if r in zone.cand[k].tolist()]
        users = [(k, c) for k, c in users if zone.owned[k][c]]
        if not users:
            continue
        for k, c in users:
            if c == 0:
                qc.x(code[k])
            _increment(qc, code[k], cnt)
            if c == 0:
                qc.x(code[k])
        for v in range(zone.w_rb + 1, n + 1):
            _pattern_mcx(qc, cnt, v, rb_flag[fi])
        for k, c in reversed(users):
            if c == 0:
                qc.x(code[k])
            _decrement(qc, code[k], cnt)
            if c == 0:
                qc.x(code[k])
    qc.barrier()

    # ---- AP admission limit ---------------------------------------------
    for fi, a in enumerate(zone.aps):
        users = []
        for k in range(n):
            for c in range(len(zone.cand[k])):
                if zone.owned[k][c] and zone.rb_owner[zone.cand[k][c]] == a:
                    users.append((k, c))
        if not users:
            continue
        for k, c in users:
            if c == 0:
                qc.x(code[k])
            _increment(qc, code[k], cnt)
            if c == 0:
                qc.x(code[k])
        for v in range(zone.w_ap + 1, n + 1):
            _pattern_mcx(qc, cnt, v, ap_flag[fi])
        for k, c in reversed(users):
            if c == 0:
                qc.x(code[k])
            _decrement(qc, code[k], cnt)
            if c == 0:
                qc.x(code[k])

    return qc, dict(n=n, w_cnt=w_cnt, n_rb=len(zone.rbs), n_ap=len(zone.aps),
                    n_val=max(1, len(novel)))


def certify(zone, lam):
    """Total-variation distance between the accepted branch and exp(lambda J_z)."""
    qc, meta = build_zone_circuit(zone, lam)
    sv = Statevector.from_instruction(qc)
    p = np.abs(np.asarray(sv.data)) ** 2

    n, w = meta["n"], meta["w_cnt"]
    n_aux = w + meta["n_rb"] + meta["n_ap"] + meta["n_val"]
    idx = np.arange(p.size, dtype=np.int64)
    code = idx & ((1 << n) - 1)                       # qubits 0..n-1
    rest = idx >> n                                   # cost, then every ancilla
    accepted = rest & ((1 << (n + n_aux)) - 1)        # all must be zero
    m = accepted == 0

    hist = np.bincount(code[m], weights=p[m], minlength=1 << n)
    mu = hist.sum()
    hist = hist / mu

    rows, pref, _ = enumerate_zone(zone, lam)
    ref = np.zeros(1 << n)
    for row, pp in zip(rows, pref):
        ref[int(sum(int(c) << k for k, c in enumerate(row)))] = pp

    tvd = 0.5 * float(np.abs(hist - ref).sum())
    spurious = int(np.count_nonzero((hist > 1e-12) & (ref == 0)))
    return dict(tvd=tvd, mu=float(mu), n_accepted=int((hist > 1e-12).sum()),
                n_feasible=len(rows), n_qubits=qc.num_qubits, depth=qc.depth(),
                spurious=spurious)


def certify_instance(g=4, seed=2, beta=1.5, max_ue=7, n_zones=6,
                     max_qubits=24, verbose=True):
    inst = make_instance(g=g, seed=seed)
    part = partition_aps(inst)
    lam = beta / utility_scale(inst)
    out = []
    for z in range(part.n_zones):
        zone = build_zone(inst, part, z)
        if zone.n_ue == 0 or zone.n_ue > max_ue:
            continue
        if any(w != 1 for w in zone.code_widths()):
            continue
        w = max(1, int(np.ceil(np.log2(zone.n_ue + 1))))
        if 2 * zone.n_ue + w + len(zone.rbs) + len(zone.aps) + 1 > max_qubits:
            continue
        r = certify(zone, lam)
        r["zone"] = z
        r["n_ue"] = zone.n_ue
        out.append(r)
        if verbose:
            print(f"  zone {z:2d}  N_z={zone.n_ue}  qubits={r['n_qubits']:2d}  "
                  f"|F_z|={r['n_feasible']:3d}  mu={r['mu']:.4f}  "
                  f"TVD={r['tvd']:.2e}  spurious={r['spurious']}")
        if len(out) >= n_zones:
            break
    return out


if __name__ == "__main__":
    print("Statevector certification: accepted branch vs exp(lambda J_z) on F_z")
    res = certify_instance()
    print(f"\nmax TVD over {len(res)} zones = {max(r['tvd'] for r in res):.3e}")
