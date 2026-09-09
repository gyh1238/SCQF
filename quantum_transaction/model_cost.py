"""
Two-qubit gate cost model for the zone-local circuit (Sec. IV, Sec. V-C).
=========================================================================
`HARDWARE_LIMITS.md` established that the executable size of these circuits
is set by the two-qubit gate count -- not by the qubit count -- because the
measured decay is p_signal = exp(-lambda * n_2q) with lambda anchored on
ibm_kingston (Heron2).  The partitioner therefore takes a *two-qubit gate
budget* as its constraint, and the scaling figure plots the same quantity.

The model counts one compute-mark-uncompute pass of the oracle:

    compute   : utility rotations, validity flags, RB counters + comparators,
                AP counters + comparators
    mark      : one superflag MCX over every constraint-outcome qubit
    uncompute : mirror of the compute phase

Amplification multiplies the whole pass by (2k+1) for k rounds; the figures
compare centralized and zone-local circuits at one pass, which is the
conservative choice (a centralized circuit has a smaller accepted mass and
would need *more* rounds than a zone).

`calibrate.py`-style validation against real Qiskit transpilation lives in
`certify_sampler.py`.
"""

from math import ceil, log2

import numpy as np

# --- measured hardware anchors, from HARDWARE_LIMITS.md ------------------
# p_signal(n_2q) = exp(-LAMBDA * n_2q);  the ceiling quoted below is the
# 50%-signal point, n = ln2 / lambda.
LAMBDA_KINGSTON = 0.00118        # measured: p ~ 0.53 at 640 heavy-hex 2q
LAMBDA_BOSTON = 0.00068          # scaled from 2q error ratio (estimate)

CEIL_KINGSTON = int(round(np.log(2) / LAMBDA_KINGSTON))   # ~587 2q at 50% signal
CEIL_BOSTON = int(round(np.log(2) / LAMBDA_BOSTON))       # ~1019 2q at 50% signal

# Operating budget used by the partitioner.  It must exceed the cost of the
# most expensive single AP, since an AP is the atom of the partition; the
# measured ceilings above are reported alongside it as device references.
BUDGET_DEFAULT = 2000


def mcx_2q(n_controls):
    """Two-qubit gate count of an n-controlled X.

    Toffoli decomposes into 6 CX; an n-control MCX with a borrowed-ancilla
    v-chain costs 2n-3 Toffolis for n >= 3.
    """
    if n_controls <= 0:
        return 0
    if n_controls == 1:
        return 1
    if n_controls == 2:
        return 6
    return 6 * (2 * n_controls - 3)


def _qft_2q(w):
    """Controlled-phase count of a QFT on w qubits (approximate QFT: same order)."""
    return w * (w - 1) // 2


def zone_2q_cost(inst, ap_set, detail=False):
    """
    Two-qubit gate count of one compute-mark-uncompute pass for the zone
    owning `ap_set`, together with its register widths.

    Follows the manuscript's structure: UE code registers of width l_i, one
    cost qubit and one validity flag per local UE, a reusable counter of
    width ceil(log2(N_z+1)), one violation flag per owned RB and AP, and a
    single superflag conjunction.
    """
    ap_set = set(int(a) for a in ap_set)
    rb_owner = inst.rb_owner
    in_zone_rb = np.array([o in ap_set for o in rb_owner])

    # local UE set U_z = {i : V_i intersects R_z}
    ue_idx = [i for i, c in enumerate(inst.cand) if in_zone_rb[c].any()]
    n_z = len(ue_idx)
    if n_z == 0:
        empty = dict(n2q=0, n_ue=0, q_state=0, n_qubits=0, n_rb_used=0, n_ap=len(ap_set))
        return empty if detail else 0

    ell = {i: inst.code_width(i) for i in ue_idx}
    q_state = sum(ell.values())                      # Eq. (state-width)
    w = max(1, ceil(log2(n_z + 1)))                  # counter width

    # ---- (1) utility rotations: every candidate of every local UE -------
    # a candidate owned outside the zone still needs its rotation, with
    # owner-local utility 0 (full-domain convention of Sec. III-C).
    c_rot = 0
    for i in ue_idx:
        for _ in inst.cand[i]:
            c_rot += 2 * mcx_2q(ell[i]) + 2         # match, cRy, unmatch

    # ---- (2) validity flags for unused codewords ------------------------
    c_val = 0
    for i in ue_idx:
        n_invalid = 2 ** ell[i] - len(inst.cand[i])
        c_val += n_invalid * mcx_2q(ell[i])

    # ---- (3) RB access limit: one counter pass per owned, used RB -------
    rbs_used = {}
    for i in ue_idx:
        for r in inst.cand[i]:
            if in_zone_rb[r]:
                rbs_used.setdefault(int(r), []).append(i)
    c_rb = 0
    for r, users in rbs_used.items():
        acc = sum(2 * mcx_2q(ell[i]) + w for i in users)   # match + Fourier add
        c_rb += 2 * _qft_2q(w) + 2 * acc + 2 * mcx_2q(w)   # qft/iqft, acc+unacc, compare

    # ---- (4) AP admission limit: one counter pass per owned AP ----------
    c_ap = 0
    for a in ap_set:
        users = [(i, r) for i in ue_idx for r in inst.cand[i] if rb_owner[r] == a]
        if not users:
            continue
        acc = sum(2 * mcx_2q(ell[i]) + w for i, _ in users)
        c_ap += 2 * _qft_2q(w) + 2 * acc + 2 * mcx_2q(w)

    # ---- (5) superflag conjunction --------------------------------------
    n_flags = n_z + n_z + len(rbs_used) + len(ap_set)     # cost, validity, RB, AP
    c_sf = mcx_2q(n_flags)

    n2q = 2 * (c_rot + c_val + c_rb + c_ap) + c_sf

    if not detail:
        return n2q

    n_qubits = (q_state + n_z + n_z + len(rbs_used) + len(ap_set)
                + w + max(1, w - 1) + 1)              # + counter, cmp workspace, superflag
    return dict(n2q=int(n2q), n_ue=n_z, ue_idx=ue_idx, q_state=q_state,
                n_qubits=int(n_qubits), n_rb_used=len(rbs_used), n_ap=len(ap_set),
                counter_w=w,
                breakdown=dict(rot=2 * c_rot, val=2 * c_val, rb=2 * c_rb,
                               ap=2 * c_ap, superflag=c_sf))


def centralized_2q_cost(inst):
    """One-pass cost of the single circuit that solves the whole instance."""
    return zone_2q_cost(inst, set(range(inst.n_ap)), detail=True)


if __name__ == "__main__":
    from model_instance import make_instance

    print(f"hardware 2q ceilings (50% signal):  kingston={CEIL_KINGSTON}  boston={CEIL_BOSTON}")
    for g in (3, 4, 5, 6):
        inst = make_instance(g=g, seed=1)
        c = centralized_2q_cost(inst)
        print(f"g={g}: centralized  UEs={c['n_ue']:3d}  Q_state={c['q_state']:3d}  "
              f"qubits={c['n_qubits']:3d}  2q={c['n2q']:7d}")
