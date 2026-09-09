"""
Does the released oracle realize the exponential weight of Eq. (19)-(21)?
========================================================================
Section IV-C claims one thing about the circuit that the rest of the paper
leans on: conditioned on acceptance, the state register follows

    p(z) proportional to exp[lambda J(z)]   on the feasible set,      Eq. (30)

and *not merely* some monotone function of J(z).  The distinction is the whole
of Sec. III-B: only an exponential turns the additive utility of Eq. (13) into
a product across zones, Eq. (15), so only an exponential lets zone-local laws
be composed.  A weight that ranks assignments correctly but bends the spacing
between them still ranks correctly inside one zone and still fails to compose.

That is exactly the failure mode of the rotation the enumeration oracle of
[jang2025cqf] shipped,

    theta = arccos sqrt(u / u_max),                                (LINEAR)

which is linear in the utility and carries half the angle Eq. (20) asks for.
This module measures both forms on the same circuits:

  * `law_from_statevector` extracts the accepted-branch law exactly, by
    projecting the statevector onto (cost register all |0>) and the feasible
    codes, then renormalizing -- the same conditioning the hardware performs
    by discarding rejected shots.
  * `tvd_to_gibbs` compares it with exp[lambda J] / Z on the feasible set.
  * `spread` reports max/min of the ratio law(z) / exp[lambda J(z)], which is
    1.0 exactly when the weight is the intended one and grows without bound as
    the weight bends.  This is the number that decides composability; a TVD can
    be small for a wrong-shaped weight on a nearly-flat instance.

Run:  python certify_rotation.py
"""

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import Statevector

import circuit_inter as inter
import circuit_intra as intra


# --------------------------------------------------------------- exact laws

def prep_law(build_prep, n_assign, n_cost, feasible):
    """
    The accepted-branch law of Eq. (29)-(30), exactly.

    `build_prep` returns the state-preparation circuit P_z = Theta_z H^(x)Q of
    Eq. (26) -- the Hadamard layer and the utility rotations, with the
    rotations *left standing*.  The oracle of Eq. (32) computes the constraint
    flags, marks, and mirrors the rotations away again, so the accepted branch
    is fixed here and amplification only rescales its total mass, Eq. (33)-(34);
    reading the law off P_z is therefore reading the law the sampler emits.

    Conditioning is the one Sec. IV-F performs on hardware: keep an outcome
    only when every cost qubit reads |0> and the decoded assignment is
    feasible, then renormalize.  Registers are laid out assign-then-cost, so
    the assignment code is the last `n_assign` characters of Qiskit's
    q[n-1]..q[0] string and the cost register sits directly above it.
    """
    qc = build_prep()
    probs = Statevector.from_instruction(qc).probabilities_dict()
    out = {}
    for bits, p in probs.items():
        if p <= 0.0:
            continue
        rev = bits[::-1]                              # rev[k] = q[k]
        if "1" in rev[n_assign:n_assign + n_cost]:    # cost not all |0>
            continue
        code = bits[-n_assign:]
        if code not in feasible:                      # rejected by the superflag
            continue
        out[code] = out.get(code, 0.0) + p
    total = sum(out.values())
    return {k: v / total for k, v in out.items()}, total


def gibbs(feasible_scores, lam):
    """exp[lambda J] / Z on the feasible set -- the law Eq. (30) asks for."""
    w = {k: np.exp(lam * j) for k, j in feasible_scores.items()}
    z = sum(w.values())
    return {k: v / z for k, v in w.items()}


def tvd(p, q):
    keys = set(p) | set(q)
    return 0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in keys)


def spread(p, q):
    """max/min of p/q over the shared support: 1.0 iff the shapes agree."""
    r = [p[k] / q[k] for k in p if q.get(k, 0.0) > 0 and p[k] > 0]
    return max(r) / min(r) if r else float("inf")


# ------------------------------------------------------- the two rotations

def _exp_angles(u_row, lam):
    """Eq. (19)-(21): per-UE exponential weight, full angle."""
    u_max = float(np.max(u_row))
    g = np.exp(-lam * (u_max - np.asarray(u_row, dtype=float)))
    return 2.0 * np.arccos(np.sqrt(np.clip(g, 0.0, 1.0)))


def _lin_angles(u_row, u_global_max, lam):
    """The rotation [jang2025cqf] released: linear weight, half the angle."""
    r = np.asarray(u_row, dtype=float) / u_global_max
    return np.arccos(np.sqrt(np.clip(r, 0.0, 1.0)))


# ------------------------------------------------------------- the two cases

def _prep_inter(utility, angles):
    """P_z for the 1-bit-per-UE inter-cell encoding: assign, then cost."""
    N = utility.shape[0]
    assign = QuantumRegister(N, "assign")
    cost = QuantumRegister(N, "cost")
    qc = QuantumCircuit(assign, cost)
    qc.h(assign)
    for i in range(N):
        for j in range(2):
            th = angles[i][j]
            if th == 0.0:
                continue
            if j == 0:
                qc.x(assign[i])
            qc.cry(th, assign[i], cost[i])
            if j == 0:
                qc.x(assign[i])
    return qc


def _prep_intra(throughput, angles):
    """P_z for the 2-bit-per-node intra-cell encoding: assign, cost, match."""
    N, M = throughput.shape
    assign = QuantumRegister(2 * N, "assign")
    cost = QuantumRegister(N, "cost")
    match = QuantumRegister(1, "match")
    qc = QuantumCircuit(assign, cost, match)
    qc.h(assign)
    for i in range(N):
        lsb, msb = intra._pair(assign, i)
        for r in range(M):
            th = angles[i][r]
            if th == 0.0:
                continue
            intra._mask_pair(qc, lsb, msb, r)
            qc.mcx([lsb, msb], match[0])
            qc.cry(th, match[0], cost[i])
            qc.mcx([lsb, msb], match[0])
            intra._mask_pair(qc, lsb, msb, r)
    return qc


def check_inter(utility, weights, cap, beta=inter.BETA):
    lam = beta / inter.utility_scale(utility)
    N = utility.shape[0]
    feas = {"".join(str(b) for b in reversed(bits)): score
            for bits, score in inter.enumerate_feasible(utility, weights, cap)}
    ref = gibbs(feas, lam)
    gmax = float(np.max(utility))

    rows = []
    for name, ang in (
            ("exponential, Eq. (19)-(21)",
             [_exp_angles(utility[i], lam) for i in range(N)]),
            ("linear, [jang2025cqf]",
             [_lin_angles(utility[i], gmax, lam) for i in range(N)])):
        law, mass = prep_law(lambda a=ang: _prep_inter(utility, a), N, N, feas)
        rows.append((name, tvd(law, ref), spread(law, ref), mass))
    return lam, len(feas), rows


def check_intra(throughput, W, beta=intra.BETA):
    lam = beta / intra.utility_scale(throughput)
    N, M = throughput.shape
    feas = {}
    for codes, score in intra.enumerate_feasible(throughput, W):
        bits = "".join(f"{c & 1}{(c >> 1) & 1}" for c in codes)   # q0 q1 q2 ...
        feas[bits[::-1]] = score
    ref = gibbs(feas, lam)
    gmax = float(np.max(throughput))

    rows = []
    for name, ang in (
            ("exponential, Eq. (19)-(21)",
             [_exp_angles(throughput[i], lam) for i in range(N)]),
            ("linear, [jang2025cqf]",
             [_lin_angles(throughput[i], gmax, lam) for i in range(N)])):
        law, mass = prep_law(lambda a=ang: _prep_intra(throughput, a),
                             2 * N, N, feas)
        rows.append((name, tvd(law, ref), spread(law, ref), mass))
    return lam, len(feas), rows


def _report(title, lam, n_feas, rows):
    print()
    print(title)
    print(f"  lambda = {lam:.4f}   |F| = {n_feas}")
    print(f"  {'rotation':<28} {'TVD to exp[lam J]':>18} {'ratio spread':>13} "
          f"{'accepted mass':>14}")
    print("  " + "-" * 76)
    for name, t, s, m in rows:
        print(f"  {name:<28} {t:>18.2e} {s:>13.2f} {m:>14.4f}")


if __name__ == "__main__":
    print("=" * 80)
    print("  Does the accepted branch follow exp[lambda J] on the feasible set?")
    print("  ratio spread = max/min of law / exp[lambda J];  1.00 iff Eq. (30) holds")
    print("=" * 80)

    utility = np.array([[6.01, 3.97],
                        [1.10, 8.36],
                        [6.88, 2.21],
                        [2.85, 9.11]])
    _report("inter-cell, 4 UE x 2 AP, cap = 2",
            *check_inter(utility, [1, 1, 1, 1], 2))

    throughput = np.array([[0.42, 0.90, 0.31],
                           [0.71, 0.62, 0.18],
                           [0.83, 0.11, 0.49]])
    _report("intra-cell, 3 node x 3 RB, W = 1",
            *check_intra(throughput, [1, 1, 1]))
