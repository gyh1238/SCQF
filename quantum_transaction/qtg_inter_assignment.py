"""
CQF Inter-cell Assignment (N UE x 2 AP) -- QTG-style feasibility + Grover amplification
=======================================================================================
QTG-style completion of `inter_assignment.py`.

The original file enforces the AP access-limit by *hardcoding* the feasible
patterns and phase-flipping states that match one of them (`oracle_access_limit`).
This version replaces that hardcoded feasibility with a **QTG-style construction**:
the per-AP load is accumulated into a capacity register with QFT arithmetic
(`controlled_add`, identical primitive to `QTG/qtg_knapsack.py`) and the access
limit is checked with an `IntegerComparator`.  Feasibility therefore scales with
the problem (no explicit pattern enumeration) while the CQF objective (link-cost
phase rotation) and Grover amplification are kept intact.

Encoding    : 1 bit/UE.  assign[i] = 0 -> AP0, 1 -> AP1.
Feasibility : load(AP1) = sum_i w_i*assign[i],  load(AP0) = total_w - load(AP1).
              feasible  <=>  load(AP0) <= cap  AND  load(AP1) <= cap
                       <=>  (total_w - cap) <= S <= cap          (two-sided)
              With unit demands (w_i=1) and cap=Ua this reproduces the paper's
              access limit "each AP serves at most Ua UEs".
Objective   : CQF-style phase rotation of link utility into cost qubits;
              superflag marks feasible AND high-utility states (cost all |0>).

Qiskit 1.1 / 2.x compatible.

Author style follows CQF/inter_assignment.py and hybrid_nego/dqna_scale.py.
"""

import numpy as np
from math import ceil, log2
from itertools import product

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit.library import IntegerComparator
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator


# ── QFT arithmetic primitives (shared with QTG/qtg_knapsack.py) ─────

def qft_on(circ, qubits):
    n = len(qubits)
    for j in range(n):
        circ.h(qubits[j])
        for k in range(j + 1, n):
            circ.cp(np.pi / (2 ** (k - j)), qubits[k], qubits[j])
    for j in range(n // 2):
        circ.swap(qubits[j], qubits[n - 1 - j])


def iqft_on(circ, qubits):
    n = len(qubits)
    for j in range(n // 2):
        circ.swap(qubits[j], qubits[n - 1 - j])
    for j in reversed(range(n)):
        for k in reversed(range(j + 1, n)):
            circ.cp(-np.pi / (2 ** (k - j)), qubits[k], qubits[j])
        circ.h(qubits[j])


def controlled_add(circ, ctrl, qubits, a):
    """ctrl=1 -> |x> -> |x + a (mod 2^n)>. qubits: MSB -> LSB."""
    if a == 0:
        return
    qft_on(circ, qubits)
    for j, qj in enumerate(qubits):
        circ.cp(+2 * np.pi * a / (2 ** (j + 1)), ctrl, qj)
    iqft_on(circ, qubits)


def controlled_sub(circ, ctrl, qubits, a):
    """ctrl=1 -> |x> -> |x - a (mod 2^n)>. qubits: MSB -> LSB."""
    if a == 0:
        return
    qft_on(circ, qubits)
    for j, qj in enumerate(qubits):
        circ.cp(-2 * np.pi * a / (2 ** (j + 1)), ctrl, qj)
    iqft_on(circ, qubits)


# ── Feasibility register layout ───────────────────────────────────

def _feasibility_layout(weights, cap):
    """
    Two-sided QTG capacity check on AP1 load S = sum w_i*assign_i.

        AP1 ok  <=>  S <= cap            <=>  NOT (S >= cap + 1)      (hi check)
        AP0 ok  <=>  total_w - S <= cap  <=>  S >= total_w - cap      (lo check)

    Returns register sizes and which checks are non-trivial.
    """
    total_w = int(sum(weights))
    hi_violation = cap + 1              # S >= this  =>  AP1 overflow
    lo_threshold = total_w - cap        # S >= this  =>  AP0 ok

    need_hi = hi_violation <= total_w   # otherwise AP1 can never overflow
    need_lo = lo_threshold > 0          # otherwise AP0 always ok

    n_cap = max(2, ceil(log2(total_w + 1)))

    n_anc_ic = 0
    for need, val in ((need_hi, hi_violation), (need_lo, lo_threshold)):
        if need:
            cmp = IntegerComparator(num_state_qubits=n_cap, value=val, geq=True)
            n_anc_ic = max(n_anc_ic, cmp.num_qubits - (n_cap + 1))

    return {
        "total_w": total_w, "cap": cap,
        "hi_violation": hi_violation, "lo_threshold": lo_threshold,
        "need_hi": need_hi, "need_lo": need_lo,
        "n_cap": n_cap, "n_anc_ic": n_anc_ic,
    }


# ── Oracle: objective phase (CQF link-utility rotation) ────────────

def oracle_objective_phase(qc, assign, cost, utility, sign=+1):
    """
    Per-UE per-AP utility encoding (CQF identity, cf. inter_assignment.py):
        theta = arccos(sqrt(U/U_max));  high utility -> small theta -> cost |0>.
    `sign=-1` applies the inverse rotation (uncompute).
    """
    N, M = utility.shape
    max_u = float(np.max(utility))
    for i in range(N):
        for j in range(M):
            u_norm = utility[i, j] / max_u
            theta = np.arccos(np.sqrt(np.clip(u_norm, 0.0, 1.0)))
            if j == 0:
                qc.x(assign[i])
            qc.cry(sign * theta, assign[i], cost[i])
            if j == 0:
                qc.x(assign[i])


# ── QTG-style feasibility + CQF superflag marking ─────────────────

def apply_qtg_oracle(qc, assign, cost, cap_reg, ge_hi, ge_lo, anc_ic, sf,
                     weights, utility, lay):
    """
    One Grover oracle iteration:
      (1) CQF objective phase rotation into `cost`.
      (2) QTG capacity accumulation  S = sum w_i*assign_i  into `cap_reg`.
      (3) IntegerComparator feasibility flags  ge_hi / ge_lo.
      (4) Superflag phase kickback iff feasible AND high-utility (cost all |0>).
      (5) Uncompute (3),(2),(1) to disentangle all auxiliaries.
    """
    n_cap = lay["n_cap"]

    # (1) objective phase
    oracle_objective_phase(qc, assign, cost, utility, sign=+1)
    qc.barrier()

    # (2) QTG weighted-sum accumulation:  S <- sum w_i on AP1
    for i in range(len(assign)):
        controlled_add(qc, assign[i], cap_reg, int(weights[i]))
    qc.barrier()

    # (3) feasibility comparators (IntegerComparator restores its own ancilla,
    #     so ge_hi and ge_lo can both stay set simultaneously)
    cap_lsb = cap_reg[::-1]                      # LSB-first for IntegerComparator
    if lay["need_hi"]:
        cmp_hi = IntegerComparator(n_cap, lay["hi_violation"], geq=True)
        qc.compose(cmp_hi, cap_lsb + [ge_hi] + list(anc_ic[:lay["n_anc_ic"]]),
                   inplace=True)
    if lay["need_lo"]:
        cmp_lo = IntegerComparator(n_cap, lay["lo_threshold"], geq=True)
        qc.compose(cmp_lo, cap_lsb + [ge_lo] + list(anc_ic[:lay["n_anc_ic"]]),
                   inplace=True)
    qc.barrier()

    # (4) superflag marking.  Fire iff:
    #        cost all |0>  (high utility)  AND  ge_hi == 0 (AP1 ok)  AND  ge_lo == 1 (AP0 ok)
    controls = list(cost)
    qc.x(cost)                                   # cost==0 -> control ready
    if lay["need_hi"]:
        qc.x(ge_hi)                              # want ge_hi==0
        controls = controls + [ge_hi]
    if lay["need_lo"]:
        controls = controls + [ge_lo]            # want ge_lo==1
    qc.mcx(controls, sf[0])
    if lay["need_hi"]:
        qc.x(ge_hi)
    qc.x(cost)
    qc.barrier()

    # (5) uncompute comparators
    if lay["need_lo"]:
        qc.compose(cmp_lo.inverse(), cap_lsb + [ge_lo] + list(anc_ic[:lay["n_anc_ic"]]),
                   inplace=True)
    if lay["need_hi"]:
        qc.compose(cmp_hi.inverse(), cap_lsb + [ge_hi] + list(anc_ic[:lay["n_anc_ic"]]),
                   inplace=True)
    qc.barrier()

    # (5) uncompute weighted sum
    for i in reversed(range(len(assign))):
        controlled_sub(qc, assign[i], cap_reg, int(weights[i]))
    qc.barrier()

    # (5) uncompute objective phase
    oracle_objective_phase(qc, assign, cost, utility, sign=-1)
    qc.barrier()


# ── Diffusion ─────────────────────────────────────────────────────

def apply_diffusion(qc, assign):
    qc.h(assign)
    qc.x(assign)
    qc.h(assign[-1])
    if len(assign) > 1:
        qc.mcx(assign[:-1], assign[-1])
    qc.h(assign[-1])
    qc.x(assign)
    qc.h(assign)
    qc.barrier()


# ── Full circuit builder ──────────────────────────────────────────

def build_circuit(utility, weights, cap, iterations=1):
    utility = np.asarray(utility, dtype=float)
    N, M = utility.shape
    assert M == 2, "binary (2-AP) inter-cell encoding"
    weights = [int(w) for w in weights]

    lay = _feasibility_layout(weights, cap)

    assign = QuantumRegister(N, "assign")
    cost   = QuantumRegister(N, "cost")
    cap_r  = QuantumRegister(lay["n_cap"], "cap")
    ge_hi  = QuantumRegister(1, "ge_hi")
    ge_lo  = QuantumRegister(1, "ge_lo")
    anc_ic = QuantumRegister(max(1, lay["n_anc_ic"]), "anc_ic")
    sf     = QuantumRegister(1, "sf")
    cl     = ClassicalRegister(N, "c")

    qc = QuantumCircuit(assign, cost, cap_r, ge_hi, ge_lo, anc_ic, sf, cl,
                        name="QTG_Inter_Grover")

    # State preparation: uniform superposition + |-> superflag
    qc.h(assign)
    qc.x(sf); qc.h(sf)
    qc.barrier()

    for _ in range(iterations):
        apply_qtg_oracle(qc, assign, cost, cap_r, ge_hi[0], ge_lo[0], anc_ic, sf,
                         weights, utility, lay)
        apply_diffusion(qc, assign)

    qc.measure(assign, cl)
    return qc, lay


# ── Classical reference ───────────────────────────────────────────

def enumerate_feasible(utility, weights, cap):
    N = utility.shape[0]
    total_w = sum(weights)
    out = []
    for bits in product([0, 1], repeat=N):
        s = sum(w for b, w in zip(bits, weights) if b == 1)
        if (total_w - cap) <= s <= cap:
            score = sum(utility[i, bits[i]] for i in range(N))
            out.append((bits, score))
    out.sort(key=lambda x: -x[1])
    return out


# ── Demo: 4 UE x 2 AP ─────────────────────────────────────────────

if __name__ == "__main__":
    # Link utility (higher = better), same numbers as inter_assignment.py.
    utility = np.array([
        [6.01, 3.97],
        [1.10, 8.36],
        [6.88, 2.21],
        [2.85, 9.11],
    ])
    N = utility.shape[0]

    weights = [1, 1, 1, 1]   # unit demand  ->  access limit = "cap UEs per AP"
    cap = 2                  # Ua = 2  ->  feasible = balanced 2-2 splits (6 states)
    iterations = 1

    qc, lay = build_circuit(utility, weights, cap, iterations=iterations)

    print("=" * 60)
    print("  QTG-style Inter-cell Assignment  (4 UE x 2 AP)")
    print(f"  weights={weights}, cap={cap}, total_w={lay['total_w']}")
    print(f"  feasible S range: [{lay['total_w']-cap}, {cap}]  "
          f"(need_hi={lay['need_hi']}, need_lo={lay['need_lo']})")
    print(f"  qubits={qc.num_qubits}, depth={qc.depth()}, iterations={iterations}")
    print("=" * 60)

    # Classical reference
    feas = enumerate_feasible(utility, weights, cap)
    print(f"\n[classical] {len(feas)} feasible states; best = "
          f"{feas[0][0]} score={feas[0][1]:.2f}")

    # Statevector (exact) marginal over assignment qubits
    qc_nm = qc.remove_final_measurements(inplace=False)
    sv = Statevector.from_instruction(qc_nm)
    probs = sv.probabilities_dict()
    assign_p = {}
    for bs, p in probs.items():
        key = bs[-N:]                       # last N chars = assign register
        assign_p[key] = assign_p.get(key, 0.0) + p

    ranked = sorted(assign_p.items(), key=lambda kv: -kv[1])
    print("\n[quantum] top assignments (bitstring q[N-1]..q[0]):")
    print(f"  {'bits':>8}  {'prob':>8}  {'assign(UE0..)':>14}  {'score':>7}  feas")
    print("  " + "-" * 52)
    for key, p in ranked[:8]:
        ue_bits = [int(b) for b in reversed(key)]   # UE0..UE_{N-1}
        s = sum(weights[i] for i in range(N) if ue_bits[i] == 1)
        feasible = (lay["total_w"] - cap) <= s <= cap
        score = sum(utility[i, ue_bits[i]] for i in range(N))
        print(f"  {key:>8}  {p:>8.4f}  {str(ue_bits):>14}  {score:>7.2f}  "
              f"{'ok' if feasible else 'X'}")

    feasible_mass = sum(
        p for key, p in assign_p.items()
        if (lay["total_w"] - cap)
        <= sum(weights[i] for i in range(N) if int(key[::-1][i]) == 1)
        <= cap)
    print(f"\n  total feasible probability mass = {feasible_mass:.4f}")

    # Shot-based confirmation
    backend = AerSimulator()
    tqc = transpile(qc, backend)
    counts = backend.run(tqc, shots=10000).result().get_counts()
    top = max(counts.items(), key=lambda kv: kv[1])
    print(f"  [shots] most frequent = {top[0]} ({top[1]/10000:.3f})")
