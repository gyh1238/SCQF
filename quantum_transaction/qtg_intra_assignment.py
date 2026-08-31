"""
CQF Intra-cell RB Assignment (N node x M RB, 2-bit encoding)
    -- QTG-style feasibility + Grover amplification
================================================================================
QTG-style completion of `intra_assignment.py`.

The original file enforces the resource-block limit with an O(N^2) pairwise
comparison oracle (`compare_nodes` over every node pair) plus a state-violation
oracle.  This version replaces the pairwise RB limit with a **QTG-style
occupancy count**: for each RB the number of nodes assigned to it is accumulated
into a counter with QFT arithmetic (`controlled_add`) and checked against W_r
with an `IntegerComparator` -- the same primitives as `QTG/qtg_knapsack.py`.
This is O(M * N) and generalises to arbitrary per-RB capacities W_r, whereas the
pairwise oracle only expresses W_r = 1.

Because two constraints (state-validity and RB-limit) can be violated *at the
same time*, a parity-based phase kickback is unsafe.  Following the CQF identity
(cf. `intra_assignment.py`), all constraint flags are computed into ancillae and
a **single superflag MCX** fires only when every constraint is satisfied AND the
objective is high -- correct regardless of how many constraints are violated.

Encoding    : 2 bits/node.  code 00->RB0, 01->RB1, 10->RB2, 11->invalid.
Feasibility : (a) state-validity  -- no node holds the invalid code 11
              (b) RB-limit        -- occupancy(RB r) <= W_r  for every r
              With W_r = 1 and N = M = 3 the feasible set is the 6 permutations.
Objective   : exponential utility weight of Eq. (19)-(21) rotated into the
              cost qubits, g = exp[-lambda(w_i^max - w_ir)], theta = 2 arccos sqrt(g),
              so the accepted branch carries exp[lambda J] over the feasible
              set, Eq. (30); superflag marks feasible AND cost all |0>.

Qiskit 1.1 / 2.x compatible.

Author style follows CQF/intra_assignment.py and hybrid_nego/dqna_scale.py.
"""

import numpy as np
from math import ceil, log2
from itertools import product

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit.library import IntegerComparator
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator

INVALID_CODE = 3   # 11 -> no such RB

# Target exponent of Sec. V-A; lambda = BETA / utility_scale(w).
BETA = 1.5


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


# ── 2-bit code helpers ────────────────────────────────────────────

def _pair(assign, i):
    """(LSB, MSB) qubits of node i.  value = LSB + 2*MSB."""
    return assign[2 * i], assign[2 * i + 1]


def _mask_pair(qc, lsb, msb, r):
    """X the pair bits that must be 0 so that MCX(pair) fires iff code == r."""
    if (r & 1) == 0:
        qc.x(lsb)
    if ((r >> 1) & 1) == 0:
        qc.x(msb)


# ── Oracle: objective phase (CQF throughput rotation, 2-bit) ───────

def utility_scale(throughput):
    """
    The utility model's scale `ubar`, so that `lambda = beta / ubar` is
    dimensionless and one `beta` means the same thing in every zone.  Defined
    from the model's dynamic range -- the mean of the per-node best throughput
    -- rather than from an instance, as `haiq_instance.utility_scale`.
    """
    return float(np.mean(np.max(np.asarray(throughput, dtype=float), axis=1)))


def oracle_objective_phase(qc, assign, cost, match_anc, throughput, lam, sign=+1):
    """
    Per-node per-RB throughput encoding, Eq. (19)-(21):

        g_{i,r} = exp[-lambda (w_i^max - w_{i,r})],
        theta_{i,r} = 2 arccos sqrt(g_{i,r}),

    so the all-zero cost event of a complete assignment carries weight
    kappa * exp[lambda J], Eq. (24).  Only an exponential weight turns the
    additive utility of Eq. (13) into a product across zones, Eq. (15); the
    linear form `arccos sqrt(w/w_max)` used by the enumeration oracle of
    [jang2025cqf] preserves the per-node ranking but not that identity.

    A node holding the invalid code 11 matches no RB and leaves cost |0>; the
    state-validity flag (not the objective) is what rejects it.
    """
    N, M = throughput.shape
    node_iter = range(N) if sign > 0 else reversed(range(N))
    for i in node_iter:
        lsb, msb = _pair(assign, i)
        w_max = float(np.max(throughput[i]))
        rb_iter = range(M) if sign > 0 else reversed(range(M))
        for r in rb_iter:
            g = np.exp(-lam * (w_max - throughput[i, r]))
            theta = 2.0 * np.arccos(np.sqrt(np.clip(g, 0.0, 1.0)))
            if theta == 0.0:       # best candidate: unit weight, no rotation
                continue
            _mask_pair(qc, lsb, msb, r)
            qc.mcx([lsb, msb], match_anc)
            qc.cry(sign * theta, match_anc, cost[i])
            qc.mcx([lsb, msb], match_anc)
            _mask_pair(qc, lsb, msb, r)


# ── QTG-style feasibility flags (self-inverse) ────────────────────

def feasibility_flags(qc, assign, anc_state, rb_viol, cnt, match_anc, ge, anc_ic,
                      N, M, W):
    """
    Compute constraint-violation flags.  Self-inverse: calling twice restores
    every ancilla to |0>, so it is used both to compute the flags (before the
    superflag) and to uncompute them (after).

        anc_state[i] ^= 1   iff node i holds the invalid code 11
        rb_viol[r]   ^= 1   iff occupancy(RB r) >= W_r + 1   (over capacity)

    RB occupancy is accumulated QTG-style into `cnt` and compared with an
    IntegerComparator; `cnt`, `ge`, `anc_ic`, `match_anc` are all restored to
    |0> internally so the counter register is reused across RBs.
    """
    # (a) state-validity: node == 11  <=>  both pair bits set
    for i in range(N):
        lsb, msb = _pair(assign, i)
        qc.mcx([lsb, msb], anc_state[i])
    qc.barrier()

    # (b) RB-limit: occupancy(RB r) counted into `cnt`, checked against W_r
    n_cnt = len(cnt)
    n_anc = len(anc_ic)
    for r in range(M):
        # count nodes on RB r
        for i in range(N):
            lsb, msb = _pair(assign, i)
            _mask_pair(qc, lsb, msb, r)
            qc.mcx([lsb, msb], match_anc)          # match_anc = (node i on RB r)
            controlled_add(qc, match_anc, cnt, 1)  # cnt += 1
            qc.mcx([lsb, msb], match_anc)          # uncompute match_anc
            _mask_pair(qc, lsb, msb, r)

        # overflow?  cnt >= W_r + 1  ->  ge = 1  ->  record in rb_viol[r]
        # `cnt` is stored MSB-first (controlled_add convention); IntegerComparator
        # expects LSB-first, so reverse it.
        cmp = IntegerComparator(n_cnt, W[r] + 1, geq=True)
        cmp_qubits = list(cnt)[::-1] + [ge] + list(anc_ic[:cmp.num_qubits - n_cnt - 1])
        qc.compose(cmp, cmp_qubits, inplace=True)
        qc.cx(ge, rb_viol[r])
        qc.compose(cmp.inverse(), cmp_qubits, inplace=True)

        # uncompute the counter so it is |0> for the next RB
        for i in reversed(range(N)):
            lsb, msb = _pair(assign, i)
            _mask_pair(qc, lsb, msb, r)
            qc.mcx([lsb, msb], match_anc)
            controlled_sub(qc, match_anc, cnt, 1)
            qc.mcx([lsb, msb], match_anc)
            _mask_pair(qc, lsb, msb, r)
    qc.barrier()


# ── Combined oracle ───────────────────────────────────────────────

def apply_qtg_oracle(qc, assign, cost, anc_state, rb_viol, cnt, match_anc,
                     ge, anc_ic, sf, throughput, N, M, W, lam):
    # (1) objective phase into cost
    oracle_objective_phase(qc, assign, cost, match_anc, throughput, lam, sign=+1)
    qc.barrier()

    # (2) constraint-violation flags (QTG occupancy counts)
    feasibility_flags(qc, assign, anc_state, rb_viol, cnt, match_anc, ge, anc_ic,
                      N, M, W)

    # (3) superflag: fire iff cost all |0> (high throughput)
    #                  AND anc_state all |0> (all codes valid)
    #                  AND rb_viol   all |0> (all RBs within W_r)
    conds = list(cost) + list(anc_state) + list(rb_viol)
    qc.x(conds)                       # every condition qubit must be |0>
    qc.mcx(conds, sf[0])
    qc.x(conds)
    qc.barrier()

    # (4) uncompute flags (self-inverse call) and objective phase
    feasibility_flags(qc, assign, anc_state, rb_viol, cnt, match_anc, ge, anc_ic,
                      N, M, W)
    oracle_objective_phase(qc, assign, cost, match_anc, throughput, lam, sign=-1)
    qc.barrier()


# ── Diffusion ─────────────────────────────────────────────────────

def apply_diffusion(qc, assign):
    qc.h(assign)
    qc.x(assign)
    qc.h(assign[-1])
    qc.mcx(assign[:-1], assign[-1])
    qc.h(assign[-1])
    qc.x(assign)
    qc.h(assign)
    qc.barrier()


# ── Full circuit builder ──────────────────────────────────────────

def build_circuit(throughput, W=None, iterations=1, beta=BETA, lam=None):
    throughput = np.asarray(throughput, dtype=float)
    if lam is None:
        lam = beta / utility_scale(throughput)
    N, M = throughput.shape
    assert M <= 3, "2-bit encoding supports up to 3 RBs (codes 00/01/10)"
    if W is None:
        W = [1] * M

    n_state = 2 * N
    n_cnt = max(1, ceil(log2(N + 1)))
    n_anc_ic = IntegerComparator(n_cnt, 1, geq=True).num_qubits - n_cnt - 1

    assign    = QuantumRegister(n_state, "assign")
    cost      = QuantumRegister(N, "cost")
    anc_state = QuantumRegister(N, "anc_state")
    rb_viol   = QuantumRegister(M, "rb_viol")
    cnt       = QuantumRegister(n_cnt, "cnt")
    match_anc = QuantumRegister(1, "match_anc")
    ge        = QuantumRegister(1, "ge")
    anc_ic    = QuantumRegister(max(1, n_anc_ic), "anc_ic")
    sf        = QuantumRegister(1, "sf")
    cl        = ClassicalRegister(n_state, "c")

    qc = QuantumCircuit(assign, cost, anc_state, rb_viol, cnt, match_anc,
                        ge, anc_ic, sf, cl, name="QTG_Intra_Grover")

    # State preparation: uniform superposition + |-> superflag
    qc.h(assign)
    qc.x(sf); qc.h(sf)
    qc.barrier()

    for _ in range(iterations):
        apply_qtg_oracle(qc, assign, cost, anc_state, rb_viol, cnt, match_anc,
                         ge[0], anc_ic, sf, throughput, N, M, W, lam)
        apply_diffusion(qc, assign)

    qc.measure(assign, cl)
    return qc


# ── Sampler of Sec. IV-F (the circuit the evaluation measures) ────

def build_sampler(throughput, W=None, k=1, beta=BETA, lam=None, measure=True):
    """
    The zone sampler of Eq. (26)-(36), as opposed to the Grover demo above.

    `build_circuit` keeps the historical arrangement: the utility rotations
    live inside the oracle and are mirrored away again, and the reflection is
    the textbook one about the uniform state.  Section IV places them
    differently:

        P_z = Theta_z H^(x)Q_z                                   Eq. (26)
        O_z = E_z^dagger Lambda_z E_z                            Eq. (32)
        D_z = P_z (2|0><0| - I) P_z^dagger                       Eq. (33)
        C_z = G_z^k P_z,   G_z = D_z O_z                         Eq. (35)

    The rotations belong to preparation and stay standing through the
    measurement, because the acceptance event of Sec. IV-F is *all utility
    qubits in |0>*.  The reflection is then about P_z|0>, not about the
    uniform state.

    Returns (qc, meta) with the offsets a measurement needs to decode.
    """
    throughput = np.asarray(throughput, dtype=float)
    if lam is None:
        lam = beta / utility_scale(throughput)
    N, M = throughput.shape
    assert M <= 3, "2-bit encoding supports up to 3 RBs (codes 00/01/10)"
    if W is None:
        W = [1] * M

    n_state = 2 * N
    n_cnt = max(1, ceil(log2(N + 1)))
    n_anc_ic = IntegerComparator(n_cnt, 1, geq=True).num_qubits - n_cnt - 1

    assign    = QuantumRegister(n_state, "assign")
    cost      = QuantumRegister(N, "cost")
    anc_state = QuantumRegister(N, "anc_state")
    rb_viol   = QuantumRegister(M, "rb_viol")
    cnt       = QuantumRegister(n_cnt, "cnt")
    match_anc = QuantumRegister(1, "match_anc")
    ge        = QuantumRegister(1, "ge")
    anc_ic    = QuantumRegister(max(1, n_anc_ic), "anc_ic")
    sf        = QuantumRegister(1, "sf")
    regs = [assign, cost, anc_state, rb_viol, cnt, match_anc, ge, anc_ic, sf]
    if measure:
        regs.append(ClassicalRegister(n_state + N, "c"))   # assign, then cost
    qc = QuantumCircuit(*regs, name="CQF_Intra_Sampler")

    def prep(inverse=False):
        if inverse:
            oracle_objective_phase(qc, assign, cost, match_anc[0], throughput,
                                   lam, sign=-1)
            qc.h(assign)
        else:
            qc.h(assign)
            oracle_objective_phase(qc, assign, cost, match_anc[0], throughput,
                                   lam, sign=+1)

    def mark():
        conds = list(cost) + list(anc_state) + list(rb_viol)
        qc.x(conds)
        qc.mcx(conds, sf[0])
        qc.x(conds)

    qc.x(sf); qc.h(sf)
    prep()
    qc.barrier()
    for _ in range(k):
        # O_z = E_z^dagger Lambda_z E_z; feasibility_flags is self-inverse
        feasibility_flags(qc, assign, anc_state, rb_viol, cnt, match_anc[0],
                          ge[0], anc_ic, N, M, W)
        mark()
        feasibility_flags(qc, assign, anc_state, rb_viol, cnt, match_anc[0],
                          ge[0], anc_ic, N, M, W)
        qc.barrier()
        prep(inverse=True)                                # D_z
        qc.x(list(assign) + list(cost))
        qc.h(cost[-1])
        qc.mcx(list(assign) + list(cost)[:-1], cost[-1])
        qc.h(cost[-1])
        qc.x(list(assign) + list(cost))
        prep()
        qc.barrier()
    if measure:
        qc.measure(list(assign) + list(cost), range(n_state + N))

    meta = dict(n_assign=n_state, n_cost=N, lam=float(lam), k=int(k),
                N=N, M=M, W=list(W))
    return qc, meta


# ── Classical reference ───────────────────────────────────────────

def decode(key, N):
    """key = qiskit bitstring (q[-1]..q[0]) of assign register -> node codes."""
    bits = [int(b) for b in reversed(key)]        # bits[k] = q[k]
    return [bits[2 * i] + 2 * bits[2 * i + 1] for i in range(N)]  # node0..node_{N-1}


def is_feasible(codes, M, W):
    if any(c == INVALID_CODE for c in codes):
        return False
    occ = [codes.count(r) for r in range(M)]
    return all(occ[r] <= W[r] for r in range(M))


def enumerate_feasible(throughput, W):
    N, M = throughput.shape
    out = []
    for codes in product(range(4), repeat=N):
        if is_feasible(codes, M, W):
            score = sum(throughput[i, codes[i]] for i in range(N))
            out.append((codes, score))
    out.sort(key=lambda x: -x[1])
    return out


# ── Demo: 3 nodes x 3 RBs ─────────────────────────────────────────

if __name__ == "__main__":
    # Link throughput (higher = better), same numbers as intra_assignment.py.
    throughput = np.array([
        [0.42, 0.90, 0.31],
        [0.71, 0.62, 0.18],
        [0.83, 0.11, 0.49],
    ])
    N, M = throughput.shape
    W = [1, 1, 1]            # each RB serves at most one node  ->  permutations
    iterations = 1

    qc = build_circuit(throughput, W=W, iterations=iterations)

    print("=" * 60)
    print("  QTG-style Intra-cell RB Assignment  (3 node x 3 RB)")
    print(f"  W (per-RB limit) = {W}")
    print(f"  qubits={qc.num_qubits}, depth={qc.depth()}, iterations={iterations}")
    print("=" * 60)

    feas = enumerate_feasible(throughput, W)
    print(f"\n[classical] {len(feas)} feasible states; best = "
          f"{feas[0][0]} score={feas[0][1]:.2f}")

    qc_nm = qc.remove_final_measurements(inplace=False)
    sv = Statevector.from_instruction(qc_nm)
    probs = sv.probabilities_dict()
    assign_p = {}
    for bs, p in probs.items():
        key = bs[-2 * N:]
        assign_p[key] = assign_p.get(key, 0.0) + p

    ranked = sorted(assign_p.items(), key=lambda kv: -kv[1])
    print("\n[quantum] top assignments:")
    print(f"  {'bits':>8}  {'prob':>8}  {'codes(node0..)':>16}  {'score':>7}  feas")
    print("  " + "-" * 54)
    for key, p in ranked[:10]:
        codes = decode(key, N)
        feasible = is_feasible(codes, M, W)
        score = sum(throughput[i, codes[i]] for i in range(N)) if all(
            c != INVALID_CODE for c in codes) else float("nan")
        rb_names = [("RB0", "RB1", "RB2", "X!")[c] for c in codes]
        print(f"  {key:>8}  {p:>8.4f}  {str(rb_names):>16}  {score:>7.2f}  "
              f"{'ok' if feasible else 'X'}")

    feasible_mass = sum(p for key, p in assign_p.items()
                        if is_feasible(decode(key, N), M, W))
    print(f"\n  total feasible probability mass = {feasible_mass:.4f}")

    backend = AerSimulator()
    tqc = transpile(qc, backend)
    counts = backend.run(tqc, shots=10000).result().get_counts()
    top = max(counts.items(), key=lambda kv: kv[1])
    print(f"  [shots] most frequent = {top[0]} -> "
          f"codes {decode(top[0], N)} ({top[1]/10000:.3f})")
