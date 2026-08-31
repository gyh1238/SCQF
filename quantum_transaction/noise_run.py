"""
The zone sampler on a noisy device (Sec. V-B).
==============================================
`certify_rotation.py` shows the accepted branch follows exp[lambda J_z] on
F_z to statevector precision.  That is a statement about a perfect device.
This module asks what a real one costs.

Two questions are kept apart, because they have different answers and only
one of them is a simulation:

  * **What does the lattice cost?**  Routing this circuit onto a heavy-hex
    coupling map inserts swaps, and the two-qubit count is what the executable
    size is set by.  This is settled by transpiling and counting -- there is
    nothing to simulate, and a 127-qubit device could not be simulated anyway.
    `routing_cost` reports it.
  * **What do noisy gates cost?**  This one is a simulation, and it is run on
    the circuit at all-to-all connectivity, where the simulated qubit count is
    the circuit's own -- sixteen, not the device's hundred and twenty-seven.
    `noisy_law` reports it.

Mixing them produces a number that answers neither: a figure measured on the
routed circuit cannot say whether the damage came from the gates being noisy
or from the lattice being sparse, and the routed circuit is exactly the one
too large to simulate.  Reported separately the two still compose, because
the depolarizing sweep is parameterized by the *product* of two-qubit error
and gate count: the routed cost is read off the sweep at the routed count.

What noise can and cannot break here is worth stating precisely, because the
paper's headline claim is *strict* feasibility:

  * It cannot put an infeasible assignment into a report.  The acceptance
    test of Sec. IV-F is applied to the **decoded assignment**, classically,
    after the shot comes back -- not to the internal flags.  A corrupted flag
    changes which branch the superflag marked, and therefore how often a shot
    survives, but a decoded assignment that violates an owned limit is thrown
    away whatever the flags did.  Leakage into a report is therefore zero by
    construction, and the quantity worth measuring is what that costs:
    `infeas` below is the share of shots that pass the cost test and are then
    discarded for decoding to an infeasible assignment.  It is *not* zero
    without noise -- an infeasible assignment can perfectly well show all its
    utility qubits in |0> -- so the noise-free value is the baseline the noisy
    one is read against.
  * It can make acceptance rarer, which is a cost in executions.
  * It can bend the accepted law away from exp[lambda J_z], which is a cost in
    the fidelity of the report the coordination stage consumes.

The noise-free distance is measured at the same shot count and printed beside
the noisy one: they are two measurements of the same quantity under two
devices, not an additive decomposition -- finite sampling and noise do not add.

Run:  python noise_run.py [--shots 10000] [--backend FakeSherbrooke] [--sweep]
"""

import argparse
import time

import numpy as np
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

import qtg_inter_assignment as inter
import qtg_intra_assignment as intra
from certify_rotation import gibbs, tvd

SHOTS = 10_000
BACKEND = "FakeSherbrooke"
BASIS = ["cx", "rz", "sx", "x"]


# ----------------------------------------------------- what the lattice costs

def routing_cost(qc, backend_name=BACKEND):
    """Two-qubit gates before and after the device's coupling map is imposed.

    Transpiling is cheap and exact; this is a gate count, not a simulation.
    The ratio is what the partition budget has to absorb.
    """
    from qiskit_ibm_runtime import fake_provider
    backend = getattr(fake_provider, backend_name)()

    flat = transpile(qc, basis_gates=BASIS, optimization_level=1)
    routed = transpile(qc, backend=backend, optimization_level=1,
                       seed_transpiler=7)

    def n2q(c):
        return sum(v for g, v in c.count_ops().items()
                   if g in ("cx", "ecr", "cz"))

    return dict(flat_2q=n2q(flat), routed_2q=n2q(routed),
                flat_depth=flat.depth(), routed_depth=routed.depth(),
                factor=n2q(routed) / max(n2q(flat), 1),
                device=backend_name,
                device_qubits=getattr(backend, "num_qubits", None))


# ----------------------------------------------------- what noisy gates cost

def depolarizing(p2, p1=None, p_ro=None):
    """A one-knob model: `p2` on two-qubit gates, the rest following it.

    One-qubit gates at `p2/10` and readout at `p2/5` is the ordering
    superconducting hardware actually has, so a sweep in `p2` moves the whole
    device along a plausible line rather than an arbitrary one.  Same
    convention as `haiq_noise.depolarizing`.
    """
    p1 = p2 / 10.0 if p1 is None else p1
    p_ro = p2 / 5.0 if p_ro is None else p_ro
    m = NoiseModel(basis_gates=BASIS)
    if p1 > 0:
        m.add_all_qubit_quantum_error(depolarizing_error(p1, 1), ["sx", "x"])
    if p2 > 0:
        m.add_all_qubit_quantum_error(depolarizing_error(p2, 2), ["cx"])
    if p_ro > 0:
        m.add_all_qubit_readout_error([[1 - p_ro, p_ro], [p_ro, 1 - p_ro]])
    return m


def device_rates(backend_name=BACKEND):
    """Median measured error rates of a calibrated backend.

    Taken as an all-qubit model rather than per-qubit: the circuit is
    simulated at all-to-all connectivity, so there is no physical qubit for a
    per-qubit rate to attach to.  What survives is the magnitude of the
    device's errors, which is the part this question needs; the spread across
    the chip is the part `routing_cost` deliberately does not mix in.
    """
    from qiskit_ibm_runtime import fake_provider
    target = getattr(fake_provider, backend_name)().target

    def median(names):
        vals = []
        for n in names:
            if n not in target:
                continue
            vals += [p.error for p in target[n].values()
                     if p is not None and p.error is not None and p.error > 0]
        return float(np.median(vals)) if vals else 0.0

    return median(["ecr", "cz", "cx"]), median(["sx", "x"]), median(["measure"])


def noisy_law(qc, meta, feasible, noise, shots, seed=11):
    """Accepted law, acceptance, discard rate and leak, at `shots` shots."""
    tqc = transpile(qc, basis_gates=BASIS, optimization_level=1)
    sim = AerSimulator(noise_model=noise)
    counts = sim.run(tqc, shots=shots, seed_simulator=seed).result().get_counts()

    n_assign, n_cost = meta["n_assign"], meta["n_cost"]
    law, n_acc, n_leak, n_tot = {}, 0, 0, 0
    for bits, c in counts.items():
        bits = bits.replace(" ", "")
        n_tot += c
        if "1" in bits[::-1][n_assign:n_assign + n_cost]:   # cost not all |0>
            continue
        code = bits[-n_assign:]
        if code not in feasible:            # decoded assignment breaks a limit
            n_leak += c                     # discarded, never reaches a report
            continue
        law[code] = law.get(code, 0) + c
        n_acc += c
    if n_acc:
        law = {k: v / n_acc for k, v in law.items()}
    ref = gibbs(feasible, meta["lam"])
    return dict(accept=n_acc / n_tot, discard=1 - n_acc / n_tot,
                infeas=n_leak / n_tot,
                tvd=tvd(law, ref) if law else float("nan"))


# --------------------------------------------------------------- the problems

def cases(k=1):
    """The two problems the version-fixed files carry, as zone samplers."""
    utility = np.array([[6.01, 3.97],
                        [1.10, 8.36],
                        [6.88, 2.21],
                        [2.85, 9.11]])
    w, cap = [1, 1, 1, 1], 2
    qc_i, meta_i = inter.build_sampler(utility, w, cap, k=k)
    feas_i = {"".join(str(b) for b in reversed(bits)): s
              for bits, s in inter.enumerate_feasible(utility, w, cap)}

    throughput = np.array([[0.42, 0.90, 0.31],
                           [0.71, 0.62, 0.18],
                           [0.83, 0.11, 0.49]])
    W = [1, 1, 1]
    qc_r, meta_r = intra.build_sampler(throughput, W=W, k=k)
    feas_r = {}
    for codes, s in intra.enumerate_feasible(throughput, W):
        b = "".join(f"{c & 1}{(c >> 1) & 1}" for c in codes)
        feas_r[b[::-1]] = s

    return [("inter-cell, 4 UE x 2 AP", qc_i, meta_i, feas_i),
            ("intra-cell, 3 node x 3 RB", qc_r, meta_r, feas_r)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shots", type=int, default=SHOTS)
    ap.add_argument("--backend", default=BACKEND)
    ap.add_argument("--k", type=int, default=1)
    ap.add_argument("--sweep", action="store_true",
                    help="also sweep the two-qubit error rate")
    ap.add_argument("--only", choices=("inter", "intra"),
                    help="run one case (the intra circuit is much the slower)")
    a = ap.parse_args()

    p2, p1, ro = device_rates(a.backend)
    print("=" * 78)
    print(f"  Zone sampler, {a.shots:,} shots, k={a.k}")
    print(f"  {a.backend} median rates: 2q {p2:.2e}, 1q {p1:.2e}, "
          f"readout {ro:.2e}")
    print("=" * 78, flush=True)

    todo = [c for c in cases(a.k)
            if not a.only or c[0].startswith(a.only)]
    for name, qc, meta, feas in todo:
        rc = routing_cost(qc, a.backend)
        print()
        print(f"{name}   |F| = {len(feas)}   lambda = {meta['lam']:.4f}   "
              f"{qc.num_qubits} qubits")
        print(f"  lattice cost, transpile only: {rc['flat_2q']} -> "
              f"{rc['routed_2q']} two-qubit gates ({rc['factor']:.1f}x) on "
              f"{rc['device']} ({rc['device_qubits']} qubits), depth "
              f"{rc['flat_depth']} -> {rc['routed_depth']}")
        print(f"  {'':<24}{'accept':>9}{'discard':>9}"
              f"{'TVD to exp[lam J]':>20}{'infeas':>8}", flush=True)

        rows = [("noise-free", None),
                (f"{a.backend} rates", depolarizing(p2, p1, ro))]
        if a.sweep:
            rows += [(f"depolarizing p2={q:.0e}", depolarizing(q))
                     for q in (1e-4, 1e-3, 1e-2)]

        base = None
        for tag, model in rows:
            t = time.time()
            r = noisy_law(qc, meta, feas, model, a.shots)
            share = "" if base is None else \
                f"   {100 * r['accept'] / base:.0f}% of noise-free"
            if base is None:
                base = r["accept"]
            print(f"  {tag:<24}{r['accept']:>9.3f}{r['discard']:>9.3f}"
                  f"{r['tvd']:>20.3f}{r['infeas']:>8.4f}"
                  f"   [{time.time() - t:.0f}s]{share}", flush=True)


if __name__ == "__main__":
    main()
