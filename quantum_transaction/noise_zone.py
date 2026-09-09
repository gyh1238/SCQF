"""
The zone circuit on a noisy device.
===================================
Which of the two noise modules is this?  `noise_oracle.py` asks the same question
of the released `circuit_inter`/`circuit_intra` oracles and is what Sec. V-B
quotes.
This one asks it of the `certify_sampler` zone circuit, one stack earlier, and is
the only place the leak metric is posed against that circuit directly.  It is
imported by nothing; run it on its own.

`certify_sampler` proves the classical sampler reproduces the circuit's accepted
branch exactly, on a noiseless statevector.  That is what licenses the figures,
and it stays: an exact identity is not something a noisy run can establish.

This module asks the other question.  The protocol's selling point is *strict*
feasibility -- an assignment outside F_z carries exactly zero accepted
probability, because the circuit measures the limit flags and rejects any
branch that violates one.  That guarantee is a statement about a perfect
device.  On a real one a flag can be flipped by noise after the violation was
recorded, and an infeasible assignment then walks through the acceptance test.
So the quantity to measure is not fidelity, it is **leak**: how much accepted
probability lands on assignments that break an owned limit.

Three numbers come back, all exact:

    mu     acceptance probability -- what fraction of runs survive the test
    tvd    total-variation distance between the accepted law and exp(lambda J_z)
    leak   accepted probability mass on infeasible assignments; zero by
           construction without noise, and the thing strict feasibility loses

They are exact because nothing here is sampled.  The circuit is evolved as a
density matrix, which carries the full noisy state, and readout error is then
applied to its diagonal as the confusion matrix it is.  Shot noise would only
blur what the noise model already says.  The price is the density matrix
itself: 2^(2n) amplitudes, so this reaches about 12 qubits, which is zones of
two or three UEs.  `MAX_QUBITS` is where that cap lives.

Two noise sources:

  * `depolarizing(p2)` -- a sweep knob.  One number, with the one-qubit rate
    at a tenth of it and readout at a fifth, which is the usual ordering on
    superconducting hardware.  This is what answers "how good does the device
    have to be".
  * `device(name)` -- a real backend's calibration through
    `qiskit_ibm_runtime.fake_provider`: measured T1, T2, gate and readout
    errors, no connectivity constraint.  This answers "what would today's
    hardware give if it were fully connected".  Routing a circuit this shape
    onto a heavy-hex lattice costs more two-qubit gates than the circuit
    contains, so keeping the two effects apart is the only way either number
    means anything.

Usage:  python noise_zone.py
"""

import numpy as np

MAX_QUBITS = 12          # density matrix is 2^(2n); 12 qubits is 268 MB
BASIS = ["u", "cx"]      # the sweep's basis; a device brings its own


def depolarizing(p2, p1=None, p_ro=None):
    """A one-knob noise model: `p2` on two-qubit gates, the rest follows it.

    One-qubit gates at `p2/10` and readout at `p2/5` is the ordering
    superconducting hardware actually has, so a sweep in `p2` moves the whole
    device along a plausible line rather than an arbitrary one.
    """
    from qiskit_aer.noise import NoiseModel, depolarizing_error
    p1 = p2 / 10.0 if p1 is None else p1
    model = NoiseModel(basis_gates=BASIS)
    model.add_all_qubit_quantum_error(depolarizing_error(p1, 1), ["u"])
    model.add_all_qubit_quantum_error(depolarizing_error(p2, 2), ["cx"])
    model.readout = p2 / 5.0 if p_ro is None else p_ro   # applied separately
    return model


def device(name="FakeSherbrooke"):
    """A real backend's calibration, minus its connectivity.

    The errors are measured ones -- T1, T2, per-gate and per-qubit readout --
    so the magnitudes are what a run today would meet.  The coupling map is
    deliberately not imposed: routing this circuit onto a heavy-hex lattice
    multiplies its two-qubit count several times over, and a number that mixes
    "the gates are noisy" with "the lattice is sparse" answers neither
    question.  The second effect is reported separately by `routing_cost`.
    """
    from qiskit_aer.noise import NoiseModel
    from qiskit_ibm_runtime import fake_provider
    backend = getattr(fake_provider, name)()
    model = NoiseModel.from_backend(backend)
    model.readout = None            # already carried by the backend model
    model.backend_name = name
    return model


def _readout_confusion(diag, n_qubits, p_ro):
    """Apply a symmetric readout flip to an exact outcome distribution.

    The density matrix carries everything the gates did; readout error happens
    after it, on the classical bits.  Sampling shots to model that would add
    noise to a calculation that has none, so the confusion matrix is applied
    where it belongs -- as the tensor product of 2x2 flips it is, one qubit at
    a time, which costs `n` passes over the distribution rather than a matrix
    of size 4^n.
    """
    if not p_ro:
        return diag
    out = diag.reshape([2] * n_qubits)
    for axis in range(n_qubits):
        out = np.moveaxis(out, axis, 0)
        stay, flip = (1.0 - p_ro) * out, p_ro * out[::-1]
        out = np.moveaxis(stay + flip, 0, axis)
    return out.ravel()


def noisy_zone_law(zone, lam, noise, max_qubits=MAX_QUBITS):
    """`mu`, `tvd` and `leak` for one zone under one noise model, exactly.

    Returns None if the zone's circuit is past `max_qubits`, since the density
    matrix would not fit; the caller is expected to say so rather than to
    quietly show a smaller zone than it claimed.
    """
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    from certify_sampler import build_zone_circuit
    from proto_zone import enumerate_zone

    qc, meta = build_zone_circuit(zone, lam)
    if qc.num_qubits > max_qubits:
        return None

    basis = getattr(noise, "basis_gates", BASIS)
    circ = transpile(qc, basis_gates=list(basis), optimization_level=1)
    n_cx = sum(c for g, c in circ.count_ops().items() if g in ("cx", "ecr", "cz"))
    circ.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=noise)
    rho = np.asarray(sim.run(circ, shots=1).result().data()["density_matrix"])
    diag = np.clip(np.real(np.diag(rho)), 0.0, None)
    diag = diag / diag.sum()
    diag = _readout_confusion(diag, qc.num_qubits, getattr(noise, "readout", None))

    # the same split `certify_sampler` makes: code register low, everything the
    # acceptance test looks at above it, and all of that has to be zero
    n = meta["n"]
    n_aux = meta["w_cnt"] + meta["n_rb"] + meta["n_ap"] + meta["n_val"]
    idx = np.arange(diag.size, dtype=np.int64)
    accepted = (idx >> n) & ((1 << (n + n_aux)) - 1)
    keep = accepted == 0

    hist = np.bincount(idx[keep] & ((1 << n) - 1), weights=diag[keep],
                       minlength=1 << n)
    mu = float(hist.sum())
    law = hist / mu if mu > 0 else hist

    rows, pref, _ = enumerate_zone(zone, lam)
    ref = np.zeros(1 << n)
    for row, p in zip(rows, pref):
        ref[int(sum(int(c) << k for k, c in enumerate(row)))] = p

    return dict(mu=mu, tvd=0.5 * float(np.abs(law - ref).sum()),
                leak=float(law[ref == 0].sum()), n_qubits=qc.num_qubits,
                n_cx=int(n_cx), n_ue=zone.n_ue, n_feasible=len(rows))


def routing_cost(zone, lam, name="FakeSherbrooke"):
    """What the device's own lattice adds, counted separately from its noise.

    The zone circuit is written for a machine that can entangle any pair.  A
    heavy-hex processor cannot, and the router pays for that in swaps.  This
    is the multiplier, and it belongs next to any claim about running the
    protocol on hardware that exists.
    """
    from qiskit import transpile
    from qiskit_ibm_runtime import fake_provider
    from certify_sampler import build_zone_circuit
    backend = getattr(fake_provider, name)()
    qc, _ = build_zone_circuit(zone, lam)
    flat = transpile(qc, basis_gates=["u", "cx"], optimization_level=1)
    routed = transpile(qc, backend=backend, optimization_level=1, seed_transpiler=7)
    two = lambda c: sum(n for g, n in c.count_ops().items()
                        if g in ("cx", "ecr", "cz"))
    return dict(all_to_all=two(flat), routed=two(routed),
                factor=two(routed) / max(two(flat), 1))


def zones_that_fit(g=3, seed=5, beta=1.5, max_qubits=MAX_QUBITS):
    """The zones a density matrix can hold, the ones with something to lose first.

    Ordered by how many codewords are *infeasible*, because that is what the
    leak has to land on.  A zone whose every codeword satisfies its limits
    reports leak zero at any noise level and has proved nothing -- there was
    nowhere for the probability to leak to.
    """
    from model_instance import make_instance, utility_scale
    from model_partition import partition_aps
    from proto_zone import build_zone, enumerate_zone
    from certify_sampler import build_zone_circuit

    inst = make_instance(g=g, seed=seed)
    part = partition_aps(inst)
    lam = beta / utility_scale(inst)
    out = []
    for z in range(part.n_zones):
        zone = build_zone(inst, part, z)
        if zone.n_ue < 2 or any(w != 1 for w in zone.code_widths()):
            continue
        qc, _ = build_zone_circuit(zone, lam)
        if qc.num_qubits > max_qubits:
            continue
        rows, _, _ = enumerate_zone(zone, lam)
        out.append((zone, z, (1 << zone.n_ue) - len(rows), len(rows)))
    out.sort(key=lambda t: (-t[2], -t[0].n_ue))
    return out, lam


def leak_per_gate(zone, lam, rates=(3e-4, 1e-3, 3e-3)):
    """How fast the leak grows with error rate and circuit length.

    Over the range a device might plausibly reach, the leak is close to
    linear in `p2 * n_cx` -- one flipped flag is one leaked branch, and at
    these rates two flips are rare enough not to matter.  The slope is what
    lets a zone of ninety gates say something about a zone of two thousand,
    which is the only way to reach the sizes the manuscript actually uses: a
    density matrix stops at twelve qubits and those zones need thirty.

    It is an extrapolation and should be read as one.  It overstates the
    damage where the leak saturates, and saturation is well past the point
    where the protocol has already failed, so the error is on the safe side.
    """
    xs, ys = [], []
    for p2 in rates:
        r = noisy_zone_law(zone, lam, depolarizing(p2))
        xs.append(p2 * r["n_cx"])
        ys.append(r["leak"])
    return float(np.polyfit(xs, ys, 1)[0]), xs, ys


if __name__ == "__main__":
    zones, lam = zones_that_fit()
    if not zones:
        raise SystemExit("no zone small enough for a density matrix")
    zone, z, n_bad, n_ok = zones[0]

    print(f"zone Z{z}: {zone.n_ue} UEs, {len(zone.rbs)} RBs, {len(zone.aps)} APs, "
          f"|F_z|={n_ok} of {1 << zone.n_ue} ({n_bad} infeasible to leak onto)")
    r = routing_cost(zone, lam)
    print(f"two-qubit gates: {r['all_to_all']} all-to-all, {r['routed']} routed "
          f"onto a heavy-hex lattice ({r['factor']:.1f}x)\n")

    print(f"{'two-qubit error':>15} {'accept mu':>10} {'TVD':>9} "
          f"{'infeasible leak':>16}")
    ideal = noisy_zone_law(zone, lam, depolarizing(0.0))
    print(f"{'0 (ideal)':>15} {ideal['mu']:>10.4f} {ideal['tvd']:>9.1e} "
          f"{ideal['leak']:>16.1e}")
    for p2 in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
        r = noisy_zone_law(zone, lam, depolarizing(p2))
        print(f"{p2:>15.0e} {r['mu']:>10.4f} {r['tvd']:>9.3f} "
              f"{r['leak']:>16.3f}")

    d = device()
    r = noisy_zone_law(zone, lam, d)
    print(f"\n{d.backend_name} calibration, all-to-all: "
          f"mu={r['mu']:.4f}  TVD={r['tvd']:.3f}  leak={r['leak']:.3f}")

    from model_cost import BUDGET_DEFAULT
    slope, _, _ = leak_per_gate(zone, lam)
    print(f"\nleak grows at {slope:.3f} per unit of (two-qubit error x gates), "
          f"so in a {BUDGET_DEFAULT}-gate zone:")
    for target in (0.01, 0.05):
        print(f"  under {target:>4.0%} leak needs two-qubit error "
              f"{target / (slope * BUDGET_DEFAULT):.0e}")
