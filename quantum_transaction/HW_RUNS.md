# Hardware runs behind Sec. V-B

Every quantum-processor execution the paper rests on, with the job identifier
that can be looked up on the provider's dashboard, the settings it ran under,
and the numbers computed from its raw counts.

Each run is stored as one JSON file in `hw_runs/`, carrying the job ID, the
provider's own record of the job, the circuit as sent (gate counts, depth,
layout, basis), the measurement map, the raw counts, and the scored metrics.
Every run is scored by `noise_oracle.law_from_counts`, the same function that
scores the simulations, so a hardware row and a simulated row are comparable.

Reproduce:

```bash
python hw_oracle.py  --dry-run                       # transpile only, no QPU
python hw_oracle.py  --backend ibm_kingston --only inter --k 1
python hw_oracle.py  --job <job_id>                  # re-score a finished job
python ionq_oracle.py --noise-model forte-enterprise-1 --shots 2000 \
    --only inter --demands 1,2,1,2 --cap 4 --k 1     # free cloud simulator
python ionq_oracle.py --backend qpu.forte-enterprise-1 --confirm-qpu \
    --no-wait --shots 500 --only inter --demands 1,2,1,2 --cap 4 --k 1
```

---

## The circuit that was run

The inter-cell test circuit of Sec. V-B: four UEs, two APs, one amplification
round, that is $\mathcal{C}_z = \mathcal{G}_z\mathcal{P}_z$ of Eq. (37). It
exercises index encoding (Fig. 4), the utility rotations (Fig. 5), the QFT
adder, the comparator and the superflag (Fig. 6), one amplification round and
the decoded acceptance test (Fig. 7).

Six demand settings were used, from three to six UEs, chosen so that the
routed two-qubit count spans 364 to 1944 gates:

| demands | cap | UEs | feasible | logical 2q | routed on boston | depth |
|---|---:|---:|---:|---:|---:|---:|
| (1,1,1) | 2 | 3 | 6 | 196 | 364 | 916 |
| (1,2,1) | 3 | 3 | 6 | 306 | 580 | 1274 |
| (1,1,1,1) | 2 | 4 | 6 | 384 | 767 | 1470 |
| (1,2,1,2) | 4 | 4 | 10 | 384 | 787 | 1595 |
| (3,5,2,6,4) | 11 | 5 | 4 | 784 | 1380 | 2757 |
| (7,11,5,13,9,6) | 31 | 6 | 6 | 1052 | 1868 | 3517 |

Unit demands make the QFT adder count UEs; the others make it accumulate
demand, and the larger ones widen the counter, which is what carries the
circuit to zone scale. The largest is half the 2000-gate partition budget in
the counting of Eq. (39).

**What it is not.** This is not a zone of the campus instance of Sec. V-C, and
it does not use the general per-AP accumulation of Fig. 6: with two APs the
load of one AP determines the other, so a single counter and a two-sided
comparison decide both capacity constraints.

---

## Runs

`k=0` carries no oracle at all (6-8 two-qubit gates): the state preparation
$\mathcal{P}_z$ and the decoded acceptance test, with no amplification. It is
both the baseline the amplified circuit has to beat and, on its own, a test of
whether a device reproduces the accepted law of Eq. (36).

Ideal values are exact, not simulated. `P(opt)` is the probability that an
accepted sample is the optimal assignment, `J/J*` the mean utility of an
accepted sample. All IBM runs use 10 000 shots, optimization level 3, and
dynamical decoupling unless the row says otherwise.

### Without amplification: the accepted law is reproduced

| demands | cap | 2q | acceptance (exact) | TVD | P(opt) (exact) | J/J* |
|---|---:|---:|---|---:|---|---:|
| (1,1,1) | 2 | 6 | 0.298 (0.284) | 0.038 | 0.425 (0.436) | 0.865 |
| (1,2,1) | 3 | 6 | 0.290 (0.284) | 0.024 | 0.422 (0.436) | 0.868 |
| (1,1,1,1) | 2 | 8 | 0.100 (0.098) | 0.030 | 0.609 (0.632) | 0.870 |
| (1,2,1,2) | 4 | 8 | 0.138 (0.141) | 0.049 | 0.436 (0.444) | 0.833 |
| (3,5,2,6,4) | 11 | 8 | 0.087 (0.081) | 0.029 | 0.741 (0.770) | 0.904 |
| (7,11,5,13,9,6) | 31 | 8 | 0.139 (0.149) | 0.028 | 0.431 (0.423) | 0.903 |

Acceptance lands within 0.92-1.08 of its exact value in all six, and the
accepted law within 0.024-0.049 of Eq. (36). The utility rotations and the
decoded acceptance test therefore work on a device. This does **not** certify
the oracle: at k=0 feasibility is screened classically after measurement.

### With one amplification round: where the round stops paying

| device | demands | 2q routed | acceptance | k=0 baseline | signal | decay/gate | P(opt) |
|---|---|---:|---:|---:|---:|---:|---:|
| boston | (1,1,1) | 364 | **0.337** | 0.298 | 0.272 | 0.00357 | 0.373 |
| boston | (1,2,1) | 580 | 0.219 | 0.290 | 0.141 | 0.00338 | 0.322 |
| kingston | (1,1,1,1) | 767 | 0.050 | 0.100 | 0.041 | 0.00416 | 0.331 |
| boston | (1,1,1,1) | 767 | 0.055 | 0.100 | 0.048 | 0.00395 | 0.364 |
| boston | (1,2,1,2) | 787 | 0.125 | 0.138 | 0.107 | 0.00284 | 0.299 |
| boston | (1,2,1,2), **DD off** | 787 | 0.034 | 0.138 | ~0 | 0.02633 | 0.079 |
| boston | (3,5,2,6,4) | 1380 | 0.031 | 0.087 | 0.027 | 0.00261 | 0.418 |
| boston | (7,11,5,13,9,6) | 1868 | 0.036 | 0.139 | 0.015 | 0.00224 | 0.222 |

`signal` is the surviving fraction of the marked branch, measured between the
fully depolarized acceptance and the exact one; `decay/gate` is its logarithm
per routed two-qubit gate.

* **The round pays off below about 400 routed gates.** At 364 gates
  acceptance rises above the k=0 baseline (0.337 against 0.298); at 580 and
  above it does not. In the counting of Eq. (39), which is unrouted, the
  break-even is about 355 logical gates.
* **The decay is stable, and mildly favourable at size**: 0.0022-0.0042 per
  routed gate over a fivefold range of circuit size, decreasing as circuits
  grow, so a simple exponential model is conservative.
* **Idle decoherence dominates.** The same 787-gate circuit with dynamical
  decoupling disabled loses its signal entirely (0.034, at the depolarized
  floor). Gate errors are unchanged between those two rows, so what separates
  them is qubits waiting through a depth-1595 circuit. Circuit duration, not
  gate count, is the binding resource - which is what the parallel-counter
  schedule of Sec. V-D reduces.
* **The optimum survives well past the acceptance break-even.** At 1868
  routed gates (1052 logical, half the partition budget) the signal is 1.5%,
  yet the optimal assignment still holds 0.222 of the accepted mass against
  0.016 for a random output, a factor of 14.

### Layout variants, and why majority voting adds nothing

The 1868-gate circuit was run under five transpiler seeds, giving five
different physical layouts (1808-1944 routed gates), and aggregated three ways.

| aggregation | TVD | P(opt) | J/J* |
|---|---:|---:|---:|
| the five layouts, separately | 0.242-0.304 | 0.209-0.251 | 0.795-0.816 |
| averaged | 0.273 | 0.229 | 0.804 |
| majority vote, median over layouts | 0.269 | 0.234 | 0.807 |
| majority vote, min over layouts | 0.279 | 0.237 | 0.804 |
| exact | 0 | 0.423 | 0.908 |

Voting moves P(opt) from 0.229 to 0.234, inside the spread of the individual
layouts. The reason is structural: plurality voting earns its gains by
deleting outcomes that noise scattered, and the decoded acceptance test has
already deleted them. What remains is a distortion **among the six feasible
assignments**, which is common to all layouts and therefore survives the vote.
The construction's own acceptance test does most of what an error-mitigated
aggregation would do. The spread across layouts also shows the measurement is
reproducible rather than one lucky run.

## Averaged against sharpened aggregation (IonQ)

IonQ's debiasing runs the circuit as several variants and aggregates them.
Averaging keeps the shape of the distribution; sharpening takes a plurality
vote across variants, which suppresses outcomes that noise scattered and
concentrates the report on the branch the variants agree on.

Same job `01a0b994-…`, the two aggregations of its results:

| | distinct outcomes | acceptance | TVD | P(optimum \| accepted) | mean J/J* of an accepted sample |
|---|---:|---:|---:|---:|---:|
| noise-free reference | 81 | 0.836 | 0.011 | 0.444 | 0.835 |
| averaged (`sharpen=false`) | 180 | 0.099 | 0.257 | 0.240 | 0.746 |
| sharpened (`sharpen=true`) | 3 | 0.927 | 0.504 | 0.460 | 0.806 |

**How to read this.** The two aggregations answer different questions.

* *Which assignment should the zone report?* Sharpened output answers this at
  essentially noise-free quality: the optimum is returned with probability
  0.46 against 0.444 noise-free, and an accepted sample carries 0.806 J\*
  against 0.835 J\*.
* *What is the accepted law?* Only the averaged output can answer it. The
  sharpened distribution is not the circuit's law — it is three outcomes —
  and its TVD of 0.504 is larger than the averaged one.

**Two numbers that must not be quoted as performance.** Sharpened acceptance
of 0.927 exceeds the noise-free 0.836, which no amplification can do; it is an
artefact of discarding outcomes that fail the cost condition, so it cannot be
set against the k=0 baseline as evidence of an amplification gain. And the
sharpened TVD is not a measurement of the accepted law.

**What this costs the distributed stage.** The boundary coordination of
Sec. IV-F consumes joint samples, not one assignment: it filters the sample
set after each commitment and recomputes the remaining marginals. A sharpened
report holds three distinct joint outcomes, so filtering exhausts it almost
immediately. Sharpening is therefore usable for selecting a zone's assignment
and not for supplying the sample set that coordination re-conditions.

---

## Why the circuit is near the limit of present devices

IonQ Forte Enterprise publishes a median two-qubit fidelity of 0.9911 and a
SPAM fidelity of 0.9964 (characterization `b2ddac47-…`, 2026-09-19). For the
executed circuit,

    0.9911^196 * 0.9964^16 = 0.16,

against a measured normalized fidelity of 0.20 (averaged aggregation, Hellinger
against the exact 8-bit output distribution). The device performs as
specified; the circuit is simply longer than the error budget allows. A vendor
figure such as "a thousand two-qubit gates" refers to benchmark circuits with a
single correct answer evaluated under mitigation and plurality voting, where
the metric tolerates exactly the distortion measured here.

The same reading applies to IBM: `model_cost.py` takes a decay of 0.00118 per
two-qubit gate on ibm_kingston, which predicts 0.40 of the signal surviving
767 gates, whereas the run leaves about 0.04. The measured decay for this
circuit is roughly 0.0042 per gate, so the 2000-gate partition budget is
optimistic for a NISQ device and should be read as a fault-tolerant-era target.

## What the noise models predict, and why they cannot be used here

The same five amplified circuits were simulated with `NoiseModel.from_backend`
for ibm_boston, which carries the device's measured gate and readout errors.

| routed 2q | logical | simulated signal | measured signal | optimism |
|---:|---:|---:|---:|---:|
| 364 | 196 | 0.916 | 0.272 | 3.4x |
| 580 | 306 | 0.925 | 0.141 | 6.6x |
| 787 | 384 | 0.845 | 0.107 | 7.9x |
| 1380 | 784 | 0.792 | 0.027 | 29x |
| 1868 | 1052 | 0.818 | 0.015 | 55x |

The simulation barely responds to circuit size: a fivefold increase in gates
costs it 0.92 -> 0.82, while the device goes 0.27 -> 0.015. The optimism is
therefore not a constant that could be divided out; it grows with exactly the
size that matters for the partition budget. The cause is the DD row above:
the dominant loss is qubits idling through a deep circuit, and a gate-attached
noise model has no such term. IonQ's cloud simulator behaves the same way -
196 and 384 native gates gave the same acceptance to within sampling error.

Cost points the same way. The 1868-gate circuit took 51 minutes to simulate at
500 shots on 22 qubits, against 10 seconds on hardware for 10 000 shots.

The practical order is therefore: size circuits by transpiling (free and
exact), predict performance from the measured decay above, use noise
simulation only for what it does model - the effect of gate errors - and
reserve metered hardware, such as IonQ, for questions the other three cannot
answer.

**One caveat on IonQ results.** `qiskit-ionq` rebuilds simulator counts by
resampling the returned histogram on the client, so repeated fetches of one
job disagree. `ionq_oracle.raw_counts` bypasses it and marginalizes the
provider's own histogram onto the measured clbits; every IonQ number here is
computed that way and is reproducible from the stored record.
