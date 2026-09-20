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

Two demand settings were used:

| setting | demands | AP capacity | feasible | ideal acceptance at k=1 | ideal at k=0 |
|---|---|---:|---:|---:|---:|
| unit | (1,1,1,1) | 2 | 6 | 0.666 | 0.098 |
| unequal | (1,2,1,2) | 4 | 10 | 0.836 | 0.140 |

The unequal setting is the one that makes the QFT adder accumulate demand
rather than count UEs.

**What it is not.** This is not a zone of the campus instance of Sec. V-C, and
it does not use the general per-AP accumulation of Fig. 6: with two APs the
load of one AP determines the other, so a single counter and a two-sided
comparison decide both capacity constraints.

---

## Runs

`k=0` carries no oracle at all (8 two-qubit gates): the state preparation
$\mathcal{P}_z$ and the decoded acceptance test, with no amplification. It is
both the baseline the amplified circuit has to beat and, on its own, a test of
whether a device reproduces the accepted law of Eq. (36).

Ideal values are exact, not simulated. Exactness of the accepted law is the
statevector check of Sec. V-B; these are the same quantities measured on a
device. `P(opt)` is the probability that an accepted sample is the optimal
assignment, and `J/J*` the mean utility of an accepted sample.

| device | setting | k | 2q gates | shots | acceptance | TVD | P(opt) | J/J* | job |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| ibm_boston (Heron r3) | unit | 0 | 8 | 10 000 | 0.100 | 0.030 | 0.609 | 0.870 | `dan7dbg2fm4c73f3sd00` |
| *exact* | unit | 0 | | | *0.098* | *0* | *0.632* | *0.878* | |
| ibm_boston | unequal | 0 | 8 | 10 000 | 0.138 | 0.049 | 0.436 | 0.833 | `dan7didr85ps73fdn6j0` |
| *exact* | unequal | 0 | | | *0.140* | *0* | *0.444* | *0.835* | |
| ibm_boston | unequal | 1 | 787 routed | 10 000 | 0.125 | 0.159 | 0.299 | 0.773 | `dan7df5r85ps73fdn6g0` |
| ibm_boston | unit | 1 | 767 routed | 10 000 | 0.055 | 0.269 | 0.364 | 0.768 | `dan7adgpqrnc7398dp0g` |
| ibm_kingston (Heron r2) | unit | 1 | 767 routed | 10 000 | 0.050 | 0.302 | 0.331 | 0.753 | `dan775n8gn2s739ma17g` |
| IonQ Forte Enterprise | unequal | 1 | 196 native | 500 | 0.099 | 0.257 | 0.240 | 0.746 | `01a0b994-aabd-715e-b452-a7655253d170` |
| *exact* | unequal | 1 | | | *0.836* | *0* | *0.444* | *0.835* | |
| *exact* | unit | 1 | | | *0.666* | *0* | *0.632* | *0.878* | |

IBM runs carry no error mitigation. The IonQ run was submitted with debiasing
enabled; the row above is its **averaged** aggregation (`sharpen=false`), and
the sharpened aggregation of the same job is discussed below.

References for reading the rows: a fully depolarized output would accept with
probability 0.023 (unit) or 0.039 (unequal), and a uniform distribution over
the feasible set sits at TVD 0.466 (unit) or 0.379 (unequal).

### What the runs show

* **The accepted law of Eq. (36) is reproduced on hardware at k=0.** On
  ibm_boston the two settings land at TVD 0.030 and 0.049 over 10 000 shots,
  with acceptance, P(opt) and mean utility matching the exact values to the
  second decimal. The utility rotations of Fig. 5 and the decoded acceptance
  test therefore work on a device, not only in simulation. What this does
  *not* certify is the oracle: at k=0 feasibility is screened classically
  after measurement.
* **Every reported assignment satisfies the constraints, on every device.**
  This holds by construction: the acceptance test of Sec. IV-F is applied to
  the decoded assignment, so a corrupted flag costs a shot but cannot put an
  infeasible assignment into a report.
* **The optimum remains the most frequent report everywhere**, including the
  amplified runs.
* **Amplification does not pay off at this circuit size, on any device
  tested.** Acceptance after one round is below what the same sampler reaches
  at k=0: 0.125 against 0.138 on ibm_boston, 0.055 against 0.100 on the unit
  setting, 0.099 against an ideal 0.140 on IonQ. One amplification round is a
  loss on present hardware; the crossover is a matter of gate fidelity, not of
  the construction. The k=0 rows are what isolates this, since they show the
  rest of the pipeline reaching its exact values on the same device.
* **The accepted law is flattened once the oracle is in the circuit.** Noise
  mostly removes weight from the optimum: at IonQ the optimum holds 0.240 of
  the accepted mass against 0.444 exactly, while the next three assignments
  stay near their exact values.

---

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

## Cross-checks on the simulators

| | acceptance | TVD |
|---|---:|---:|
| IonQ cloud simulator, `forte-enterprise-1` noise, 2000 shots | 0.110 | 0.215 |
| IonQ Forte Enterprise hardware, 500 shots, averaged | 0.099 | 0.257 |

The vendor noise model predicts the device within the sampling error of these
runs, so the hardware result is not an unlucky draw.

**One caveat on IonQ results.** `qiskit-ionq` rebuilds simulator counts by
resampling the returned histogram on the client, so repeated fetches of one
job disagree. `ionq_oracle.raw_counts` bypasses it and marginalizes the
provider's own histogram onto the measured clbits; every IonQ number here is
computed that way and is reproducible from the stored record.
