# Evaluation artifacts for `manuscript_260901.tex` (Sec. V)

Two new figures, independent of Fig. 8 (which stays as the single-zone
complexity model). These cover the multi-zone case: a snapshot of one
partitioned region, and the average behaviour as the region grows.

## Figures

| file | what it shows |
|---|---|
| `fig/fig_snapshot.{pdf,png}` | one region: the partition the budget produced (zones as territories, each boundary UE drawn joined to the APs that offer it a candidate), the accepted law of four zones against the exact Gibbs reference, boundary coordination on the three most contested UEs, and the decimation order |
| `fig/fig_scaling.{pdf,png}` | growing the region: per-zone circuit cost stays bounded while the centralized circuit grows, at unchanged utility ratio |

## Reproducing

```bash
python certify_rotation.py               # V-B: does the released oracle realize Eq. (19)-(21)?
python certify_sampler.py                   # V-B: statevector check of the zone law
python noise_oracle.py                      # V-B: the same sampler on a calibrated device
python fig_scaling.py --recollect    # V-E: collect + plot
python fig_snapshot.py               # V-C: plot
python backoff_study.py                  # V-E: execution budget vs utility
python eval_numbers.py                   # every number Sec. V quotes, printed
```

`fig_scaling.py` caches its sweep in `scaling_data.npz`; run without
`--recollect` to re-plot from the cache. `fig_snapshot.py` caches its
ablation in `ablation_data.npz`; delete that file to recollect. Requires
numpy, scipy, matplotlib, qiskit, qiskit-aer, and — for `noise_oracle.py` —
qiskit-ibm-runtime for the calibrated fake backends.

**Amplification rounds.** `proto_zone.K_ROUNDS` is the `k` of Eq. (35) and is
the single knob that sets what a report costs: a zone pays `K_z / P_z(k)`
executions, so raising `k` relaxes the exponent backoff and lifts every
`beta_z` and effective sample size with it. It is `1`, matching Sec. V-A.
Changing it invalidates both caches.

**Sec. V numbers.** Several sentences in Sec. V quote quantities no panel is
annotated with. `eval_numbers.py` prints all of them from the same run that
produced the figures, so a sentence can be checked against a printed line
rather than read off an axis.

## Modules

Every `.py` in this directory, grouped by which question it answers. The
groups are a ladder: each one covers the range the one above it cannot reach.
A 25-AP zone circuit is not simulable and not runnable on today's hardware, so
the circuit is certified exactly where that *is* possible, the protocol is
then evaluated at scale on the law that certification fixes, and the cost of
sizes beyond either is read off the model.

    (A) circuit         does the released oracle produce exp[lambda J_z]?
                        exact, qiskit, small instances only
              |  the certified law is the sampler's law
              v
    (B) protocol        what does the whole multi-zone protocol do?
                        classical, at the sizes the figures plot
              |  beyond what sampling reaches
              v
    (C) model           where is the crossover?
                        closed-form cost, no run

The first word of a filename names the group, so `ls` sorts the directory into
it:

| prefix | group | what the file is about |
|---|---|---|
| `model_` | shared | the instance, its geography, its cost, and the optima it is scored against. No circuits. |
| `circuit_` | (A) | the two released oracles themselves. |
| `certify_` | (A) | an exactness check against one of them, suffix naming what is certified. |
| `noise_` | (A) | the same question under a device noise model, suffix naming the target circuit. |
| `proto_` | (B) | the multi-zone protocol: the zone-local law and the coordination over it. |
| `fig_` | (B) | a figure the paper prints, one script per figure. |

The remaining five carry their own names because they answer one-off
questions: `backoff_study`, `eval_numbers`, `preview_seeds`,
`hybrid_assignment_example`, and `fig/complexity/plot_qtg_assignment_runtime`.

`lib` below means the file is imported by other modules; `run` means it is a
script that produces something. Everything also has a `__main__`, so any file
can be run for a self-check.

### Shared model — no circuits, imported by both (A) and (B)

| file | kind | role |
|---|---|---|
| `model_instance.py` | lib | coverage instance generator: APs on a jittered grid, RB pools, coverage-dependent candidate sets `V_i`, utilities. Growing `g` at fixed densities is the growth axis. |
| `model_geo.py` | lib | campus geography for `geo=True`: registered base stations, open-ground masks. Lazily imported at `model_instance.py:99`, so it does not appear in a plain import graph. |
| `model_cost.py` | lib | two-qubit gate model of one compute-mark-uncompute pass, plus the measured Heron2/Heron3 anchors from `HARDWARE_LIMITS.md`. |
| `model_partition.py` | lib | **the missing upstream of Sec. III-C**: builds the AP co-coverage graph and partitions it to minimize boundary UEs subject to a two-qubit budget per zone. |
| `model_reference.py` | lib | centralized strict optimum (MILP, HiGHS) — the denominator, not a competing protocol. |

### (A) Circuit — exact, qiskit

Two circuits, and they stay two: **inter-cell picks the AP, intra-cell picks
the RB within it.** They are the paper's two released oracles, each meant to
be readable on its own, so the QFT adders they share are deliberately carried
in both files rather than factored into a common module.

| file | kind | role |
|---|---|---|
| `circuit_inter.py` | lib + run | **inter-cell circuit** (N UE x 2 AP): 1 bit/UE, QFT rate accumulator, `IntegerComparator`, one superflag. `build_sampler` assembles Eq. (26)-(36); `build_circuit` keeps the older Grover demo, called only by this file's own `__main__`. |
| `circuit_intra.py` | lib + run | **intra-cell circuit** (N node x M RB): 2 bits/node, QFT occupancy counter, validity flags, one superflag. Same two entry points, same deliberate duplication of `qft_on`/`controlled_add`/`controlled_sub`. |
| `certify_rotation.py` | run | does the accepted branch follow `exp[lambda J]`? Measures the corrected exponential rotation against the linear one the earlier oracle shipped. Drives `circuit_inter`/`circuit_intra`. |
| `certify_sampler.py` | run | statevector check that the **classical sampler** of `proto_zone.py` is the circuit's accepted branch. This is the join between (A) and (B): it is what licenses evaluating the protocol classically. |
| `noise_oracle.py` | run | the released sampler under a calibrated backend's noise model, transpiled to its basis **and coupling map**: acceptance, distance to the exact law, discard rate, measured leak. Built on the `certify_rotation` / `circuit_*` stack. |
| `noise_zone.py` | run | the same noise question asked one stack earlier, against the `certify_sampler` zone circuit instead of the `circuit_*` oracles. Superseded by `noise_oracle.py` for the Sec. V-B numbers and imported by nothing; kept because it is the only place the leak metric is posed against the zone circuit directly. |
| `hybrid_assignment_example.py` | run | self-contained five-step decompose / solve / exchange / reconcile / benchmark walkthrough. Shares no code with the modules above — it predates them and reads as a standalone illustration. |

### (B) Protocol — classical, at figure scale

Deliberately free of qiskit: `sample_zone` reproduces the rejection process
the accepted branch implements, on the law (A) certified.

| file | kind | role |
|---|---|---|
| `proto_zone.py` | lib | the zone-local law, its sampler, and the execution-exponent choice under a shot budget. `K_ROUNDS` lives here. |
| `proto_coordination.py` | lib | reconstruction to the common exponent, belief product over boundary marginals, confidence-ordered decimation. |
| `fig_snapshot.py` | run | one region opened up. Writes `fig/fig_snapshot.{pdf,png}` and every `fig/panels/snapshot_*`; caches its ablation in `ablation_data.npz`. |
| `fig_scaling.py` | run | growth at fixed density. Writes `fig/fig_scaling.{pdf,png}` and `fig/panels/scaling_*`; caches its sweep in `scaling_data.npz` (`--recollect` to redo it). |
| `backoff_study.py` | run | executions demanded at the target exponent, uncapped vs capped, and what the cap costs in utility. |
| `eval_numbers.py` | run | every quantity Sec. V quotes, recomputed from the figures' own run. |
| `preview_seeds.py` | run | contact sheet of candidate snapshot regions, `fig/seed_preview.png`. A selection aid for choosing `SEED`, not a paper figure. |

### (C) Model — closed form, nothing is run

| file | kind | role |
|---|---|---|
| `fig/complexity/plot_qtg_assignment_runtime.py` | run | the zone-local classical-quantum crossover. Quantum curves are modeled one-pass logical latencies; the classical curve is `0.400 N_z^3` ns, an interior coefficient of the 0.058-0.740 interval. Writes `fig/complexity/qtg_assignment_runtime_comparison.{pdf,png}`. Independent of everything above — it imports no project module. |

## Three things worth knowing before writing the text

**1. The classical sampler is not an approximation.** From Eq. (theta)–(gibbs-local),
the accepted branch is `q_0 * prod_i g_{z,i}(r_i) * 1[feasible]`, which *is*
independent per-UE sampling with weight `g` followed by rejection on strict
feasibility. `certify_sampler.py` confirms this against a real circuit:
max TVD **2.7e-15** with zero spurious states over the zones checked. The
quantum and classical versions differ only in cost, so scaling the protocol
experiment past what any current device can run is licensed, not assumed.

**2. The partition budget is a two-qubit count, not a qubit count.**
`HARDWARE_LIMITS.md` measured `lambda = 0.00118` on ibm_kingston, i.e. a
half-signal point near 590 two-qubit gates; qubit width is never the binding
constraint at these zone sizes. The partitioner therefore takes a 2q budget.
The default is 2000, which must exceed the most expensive *single* AP (1906
over the instances used) because an AP is the atom of the partition. Zones then
hold 4–9 UEs. The measured Heron2/Heron3 ceilings are drawn on the scaling
figure for reference, so the figure does not overstate what runs today.

**3. A shot budget is what makes the protocol executable, and it costs little.**
At the target exponent the tightest zones need 10^6-10^8 shots for their
reports, which no run can pay. Each zone therefore executes at the largest
`beta_z <= beta` its budget allows and is reweighted to the target afterwards
(Sec. IV-C), paying in effective sample size rather than bias. Capping every
zone at 10,000 shots costs 0.05 percentage points of utility on average and
under 1 point on the worst of 16 instances, against a 25,000-fold reduction in
the peak zone shot count. `fig/README.md` carries the comparison table.

**4. Decimation can dead-end, and the figure reports it.** Confidence-ordered
commitment is irreversible, so on rare instances it strands a zone. The guard
only commits to a value that leaves every zone completable, and looks one step
ahead at the neighbouring boundary UEs; when that still fails the zones are
re-sampled and coordination is retried (`attempts`), which is the manuscript's
own remedy. Over the 40 instances in the sweep, every run closed on a strictly
feasible global assignment. Two generated
instances were *globally* infeasible and are excluded — that is a property of
the random draw, not of the method, and the count is stated on the figure.

## Numbers currently on the figures

Scaling sweep (40 instances, 6 seeds per size, target beta=1.5, K_z=200,
shot budget 10 000 per zone):

| zones | UEs | utility / centralized optimum | max zone 2q | centralized 2q | boundary density |
|---:|---:|---:|---:|---:|---:|
| 5.3 | 18 | 99.87 ± 0.25 % | 1693 | 10 801 | 61 % |
| 8.7 | 32 | 98.81 ± 0.99 % | 1786 | 24 227 | 54 % |
| 13.7 | 50 | 99.20 ± 0.72 % | 1919 | 37 755 | 55 % |
| 20.3 | 72 | 99.09 ± 0.63 % | 1871 | 64 120 | 66 % |
| 27.5 | 98 | 99.15 ± 0.42 % | 1894 | 88 647 | 63 % |
| 36.2 | 128 | 99.01 ± 0.28 % | 1926 | 138 690 | 65 % |
| 45.0 | 162 | 99.26 ± 0.28 % | 1950 | 173 983 | 64 % |

Per-zone circuit cost is flat to within 15 % across a 9× growth in problem
size while the centralized circuit grows 16×, and the utility ratio does not
degrade.

## Open points for the manuscript text

- Sec. III-C currently *assumes* the partition ("Partition the APs into zones")
  while Sec. V cites the bounded-density regime as evidence. `model_partition.py`
  produces that regime and should be described, otherwise the regime has no
  source.
- Sec. V-D's protocol-resource discussion is still prose. The sweep records
  `comm_bits` and `exceptions` per run, so a coordination-cost panel can be
  added from the existing cache without new experiments.
- The utility ratio sits at 98.8–99.9 % everywhere, i.e. the decomposition and
  the shot budget together cost about a point. Claiming exactness would be
  wrong; "within about 1 % of the centralized strict optimum" is what the data
  supports.
