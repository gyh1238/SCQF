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
python haiq_certify.py                   # V-B: statevector check of the zone law
python noise_run.py                      # V-B: the same sampler on a calibrated device
python make_fig_scaling.py --recollect    # V-E: collect + plot
python make_fig_snapshot.py               # V-C: plot
python backoff_study.py                  # V-E: execution budget vs utility
python eval_numbers.py                   # every number Sec. V quotes, printed
```

`make_fig_scaling.py` caches its sweep in `scaling_data.npz`; run without
`--recollect` to re-plot from the cache. `make_fig_snapshot.py` caches its
ablation in `ablation_data.npz`; delete that file to recollect. Requires
numpy, scipy, matplotlib, qiskit, qiskit-aer, and — for `noise_run.py` —
qiskit-ibm-runtime for the calibrated fake backends.

**Amplification rounds.** `haiq_zone.K_ROUNDS` is the `k` of Eq. (35) and is
the single knob that sets what a report costs: a zone pays `K_z / P_z(k)`
executions, so raising `k` relaxes the exponent backoff and lifts every
`beta_z` and effective sample size with it. It is `1`, matching Sec. V-A.
Changing it invalidates both caches.

**Sec. V numbers.** Several sentences in Sec. V quote quantities no panel is
annotated with. `eval_numbers.py` prints all of them from the same run that
produced the figures, so a sentence can be checked against a printed line
rather than read off an axis.

## Modules

| file | role |
|---|---|
| `haiq_instance.py` | coverage instance generator: APs on a jittered grid, RB pools, coverage-dependent candidate sets `V_i`, utilities. Growing `g` at fixed densities is the growth axis. |
| `haiq_cost.py` | two-qubit gate model of one compute-mark-uncompute pass, plus the measured Heron2/Heron3 anchors from `HARDWARE_LIMITS.md`. |
| `haiq_partition.py` | **the missing upstream of Sec. III-C**: builds the AP co-coverage graph and partitions it to minimize boundary UEs subject to a two-qubit budget per zone. |
| `haiq_zone.py` | the zone-local law, its sampler, and the execution-exponent choice under a shot budget. |
| `haiq_protocol.py` | reconstruction to the common exponent, belief product over boundary marginals, confidence-ordered decimation. |
| `haiq_reference.py` | centralized strict optimum (MILP, HiGHS) — the denominator, not a competing protocol. |
| `haiq_certify.py` | Qiskit statevector check that the sampler is the circuit's accepted branch. |
| `qtg_inter_assignment.py` | **the paper's inter-cell circuit**: 1 bit/UE, QFT rate accumulator, `IntegerComparator`, one superflag. `build_sampler` assembles Eq. (26)-(36); `build_circuit` keeps the older Grover demo. |
| `qtg_intra_assignment.py` | **the paper's intra-cell circuit**: 2 bits/node, QFT occupancy counter, validity flags, one superflag. Same two entry points. |
| `hybrid_assignment_example.py` | the five-step decompose / solve / exchange / reconcile / benchmark walkthrough. |
| `certify_rotation.py` | does the accepted branch follow `exp[lambda J]`? Measures the corrected exponential rotation against the linear one the earlier oracle shipped. |
| `noise_run.py` | the same sampler under a calibrated backend's noise model, transpiled to its basis **and coupling map**: acceptance, distance to the exact law, discard rate, and measured leak. |
| `backoff_study.py` | executions demanded at the target exponent, uncapped vs capped, and what the cap costs in utility. |
| `eval_numbers.py` | every quantity Sec. V quotes, recomputed from the figures' own run. |

## Three things worth knowing before writing the text

**1. The classical sampler is not an approximation.** From Eq. (theta)–(gibbs-local),
the accepted branch is `q_0 * prod_i g_{z,i}(r_i) * 1[feasible]`, which *is*
independent per-UE sampling with weight `g` followed by rejection on strict
feasibility. `haiq_certify.py` confirms this against a real circuit:
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
  while Sec. V cites the bounded-density regime as evidence. `haiq_partition.py`
  produces that regime and should be described, otherwise the regime has no
  source.
- Sec. V-D's protocol-resource discussion is still prose. The sweep records
  `comm_bits` and `exceptions` per run, so a coordination-cost panel can be
  added from the existing cache without new experiments.
- The utility ratio sits at 98.8–99.9 % everywhere, i.e. the decomposition and
  the shot budget together cost about a point. Claiming exactness would be
  wrong; "within about 1 % of the centralized strict optimum" is what the data
  supports.
