# Scalability figures for `manuscript.tex`

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
python haiq_certify.py          # statevector check of the zone law   (~2 min)
python make_fig_scaling.py --recollect   # collect + plot              (~2 min)
python make_fig_snapshot.py              # plot                        (~1 min)
```

`make_fig_scaling.py` caches its sweep in `scaling_data.npz`; run without
`--recollect` to re-plot from the cache. Requires numpy, scipy, matplotlib,
qiskit (statevector only — no Aer, no networkx).

## Modules

| file | role |
|---|---|
| `haiq_instance.py` | coverage instance generator: APs on a jittered grid, RB pools, coverage-dependent candidate sets `V_i`, utilities. Growing `g` at fixed densities is the growth axis. |
| `haiq_cost.py` | two-qubit gate model of one compute-mark-uncompute pass, plus the measured Heron2/Heron3 anchors from `HARDWARE_LIMITS.md`. |
| `haiq_partition.py` | **the missing upstream of Sec. III-C**: builds the AP co-coverage graph and partitions it to minimize boundary UEs subject to a two-qubit budget per zone. |
| `haiq_zone.py` | the zone-local law and its sampler. |
| `haiq_protocol.py` | belief product over boundary marginals + confidence-ordered decimation. |
| `haiq_reference.py` | centralized strict optimum (MILP, HiGHS) — the denominator, not a competing protocol. |
| `haiq_certify.py` | Qiskit statevector check that the sampler is the circuit's accepted branch. |

## Three things worth knowing before writing the text

**1. The classical sampler is not an approximation.** From Eq. (theta)–(gibbs-local),
the accepted branch is `q_0 * prod_i g_{z,i}(r_i) * 1[feasible]`, which *is*
independent per-UE sampling with weight `g` followed by rejection on strict
feasibility. `haiq_certify.py` confirms this against a real circuit:
max TVD **2.7e-15** with zero spurious states over the zones checked. The
quantum and classical versions differ only in cost — `~K/mu_z` draws classically
against `O(mu_z^{-1/2})` amplification rounds — so scaling the protocol
experiment past what any current device can run is licensed, not assumed.

**2. The partition budget is a two-qubit count, not a qubit count.**
`HARDWARE_LIMITS.md` measured `lambda = 0.00118` on ibm_kingston, i.e. a
half-signal point near 590 two-qubit gates; qubit width is never the binding
constraint at these zone sizes. The partitioner therefore takes a 2q budget.
The default is 2000, which must exceed the most expensive *single* AP (1906
over the instances used) because an AP is the atom of the partition. Zones then
hold 4–9 UEs. The measured Heron2/Heron3 ceilings are drawn on the scaling
figure for reference, so the figure does not overstate what runs today.

**3. Decimation can dead-end, and the figure reports it.** Confidence-ordered
commitment is irreversible, so on rare instances it strands a zone. The guard
only commits to a value that leaves every zone completable, and looks one step
ahead at the neighbouring boundary UEs; when that still fails the zones are
re-sampled and coordination is retried (`attempts`), which is the manuscript's
own remedy. Over the 40 instances in the sweep, every run closed on a strictly
feasible global assignment; one instance needed 3 attempts. Two generated
instances were *globally* infeasible and are excluded — that is a property of
the random draw, not of the method, and the count is stated on the figure.

## Numbers currently on the figures

Scaling sweep (40 instances, 6 seeds per size, beta=1.5, K_z=400):

| zones | UEs | utility / centralized optimum | max zone 2q | centralized 2q | boundary density |
|---:|---:|---:|---:|---:|---:|
| 5.3 | 18 | 99.63 ± 0.54 % | 1693 | 10 801 | 61 % |
| 8.7 | 32 | 99.45 ± 0.46 % | 1786 | 24 227 | 54 % |
| 13.7 | 50 | 99.42 ± 0.28 % | 1919 | 37 755 | 55 % |
| 20.3 | 72 | 99.11 ± 0.46 % | 1871 | 64 120 | 66 % |
| 27.5 | 98 | 99.34 ± 0.25 % | 1894 | 88 647 | 63 % |
| 36.2 | 128 | 99.66 ± 0.20 % | 1926 | 138 690 | 65 % |
| 45.0 | 162 | 99.28 ± 0.29 % | 1950 | 173 983 | 64 % |

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
- The utility ratio sits at 99.1–99.7 % everywhere, i.e. the decomposition
  loss is small but non-zero. Claiming exactness would be wrong; "within 1 % of
  the centralized strict optimum" is what the data supports.
