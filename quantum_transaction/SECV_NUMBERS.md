# Sec. V: every parameter and measured number

The manuscript's Section V was trimmed to read as prose, so the parameter
settings and the secondary measurements no longer appear in the text. They are
here instead, with the command that produces each one, so a sentence in Sec. V
can be checked against a printed line rather than reconstructed from the code.

Nothing here is a new experiment. Every figure and every number in Sec. V comes
from these settings.

Regenerate everything:

```bash
python certify_sampler.py               # V-B  exact statevector check
python certify_rotation.py              # V-B  the released oracle's rotation
python noise_oracle.py --sweep          # V-B  calibrated backend + error sweep
python fig_snapshot.py                  # V-C  the campus instance and its panels
python fig_scaling.py --recollect       # V-E  the 7-size sweep
python backoff_study.py                 # V-A  execution budget vs utility
python eval_numbers.py                  # every number Sec. V quotes, printed
```

`fig_scaling.py` caches in `scaling_data.npz`, `fig_snapshot.py` caches its
ablation in `ablation_data.npz`. Delete a cache to recollect. Changing
`proto_zone.K_ROUNDS` invalidates both.

---

## 1. Instance generator

`model_instance.make_instance`, defaults at `model_instance.py:79-80`:

| parameter | value | meaning |
|---|---|---|
| `g` | varies | AP grid side, in AP spacings. The growth axis. |
| `ue_per_ap` | 2.0 | UEs per AP, held fixed as `g` grows |
| `n_rb_per_ap` | 4 | RBs owned by each AP |
| `radius` | 1.2 | coverage radius, in AP spacings |
| `max_deg` | 2 | candidates kept per UE, so `|C_i| <= 2` and `l_i = 1` |
| `w_rb` | 1 | RB occupancy limit `W_r` |
| `w_ap` | 4 | AP admission budget `W_a` |
| `jitter` | 0.18 | grid perturbation |
| `geo` | `False` / `True` | bare square, or campus rasters via `model_geo.py` |

An AP covering a UE offers it exactly one RB, allocated round robin by
proximity (`model_instance.py:131-136`), so a UE's candidate count equals the
number of APs covering it, before `max_deg` truncation. The truncation binds:
at `g=5, seed=1` the covering-AP count runs 1 to 6 with mean 3.66.

`d_i = 1` for every UE. Utilities decay with UE-AP distance with a per-RB
perturbation; `model_instance.py:157-164` defines the scale

```
ubar = mean over UEs of ( max utility among that UE's candidates )
```

which is the `u_bar` of `lambda = beta / u_bar`. It is a property of the
utility model, not of an instance, so every zone shares it.

**Which instance is which.** The three experiment families use different
instances, and this is easy to lose:

| experiment | instance family | source |
|---|---|---|
| V-B exact check | `g=4`, seed 2 | `certify_sampler.py:177` |
| V-B noise | two hand-sized circuits, below | `noise_oracle.py:172` |
| V-C snapshot | `g=5`, seed 0, `geo=True`, peers = seeds 0-15 | `fig_snapshot.py:63-64,169` |
| V-C ablation | `g=4..7` x seeds 0-3 = 16, 14 feasible | `fig_snapshot.py:80-81` |
| V-D | no instance; closed-form model | `fig/complexity/` |
| V-E sweep | `g=3..9` x seeds 0-5 = 42, 40 feasible | `fig_scaling.py:43-44` |
| V-A backoff | `g=4..7` x seeds 0-3 = 16, 281 zones | `backoff_study.py:35-36` |

---

## 2. Protocol

| parameter | value | where |
|---|---|---|
| `beta` (target common exponent) | 1.5 | `fig_scaling.py:45`, `fig_snapshot.py:65` |
| `lambda` | `beta / ubar` | `model_instance.py:157` |
| `K_ROUNDS` (`k` of Eq. 35) | 1 | `proto_zone.py:222` |
| `K_z` accepted joint samples | 200 | `K_ACCEPT` |
| execution cap per zone | 10 000 | `SHOT_BUDGET` |
| `K_min`, re-sample trigger | 25 | `proto_coordination.py:64` |
| partition budget | 2000 two-qubit gates | `model_cost.py:42` |

Partition: greedy co-coverage region growing with single-AP boundary
refinement, minimizing `|B|` subject to `n2q(z) <= budget` for every zone
(`model_partition.py:89-142`). The constraint is `<=`, so a zone may sit
exactly on 2000; one instance of the 40 does.

---

## 3. Hardware anchors and the cost model

`model_cost.py:33-42`. Two-qubit gates set the executable size, not qubit
width:

| constant | value | meaning |
|---|---|---|
| `LAMBDA_KINGSTON` | 0.00118 | measured, `p ~ 0.53` at 640 heavy-hex 2q |
| `CEIL_KINGSTON` | 587 | 2q at the 50 % signal point, `ln2 / lambda` |
| `LAMBDA_BOSTON` | 0.00068 | scaled from the 2q error ratio (estimate) |
| `CEIL_BOSTON` | 1019 | same, for Boston |
| `BUDGET_DEFAULT` | 2000 | partition budget |

Signal model `p_signal = exp(-lambda * n_2q)`. The 2000-gate budget therefore
sits above the measured 50 % ceiling; it marks roughly the largest oracle a
present device runs with a useful fraction of shots error-free, not a
noise-free circuit.

---

## 4. V-B, exact validation

`python certify_sampler.py`, on `g=4` seed 2, under a cap of seven UEs and 24
qubits. Five zones of that partition qualify:

| zone | `N_z` | qubits | `|F_z|` | `mu_z` | TVD to `exp[lambda J_z]` | accepted outside `F_z` |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 6 | 1 | 0.5000 | 0.00e+00 | 0 |
| 3 | 2 | 10 | 4 | 0.3934 | 1.14e-16 | 0 |
| 7 | 6 | 21 | 36 | 0.0316 | 1.92e-15 | 0 |
| 8 | 6 | 21 | 36 | 0.0214 | 1.73e-15 | 0 |
| 9 | 6 | 21 | 36 | 0.0168 | 2.75e-15 | 0 |

Largest TVD **2.75e-15**, largest circuit **six UEs on 21 qubits**, no accepted
probability outside `F_z` anywhere. This is the exactness result the classical
sampler of Sec. V-C and V-E rests on.

---

## 5. V-B, noise

`python noise_oracle.py --sweep`. Backend `FakeSherbrooke`, 127 qubits, median
rates **2q 7.79e-03, 1q 2.36e-04, readout 1.98e-02**, 10 000 shots per point.

Two questions are kept apart. **Routing is a transpile, not a simulation** —
there is nothing to simulate and a 127-qubit device could not be simulated
anyway. **Noise is simulated at all-to-all connectivity**, on the circuit's own
qubit count, not the device's.

### Lattice cost (transpile only)

| circuit | qubits | 2q before | 2q after | factor | depth before | depth after |
|---|---:|---:|---:|---:|---:|---:|
| inter-cell, 4 UE x 2 AP | 16 | 384 | 969 | 2.5x | 773 | 2717 |
| intra-cell, 3 node x 3 RB | 21 | 1354 | 3373 | 2.5x | 2036 | 8697 |

### Noise, inter-cell circuit (`|F| = 6`, `lambda = 0.1976`)

| condition | accept | discard | TVD to `exp[lambda J]` | decoded infeasible | accepted mass vs noise-free |
|---|---:|---:|---:|---:|---:|
| noise-free | 0.666 | 0.334 | 0.005 | 0.0496 | — |
| FakeSherbrooke rates | 0.123 | 0.877 | 0.134 | 0.1228 | 18 % |
| depolarizing `p2=1e-04` | 0.650 | 0.350 | 0.011 | 0.0487 | 98 % |
| depolarizing `p2=1e-03` | 0.513 | 0.487 | 0.016 | 0.0733 | 77 % |
| depolarizing `p2=1e-02` | 0.104 | 0.896 | 0.137 | 0.1115 | 16 % |

Intra-cell, noise-free reference: accept 0.209, discard 0.791, TVD 0.011,
decoded infeasible 0.2972. Its noisy rows take about 50 s each; run
`noise_oracle.py --sweep` without a timeout to complete them.

**Reading this.** Noise thins the accepted branch and tilts the accepted
distribution; it never admits an infeasible assignment. The acceptance test
applies to the *decoded* assignment, not to the internal flags, so a corrupted
outcome is discarded. That is why the "decoded infeasible" column is nonzero
even noise-free (0.0496): an infeasible assignment can leave every utility
qubit in `|0>` and still fail the decode check. `p2 = 1e-03`, roughly current
hardware, holds the distance at 0.016 with 77 % of the accepted mass — which is
the regime the 2000-gate budget targets.

---

## 6. V-C, the campus snapshot

`python eval_numbers.py`. Instance `g=5`, seed 0, `geo=True`.

| quantity | value |
|---|---|
| APs / UEs | 25 / 50 |
| zones `|Z|` | 14 |
| boundary UEs `|B|` | 23 (46 % of UEs) |
| utility | 98.66 % of `J*` |
| peers | 16 seeds, 13 globally feasible, 97.34-99.97 %, mean 99.19 % |
| largest zone circuit | 1910 two-qubit gates |
| sampling distance, 4 displayed zones | max 0.184, mean 0.157 |
| sampling distance, all 14 zones | max 0.184, mean 0.115 |
| `beta_z` range | 1.36 to 1.50 |
| effective sample size | min 198 of 200 |
| exception re-samples | 5 |
| executions | 19 224 initial, +11 973 from re-samples (+62 %) |
| executions per zone | 201 to 9 997 (cap 10 000) |
| `P_z(k) / mu_z` | 5.9x at `k=1`; the `mu->0` bound `(2k+1)^2` is 9 |

**`beta_z` and ESS are per-instance.** For this instance they are 1.36-1.50 and
198/200. The values 0.38 and 88/200 that appear in some drafts are the minima
over all zones of the 40-instance V-E sweep — correct there, not here.

### The two most dependent boundary pairs

`strongest_pair` in `fig_snapshot.py:407`, ranked by TV between the empirical
joint and the outer product of its own marginals.

**Zone Z13, UEs 12 & 14** — TV 0.385. Candidates `{RB52, RB38}` and
`{RB52, RB35}`, sharing RB52 under `W_r = 1`.

|  | UE14 = RB52 | UE14 = RB35 |
|---|---:|---:|
| **UE12 = RB52** | 0.000 (infeasible) | 0.307 |
| **UE12 = RB38** | 0.629 | 0.065 |

Marginals: UE12 `[0.307, 0.693]`, UE14 `[0.629, 0.371]`. Their product puts
**19.3 %** on the cell the joint report excludes.

**Zone Z0, UEs 30 & 38** — TV 0.377. Candidates `{RB85, RB81}` and
`{RB81, RB28}`, sharing RB81. Joint `[[0.555, 0.105], [0.000, 0.340]]`;
product puts **19 %** on the excluded cell.

### Commitment order for that pair

From the decimation trace of the snapshot run:

```
step  7:  UE12 -> index 1 = RB38   confidence 0.885   belief [0.115, 0.885]
step 19:  UE14 -> index 1 = RB35   confidence 0.533   belief [0.467, 0.533]
final:    UE12 -> RB38,  UE14 -> RB35
```

UE12 commits **first**. The belief `b_i` is a product over every zone holding
that UE, so it is not the zone-local marginal above, and twelve other
commitments fall between the two steps. The point survives either way: once
UE12 holds RB38, filtering removes every draw in which UE12 holds RB52, so the
(RB52, RB52) cell the one-pass product favours leaves the report before UE14 is
considered.

### Joint report vs one-pass marginals

`ablation_data.npz`, 16 instances at `g=4..7` seeds 0-3, 14 globally feasible.
Same sampling, same commitment rule; only the report format differs.

| | value |
|---|---|
| re-conditioning wins | 11 of 14 |
| ties / marginal wins | 0 / 3 |
| mean improvement | +0.87 percentage points of `J*` |
| median | +0.53 pp |
| range | -2.36 to +3.42 pp |

---

## 7. V-D, the complexity model

`fig/complexity/plot_qtg_assignment_runtime.py`. **Closed form, nothing is
run.** No project module is imported.

| constant | value |
|---|---|
| `TAU_NS` (`t_layer`) | 12.5 ns per logical layer |
| `CLASSICAL_COEFFICIENT_NS` (`rho`) | **0.400** |
| `CANDIDATES_PER_UE` (`D_z = 2 N_z`) | 2.0 |
| `CAPACITY_CONDITIONS_PER_UE` | 1.0 |
| `RESOURCE_COUNTER_BITS` (`Q_z^cnt`) | 3 |
| `PARALLEL_COUNTERS` (`R_cnt`, enhanced) | 4 |
| `ARITHMETIC_DEPTH_FACTOR` (enhanced) | 0.75 |

`rho = 0.400` is an interior point of the 0.058-0.740 interval measured for the
Hungarian method in the prior work — a 13x-wide interval, so the crossover
location inherits that width. The figure legend prints the literal string
`$0.400N_z^3$ ns`; the text must not say 0.399.

### Modeled latency per zone UE

| method | ns per `N_z`, at `N_z = 10 / 40 / 160` | asymptotic form |
|---|---|---|
| proposed, enhanced | 133.5 / 119.0 / 114.4 | ~113 `N_z` |
| proposed, serial | 392.2 / 385.2 / 382.6 | ~381 `N_z` |
| GAS | 477.6 / 537.5 / 607.8 | ~331 `N_z` + 37.5 `N_z log2 N_z` |
| QTG | 640.7 / 704.9 / 776.2 | ~500 `N_z` + 37.5 `N_z log2 N_z` |

The `37.5 N_z log2 N_z` term is the register accumulating a zone-wide numeric
objective. The proposed oracle carries no such register: the real-valued
utilities live in utility-qubit rotations.

### Crossover against `0.400 N_z^3`

| method | `N_z` |
|---|---:|
| proposed, enhanced | 17.71 |
| proposed, serial | 31.07 |
| GAS | 36.51 |
| QTG | 42.05 |

One oracle pass only. Excluded: repeated search rounds, accepted-draw
repetitions, readout, queueing, boundary coordination. This is oracle-level
scaling, not an end-to-end speedup.

---

## 8. V-E, the region sweep

`scaling_data.npz`, `g = 3..9` x seeds 0-5. **42 instances, 2 globally
infeasible and dropped, 40 remain. All 40 finish globally feasible.**

| `g` | UEs | zones | utility mean | min | max | largest zone 2q (mean) | centralized 2q | traffic KB | boundary |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 18 | 5.3 | 99.32 % | 96.79 | 100.00 | 1693 | 10 801 | 7.3 | 61 % |
| 4 | 32 | 8.7 | 99.13 % | 98.17 | 100.00 | 1786 | 24 227 | 11.2 | 54 % |
| 5 | 50 | 13.7 | 99.09 % | 97.74 | 99.69 | 1919 | 37 755 | 18.2 | 55 % |
| 6 | 72 | 20.3 | 98.77 % | 98.05 | 99.50 | 1871 | 64 120 | 28.8 | 66 % |
| 7 | 98 | 27.5 | 99.17 % | 98.75 | 99.51 | 1894 | 88 647 | 37.5 | 63 % |
| 8 | 128 | 36.2 | 98.76 % | 97.83 | 99.46 | 1926 | 138 690 | 52.6 | 65 % |
| 9 | 162 | 45.0 | 99.15 % | 98.80 | 99.59 | 1950 | 173 983 | 65.5 | 64 % |

The plotted "largest zone circuit" is the **mean over seeds of each instance's
largest zone**, which is what makes it 1693 -> 1950. The worst single zone
anywhere in the sweep is **2000**, in one instance — exactly on the budget, not
over it.

| claim | measured |
|---|---|
| centralized growth | 10 801 -> 173 983, **x16.11** |
| largest zone growth | 1693 -> 1950, **x1.152** |
| zone count growth | 5.3 -> 45.0, **x8.44** |
| regional traffic growth | **x8.97** |
| traffic per zone | 11 206 -> 11 916 bits, **+6.3 %** |
| utility, mean per size | **98.76 - 99.32 %** of `J*` |
| utility, individual runs | **96.79 - 100.00 %** |
| globally feasible | 40 of 40 |

Traffic model: `(1 + n_re) * K_z * ( sum_{i in B_z} ceil(log2 |C_i|) + 32 )`
bits per zone. The 32-bit field carries `J_z` for the reweighting.

---

## 9. V-A, exponent backoff

`python backoff_study.py`. 16 instances at `g=4..7` seeds 0-3, **281 zones**.

| quantity | value |
|---|---|
| executions at the target exponent, per zone | 1 386 to 28 358 487 |
| peak single zone | 2.84e+07 |
| capped at 10 000 | reduction **x2 836** |
| rounds maximizing `P_z(k)` | 0 to **887** |
| mean utility, uncapped | 99.19 % of `J*` |
| mean utility, capped | 99.19 % of `J*` |
| loss from the cap | 0.00 pp mean, **0.74 pp worst** |

The spread over zones is four orders of magnitude, which is why the cap exists.
Reconstruction by `w_z^(n) = exp[(lambda - lambda_z) J_z]` is unbiased, so the
cap trades effective sample size for execution cost, not accuracy — the mean
utility is unchanged to two decimals. The round count maximizing `P_z(k)`
reaching 887 for the most tightly constrained zones, against 0 for the
sparsest, is what motivates `k = 1` as the practical setting rather than a
per-zone optimum.

---

## 10. Still open

- **`1,906` gates.** `manuscript_260901.tex` Sec. V-D says "the largest
  measured circuit costing 1,906 gates". Not reproduced. The snapshot
  partition's largest zone is 1910; the largest anywhere in the V-E sweep is
  2000. Re-derive before submission.
- **`rho` in the text.** Some drafts say 0.399; the code and the figure legend
  say 0.400.
- **Peer upper bound.** `97.3-100.0 %` in V-C rounds 99.97 % into reaching the
  optimum. No instance of that family did. V-E's `96.79-100 %` is exact — that
  sweep does contain a 100.00 %.
- **RB selection is not searched.** Sec. V evaluates joint AP-RB candidates,
  but which RB an AP offers is fixed at instance-build time by the round-robin
  above, not chosen by the circuit. `W_r = 1` still binds because the offer
  wraps: one AP covers up to 14 UEs while owning 4 RBs. See
  `ENCODING_NOTES.md`.
