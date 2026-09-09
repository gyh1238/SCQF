# Section V supplement

Everything Section V rests on: what it states, what it does not state but the
code measures, what it states that the code contradicts, and what nothing in
this repository measures at all.

Section V was trimmed to read as prose, so most of the settings and the
secondary measurements left the text. They are here, each with the command
that produces it, so a sentence can be checked against a printed line rather
than reconstructed from the code.

Every item carries one mark:

| mark | meaning |
|---|---|
| **P** | stated in Sec. V, and the code agrees |
| **P!** | stated in Sec. V, and the code disagrees — see [Corrections](#corrections-to-apply) |
| **R** | measured here, not stated in Sec. V |
| **—** | not measured anywhere in this repository |

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

## 0. At a glance

| | **P** agrees | **P!** disagrees | **R** unstated | **—** unmeasured |
|---|---:|---:|---|---:|
| V-A Setup | 8 | 1 | 6 | 0 |
| V-B Validation | 5 | 0 | the whole noise study | 2 |
| V-C Coordination | 10 | 4 | 7 | 1 |
| V-D Cost model | 6 | 1 | 5 | 2 |
| V-E Region scale | 10 | 2 | 4 | 1 |
| Exponent backoff | 0 | 0 | all of it | 0 |

The two gaps that matter most are in [Gaps](#gaps): **no real-QPU execution
exists in this repository**, and **`HARDWARE_LIMITS.md`, the source of the only
claimed hardware measurement, is not in it either.**

---

## V-A. Setup

### Stated, and correct

| claim | measured | |
|---|---|---|
| `\|C_i\| <= 2` | `max_deg=2`, so `l_i = 1` | **P** |
| 2000 two-qubit gate budget per zone oracle | `BUDGET_DEFAULT = 2000` | **P** |
| `beta = 1.5` | `BETA = 1.5` | **P** |
| `lambda = beta / ubar`, `ubar` = mean over UEs of the largest utility among their candidates | `model_instance.py:157-164` | **P** |
| `k = 1` amplification round | `K_ROUNDS = 1` | **P** |
| `K_z = 200` accepted joint samples | `K_ACCEPT = 200` | **P** |
| execution cap `10^4` per zone | `SHOT_BUDGET = 10_000`; measured 201 to 9 997 per zone | **P** |
| `K_min = 25` re-sample trigger | `proto_coordination.py:64` | **P** |

### Stated, and wrong

| claim | measured | |
|---|---|---|
| one round gives "approximately a ninefold reduction" | `P_z(k)/mu_z = 5.9`. `(2k+1)^2 = 9` is the `mu->0` bound, and these zones are not in that limit | **P!** |

### Measured, not stated

| | value | |
|---|---|---|
| UEs per AP | 2.0, held fixed as `g` grows | **R** |
| RBs per AP | 4 | **R** |
| coverage radius | 1.2 AP spacings | **R** |
| `W_r`, `W_a`, `d_i` | 1, 4, 1 | **R** |
| grid jitter | 0.18 | **R** |
| partition method | greedy co-coverage region growing + single-AP boundary refinement, minimizing `\|B\|` subject to `n2q(z) <= budget` (`model_partition.py:89-142`) | **R** |

The budget constraint is `<=`, not `<`, so a zone may sit exactly on 2000. One
instance of the 40 does — see [V-E](#v-e-region-scale).

**Why `\|C_i\| <= 2` is an imposed cut, not a property of the deployment.** An
AP covering a UE offers it exactly one RB, allocated round robin by proximity
(`model_instance.py:131-136`), so a UE's candidate count equals its covering-AP
count before truncation. At `g=5, seed=1` that count runs 1 to 6, mean 3.66;
`max_deg=2` then leaves 48 UEs with two candidates and 2 with one. The cut buys
`l_i = 1`, no invalid codewords, and a statevector-certifiable circuit.

---

## V-B. Zone-local distribution validation

### Stated, and correct

| claim | measured | |
|---|---|---|
| five zones admit statevector checks | 5 zones of the `g=4` seed-2 partition | **P** |
| largest reaching six UEs and 21 qubits | zones 7, 8, 9 at `N_z=6`, 21 qubits | **P** |
| total-variation distance below `3e-15` | max **2.75e-15** | **P** |
| no infeasible assignment among accepted outcomes | 0 in every zone | **P** |
| the classical sampler reproduces the accepted law | `P_z` factorizes over UEs before screening; screening then gives `p_z` exactly | **P** |

Full result, `python certify_sampler.py`:

| zone | `N_z` | qubits | `\|F_z\|` | `mu_z` | TVD to `exp[lambda J_z]` | accepted outside `F_z` |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 6 | 1 | 0.5000 | 0.00e+00 | 0 |
| 3 | 2 | 10 | 4 | 0.3934 | 1.14e-16 | 0 |
| 7 | 6 | 21 | 36 | 0.0316 | 1.92e-15 | 0 |
| 8 | 6 | 21 | 36 | 0.0214 | 1.73e-15 | 0 |
| 9 | 6 | 21 | 36 | 0.0168 | 2.75e-15 | 0 |

### Measured, not stated

Sec. V says only that "transpiled gate counts and noise sweeps quantifying
these effects are provided in the repository". This is what is there.

`python noise_oracle.py --sweep`. Backend `FakeSherbrooke`, 127 qubits, median
rates **2q 7.79e-03, 1q 2.36e-04, readout 1.98e-02**, 10 000 shots per point.

Two questions are kept apart. **Routing is a transpile, not a simulation** —
there is nothing to simulate, and a 127-qubit device could not be simulated
anyway. **Noise is simulated at all-to-all connectivity**, on the circuit's own
qubit count, not the device's.

**Lattice cost (transpile only)** — **R**

| circuit | qubits | 2q before | 2q after | factor | depth before | depth after |
|---|---:|---:|---:|---:|---:|---:|
| inter-cell, 4 UE x 2 AP | 16 | 384 | 969 | 2.5x | 773 | 2717 |
| intra-cell, 3 node x 3 RB | 21 | 1354 | 3373 | 2.5x | 2036 | 8697 |

**Noise, inter-cell circuit** (`\|F\| = 6`, `lambda = 0.1976`) — **R**

| condition | accept | discard | TVD to `exp[lambda J]` | decoded infeasible | accepted mass vs noise-free |
|---|---:|---:|---:|---:|---:|
| noise-free | 0.666 | 0.334 | 0.005 | 0.0496 | — |
| FakeSherbrooke rates | 0.123 | 0.877 | 0.134 | 0.1228 | 18 % |
| depolarizing `p2=1e-04` | 0.650 | 0.350 | 0.011 | 0.0487 | 98 % |
| depolarizing `p2=1e-03` | 0.513 | 0.487 | 0.016 | 0.0733 | 77 % |
| depolarizing `p2=1e-02` | 0.104 | 0.896 | 0.137 | 0.1115 | 16 % |

Intra-cell noise-free reference: accept 0.209, discard 0.791, TVD 0.011,
decoded infeasible 0.2972. Its noisy rows take about 50 s each.

**How to read the table.** Noise thins the accepted branch and tilts the
accepted distribution; it never admits an infeasible assignment, because the
acceptance test applies to the *decoded* assignment rather than to the internal
flags. That is also why "decoded infeasible" is nonzero even noise-free
(0.0496): an infeasible assignment can leave every utility qubit in `|0>` and
still fail the decode check. `p2 = 1e-03`, roughly current hardware, holds the
distance at 0.016 with 77 % of the accepted mass — the regime the 2000-gate
budget targets.

### Not measured

| | |
|---|---|
| **Real-QPU execution.** Sec. V says "selected zone circuits are further executed on a quantum processor"; the abstract and Sec. I say the algorithm "is executed on a real NISQ processor". No such run is in this repository: every execution path is `AerSimulator` or a `qiskit_ibm_runtime.fake_provider` backend, and there is no `QiskitRuntimeService` call anywhere. If the run happened, its job IDs and counts are not in the artifact. | **—** |
| **Certification above six UEs.** The statevector check stops at 21 qubits, but the partition produces zones of up to 13 UEs. Nothing certifies the accepted law for a zone of the size actually used. | **—** |

---

## V-C. Boundary coordination in a partitioned region

Instance `g=5`, seed 0, `geo=True`, peers = seeds 0-15.

### Stated, and correct

| claim | measured | |
|---|---|---|
| `M=25` APs, `N=50` UEs | 25 / 50 | **P** |
| `\|Z\| = 14` zones | 14 | **P** |
| `\|B\| = 23`, 46 % of UEs | 23, 46 % | **P** |
| largest sampling distance 0.184 at `K_z=200` | 0.184 (4 displayed zones and all 14) | **P** |
| 19 % of the product's belief on excluded joint choices | 19.3 % and 19 % on the two pairs | **P** |
| re-conditioning wins 11 of 14 | 11 of 14 | **P** |
| improvement 0.9 percentage points | +0.87 pp mean | **P** |
| five selective re-sampling events | 5 | **P** |
| 98.7 % of `J*` | 98.66 % | **P** |
| mean of peers 99.2 % | 99.19 % | **P** |

### Stated, and wrong

| claim | measured | |
|---|---|---|
| `beta_z` ranges from 0.38 to 1.50 | **1.36 to 1.50** for this instance. 0.377 is the minimum over all zones of the 40-instance V-E sweep | **P!** |
| at least 88 of 200 samples effective in every zone | **198 of 200** for this instance. 87.7 is the same sweep-wide minimum | **P!** |
| peers range 97.3 % to 100.0 % | 97.34 % to **99.97 %**. No instance of this family reached `J*` | **P!** |

Both misattributions take a sweep-wide minimum into a paragraph describing one
instance. In `manuscript_260901.tex` the same pair sits in the Setup under
"across zones", where it is correct for the evaluation as a whole.

### Measured, not stated

| | value | |
|---|---|---|
| peer family size | 16 seeds, **13 globally feasible**; the other 3 are dropped | **R** |
| ablation family | `g=4..7` x seeds 0-3 = 16, **14 feasible** | **R** |
| largest zone circuit, this instance | 1910 two-qubit gates | **R** |
| sampling distance, mean | 0.157 over the 4 displayed zones, 0.115 over all 14 | **R** |
| executions | 19 224 initial, +11 973 from re-samples (+62 %) | **R** |
| executions per zone | 201 to 9 997 (cap 10 000) | **R** |
| ablation spread | median +0.53 pp, range −2.36 to +3.42 pp, 3 losses | **R** |

### The two most dependent boundary pairs

`strongest_pair` (`fig_snapshot.py:407`) ranks by total variation between the
empirical joint and the outer product of its own marginals.

**Zone Z13, UEs 12 & 14** — TV 0.385. Candidates `{RB52, RB38}` and
`{RB52, RB35}`, sharing RB52 under `W_r = 1`.

|  | UE14 = RB52 | UE14 = RB35 |
|---|---:|---:|
| **UE12 = RB52** | 0.000 (infeasible) | 0.307 |
| **UE12 = RB38** | 0.629 | 0.065 |

Marginals UE12 `[0.307, 0.693]`, UE14 `[0.629, 0.371]`; their product puts
**19.3 %** on the cell the joint report excludes.

**Zone Z0, UEs 30 & 38** — TV 0.377. Candidates `{RB85, RB81}` and
`{RB81, RB28}`, sharing RB81. Joint `[[0.555, 0.105], [0.000, 0.340]]`; the
product puts **19 %** on the excluded cell.

### Commitment order — **P!**

The appendix walkthrough has this backwards. From the decimation trace:

```
step  7:  UE12 -> index 1 = RB38   confidence 0.885   belief [0.115, 0.885]
step 19:  UE14 -> index 1 = RB35   confidence 0.533   belief [0.467, 0.533]
final:    UE12 -> RB38,  UE14 -> RB35
```

UE12 commits **first**, and UE14 ends on RB35, not RB52. The belief `b_i` is a
product over every zone holding that UE, so it is not the zone-local marginal
above, and twelve other commitments fall between the two steps. The point
survives either way: once UE12 holds RB38, filtering removes every draw in
which UE12 holds RB52, so the (RB52, RB52) cell the one-pass product favours
leaves the report before UE14 is considered.

### Not measured

| | |
|---|---|
| **A classical matching baseline.** Once AP admission is fixed, RB assignment inside an AP is a small bipartite matching, solvable exactly in classical polynomial time. Nothing here measures "quantum picks APs, classical picks RBs". See `ENCODING_NOTES.md`. | **—** |

---

## V-D. Zone-local circuit cost and latency

**Closed form. Nothing is run.** `fig/complexity/plot_qtg_assignment_runtime.py`
imports no project module.

### Stated, and correct

| claim | measured | |
|---|---|---|
| `t_layer = 12.5` ns per logical layer | `TAU_NS = 12.5` | **P** |
| `D_z = 2N_z`, `Q_z^cnt = 3`, `\|M_z\| = N_z` | matches the module constants | **P** |
| proposed enhanced ~ `113 N_z` ns | 114.4 at `N_z=160` | **P** |
| proposed serial ~ `381 N_z` ns | 382.6 at `N_z=160` | **P** |
| GAS `331N_z` and QTG `500N_z`, plus a common `37.5 N_z log2 N_z` | 607.8 against 605.6; 776.2 against 774.6 | **P** |
| crossovers 18 / 31 / 37 / 42 | 17.71 / 31.07 / 36.51 / 42.05 | **P** |

### Stated, and wrong

| claim | measured | |
|---|---|---|
| `rho = 0.399` | `CLASSICAL_COEFFICIENT_NS = 0.400`, and the figure legend prints the literal string `$0.400N_z^3$ ns` | **P!** |

### Measured, not stated

| | value | |
|---|---|---|
| `R_cnt` parallel counters, enhanced schedule | 4 | **R** |
| arithmetic depth factor, enhanced schedule | 0.75 | **R** |
| slope at small `N_z` | enhanced 133.5 at `N_z=10`, falling to 114.4 at 160 — the quoted `113N_z` is asymptotic | **R** |
| GAS / QTG slope at `N_z=10` | 477.6 / 640.7 ns per UE | **R** |
| what the `37.5 N_z log2 N_z` term is | the register accumulating a zone-wide numeric objective, which the proposed oracle does not carry | **R** |

### Not measured

| | |
|---|---|
| **GAS and QTG are not implemented.** Their curves come from depth formulas in the same module, not from built circuits. No transpiled gate count exists for either. | **—** |
| **`rho` is a point in a 13x-wide interval.** 0.400 sits inside 0.058-0.740, measured for the Hungarian method in the prior work. Every crossover in the figure inherits that width, and no sensitivity band is plotted. | **—** |

---

## V-E. Region-scale scalability

`scaling_data.npz`, `g = 3..9` x seeds 0-5.

### Stated, and correct

| claim | measured | |
|---|---|---|
| seven sizes, 40 instances | 42 built, 2 infeasible and dropped, 40 remain | **P** |
| `N = 18` to `162` | 18 to 162 | **P** |
| 5.3 to 45.0 zones | 5.3 to 45.0 | **P** |
| centralized 10 801 -> 173 983, factor 16 | x16.11 | **P** |
| largest zone 1 693 -> 1 950, factor 1.15 | x1.152 | **P** |
| traffic factor 9.0, zone count 8.4 | x8.97, x8.44 | **P** |
| per-zone traffic +6 % | 11 206 -> 11 916 bits, +6.3 % | **P** |
| mean utility 98.76 % to 99.32 % | exact | **P** |
| individual runs 96.79 % to 100 % | exact — the sweep does contain a 100.00 % | **P** |
| all 40 globally feasible | 40 of 40 | **P** |

### Stated, and wrong

| claim | measured | |
|---|---|---|
| "remains **below** the 2000-gate budget for every region size" | the plotted per-size mean does, but one instance's worst zone is **exactly 2000**. The constraint is `<=`, so this is legal; "never exceeds" is the exact wording | **P!** |
| "several instances at each size, giving 40 in total" | six seeds per size, 42 built, **2 dropped as globally infeasible** — the drop is unstated | **P!** |

### Measured, not stated

| | value | |
|---|---|---|
| the plotted "largest zone circuit" is the **mean over seeds** of each instance's largest zone, not the max | which is what makes it 1693 -> 1950 | **R** |
| worst single zone anywhere in the sweep | 2000 | **R** |
| boundary density per size | 61 / 54 / 55 / 66 / 63 / 65 / 64 % | **R** |
| per-size utility spread | see table | **R** |

Full sweep:

| `g` | UEs | zones | utility mean | min | max | largest zone 2q (mean) | centralized 2q | traffic KB | boundary |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 18 | 5.3 | 99.32 % | 96.79 | 100.00 | 1693 | 10 801 | 7.3 | 61 % |
| 4 | 32 | 8.7 | 99.13 % | 98.17 | 100.00 | 1786 | 24 227 | 11.2 | 54 % |
| 5 | 50 | 13.7 | 99.09 % | 97.74 | 99.69 | 1919 | 37 755 | 18.2 | 55 % |
| 6 | 72 | 20.3 | 98.77 % | 98.05 | 99.50 | 1871 | 64 120 | 28.8 | 66 % |
| 7 | 98 | 27.5 | 99.17 % | 98.75 | 99.51 | 1894 | 88 647 | 37.5 | 63 % |
| 8 | 128 | 36.2 | 98.76 % | 97.83 | 99.46 | 1926 | 138 690 | 52.6 | 65 % |
| 9 | 162 | 45.0 | 99.15 % | 98.80 | 99.59 | 1950 | 173 983 | 65.5 | 64 % |

Traffic model: `(1 + n_re) * K_z * ( sum_{i in B_z} ceil(log2 |C_i|) + 32 )`
bits per zone; the 32-bit field carries `J_z` for the reweighting.

### Not measured

| | |
|---|---|
| **A coordination-cost panel.** The sweep records `comm_bits` and `exceptions` per run, so the panel can be built from the existing cache without new experiments; V-E currently carries this in prose. | **—** |

---

## Exponent backoff (Sec. V-A, entirely unstated in the trimmed draft)

`python backoff_study.py`. 16 instances at `g=4..7` seeds 0-3, **281 zones**.
All **R**.

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
cap trades effective sample size for execution cost, not accuracy — mean
utility is unchanged to two decimals. The optimal round count reaching 887 for
the most tightly constrained zones against 0 for the sparsest is what motivates
`k = 1` as a practical setting rather than a per-zone optimum.

---

## Reference: which instance is which

The three experiment families use different instances, and the trimmed Sec. V
no longer says so.

| experiment | instance family | source |
|---|---|---|
| V-B exact check | `g=4`, seed 2 | `certify_sampler.py:177` |
| V-B noise | two hand-sized circuits | `noise_oracle.py:172` |
| V-C snapshot | `g=5`, seed 0, `geo=True`; peers seeds 0-15 | `fig_snapshot.py:63-64,169` |
| V-C ablation | `g=4..7` x seeds 0-3 = 16, 14 feasible | `fig_snapshot.py:80-81` |
| V-D | none; closed-form model | `fig/complexity/` |
| V-E sweep | `g=3..9` x seeds 0-5 = 42, 40 feasible | `fig_scaling.py:43-44` |
| V-A backoff | `g=4..7` x seeds 0-3 = 16, 281 zones | `backoff_study.py:35-36` |

"Five zones of the partition" in V-B is not five zones of the snapshot
partition.

## Reference: hardware anchors

`model_cost.py:33-42`. Two-qubit gates set the executable size, not qubit width;
`p_signal = exp(-lambda * n_2q)`.

| constant | value | meaning |
|---|---|---|
| `LAMBDA_KINGSTON` | 0.00118 | measured, `p ~ 0.53` at 640 heavy-hex 2q |
| `CEIL_KINGSTON` | 587 | 2q at the 50 % signal point, `ln2 / lambda` |
| `LAMBDA_BOSTON` | 0.00068 | scaled from the 2q error ratio (estimate) |
| `CEIL_BOSTON` | 1019 | same, for Boston |
| `BUDGET_DEFAULT` | 2000 | partition budget |

The 2000-gate budget sits **above** the measured 50 % ceiling. It marks roughly
the largest oracle a present device runs with a useful fraction of shots
error-free, not a noise-free circuit. See [Gaps](#gaps) for the sourcing
problem with these two constants.

---

## Corrections to apply

Replacement wording for the six passages marked **P!** above. Line references
are to the trimmed working draft; `manuscript_260901.tex` already carries C2,
C4 and C5 correctly, and has no C6.

**C1 — V-C, `beta_z` and effective sample size**

> Was: The local sampling exponents $\beta_z$ range from $0.38$ to the common
> target of $1.50$ and are reweighted before coordination. At least $88$ of the
> $200$ samples remain effective in every zone after reweighting…

> Replace: The local sampling exponents $\beta_z$ range from $1.36$ to the
> common target of $1.50$ and are reweighted before coordination. At least
> $198$ of the $200$ samples remain effective in every zone after reweighting…

**C2 — V-D, classical coefficient**

> Was: … with $\rho=0.399$ …  →  Replace: … with $\rho=0.400$ …

**C3 — V-C, peer range and excluded seeds**

> Was: Across 13 instances on the same map, the achieved utility ranges from
> $97.3\%$ to $100.0\%$, with a mean of $99.2\%$.

> Replace: The instance is drawn from a family of 16 seeds on the same map, 13
> of which are globally feasible. Across those 13, the achieved utility ranges
> from $97.3\%$ to $99.97\%$, with a mean of $99.2\%$.

Do **not** apply this to V-E: there the range really is 96.79–100 %.

**C4 — V-E, the budget is met, not undercut**

> Was: … and remains below the 2000-gate partition budget for every region size.

> Replace: … and never exceeds the 2000-gate partition budget at any region size.

and, same subsection:

> Was: The sampled region is varied over seven sizes with several instances at
> each size, giving 40 instances in total.

> Replace: The sampled region is varied over seven sizes with six seeds each;
> two of the 42 instances are globally infeasible and dropped, leaving 40.

**C5 — V-A, ninefold is a bound**

> Was: … giving approximately a ninefold reduction in the sampling cost for $k=1$.

> Replace: … so $(2k+1)^2=9$ bounds the reduction in sampling cost at $k=1$; the
> ratio $P_z(k)/\mu_z$ measured over the zones of this evaluation is $5.9$.

**C6 — Appendix C, commitment order** (and delete the `% FIXME` that follows it)

> Was: … When the consensus fixes UE $14$ to RB52, filtering by
> \eqref{eq:app-filter} retains the $0.63$ branch, the recomputed marginal of UE
> $12$ drops RB52 to zero, and the next commitment selects RB38 …

> Replace: … Decimation commits UE $12$ first, at confidence $0.89$, to RB38.
> Filtering by \eqref{eq:app-filter} then retains only the draws in which UE $12$
> holds RB38, so the pair (RB52,~RB52) favoured by the one-pass product leaves
> the zone's report before UE $14$ is considered. UE $14$ is committed later, at
> confidence $0.53$, to RB35, and the excluded pair never enters a commitment.

Everything else in that appendix paragraph is verified: the candidate sets, the
three joint probabilities 0.63 / 0.31 / 0.06, the RB52 marginals 0.31 and 0.63,
and the 19 %.

---

## Gaps

Ordered by what a reviewer would ask first.

**1. No real-QPU execution exists in this repository.** The abstract says the
framework is "evaluated through … quantum hardware execution"; Sec. I says the
algorithm "is executed on a real NISQ processor"; V-B says "selected zone
circuits are further executed on a quantum processor". Every execution path
here is `AerSimulator` (`circuit_inter.py:476`, `circuit_intra.py:461`,
`hybrid_assignment_example.py:182`) or a `qiskit_ibm_runtime.fake_provider`
backend (`noise_oracle.py:79`, `noise_zone.py:88`). There is no
`QiskitRuntimeService` call anywhere. A fake backend is a calibrated noise
model, not a device run. If the hardware run happened, its job IDs, backend
name, date and counts are not in the artifact, and the repository the paper
points at cannot support the sentence.

**2. `HARDWARE_LIMITS.md` is not in the repository.** `model_cost.py` cites it
three times as the source of `LAMBDA_KINGSTON = 0.00118` ("measured: p ~ 0.53
at 640 heavy-hex 2q" on ibm_kingston). That constant is what sets
`CEIL_KINGSTON = 587` and justifies the 2000-gate budget — the only
real-hardware measurement the whole cost argument rests on. The file it comes
from does not exist here.

**3. RB selection is not searched.** Sec. V evaluates joint AP–RB candidates,
but which RB an AP offers is fixed at instance-build time by the round robin in
`model_instance.py:131-136`, not chosen by the circuit. `W_r = 1` still binds,
because the offer wraps: one AP covers up to 14 UEs while owning 4 RBs, and of
80 RBs appearing as candidates, 17 are claimed by more than one UE. The joint
candidate set is a subset of `C_i` with RB multiplicity substituted by coverage
overlap. The generator's docstring says this; the manuscript does not. See
`ENCODING_NOTES.md`.

**4. Certification stops well below the zone sizes used.** The statevector
check reaches six UEs on 21 qubits; the partition produces zones of up to 13
UEs. Nothing certifies the accepted law at the size actually run.

**5. GAS and QTG are formulas, not circuits.** Their curves in Fig. 8 come from
depth expressions in the same plotting module. No circuit was built or
transpiled for either, so the comparison is between one implemented method and
two modelled ones.

**6. The classical crossover inherits a 13x-wide interval.** `rho = 0.400` is
an interior point of 0.058–0.740. Every crossover in Fig. 8 moves with it, and
no sensitivity band is drawn.

**7. `1,906` gates could not be reproduced.** `manuscript_260901.tex` V-D says
"the largest measured circuit costing 1,906 gates". The snapshot partition's
largest zone is 1910; the largest anywhere in the V-E sweep is 2000. Re-derive
before submission.

**8. A coordination-cost panel is missing.** V-E carries coordination cost in
prose while the sweep already records `comm_bits` and `exceptions` per run — a
panel needs no new experiments.

**9. Sec. III-C assumes the partition.** It states "partition the APs into
zones" while Sec. V cites the bounded-density regime as evidence.
`model_partition.py` produces that regime and is not described anywhere in the
paper, so the regime has no stated source.
