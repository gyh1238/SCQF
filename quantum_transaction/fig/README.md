# Figures: what is plotted, and what every symbol means

Two composite figures, plus every panel again as a standalone editable file.
All output is written on a **transparent background**, so the panels drop onto
a slide or a paper column without a white box.

## 1. Files

| file | contents |
|---|---|
| `fig_snapshot.pdf` / `.png` | composite: one region, partitioned and solved |
| `fig_scaling.pdf` / `.png` | composite: what stays bounded as the region grows |
| `panels/*.pdf` | each panel alone, vector, for LaTeX `\includegraphics` |
| `panels/*.svg` | each panel alone, vector, for editing in Inkscape/Illustrator |

Panels:

| file | panel |
|---|---|
| `panels/snapshot_a_partition.*` | (a) the partition |
| `panels/snapshot_b1_zone_law_Z2.*` | (b1) zone Z2 accepted law |
| `panels/snapshot_b2_zone_law_Z4.*` | (b2) zone Z4 accepted law |
| `panels/snapshot_b3_zone_law_Z12.*` | (b3) zone Z12 accepted law |
| `panels/snapshot_b4_zone_law_Z3.*` | (b4) zone Z3 accepted law |
| `panels/snapshot_c1_boundary_UE27.*` | (c1) boundary UE 27, zones Z0+Z2 |
| `panels/snapshot_c2_boundary_UE40.*` | (c2) boundary UE 40, zones Z5+Z6 |
| `panels/snapshot_c3_boundary_UE14.*` | (c3) boundary UE 14, zones Z8+Z9 |
| `panels/snapshot_d_decimation_order.*` | (d) order the boundary UEs were fixed |
| `panels/scaling_a_circuit_cost.*` | (a) two-qubit gate count vs. zones |
| `panels/scaling_b_utility_ratio.*` | (b) utility ratio vs. zones |

The zone indices and UE indices in the panel names are those of the snapshot
instance (`g=5, seed=2`); they change if that instance changes.

Regenerate everything with:

```bash
python make_fig_scaling.py     # add --recollect to redo the sweep
python make_fig_snapshot.py
```

Both scripts write the composite and the panels in one run. Panels are drawn
by the same functions as the composite (`panel_map`, `panel_zone_law`,
`panel_boundary`, `panel_order`, `panel_cost`, `panel_quality`), so editing a
panel function changes both.

## 2. Symbols

### Network model (manuscript Sec. III)

| symbol | meaning | value in these figures |
|---|---|---|
| $\mathcal{U}$, $\vert\mathcal{U}\vert$ | UE set and its size | 18 – 162 |
| $\mathcal{A}$ | AP set; each AP $a$ owns a disjoint RB pool $\mathcal{R}_a$ | 9 – 81 |
| $\mathcal{R}$, $a(r)$ | RB set; owner AP of RB $r$ | 4 RBs per AP |
| $\mathcal{V}_i \subseteq \mathcal{R}$ | candidate RBs of UE $i$ (coverage-dependent) | $\vert\mathcal{V}_i\vert \le 2$ |
| $y_{i,r}\in\{0,1\}$ | UE $i$ selects RB $r$; exactly one per UE | decision variable |
| $u_{i,r}$ | utility of assigning UE $i$ to RB $r$ | distance-decaying + per-RB jitter |
| $W_r$ | RB access limit (per RB) | 1 |
| $W_a$ | AP admission limit (per AP) | 4 |
| $J=\sum_i\sum_r u_{i,r}y_{i,r}$ | global utility, the objective | — |

Selecting an RB also fixes the serving AP, so the RB access limit and the AP
admission limit constrain the same choice.

### Zones (Sec. III-C)

| symbol | meaning | value |
|---|---|---|
| $\mathcal{Z}$ | zone set | 5 – 46 per instance; 5.3 – 45.0 as plotted means |
| $\mathcal{A}_z$ | APs owned by zone $z$ — **more than one in general** | 1 – 4, median 2 |
| $\mathcal{R}_z=\bigcup_{a\in\mathcal{A}_z}\mathcal{R}_a$ | RBs owned by zone $z$ | — |
| $\mathcal{U}_z=\{i:\mathcal{V}_i\cap\mathcal{R}_z\neq\emptyset\}$ | UEs the zone concerns | — |
| $N_z=\vert\mathcal{U}_z\vert$ | zone size | 1 – 13, median 6 |
| $\mathcal{B}$ | boundary UEs: those appearing in $\ge 2$ zones | 44 – 82 % of UEs |
| $\mathcal{B}_z=\mathcal{U}_z\cap\mathcal{B}$ | boundary UEs held by zone $z$ | — |
| $\ell_i=\lceil\log_2\vert\mathcal{V}_i\vert\rceil$ | code width of UE $i$ | 1 |
| $Q_z=\sum_{i\in\mathcal{U}_z}\ell_i$ | state-register width of zone $z$ | 1 – 13 qubits, median 6 |
| $D_z=\sum_{i\in\mathcal{U}_z}\vert\mathcal{V}_i\vert$ | zone-local candidate edges | — |
| $J_z$ | owner-local utility: only RBs owned by $z$ contribute | — |
| $\phi_z$ | strict feasibility of a local assignment (all owned RB and AP limits, all codes valid) | — |
| $\mathcal{F}_z=\{\bm r_z:\phi_z=1\}$ | feasible set of zone $z$ | 48 – 96 in the panels shown |

### Circuit and its output law (Sec. IV)

| symbol | meaning |
|---|---|
| $\lambda$ | Gibbs exponent of the utility weight |
| $\bar u$ | utility scale of the model; $\lambda=\beta/\bar u$ |
| $\beta$ | dimensionless exponent, shared by every zone — **1.5 in these figures** |
| $g_{z,i}(r)=\exp[-\lambda(u^{\max}_{z,i}-u^{(z)}_{i,r})]$ | per-candidate rotation weight |
| $\theta_{z,i,r}=2\arccos\sqrt{g_{z,i}(r)}$ | controlled-$R_Y$ angle |
| $p_z(\bm r_z)\propto e^{\lambda J_z}$ on $\mathcal{F}_z$ | the accepted-branch law — what panels (b) show |
| $\mu_z$ | acceptance mass: probability one circuit execution lands in the accepted branch, before amplification |
| $k$ | amplitude-amplification rounds. **These figures use $k=0$**; see below |
| $P_z(k)=\sin^2[(2k+1)\arcsin\sqrt{\mu_z}]$ | acceptance probability after $k$ rounds; $P_z(0)=\mu_z$ |

A smaller $\mu_z$ means a tighter zone and a more expensive draw: collecting
$K_z$ accepted samples costs $\approx K_z/P_z(k)$ shots, which at $k=0$ is
simply $K_z/\mu_z$.

> **Why the panels report shots and not rounds.** Amplification does not change
> the output law: it moves probability between the accepted and rejected
> branches and nothing else, so a measurement landing in the accepted branch
> follows $p_z$ for every $k$, including $k=0$. It only trades circuit depth
> for shot count — and on the devices measured here that trade is not yet
> available. The optimal schedule for the tightest zone in the snapshot
> ($\mu_z=0.0009$) is $k=26$, i.e. $2k+1=53$ passes and $54{,}166$ two-qubit
> gates, against a measured half-signal point of 587 (Heron2) or 1019
> (Heron3). Quoting $k$ would suggest the protocol requires 26 Grover
> iterations; it does not require any. The panels therefore report $\mu_z$ and
> the unamplified shot count, which is what a run actually pays. The
> asymptotic $1/\mu_z \rightarrow 1/\sqrt{\mu_z}$ gain from amplification is
> real, but it is not what these figures measure.

### Coordination (Sec. IV-C)

| symbol | meaning | value |
|---|---|---|
| $K_z$ | accepted draws each zone reports | 400 (scaling), 2000 (snapshot) |
| $K_{\min}$ | retained-list floor that triggers re-sampling | 40 |
| $\pi_{z,i}(v)$ | marginal of boundary UE $i$ held by zone $z$ (Laplace-smoothed, $\alpha=0.5$) | — |
| $b_i(v)\propto\prod_z\pi_{z,i}(v)$ | consensus belief over UE $i$'s candidates | — |
| confidence | $\max_v b_i(v)$; decimation commits the highest first | — |
| exception re-sample | a zone re-drawn because conditioning emptied its list | 1 in the snapshot |

### Hardware and circuit cost (Sec. V-C)

| symbol | meaning | value |
|---|---|---|
| two-qubit gates | count for **one compute-mark-uncompute pass** of the oracle; amplification multiplies it by $2k+1$ | plotted |
| partition budget | the per-zone two-qubit ceiling the partitioner must respect | 2000 |
| Heron2 half-signal | measured on ibm_kingston: decay $0.00118$ per 2q gate $\Rightarrow$ 50 % signal at | 587 |
| Heron3 half-signal | estimated for ibm_boston: decay $0.00068$ per 2q gate $\Rightarrow$ 50 % signal at | 1019 |

> **Symbol clash worth knowing.** `HARDWARE_LIMITS.md` uses $\lambda$ for the
> per-gate decay rate, while the manuscript uses $\lambda$ for the Gibbs
> exponent. They are unrelated. This file writes the hardware one as "decay".

The budget must exceed the cost of the most expensive *single* AP (1906 over
the instances used) because an AP is the atom of the partition.

### Instance generator

| parameter | meaning | value |
|---|---|---|
| `g` | region side; APs on a jittered `g × g` unit grid | 3 – 9 (5 in the snapshot) |
| `ue_per_ap` | UE density | 2.0 |
| `n_rb_per_ap` | RBs per AP | 4 |
| `radius` | coverage radius (AP spacing is 1.0) | 1.2 |
| `max_deg` | candidates kept per UE, strongest APs first | 2 |
| `seed` | instance seed | 0 – 5 (2 in the snapshot) |

Only `g` changes along the growth axis. Every density is held fixed, so
growing `g` grows the *global* problem without changing local structure.

## 3. `fig_snapshot` — one region, opened up

Instance `g=5, seed=2`: 25 APs, 50 UEs, 13 zones, 26 boundary UEs (52 %).
It is the **median** of 8 seeds by utility ratio, so the picture is typical
rather than selected.

### (a) the partition

- **axes**: the plane. No units; AP grid spacing is 1.0. Ticks are suppressed.
- **shading**: zone territory — each point takes the colour of the zone owning
  its *nearest* AP. This is a reading aid for where a zone sits, not a model
  quantity; zone membership of a UE is set by $\mathcal{V}_i$, not by distance.
- **white lines**: territory borders.
- **squares**: APs, coloured by zone.
- **open orange circles**: boundary UEs.
- **small filled dots**: interior UEs, coloured by their single zone.
- **thin orange lines**: each boundary UE joined to the APs that offer it a
  candidate RB. A link crossing a border is exactly the coupling that the
  coordination stage has to resolve.
- **zone label** `Zn / Nu Qq / G`: zone index, $N_z$ UEs, $Q_z$ state qubits,
  $G$ two-qubit gates. Every $G$ is below the 2000 budget, so the label makes
  the panel self-checking.

### (b1)–(b4) accepted law of four zones

The four zones with the largest $\vert\mathcal{F}_z\vert$.

- **x**: local assignments of that zone, ranked by $J_z$, best first. Rank
  runs over $\mathcal{F}_z$ only.
- **y**: probability of being the accepted outcome.
- **coloured bars**: measured frequency over the $K_z=2000$ accepted draws.
- **black step**: the exact reference $p_z\propto e^{\lambda J_z}$.
- **shaded band on the right**: the infeasible region. It carries **exactly
  zero** mass — that is strict feasibility made visible, and it is the
  difference from a penalty formulation, which would place small but non-zero
  mass there.
- **title**: $N_z$, $\vert\mathcal{F}_z\vert$, $\mu_z$, and $K_z/\mu_z$ — the
  shots that collecting this zone's report costs without amplification.

Bars tracking the step curve is the claim that the sampler reproduces the
circuit's law; `haiq_certify.py` checks the same statement exactly against a
Qiskit statevector (max TVD $2.7\times10^{-15}$, zero spurious states).

### (c1)–(c3) boundary coordination

The three boundary UEs with the *lowest* confidence, i.e. the hardest calls.

- **x**: that UE's candidate RBs, labelled by RB index.
- **y**: probability.
- **coloured bars**: $\pi_{z,i}$, one per owning zone, in that zone's colour.
- **black bar**: the consensus belief $b_i\propto\prod_z\pi_{z,i}$.
- **orange band**: the value decimation committed.
- **title**: UE index, the two owning zones, and the confidence.

Where the two coloured bars disagree, the panel shows two zones pulling in
opposite directions and the product resolving it. A scalar-preference exchange
would transmit one number per UE and could not express this.

### (d) the order decimation fixed the boundary UEs

- **x**: the boundary UEs, ordered by when decimation fixed them, 1 …
  $\vert\mathcal{B}\vert$. **One point is one UE, not one iteration.** Each
  boundary UE is committed exactly once, so the axis says which overlaps were
  resolved early — it is not a count of repeated passes over anything.
  Repetition, where it occurs, is counted separately as the exception
  re-samples named in the title.
- **y**: the confidence $\max_v b_i(v)$ of that UE, taken when decimation
  opened it for commitment.

The downward trend is the procedure working as intended: the least ambiguous
boundary UEs are fixed first, and each commitment conditions the retained
lists before the next choice. The title reports how many exception re-samples
the run needed.

## 4. `fig_scaling` — growing the region

40 instances: 7 region sizes × 6 seeds, minus 2 globally infeasible draws.
Each marker is the mean over the seeds at one size; shading is one standard
deviation. Both panels share the growth axis.

**x (both panels)**: the number of zones the partitioner produced, averaged
over seeds. The lower tick row on panel (b) gives the corresponding UE count,
so the axis reads in both units: 5 zones / 18 UEs up to 45 zones / 162 UEs.

### (a) circuit cost

- **y**: two-qubit gates for one oracle pass, log scale.
- **red squares, dashed**: the centralized circuit — one circuit for the whole
  region. The plotted means grow 10.8 k → 174 k.
- **green circles, solid**: the largest zone circuit, $\max_z$, with a ±1 σ
  band. Flat at 1.7 k – 2.0 k; no single zone exceeded the 2000 budget.
- **green dashed line**: the 2000 partition budget the green curve respects.
- **grey dotted / dash-dot lines**: the measured Heron2 and Heron3 half-signal
  points, for scale.

Comparing at one oracle pass is deliberate and conservative. Neither curve
includes amplification, which the figures do not use; if it were added, the
centralized circuit would need *more* rounds than a zone because its accepted
mass is smaller, so the gap shown understates the real one.

### (b) utility ratio

- **y**: utility of the distributed result as a percentage of the centralized
  strict optimum, computed exactly by MILP (HiGHS) on the same instance.
- **grey line at 100 %**: the optimum.
- **blue triangles with band**: mean ± 1 σ over seeds, 99.1 – 99.7 %. Single
  instances span 98.1 – 100.0 %.

The centralized optimum appears only as the denominator; it is not a competing
protocol. Flat here means decomposition and boundary coordination cost about
1 % of utility, and that this does not worsen as the region grows.

Read together: the centralized circuit leaves the executable region while the
zone circuits stay inside it, at unchanged solution quality.

## 5. Editing notes

- Panels are vector. Text in the `.svg` files stays as text, so it can be
  restyled or translated without re-running anything.
- Colours: zones use `tab20`. The fixed roles are utility ratio `#1f6fb4`,
  zone cost `#2e8b57`, centralized cost `#c1440e`, hardware ceilings `#8a8a8a`.
  Boundary UEs and committed values reuse `#c1440e`.
- Transparent background means anything relying on white behind the text will
  need a background added at placement time. All text is dark, so a light
  backdrop is assumed.
- Composite layouts are set with `subplots_adjust`, not `tight_layout`, so
  panel positions are stable when labels change length.
