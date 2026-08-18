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
| `panels/*.png` | each panel alone, raster at 200 dpi, for quick viewing |
| `seed_preview.png` | candidate snapshot regions side by side (`preview_seeds.py`) |
| `campus.png` | basemap the snapshot region is cut from, and drawn over |
| `ap_allowed_mask.png` | where an AP may stand: rooftops and lots, roads removed |
| `ue_allowed_mask.png` | where a UE may stand: open ground only |
| `placement_mask_preview.png` | the two masks over the basemap, for checking the trace |

Panels:

| file | panel |
|---|---|
| `panels/snapshot_a_partition.*` | (a) the partition |
| `panels/snapshot_b1_zone_law_Z13.*` | (b1) zone Z13, 3 APs |
| `panels/snapshot_b2_zone_law_Z5.*` | (b2) zone Z5, 2 APs |
| `panels/snapshot_b3_zone_law_Z1.*` | (b3) zone Z1, 2 APs |
| `panels/snapshot_b4_zone_law_Z8.*` | (b4) zone Z8, 1 AP |
| `panels/snapshot_c1_joint_Z4_UE3_38.*` | (c1) joint vs. marginals, zone Z4 |
| `panels/snapshot_c2_joint_Z13_UE29_36.*` | (c2) joint vs. marginals, zone Z13 |
| `panels/snapshot_c3_ablation_joint_vs_marginals.*` | (c3) what dropping the joint list costs |
| `panels/snapshot_d_decimation_order.*` | (d) order the boundary UEs were fixed |
| `panels/scaling_a_circuit_cost.*` | (a) two-qubit gate count vs. zones |
| `panels/scaling_b_coordination_cost.*` | (b) classical report volume vs. zones |
| `panels/scaling_c_utility_ratio.*` | (c) utility ratio vs. zones |

The zone indices and UE indices in the panel names are those of the snapshot
instance (`g=5, seed=14`); they change if that instance changes, and each
script clears its own panels before writing so stale names cannot linger.

Regenerate everything with:

```bash
python make_fig_scaling.py               # --recollect to redo the sweep
python make_fig_snapshot.py
python preview_seeds.py 12               # only to re-choose the region
```

Both figure scripts write the composite and the panels in one run. Panels are
drawn by the same functions as the composite — `panel_map`, `panel_zone_law`,
`panel_boundary`, `panel_ablation`, `panel_order` for the snapshot,
`panel_cost`, `panel_comm`, `panel_quality` for the scaling figure — so
editing a panel function changes both.

Two caches hold the sweeps, since both take minutes to recompute:

| cache | holds | force a recompute |
|---|---|---|
| `../scaling_data.npz` | the 40-instance growth sweep behind `fig_scaling` | `--recollect` |
| `../ablation_data.npz` | the joint-vs-marginal comparison behind (c3): the 12 of 16 campus instances with a centralized optimum | delete the file |

How many panels each row shows is set by `N_ZONE_PANELS` (4) and
`N_BND_PANELS` (2) in `make_fig_snapshot.py`; the grid widths follow.

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
> available, since $2k+1$ passes of a 1–2 k gate circuit leave the coherence
> budget immediately. Quoting $k$ would suggest the protocol requires Grover
> iterations; it requires none. The panels report $\mu_z$ and the unamplified
> shot count, which is what a run actually pays. The asymptotic
> $1/\mu_z \rightarrow 1/\sqrt{\mu_z}$ gain from amplification is real, but it
> is not what these figures measure.

### Coordination (Sec. IV-C)

| symbol | meaning | value |
|---|---|---|
| $K_z$ | accepted draws each zone reports | 200 |
| $eta_z\leeta$ | exponent a zone actually executes; see Sec. 6 | 0.08 – 1.50 |
| $K_{\min}$ | effective-sample floor that triggers re-sampling | 25 |
| $\pi_{z,i}(v)$ | marginal of boundary UE $i$ held by zone $z$ (Laplace-smoothed, $\alpha=0.5$) | — |
| $b_i(v)\propto\prod_z\pi_{z,i}(v)$ | consensus belief over UE $i$'s candidates | — |
| confidence | $\max_v b_i(v)$; decimation commits the highest first | — |
| exception re-sample | a zone re-drawn because conditioning drove its ESS below $K_{\min}$ | 5 in the snapshot |

### Hardware and circuit cost (Sec. V-C)

| symbol | meaning | value |
|---|---|---|
| two-qubit gates | count for **one compute-mark-uncompute pass** of the oracle; amplification multiplies it by $2k+1$ | plotted |
| partition budget | the per-zone two-qubit ceiling the partitioner must respect | 2000 |
| Heron2 half-signal | measured on ibm_kingston: decay $0.00118$ per 2q gate $\Rightarrow$ 50 % signal at | 587 |
| Heron3 half-signal | estimated for ibm_boston: decay $0.00068$ per 2q gate $\Rightarrow$ 50 % signal at | 1019 |

The two device points are why the partition budget is a two-qubit count
rather than a qubit count. They are not drawn on the figures: the budget line
is the constraint the partitioner actually enforces, and adding device
ceilings beside it invited the figure to be read as a claim about which
processor runs this today, which is a separate question.

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
| `seed` | instance seed | 0 – 15 (14 in the snapshot) |
| `geo` | take positions from the campus rasters | on for `fig_snapshot`, off for `fig_scaling` |

Only `g` changes along the growth axis. Every density is held fixed, so
growing `g` grows the *global* problem without changing local structure.

### Campus placement (`geo=True`, `haiq_geo.py`)

The densities above are unchanged; what changes is *where* the points may
fall. One pixel unit is fixed at `PX_PER_UNIT` = 150 px per AP spacing, so a
`g × g` instance is a `g × g` unit window of the campus, centred on
`CENTER_PX` and grown about that point — `g` is still a pure size knob, the
ground under it is simply not flat.

- **APs**: the same jittered `g × g` lattice, each AP then snapped to the
  nearest pixel of `ap_allowed_mask` within `AP_SNAP_MAX` = 0.55 units, no two
  closer than `AP_MIN_SEP` = 0.35. A cell with nowhere legal to mount yields no
  AP, so `n_ap` can fall below `g²`.
- **UEs**: drawn uniformly over the *allowed area* of `ue_allowed_mask` rather
  than over the square, so they follow streets and courtyards.

Two consequences are worth knowing before reading the snapshot numbers
against the scaling figure:

- boundary density runs higher than on the square (≈ 60 % here against ≈ 50 %),
  because irregular AP spacing puts more UEs within reach of two zones;
- some seeds have **no centralized optimum at all** — UEs bunch onto open
  ground faster than the APs facing it can admit them, and `solve_centralized`
  returns infeasible. Both scripts skip those seeds and `preview_seeds.py`
  prints which ones went.

`python haiq_geo.py` re-derives `CENTER_PX`: it scans the campus for the
window holding open ground *and* rooftops at every size in use.

The scaling figure deliberately stays on the square. Its claim is about
bounded density, and a flat lattice states that assumption without borrowing
one campus's street plan; `scaling_data.npz` is therefore unaffected by any of
the above.

## 3. `fig_snapshot` — one region, opened up

Instance `g=5, seed=14`, cut from the campus: 25 APs, 50 UEs, 14 zones,
30 boundary UEs (60 %).

The seed is chosen for legibility — compact zones, none of them holding a
single UE, boundary links that can be traced. To keep that presentation choice
from turning into a quality one, the caption reports where this instance's
utility ratio falls among the candidates: 99.5 % against a candidate mean of
99.4 % and a spread of 98.4–100 % over the 9 of 16 seeds that have a
centralized optimum. `preview_seeds.py` renders the candidates side by side
and prints the numbers behind the choice.

### (a) the partition

- **axes**: the campus window the instance was cut from. No units; AP grid
  spacing is 1.0, which is 150 px of the rasters. Ticks are suppressed.
- **basemap**: the campus, blended 22 % to white. It is context, not a model
  quantity — but it is the reason the APs sit where they do, since each stands
  on ground `ap_allowed_mask` permits and each UE on ground `ue_allowed_mask`
  permits.
- **coloured outlines**: zone territory — the border of the region whose every
  point has its *nearest* AP in that zone, drawn in the zone's colour over a
  white line. Territories are outlined rather than filled here: a fill heavy
  enough to identify a zone also buries the ground under it, which is the one
  thing the map was put there to show. On the square (`GEO = False`) they are
  filled instead, there being nothing behind them to hide. Either way this is
  a reading aid for where a zone sits, not a model quantity; zone membership
  of a UE is set by $\mathcal{V}_i$, not by distance.
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
- **coloured bars**: the zone's own draws, carrying their reconstruction
  weights. A zone that executed at $eta_z<eta$ contributes draws from a
  flatter law; the weights are what put them back on the target exponent, so
  the bars and the reference curve are comparable either way.
- **black step**: the exact reference $p_z\propto e^{\lambda J_z}$.
- **shaded band on the right**: the infeasible region. It carries **exactly
  zero** mass — that is strict feasibility made visible, and it is the
  difference from a penalty formulation, which would place small but non-zero
  mass there.
- **title**: $N_z$, $\vert\mathcal{F}_z\vert$, $\mu_z$, and $K_z/\mu_z$ — the
  shots that collecting this zone's report costs without amplification.

Bars tracking the step curve is the claim that the sampler reproduces the
circuit's law. Departures are finite-sample noise at $K_z=200$, larger where
$\widetilde K_z$ is smaller; `haiq_certify.py` checks the same statement exactly against a
Qiskit statevector (max TVD $2.7\times10^{-15}$, zero spurious states).

### (c1)–(c3) why the boundary report is a joint list

The row answers one question — why a zone sends a list of retained *joint*
draws instead of one marginal per boundary UE — in two steps.

**(c1), (c2): what marginals lose.** Zones and pairs are chosen by measurement,
not by eye. For every zone, `strongest_pair` takes each pair of its boundary
UEs, forms the empirical joint over that pair and the outer product of the two
marginals, and scores the pair by the total-variation distance between them;
the zone keeps its worst pair, and the row shows the two worst zones. In the
displayed instance that ranking runs TV = 0.37 (Z4), 0.35 (Z13), 0.35 (Z1),
0.33, 0.32 … down to exactly 0.00 in the zones whose owned limits do not
couple their boundary UEs — the dependence appears precisely where the
constraints bind, which is why the panels are worth showing at all.

The TV value is not printed on the panel: it selects what to show, and the
red bar then says the same thing in a form that needs no scale.

- **x**: the four joint choices the two UEs can make, labelled by RB and
  ordered by the joint report, largest first. The impossible combination
  therefore lands last, where its bar is easiest to compare against nothing.
- **solid bars**: the zone's joint report, restricted to that pair.
- **hatched bars**: the outer product of that zone's own marginals — what a
  marginal or preference exchange would reconstruct.
- **red bar, marked "impossible"**: the product still backs a combination that
  is absent from $\mathcal{F}_z$ altogether — the two UEs would break an owned
  RB or AP limit. It carries about a fifth of the product's belief and exactly
  none of the joint's, which is why it stands alone with no solid bar beside
  it.

**(c3): what losing them costs.** 16 instances — region sizes `g` = 4, 5, 6, 7
crossed with seeds 0–3, set by `ABLATION_G` and `ABLATION_SEEDS`. This is a
smaller sweep than the one behind `fig_scaling`, because each point costs two
full protocol runs; it is not the same set, and the panel makes no claim about
how the gap varies with size. One point per instance. The horizontal axis
is the utility reached when each zone's marginals are taken once and not
revisited; the vertical axis is the utility reached when the retained joint
list is re-conditioned after every commitment. Points above the grey diagonal, labelled *equal*, are instances where keeping
the joint list helped; the corner note gives the mean gap and how many of the
16 improved. Everything else is held fixed
between the two runs — sampler, feasibility guard, commitment order — so the
gap is attributable to the report format alone.

Worth being precise here, because the panel invites a stronger reading than
it supports: **the manuscript does use marginals.** The belief
$b_i \propto \prod_z \pi_{z,i}$ is a product of marginals, and Sec. IV-C
calls $\pi_{z,i}$ exactly that. The difference the panel measures is not
marginals versus no marginals; it is whether those marginals are re-taken
from a list that has been conditioned on each commitment, or read once at the
start and left alone. The joint list is what makes re-taking them possible.

The mean gap is **+1.8 points** of utility, and all 12 instances improve. That is worth reading against the
total decomposition loss: the distributed result sits about 1 point below the
centralized optimum, so the joint list is worth roughly as much as the entire
remaining gap. This is the measured version of the manuscript's statement that
the correlations "re-enter through the conditioning performed after each
commitment".

The panels deliberately do not mark where decimation committed. Commitment
happens one UE at a time, later, after conditioning, and using a product
across the zones that hold that UE — not within one zone's pairwise joint. It
therefore need not land on the largest bar, and marking it invited exactly
that misreading.

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
deviation. All three panels share the growth axis.

**x (all panels)**: region size. The independent variable is the side of the
region; the zone count and the UE count are both consequences of growing it at
fixed density, and panel (c) labels both tick rows: 5 zones / 18 UEs up to 45
zones / 162 UEs, a nine-fold growth. Zones are the axis quantity because the
per-zone claims in (a) and (b) are indexed by them.

### (a) circuit complexity

- **y**: two-qubit gates for one oracle pass, log scale. The axis names the
  quantity and the metric on two lines, as (b) and (c) do.
- **red squares, dashed**: the centralized circuit — one circuit for the whole
  region. **×16** across the range.
- **green circles, solid**: the largest zone circuit, $\max_z$, with a ±1σ
  band. **×1.15** across the same range, i.e. flat.
- **green dashed line**: the 2000 partition budget the green curve respects.

Comparing at one oracle pass is deliberate and conservative. Neither curve
includes amplification, which the figures do not use; if it were added, the
centralized circuit would need *more* rounds than a zone because its accepted
mass is smaller, so the gap shown understates the real one.

### (b) communication volume

- **y**: the classical traffic the coordination stage needs, log scale.

A zone's report is its $K_z$ retained draws restricted to its boundary UEs:
$\sum_{i\in\mathcal{B}_z}\lceil\log_2ert\mathcal{V}_iert
ceil$ bits
of code words per draw, plus the one recorded utility $J_z$ that lets the merge
point reweight it. Every zone sends one, and a zone re-sampled after a
commitment sends its new list as well, so the count here is
$(1+	ext{re-samples}_z)\,K_z(\sum\ell_i + 32)$ bits summed over zones. This is
the whole of the inter-zone communication: the local stage exchanges nothing.

- **red squares, dashed**: the whole region. **×8.8** for **×8.4** zones — it
  tracks the zone count, not something worse.
- **green circles, solid**: what a single zone sends. **×1.05**, i.e. a zone's
  message does not notice how large the region became.

Note the asymmetry with (a): the green curve there is the **largest** zone,
$\max_z$, while the green curve here is the **mean** zone, total ÷ zones. That
is deliberate — a circuit has to fit on a device, so the worst zone is what
decides feasibility, whereas traffic is carried by a network and what matters
is the load it sees. But it does mean this panel does not bound the largest
single report, only the typical one.

This is the panel that answers the obvious objection to (a): bounding the
per-zone circuit is not interesting if the coordination it requires explodes
instead. It does not. Sec. V-D argues this in prose — that under bounded zone
density the largest $Q_z$ is independent of the number of zones while the
report grows only with the overlap — and this is that statement measured.

### (c) solution accuracy

- **y**: utility of the distributed result as a percentage of the centralized
  strict optimum, computed exactly by MILP (HiGHS) on the same instance.
- **grey line at 100 %**: the optimum. The axis runs down to 90 so that the
  distance to it is visible rather than magnified: the curve sits in the top
  tenth of the panel throughout.
- **blue triangles**: mean over seeds; the band is the observed range across
  them. A ratio cannot exceed 100 by construction, so a symmetric ±1σ band
  would have reached past the bound — the range cannot.

The centralized optimum appears only as the denominator; it is not a competing
protocol. Flat here means decomposition, the shot budget and boundary
coordination together cost about a point of utility, and that this does not
worsen as the region grows.

**Why the band is widest at 9 zones.** The loss is not a continuous quantity:
it is a small count of events, each one a UE pushed to a worse RB, and each
worth $100/ert\mathcal{U}ert$ percentage points. Measured in those units
the spread is about **one such event at every size** — 0.89, 1.08, 1.22, 1.31,
0.94, 1.23 from 9 zones upwards — so what shrinks along the axis is not the
error but the percentage one error is worth: 3.1 points at 32 UEs against 0.6
at 162. The smallest size is narrow for the opposite reason: at 18 UEs the
instance is easy enough that four of six seeds find the optimum exactly, so
the point sits at a ceiling rather than in the middle of a spread. The dip
from the first size to the second is that ceiling ending, not the method
degrading.

Read together: the work one processor does and the traffic one zone sends are
both unchanged across a nine-fold growth in problem size, while the
centralized circuit and the total traffic grow with it, at no cost in quality.

## 5. Editing notes

- Panels are vector. Text in the `.svg` files stays as text, so it can be
  restyled or translated without re-running anything.
- Colours: zones use `tab20`. The fixed roles are distributed/per-zone
  `#2e8b57`, centralized/whole-region `#c1440e`, utility ratio `#1f6fb4`.
  Boundary UEs, their AP links and the impossible bar reuse `#c1440e`, so the
  same red always means "this is the thing that costs you".
- **Fills are pre-blended, not transparent.** Every large tinted area — zone
  territories, the infeasible band, the confidence bands — is drawn opaque at
  a colour already mixed toward white by `_tint`. Reducing an alpha in the SVG
  will therefore do nothing; change the blend fraction instead. This exists
  because a partially transparent fill on a transparent background renders at
  full saturation in any viewer that flattens or ignores the alpha channel,
  which made the PNGs disagree with the PDFs. Over the campus the question
  does not arise: (a) outlines its territories instead of filling them, so
  there is no large tinted area in that panel at all. The outline colour is
  `_shade(zone, 0.25)` — darkened, not tinted, because a line has far less
  area than a fill to carry a colour with.
- Type sizes come from four constants at the top of `make_fig_snapshot.py` —
  `FS_TITLE` 8, `FS_LABEL` 7.5, `FS_TICK` 6.5, `FS_LEGEND` 6.2 — so panels
  cannot drift apart. The scaling figure sets its sizes inline, being three
  panels rather than eleven.
- Neither composite carries a title or a caption block: both belong to the
  document that places the figure. The numbers that would go in a caption are
  printed to stdout when the script runs.
- Transparent background means anything relying on white behind the text will
  need a background added at placement time. All text is dark, so a light
  backdrop is assumed.
- Composite layouts are set with `subplots_adjust`, not `tight_layout`, so
  panel positions are stable when labels change length.

## 6. Supplementary: how the reports are paid for

None of this appears on the figures; it is recorded here because the
numbers behind them depend on it. A tight zone has a small $\mu_z$, and at $k=0$ the cost of $K_z$ accepted draws
is $K_z/\mu_z$ — which for the tightest zones ran to $10^6$–$10^8$ shots, far
past any real run. The remedy is the manuscript's own: **a zone executes at the
largest exponent its shot budget allows, and the target exponent is restored
after measurement rather than before it.**

| symbol | meaning | value |
|---|---|---|
| shot budget | shots one zone may spend on its report | 10 000 |
| $\beta_z \le \beta$ | exponent zone $z$ actually executes, the largest with $K_z/\mu_z(\beta_z)\le$ budget | 0.08 – 1.50 |
| $w^{(k)}=\exp[(\lambda-\lambda_z)J_z^{(k)}]$ | reconstruction weight of draw $k$ | 1 when $\beta_z=\beta$ |
| $\widetilde K_z=(\sum w)^2/\sum w^2$ | effective sample size after reweighting | $\ge 73$ of $K_z=200$ |

A lower $\beta_z$ makes the accepted law flatter and therefore cheaper to hit;
the reweighting puts each draw back on the target exponent, and the price is
paid in effective sample size, not in bias. $\widetilde K_z$ is the honest
sample count behind a zone's report and is what the exception threshold
$K_{\min}$ is compared against.

**What the budget costs.** Over 16 instances, capping every zone at 10 000
shots against no cap at all:

| | utility vs. optimum | worst instance | peak shots in a zone |
|---|---|---|---|
| shot budget 10 000 | 99.26 % | 97.99 % | $10^4$ |
| no budget | 99.31 % | 98.72 % | $2.6\times10^8$ |

A 25 000-fold reduction in the worst zone's shot count costs 0.05 percentage
points of utility on average, and under 1 point on the worst instance. That is
the accuracy deliberately traded away to make the protocol executable, and it
is why the (b) panels are noisier than an unbudgeted run would give.
