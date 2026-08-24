# 10. Complexity of the Unified Assignment Circuit

This section analyzes the joint AP--RB assignment circuit in the paper and compares it with approach-level GAS and QTG baselines. Inter-cell association and intra-cell RB allocation are not treated as separate algorithms: a candidate

$$
j=(a,r)
$$

selects an AP and an RB simultaneously. The distributed execution protocol is likewise not a fourth quantum curve; it partitions the same local circuit and reconciles the boundary variables after local sampling.

## 10.1 Comparison unit

All quantum curves refer to one logical search pass over the same zone-local joint assignment instance. The comparison includes state preparation, constraint arithmetic, marking, and the corresponding inverse operations required by the modeled pass. It excludes the algorithm-dependent number of GAS threshold updates, QTG amplification repetitions, accepted-sample repetitions, physical shots, readout, queueing, and classical boundary coordination.

The exclusions matter because GAS, QTG, and the proposed sampler do not return an answer through the same outer loop. Consequently, the graph is a comparison of per-pass circuit resources, not an end-to-end time-to-solution claim.

## 10.2 Notation for the joint problem

For a zone $z$, define

$$
N_z=|\mathcal U_z|,
\qquad
k_i=|\mathcal V_i|,
\qquad
D_z=\sum_{i\in\mathcal U_z} k_i,
$$

where $D_z$ is the number of UE--candidate edges evaluated by the zone. The compact candidate-index encoding uses

$$
\ell_i=\left\lceil\log_2 k_i\right\rceil,
\qquad
\ell_{\max}=\max_i\ell_i,
\qquad
Q_z=\sum_{i\in\mathcal U_z}\ell_i
$$

assignment qubits. Further let

$$
C_z=|\mathcal A_z|+|\mathcal R_z|
$$

be the number of capacity conditions, and let $b_m$ be the width of the counter for resource $m$. We use

$$
b_z=\max_m b_m.
$$

Each candidate edge $(i,j)$ affects $r_{i,j}$ capacity registers. For joint AP--RB assignment, it contributes to one AP load and one RB occupancy, hence

$$
r_{i,j}=2,
\qquad
H_z=\sum_{i,j}r_{i,j}=2D_z.
$$

The symbol $p_z$ denotes the width of a numeric objective-value register when a baseline requires one. The proposed circuit encodes utility by controlled rotations and therefore does not allocate this register.

## 10.3 Primitive resource model

To keep circuit multiplicities separate from a particular compiler, define the following primitive costs.

| Primitive | $T$-count | $T$-depth | Meaning |
|---|---:|---:|---|
| $\ell$-bit candidate match | $E_T(\ell)$ | $E_D(\ell)$ | Test whether UE $i$ selected candidate $j$ |
| Controlled $b$-bit addition | $A_T(b)$ | $A_D(b)$ | Add a candidate load to a resource counter |
| $b$-bit capacity comparison | $K_T(b)$ | $K_D(b)$ | Compare accumulated load with capacity |
| Utility rotation | $Y_T(\epsilon)$ | $Y_D(\epsilon)$ | Approximate controlled $R_Y$ to precision $\epsilon$ |
| $q$-control conjunction | $C_T(q)$ | $C_D(q)$ | Form the feasible superflag or reflection |

For approximate-QFT arithmetic and standard multi-control decompositions,

$$
E_T(\ell)=O(\ell),
\qquad
A_T(b)=O(b),
\qquad
K_T(b)=O(b),
$$

up to the chosen rotation-synthesis precision. A ripple implementation changes constants and workspace but not the edge multiplicities below.

## 10.4 Proposed joint AP--RB circuit

### 10.4.1 Logical width

A direct implementation contains the compact assignment register, one utility qubit and one validity flag per UE, one violation flag per capacity condition, one reused load counter, comparator workspace, a reusable match ancilla, and one superflag. Thus,

$$
\begin{aligned}
Q_{\mathrm{prop},z}
= {}& Q_z
+N_z
+N_z
+C_z
+b_z
+a_{\mathrm{cmp}}(b_z)
+O(1) \\
= {}& Q_z+2N_z+C_z+b_z+a_{\mathrm{cmp}}(b_z)+O(1).
\end{aligned}
$$

Flags may be batched and recomputed to reduce width. This trades additional depth for fewer simultaneously live flag qubits; it does not remove any logical constraint.

### 10.4.2 Compute--mark--uncompute cost

Every candidate edge is matched against its UE code, contributes one utility-controlled rotation, and updates the two resource types attached to the joint candidate. Retaining the primitive costs gives the following one-pass expression:

$$
\begin{aligned}
T_{\mathrm{prop},z}^{(1)}
= {}& 2D_zE_T(\ell_{\max})
+2H_zA_T(b_z)
+2C_zK_T(b_z) \\
&+2D_zY_T(\epsilon)
+T_{\mathrm{flag},z},
\end{aligned}
$$

where $T_{\mathrm{flag},z}$ contains the feasible conjunction and the search reflection. The factors of two account for reversible compute and uncompute in the modeled search pass. Since $H_z=2D_z$ for joint AP--RB candidates,

$$
T_{\mathrm{prop},z}^{(1)}
=O\!\left(
D_z\ell_{\max}
+D_zb_z
+C_zb_z
+D_zY_T(\epsilon)
+T_{\mathrm{flag},z}
\right).
$$

The corresponding serial depth has the same edge multiplicities:

$$
\begin{aligned}
D_{\mathrm{prop},z}^{(1)}
= {}& 2D_zE_D(\ell_{\max})
+2H_zA_D(b_z)
+2C_zK_D(b_z) \\
&+2D_zY_D(\epsilon)
+D_{\mathrm{flag},z}.
\end{aligned}
$$

Under bounded candidate degree, bounded resource density, and $b_z=O(\log N_z)$,

$$
D_z=\Theta(N_z),
\qquad
T_{\mathrm{prop},z}^{(1)}=O(N_z\log N_z),
\qquad
D_{\mathrm{prop},z}^{(1)}=O(N_z\log N_z).
$$

This order comes from candidate-edge matching and resource arithmetic. It is not a carry-lookahead population model.

### 10.4.3 Parallel accumulation

With $s$ independent counters, only mutually independent resource-accumulation work is distributed. Candidate matching, utility rotations, comparisons, the superflag, and the reflection remain. Write the serial depth as

$$
D_{\mathrm{prop},z}^{(1)}
=D_{\mathrm{match}}
+D_{\mathrm{util}}
+D_{\mathrm{acc}}
+D_{\mathrm{cmp}}
+D_{\mathrm{flag}}.
$$

Then the scheduled model is

$$
D_{\mathrm{prop},z}^{(s)}
=D_{\mathrm{match}}
+D_{\mathrm{util}}
+\frac{D_{\mathrm{acc}}}{s}
+D_{\mathrm{cmp}}
+D_{\mathrm{flag}},
$$

with approximately $s b_z$ counter qubits. Dividing the entire Proposed curve by $s$ would therefore be incorrect.

## 10.5 GAS baseline

GAS represents the candidate objective and constraint values in quantum value registers, marks feasible assignments that improve on the current threshold, measures a candidate, and then updates the threshold. For the direct candidate-edge adaptation of the joint problem, its width is schematically

$$
Q_{\mathrm{GAS},z}
=D_z+p_z+\sum_{m=1}^{C_z}b_m+Q_{\mathrm{work},z}.
$$

One threshold-oracle pass performs the same candidate and feasibility arithmetic as the joint problem and additionally accumulates the numeric objective and compares it with the threshold:

$$
\begin{aligned}
T_{\mathrm{GAS},z}^{\mathrm{oracle}}
=O\!\left(&
D_z\ell_{\max}
+H_zb_z
+C_zb_z \\
&+D_zp_z
+p_z
+T_{\mathrm{flag},z}
\right).
\end{aligned}
$$

If $R_{\mathrm{GAS},z}$ threshold-oracle calls are made, the algorithmic cost is

$$
T_{\mathrm{GAS},z}^{\mathrm{total}}
=R_{\mathrm{GAS},z}
T_{\mathrm{GAS},z}^{\mathrm{oracle}}.
$$

$R_{\mathrm{GAS},z}$ is instance- and stopping-rule-dependent. It is deliberately excluded from the graph so that the plotted GAS curve represents one threshold-oracle pass.

## 10.6 QTG baseline

QTG organizes feasible-state preparation as a sequence of item layers. At each layer, branching is controlled by remaining capacity; the selected branch updates the relevant remaining-capacity register and the profit register. The original single-knapsack form has one decision qubit per item. A direct joint-assignment extension over binary candidate edges therefore uses

$$
Q_{\mathrm{QTG},z}
=D_z+\sum_{m=1}^{C_z}b_m+p_z+Q_{\mathrm{work},z}.
$$

The extension must also enforce one selected candidate per UE and update both remaining resources associated with $j=(a,r)$. With $r=2$, the state-preparation depth has the structural form

$$
D_{\mathrm{QTG},z}^{\mathrm{prep}}
=O\!\left(
D_z\ell_{\max}
+D_z(2r\,b_z+p_z)
\right).
$$

The $2r\,b_z$ term distinguishes remaining-capacity checking from the subsequent register update for each affected capacity dimension. A search iteration uses the preparation and its inverse, and a complete QTG solve may repeat the iteration under amplitude amplification:

$$
D_{\mathrm{QTG},z}^{\mathrm{total}}
=R_{\mathrm{QTG},z}
\left(2D_{\mathrm{QTG},z}^{\mathrm{prep}}+D_{\mathrm{mark},z}\right).
$$

The graph excludes $R_{\mathrm{QTG},z}$ and plots one joint-capacity preparation/search pass. The supplied `qtg_inter_assignment.py` and `qtg_intra_assignment.py` are restricted adapters for separate subproblems; they are useful implementation references but are not separate baselines for the unified paper.

## 10.7 Approach-level comparison

| Approach | Decision representation | Objective handling | Constraint handling | One-pass dominant depth | Workspace characteristic | Repetition excluded from graph |
|---|---|---|---|---|---|---|
| Proposed unified circuit | Compact candidate index, $Q_z=\sum_i\ell_i$ | Candidate-controlled amplitude rotations | Complete assignment followed by AP and RB accumulation/comparison | $O(D_z(\ell_{\max}+b_z)+C_zb_z)$ | One reusable counter in the width-minimal form | Accepted-draw and amplification repetitions |
| GAS | Direct binary candidate edges in the plotted adaptation | Numeric objective register and adaptive threshold | Numeric constraint registers and comparisons | $O(D_z(\ell_{\max}+b_z+p_z)+C_zb_z)$ | Objective plus constraint-value registers | Adaptive threshold-oracle calls |
| QTG | One binary decision per candidate edge in the joint extension | Profit register updated layer by layer | Remaining-capacity-controlled branching and updates | $O(D_z(\ell_{\max}+2r b_z+p_z))$ | Multiple remaining-capacity registers plus profit | Amplitude-amplification repetitions |

The asymptotic orders can coincide when candidate degree, resource incidence, and numeric register widths are bounded in the same way. The approaches nevertheless have different constants, widths, and outer-loop semantics. Those differences are why the detailed expressions should not be collapsed into one common formula.

## 10.8 Reproducible trend model used by the graph

The graph uses the paper's joint-assignment profile with two candidates per UE and a bounded-density reference with one capacity condition per UE:

$$
k=2,
\qquad
N_z=\frac{D_z}{2},
\qquad
C_z=N_z,
\qquad
r=2,
$$

The capacity-counter width is determined by the largest local load range of a resource, rather than by the total number of UEs in the zone:

$$
b_m=
\left\lceil
\log_2\!\left(1+L_{m,z}^{\max}\right)
\right\rceil,
\qquad
L_{m,z}^{\max}
=\sum_{(i,j)\in E_{m,z}}c_{m,i,j}.
$$

For the capacity-aware reference curve, the bounded local incidence is represented by

$$
b=3,
\qquad
\ell=\log_2 k=1.
$$

GAS and QTG additionally maintain a zone-wide numeric objective or profit register. Its continuous trend width is

$$
p=\log_2(N_z+1).
$$

An arithmetic-depth factor

$$
\eta=0.75
$$

is applied to controlled additions, comparisons, and numeric objective/profit updates. It represents the scheduled arithmetic decomposition used by the plotted scenario. Candidate matching, utility rotations, conjunctions, and reflections are not multiplied by $\eta$. Lower-cost addition and control decompositions that motivate this resource model include [temporary logical-AND adders](https://quantum-journal.org/papers/q-2018-06-18-74/), [relative-phase Toffoli substitutions](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.93.022311), and [conditionally clean-ancilla comparators](https://quantum-journal.org/papers/q-2025-05-21-1752/).

Define the common flag and reflection term as

$$
F_z=2\log_2(N_z+C_z+1).
$$

The width-minimal serial Proposed curve is

$$
\widetilde D_{\mathrm{prop,ser}}
=2D_z\ell
+2D_z
+\eta\left(2rD_zb+2C_zb\right)
+F_z.
$$

The first two terms count forward/inverse candidate matching and utility rotations. In the enhanced schedule, match and utility operations acting on different UE registers form parallel UE layers. With $k=2$ and $s=4$ resource counters, the depth proxy becomes

$$
\widetilde D_{\mathrm{prop,enh}}
=2k\ell
+2k
+\eta\left(
\frac{2rD_zb}{s}
+2C_zb
\right)
+F_z.
$$

Only the accumulation term is divided by $s$; the comparison and flag terms remain intact. The GAS threshold-oracle proxy is

$$
\widetilde D_{\mathrm{GAS}}
=2D_z\ell
+\eta\left(
2rD_zb
+2D_zp
+2(C_z+1)b
\right)
+F_z,
$$

and the joint-capacity QTG search-pass proxy is

$$
\widetilde D_{\mathrm{QTG}}
=2D_z\ell
+\eta\left[
2D_z(2r\,b+p)
\right]
+F_z.
$$

The display converts each quantum depth to a common logical-latency scale:


$$
\lambda_z=\tau\widetilde D_z,
\qquad
\tau=12.5\ \mathrm{ns}.
$$

The representative classical scenario is

$$
t_{\mathrm{class}}(N_z)
=0.400N_z^3\ \mathrm{ns}.
$$

The coefficient $0.400$ is a rounded representative value within the previously considered $0.058$--$0.740$ interval. It is not obtained by rescaling the quantum problem variable. The graph uses continuous $N_z$ only for smooth curves and determines each crossing by log--log interpolation.

| Curve | Modeled crossover $N_z$ | Corresponding $D_z=2N_z$ |
|---|---:|---:|
| Proposed enhanced, capacity-aware $s=4$ | $17.7$ | $35.4$ |
| Proposed serial accumulation | $31.1$ | $62.1$ |
| GAS threshold-oracle pass | $36.5$ | $73.0$ |
| QTG joint-capacity search pass | $42.1$ | $84.1$ |

GAS threshold repetitions, QTG amplification repetitions, accepted-draw repetitions, physical shots, readout, queueing, and boundary coordination are not multiplied into these one-pass curves.

## 10.9 Distributed protocol cost

The zone partition does not change the local quantum formula; it changes the largest $D_z$ that any device must execute. If $K_z$ accepted draws are retained and boundary UE $i$ has $k_i$ candidates, the assignment portion of zone $z$'s report is

$$
B_z^{\mathrm{assign}}
=K_z\sum_{i\in\mathcal B_z}\left\lceil\log_2 k_i\right\rceil
$$

bits, plus the stored local utility and protocol metadata. Confidence-ordered decimation makes at most $|\mathcal B|$ irreversible boundary decisions. Recomputing affected empirical marginals and filtering the retained lists is classical post-processing and should not be folded into any of the three quantum depth curves.

## 10.10 Interpretation for the paper

Under the capacity-aware bounded-density model, the enhanced Proposed schedule reaches the representative classical curve first, followed by the serial Proposed circuit, GAS, and QTG. The separation comes from compact candidate-index encoding, amplitude-based utility representation without a numeric objective register, UE-parallel matching and utility layers, and resource-banked accumulation. The $s=4$ factor applies only to accumulation; the remaining Proposed terms retain their explicit costs.

The reported crossings are model results for the stated one-pass schedule. An end-to-end runtime comparison additionally requires the measured oracle or amplification counts, accepted draws, transpiled depth, physical shots, and matched classical-solver runtimes.
