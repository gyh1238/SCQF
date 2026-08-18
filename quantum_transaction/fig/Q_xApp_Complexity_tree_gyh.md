# 10. Complexity

This section derives projected logical-runtime curves for the Q-xApp circuits used for traffic steering (TS), network energy saving (NES), and QoS-based resource allocation (QoS-RA). The model preserves the constraint logic of the source circuits but replaces the ripple propagation inside each controlled population-counter update with quantum carry-lookahead arithmetic. Resource-local capacity checks are scheduled in parallel, whereas updates that target the same counter remain sequential.

The resulting curves describe an optimized logical circuit schedule. They are not Aer wall-clock measurements and are not claimed to be the transpiled depth of the current Python implementation.

## 10.1 Carry-lookahead population model

Let $N$ denote the number of entities being assigned, $M$ the number of candidate O-RUs, and $R$ the number of candidate DRBs. All logarithms are base 2. The exact address and counter widths are

$$
k=\left\lceil\log_2 M\right\rceil,
\qquad
w=\left\lceil\log_2(N+1)\right\rceil.
$$

The current counter update propagates a carry through a $w$-bit register, giving $O(w)$ depth per controlled increment. A quantum carry-lookahead increment instead forms the carry information through a prefix network. It therefore has

$$
D_{\mathrm{inc}}(w)=O(\!\log w)
$$

depth and requires $O(w)$ carry ancillas. The assignment predicates are not converted into $N$ independent population registers. Within each O-RU or DRB lane, the $N$ entities still update one counter sequentially.

For the smooth runtime curves, the logarithmic carry network is represented by

$$
L(w)=\log_2w+1.
$$

Counting $N$ assignment predicates and then clearing the counter gives the population term

$$
\begin{aligned}
P_{\mathrm{CLA}}(N)
&=2N L\!\left(\log_2(N+1)\right)\\
&=2N\left[\log_2\!\left(\log_2(N+1)\right)+1\right].
\end{aligned}
$$

Consequently,

$$
P_{\mathrm{CLA}}(N)=O(N\log\log N).
$$

The factor $2$ accounts for computation and uncomputation. The additive constant in $L(w)$ represents the control and cleanup stage retained in the smoothed schedule. Exact Clifford+T constants depend on the chosen controlled-adder synthesis; the expression above is the normalized depth proxy used consistently in the plotted curves.

### Parallelism boundary and workspace

Different resource counters are independent and may be evaluated in parallel. Updates within the same counter cannot be parallelized because they share a target register. The assumed schedule is therefore

```text
resource lanes:       parallel
UE updates per lane:  sequential
carry bits per update: prefix-parallel
```

For TS, $M$ resource-local counter and carry workspaces require

$$
Q_{\mathrm{work,TS}}=O(Mw)=O(M\log N)
$$

ancillas. The corresponding QoS-RA generalization uses $O(Rw)$ workspace. This is substantially smaller than a full UE-level population tree, which would require $O(MN)$ or $O(RN)$ intermediate workspaces. Coherent distribution of shared address controls adds fan-out ancillas and logarithmic scheduling depth; those terms are included among the non-population terms below.

For comparison, the four relevant schedules are

| Counter schedule | Population depth | Workspace scale | Status in this model |
|---|---:|---:|---|
| Globally reused ripple counter | $O(MN\log N)$ | $O(\log N)$ | Current source-style schedule |
| Resource-local ripple counters | $O(N\log N)$ | $O(M\log N)$ | Feasible intermediate schedule |
| Resource-local carry-lookahead counters | $O(N\log\log N)$ | $O(M\log N)$ | Runtime-curve schedule |
| Full UE population tree | $O(\log^2N)$ | $O(MN)$ | Excluded |

The quantum runtime is obtained from

$$
t_Q(N)=\tau D(N),
$$

using the projected enhanced logical-layer latency

$$
\tau=12.5\,\mathrm{ns}.
$$

This latency is a forward-looking constant-factor assumption. It is not attributed to the carry-lookahead transformation itself.

## 10.2 Circuit-derived runtime curves

### Traffic steering

For the smooth TS curve,

$$
M=\max\!\left(2,\frac{N}{10}\right),
\qquad
k=\log_2M.
$$

Replacing the counter term in the scheduled TS oracle by $P_{\mathrm{CLA}}(N)$ gives

$$
\begin{aligned}
D_{\mathrm{TS}}(N)
={}&P_{\mathrm{CLA}}(N)
+4\log_2 k
+2\log_2N
+2k\\
&+2\log_2(M+N)
+2\log_2(NM)\\
&+2\log_2(Nk-1)
+4.
\end{aligned}
$$

The carry-lookahead population term dominates the remaining logarithmic address matching, comparison, aggregation, utility-control, and reflection terms. Therefore,

$$
D_{\mathrm{TS}}(N)=O(N\log\log N),
$$

and

$$
t_{\mathrm{TS}}(N)
=12.5D_{\mathrm{TS}}(N)\ \mathrm{ns}.
$$

### Network energy saving

In `dqna_42.py`, `quality_oracle()` calls `_flag_feasible()` before and after utility marking. Each call computes and clears the population counter. Preserving these two feasibility passes gives

$$
D_{\mathrm{NES}}(N)
=2P_{\mathrm{CLA}}(N)
+4\log_2\log_2N
+2\log_2(N+1)
+2\log_2(N-1)
+4.
$$

Thus,

$$
D_{\mathrm{NES}}(N)=O(N\log\log N),
$$

with

$$
t_{\mathrm{NES}}(N)
=12.5D_{\mathrm{NES}}(N)\ \mathrm{ns}.
$$

### QoS-based resource allocation

The present `dqna_qos.py` circuit implements a fixed two-UE distinctness test. The scaling curve generalizes this predicate to $N$ UEs and $R=N$ candidate DRBs, with DRB-local occupancy counters evaluated in parallel. Using carry-lookahead increments inside each DRB lane gives

$$
D_{\mathrm{QoS}}(N)
=P_{\mathrm{CLA}}(N)
+4\log_2\log_2N
+10\log_2N
+2\log_2\!\left(N\log_2N-1\right)
+6.
$$

Therefore,

$$
D_{\mathrm{QoS}}(N)=O(N\log\log N),
$$

and

$$
t_{\mathrm{QoS}}(N)
=12.5D_{\mathrm{QoS}}(N)\ \mathrm{ns}.
$$

## 10.3 Clifford+T interpretation

The control logic is counted using balanced relative-phase Toffoli networks. An $r$-controlled logical operation is represented by

$$
\Delta_r\simeq2\left\lceil\log_2r\right\rceil
$$

non-Clifford layers, while a reversible comparison over a $w$-bit population register contributes approximately

$$
D_{\mathrm{cmp}}(w)\simeq2w.
$$

The arithmetic change concerns the carry path inside each controlled increment:

$$
O(w)\quad\longrightarrow\quad O(\log w).
$$

Because the same lane counter is updated by $N$ predicates, the UE factor remains. Computation and uncomputation therefore produce

$$
P_{\mathrm{CLA}}(N)
=2N\left[\log_2\!\left(\log_2(N+1)\right)+1\right]
=O(N\log\log N).
$$

This model does not claim constant-depth counting. With bounded-fan-in gates, forming a carry that depends on all lower counter bits still requires logarithmic depth in the register width. Hardware connectivity, routing, magic-state throughput, and control fan-out can increase the physical depth beyond the logical schedule.

## 10.4 Runtime comparison and crossover

The three Q-xApp curves derived in Section 10.2 are compared with two references that serve different purposes. The quantized curve represents the earlier discretized assignment approach and is retained only to show its faster growth with problem size. The bounded-Hungarian curve represents a conventional assignment solver and is used as the primary reference for determining the quantum--classical crossover.

| Reference | Role in the comparison | Dominant growth |
|---|---|---:|
| Quantized assignment | Contextual comparison with the earlier discretized formulation | $O(N^2\log\log N)$ |
| Bounded Hungarian | Primary reference for the crossover calculation | $O(N^3)$ |

The empirical trendline retained for the quantized approach is

$$
t_{\mathrm{quant}}(N)
=9.98N^2
\log_2\!\left(\log_2\frac{N}{10}\right)
-27.4N+1196
\ \mathrm{ns}.
$$

This curve is plotted only as a scaling reference and is not used to obtain the reported crossover. In particular, it should not be extrapolated below the range in which the fitted runtime remains physically meaningful.

### Bounded-Hungarian reference

The classical comparison models the assignment of $N$ source entities to $N$ physical resource slots augmented by $N$ dummy slots. The resulting cost matrix has $N$ rows and $2N$ columns. Applying the representative rectangular work model $r^2c$ to the reference coefficient $0.182\,\mathrm{ns}$ gives

$$
\begin{aligned}
t_{\mathrm{class}}(N)
&=0.182N^2(2N)\\
&=0.364N^3\ \mathrm{ns}.
\end{aligned}
$$

The value $0.364\,\mathrm{ns}$ is used as a representative comparison coefficient rather than a hardware-independent constant. It lies between the measured coefficients of continuous-cost and tie-heavy dense assignment instances and also reflects the additional work introduced by the physical and dummy resource columns.

### Crossover result

The crossover is obtained by comparing each quantum runtime from Section 10.2 directly with $t_{\mathrm{class}}(N)$. With the smooth carry-lookahead schedule, $\tau=12.5\,\mathrm{ns}$, and $M=\max(2,N/10)$, the resulting intersections are

| Q-xApp curve | Continuous intersection | First integer with $t_Q<t_{\mathrm{class}}$ | Quantum runtime | Classical runtime |
|---|---:|---:|---:|---:|
| TS | $N\simeq17.10$ | $N=18$ | $1.908\,\mu\mathrm{s}$ | $2.123\,\mu\mathrm{s}$ |
| NES | $N\simeq22.00$ | $N=23$ | $4.061\,\mu\mathrm{s}$ | $4.429\,\mu\mathrm{s}$ |
| QoS-RA | $N\simeq18.46$ | $N=19$ | $2.346\,\mu\mathrm{s}$ | $2.497\,\mu\mathrm{s}$ |

The continuous curves therefore place the onset of the modeled advantage at approximately $18$--$22$ UEs. Under a strict integer comparison, the first advantageous problem sizes are 18 UEs for TS, 23 UEs for NES, and 19 UEs for QoS-RA. The apparent one-UE discrepancy for NES is only a rounding effect: at $N=22$, the two modeled runtimes differ by less than $0.1\,\mathrm{ns}$.

![Carry-lookahead Q-xApp runtime curves](qxapp_complexity_cla_trendlines.png)

The figure should thus be read as a projected crossover against the bounded-Hungarian reference; the quantized curve provides scaling context but does not determine the reported intersection.

## 10.5 Relation to the source circuits

### Traffic steering (`dqna_ts.py`, `dqna_modes.py`, `dqna_constraints.py`)

`UnitCountCapacityConstraint` currently reuses the same counter sequentially across cells, and `ctrl_increment()` uses a ripple-style sequence of multi-controlled operations. The plotted model preserves the capacity predicate, comparison, joint mark, cleanup, and assignment reflection, but applies two scheduling changes:

1. independent O-RU capacity checks receive resource-local counter lanes; and
2. each ripple increment is replaced by a controlled carry-lookahead increment.

The resulting curve is therefore an optimized realization of the implemented oracle logic, not a direct depth measurement of the current source circuit.

### Network energy saving (`dqna_42.py`)

`compute_count()` and `uncompute_count()` update the cell population sequentially. The projected realization replaces the internal carry propagation while preserving the source order over UEs and the two `_flag_feasible()` calls in `quality_oracle()`. This is why the NES equation contains $2P_{\mathrm{CLA}}(N)$.

### QoS-based resource allocation (`dqna_qos.py`)

The source uses an in-place two-UE XOR distinctness test. The scaling equation is explicitly a generalized $N=R$ resource-allocation circuit in which each DRB has a local occupancy lane. It should therefore be described as a projected extension of the implemented two-UE logic.

## 10.6 Scope of the claim

The runtime comparison supports the following statement:

> With resource-local capacity lanes and logarithmic-depth carry-lookahead increments, the projected Q-xApp logical runtime scales as $O(N\log\log N)$ and intersects the representative bounded-Hungarian reference at approximately 18--22 assigned entities under the enhanced $12.5\,\mathrm{ns}$ logical-layer assumption.

The claim does not imply that:

- the current Python circuits already implement carry-lookahead counters;
- the present source has transpiled depth $O(N\log\log N)$;
- $12.5\,\mathrm{ns}$ is a measured hardware cycle time;
- $0.364N^3$ is universal across CPUs and assignment implementations; or
- a full $O(\log^2N)$ UE population tree is used.

The carry-lookahead arithmetic basis follows T. G. Draper, S. A. Kutin, E. M. Rains, and K. M. Svore, [“A Logarithmic-Depth Quantum Carry-Lookahead Adder”](https://arxiv.org/abs/quant-ph/0406142), *Quantum Information and Computation*, vol. 6, no. 4--5, pp. 351--369, 2006. Controlled hard-coded increment/decrement use is also discussed by E. Campbell, A. Khurana, and A. Montanaro, [“Applying Quantum Algorithms to Constraint Satisfaction Problems”](https://quantum-journal.org/papers/q-2019-07-18-167/), *Quantum*, vol. 3, p. 167, 2019.
