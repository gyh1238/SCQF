# Assignment encoding: what the qubits hold, and what it costs to widen them

Working notes on how a UE's choice is represented, why the evaluation sits at
`max_deg=2`, and what it would take to give a UE more RB choices. Numbers below
were measured on `make_instance(g=5, seed=1)` unless stated otherwise.

## 1. What an assign bit holds

The **index of a candidate**, not an AP or RB identifier. Eq. (state-width),
`manuscript_260901.tex:505-522`:

> A one-hot encoding of the binary variables `z_{i,j}` would use one qubit for
> each candidate. Instead, the proposed representation encodes the index of the
> selected candidate using `ceil(log2 |C_i|)` qubits.

    Q_z = sum_{i in U_z} ceil(log2 |C_i^(z)|)

The resource is recovered by lookup, not carried in the register:

    r = zone.cand[k][ci]        # code -> RB          (proto_zone.py:85)
    a = zone.rb_owner[r]        # RB   -> owning AP   (proto_zone.py:86)

`cand[i]` and `util[i]` are index-aligned (`model_instance.py:44-46`), which is
what makes the index sufficient. `enumerate_zone` therefore iterates over code
tuples, not over RB combinations.

In the Sec. V evaluation `|C_i| <= 2`, so `l_i = 1`: **one UE = one bit**,
`0 -> cand[i][0]`, `1 -> cand[i][1]`. `certify_sampler.py:56` asserts this.

Codewords with no candidate (when `|C_i|` is not a power of two) are killed by a
validity flag, Eq. (validity-flag) at `manuscript_260901.tex:623`; the builder
collects them as `novel` in `certify_sampler.py:68-73`. At `|C_i| = 2` the set is
empty and that path is inert.

The recurring `if c == 0: qc.x(code[k])` idiom is control polarity: a code qubit
controls on 1, so selecting value 0 means conjugating by X.

## 2. Zone is the execution unit, not AP

One zone = one circuit = one oracle `O_z`. APs are what gets *grouped into* a
zone, which is exactly what `model_partition.py` does under a two-qubit budget.
`manuscript_260901.tex:402`:

> Each zone contains a subset of APs and the UEs with assignment candidates in
> that zone.

The number of APs in a zone does not change the code width. It changes only the
number of `W_a` flags -- one per owned AP (`certify_sampler.py:57-58`). Measured
partition at g=5, seed=1:

| zone | APs | RBs | UEs | Q_z | ap flags |
|---|---:|---:|---:|---:|---:|
| Z13 | 1 | 4 | 9 | 9 | 1 |
| Z4  | 3 | 10 | 7 | 7 | 3 |
| Z0  | 3 | 6 | 4 | 4 | 3 |

Z4 holds three APs but seven UEs, so its code is seven bits; Z13 holds one AP
and nine UEs, so nine. The width tracks UEs only.

The counter `cnt` is reused across APs rather than duplicated, so extra APs cost
gates, not qubits -- which is why the partition budget is a two-qubit gate count
and not a qubit count.

What crosses a zone boundary is **classical**: boundary code indices plus a
32-bit `J_z` (`manuscript_260901.tex:927`). This is independent QPUs plus
classical coordination, not entanglement distribution. The paper does not claim
one physical QPU per zone; a zone is a circuit-decomposition unit.

## 3. Why each UE ends up with two candidates

It is imposed, not emergent. `model_instance.py:139-150` keeps the top `max_deg`
candidates by utility:

    order = np.argsort(-u)[:max_deg]        # "keep the strongest APs"

The truncation is binding. Covering APs per UE *before* truncation, radius 1.2:

| covering APs | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|
| UEs | 2 | 10 | 7 | 19 | 8 | 4 |

Mean 3.66, max 6. After truncation: 48 UEs with two candidates, 2 with one.

Candidate count equals covering-AP count because an AP offers each covered UE
**exactly one** RB from its pool, allocated round robin by proximity
(`model_instance.py:131-136`). Per the module's own comment, it is coverage
overlap -- not RB multiplicity inside one AP -- that creates boundary UEs.

`max_deg=2` is an evaluation choice, not physics: it yields `l_i = 1`, no invalid
codewords, and a statevector-certifiable circuit.

## 4. The RB is not searched in Sec. V

A candidate is already the pair `(a, r)`, so one bit fixes both.
`manuscript_260901.tex:858`:

> A candidate `j=(a,r) in C_i` consumes two capacities at once ... AP `a` is
> loaded by the demand `d_i` against `W_a`, and RB `r` by one unit against `W_r`.

Which RB an AP offers is decided at instance-build time by the round robin above,
so **the RB-within-AP choice is not a search variable in the scaling and snapshot
experiments**. The quantum circuit for RB selection exists separately
(`circuit_intra.py`, 2 bits/node, `00->RB0, 01->RB1, 10->RB2, 11->invalid`).

`W_r = 1` still bites, because `rank % len(pool)` wraps: one AP covers up to 14
UEs but owns 4 RBs. Of the 80 RBs appearing as candidates, 17 are claimed by more
than one UE (16 by two, 1 by three).

**Open point.** Sec. V presents joint AP-RB assignment while RB-within-AP is a
fixed offer. If raised in review, the answer is that the joint candidate set is a
subset of `C_i` and RB multiplicity is substituted by coverage overlap; the
generator's docstring says as much, but the manuscript does not.

## 5. Cost of widening the index

RB choices go into the *same* register, multiplicatively: `|C_i| = k*m` for `k`
APs and `m` RBs offered per AP, so `l_i = ceil(log2(k*m))`. Today `m = 1`, `k = 2`.

Qubits grow logarithmically; gates grow linearly, and the gates are what binds.
Eq. (zone_one_pass_cost), `manuscript_260901.tex:901`:

    T_z, G_z = O( D_z (l_max,z + Q_z^cnt) + |M_z| Q_z^cnt ),   D_z = sum_i |C_i|

Measured with `zone_2q_cost`, worst single AP at g=5, seed=1:

| max_deg | mean \|C_i\| | l_i | worst single-AP 2q |
|---:|---:|---:|---:|
| 2 | 1.96 | 1 | 1554 |
| 3 | 2.72 | 2 | 3490 |
| 4 | 3.34 | 2 | 3930 |
| 8 | 3.66 | 3 | 7946 |

The 1554 -> 3490 jump at `l = 1 -> 2` has two causes in `model_cost.py:96-105`:
each candidate pays `2*mcx_2q(l) + 2` for its match, and unused codewords
(`2**l - |C_i|`) each add a validity flag.

Adding APs grows both AP and RB flags; adding RBs inside one AP grows only RB
flags, so RB multiplicity is slightly cheaper in the `|M_z|` term -- but `D_z`
dominates, so the saving is small.

**Partition constraint.** An AP is the atom of the partition, so the budget must
exceed the most expensive single AP. At `l = 2` that is already 3490 against the
current 2000, so raising `m` forces the budget up and further above the measured
Heron2 ceiling (~587 2q at the 50% signal point, `model_cost.py:31-38`).

## 6. Cheaper encodings

Modeled with the same primitives (`mcx_2q`, `_qft_2q`) for one UE at `k=2`,
`m=4`, `|C_i| = 8`, counter width 3:

| encoding | qubits/UE | rotations | AP counter | RB counter | total |
|---|---:|---:|---:|---:|---:|
| flat binary (today's scheme, widened) | 3 | 304 | 312 | 312 | **928** |
| factored (AP field \| RB field) | 3 | 304 | 10 | 312 | **626** |
| factored + one-hot RB field | 5 | 112 | 10 | 120 | **242** |

**(a) Factor the code into an AP field and an RB field.** The `W_a` check does
not need to know which RB was taken, so the AP counter can be controlled on the
AP field alone and pays one accumulate per `(i,a)` pair instead of per `(i,a,r)`
edge -- `k` instead of `k*m`. Free in qubits, and it matches the paper's own
inter/intra split. Today the accumulate is controlled on the full codeword
(`model_cost.py:118-124`).

**(b) Make the RB field one-hot.** The binding resource is two-qubit gates, not
qubits -- `model_cost.py` header: "the executable size of these circuits is set by
the two-qubit gate count -- not by the qubit count", and `README_figs.md`: "qubit
width is never the binding constraint at these zone sizes". Trading qubits for
gates is therefore favourable here. One-hot turns each match from a
multi-controlled MCX into a single control: at `l=3`, `2*mcx_2q(3)+2 = 38` per
candidate becomes `2*1+2 = 4`.

Caveat: one-hot introduces an exactly-one condition. Checking it with a flag eats
the gain. Prepare the field with a chain of controlled `R_Y` rotations (a
W-state-shaped preparation) instead of `H^m` and exactly-one holds by
construction, no flag needed -- which also removes the `novel`/`val_flag` path in
`certify_sampler.py:68-73`. In the manuscript only the definition of `q_0` changes;
the accepted-branch derivation is untouched.

**(c) Replace fixed top-k with a utility-gap rule.** `model_instance.py:148`
currently gives every UE the same `max_deg` candidates, but many UEs have a
dominant first choice (12 UEs are covered by at most two APs). Keeping only
candidates within `delta` of the best makes `D_z` follow the mean rather than the
worst case. No circuit change; implementable immediately.

Suggested order: (c), then (a), then (b). (c) needs only a re-run; (a) is an
oracle change consistent with the existing exposition; (b) gives the largest
reduction but touches `q_0` preparation and the certification script.

## 7. Baseline worth measuring

Given `W_r = 1` and RBs of one AP differing only by `rb_gain`, once AP admission
is fixed the RB assignment inside an AP is a small bipartite matching solvable
exactly in classical polynomial time. A "quantum picks APs, classical picks RBs"
baseline should be measured before claiming anything for quantum RB selection.
The current code is effectively a weaker form of this -- a fixed proximity offer
rather than an optimal matching -- so the comparison is close at hand.

## Reproducing the measurements here

All at `g=5, seed=1`:

- coverage degrees before truncation: distances from `make_instance`, counted
  against `radius=1.2`
- per-zone widths: `partition_aps` + `build_zone`, reading `code_widths()`,
  `len(zone.aps)`, `len(zone.rbs)`
- cost vs `max_deg`: `max(zone_2q_cost(inst, [a]) for a in range(inst.n_ap))`

The table in Sec. 6 evaluates `mcx_2q`/`_qft_2q` directly rather than going
through `zone_2q_cost`, since the factored and one-hot variants are not
implemented.
