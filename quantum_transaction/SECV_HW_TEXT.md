# Hardware results as manuscript text

What the runs in [HW_RUNS.md](HW_RUNS.md) let Sec. V say, and the edits that
say it. The LaTeX is in [SECV_HW_PATCH.tex](SECV_HW_PATCH.tex), numbered to
match the list below.

The draft claims hardware execution in four places - the abstract, Sec. I, the
third contribution and the conclusion - and Table I marks "Gate-level QPU
demo" for the proposed method. Those claims are now backed by 18 runs on three
devices. What changes is what they claim.

---

## What the runs support, and what they do not

| claim | status |
|---|---|
| The circuits were executed on real quantum processors | **yes** - ibm_boston, ibm_kingston, IonQ Forte Enterprise |
| Every reported assignment satisfies the constraints | **yes**, every device, every run |
| The accepted law of Eq. (36) is reproduced on hardware | **yes without amplification**: TVD 0.024-0.049 over six settings on ibm_boston |
| The optimal assignment is the most frequent report | **yes**, amplified or not, up to 1868 routed gates |
| One amplification round raises acceptance on hardware | **yes below ~400 routed gates** (0.337 against 0.298 at 364); no above |
| The full oracle reproduces the accepted law on hardware | **no** - TVD 0.159 at best, 0.248 at zone-half scale |
| The 2000-gate budget is executable today with an undistorted law | **no** - the measured break-even is ~355 gates in that counting |
| Error-mitigated aggregation would rescue the larger circuits | **no** - majority voting over five layouts moves P(opt) 0.229 to 0.234 |

Two numbers must never be quoted as performance: IonQ's **sharpened acceptance
of 0.927**, which exceeds the noise-free 0.836 and is an artefact of the
aggregation, and the **sharpened TVD**, which does not measure the accepted
law.

---

## The edits

1. **V-A**, the ninefold sentence. The bound is 9; the measured ratio is 5.9.
   Independent of the hardware work, but wrong as written.
2. **V-B**, new hardware paragraph (two paragraphs in the patch). The claim
   the abstract and Sec. I rest on.
3. **V-B**, routing paragraph. One clause on what the budget bounds: the
   circuit a zone must fit into, met in qubit width (37.7-43.0 against 156),
   not the size at which a present device returns an undistorted law.
4. **V-B**, after the simulated-noise table. The sweep models gate errors; it
   overestimates the surviving signal by 3.4x to 55x against the runs, because
   idling is what grows with size.
5. **V-D**, after the parallel-counter sentence. The DD result makes depth
   reduction a present-device measure, not only a latency one.
6. **Table caption**, V-B noise table: 2000 shots -> 10 000, which is what
   `noise_oracle.py` runs.
7. **Optional table** for the hardware runs, if they get their own float.
8. **Abstract**, scope of the reconstruction claim.
9. **Sec. I**, name the device families.

Unchanged: V-C, V-E and their numbers; the partition budget of 2000; Table I;
the conclusion. The brown draft already carries the V-C and V-E corrections.

## Why the budget is not restated

Lowering it does not produce smaller circuits. The partition's atom is one AP
with its UEs, which already costs 686 to 1854 gates in the counting of
Eq. (39) on the campus instance, so a budget below about 1000 returns the same
one-AP-per-zone partition. Measured on that finest partition, coordination
still reaches 98.1-98.6% of $J^\star$ against 99.6% at the 2000-gate budget,
which is the sensitivity worth reporting if a reviewer asks what the budget
buys.

---

## Answers to the questions a reviewer is likely to ask

**"Which circuit was executed?"** The inter-cell test circuit of Sec. V-B, at
three to six UEs, with and without one amplification round. It is not a zone
of the campus instance, and with two APs it decides both capacity constraints
from one counter and a two-sided comparison rather than the general per-AP
accumulation of Fig. 6. HW_RUNS.md says so in its first section.

**"Why is the amplified acceptance below the unamplified one at 787 gates?"**
Because the round costs more of the marked branch than it adds at that size.
The k=0 rows, which reach their exact values on the same device, place the
loss in the oracle rather than in the sampler or the readout, and the 364-gate
circuit shows the round paying off below the break-even.

**"Does the 2000-gate budget hold on these devices?"** Not for reproducing the
accepted law: the measured break-even is about 355 gates in that counting, and
about 519 with the parallel-counter schedule of Sec. V-D. It holds for what it
is used for, bounding the circuit each zone must fit, and selection survives
to at least half the budget.

**"Is any of this error mitigation?"** Dynamical decoupling is on, and it
matters: disabling it costs the 787-gate circuit its whole marked branch. No
other mitigation is used on IBM. The IonQ run was debiased, and the repository
publishes both of its aggregations.

**"Why 500 shots on IonQ and 10 000 on IBM?"** IonQ bills by shot. At 500
shots the sampling floor on the accepted law is about 0.05 and the acceptance
estimate carries +-0.03, both small against the effects reported.
