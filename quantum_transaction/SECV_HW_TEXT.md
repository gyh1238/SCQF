# Hardware results as manuscript text

Replacement wording for the passages that claim quantum-processor execution,
written against the runs recorded in [HW_RUNS.md](HW_RUNS.md). Each snippet is
followed by the measurement it rests on, so a reviewer's question can be
answered from a stored job rather than from memory.

The draft already claims hardware execution in four places — the abstract,
Sec. I, the third contribution and the conclusion — and Table I marks
"Gate-level QPU demo" for the proposed method. Those claims are now backed;
what changes is *what* they claim.

---

## What the runs support, and what they do not

| claim | status |
|---|---|
| The circuits were executed on real quantum processors | **yes** — ibm_kingston, ibm_boston, IonQ Forte Enterprise |
| Every reported assignment satisfies the constraints | **yes**, on every device |
| The accepted law of \eqref{eq:conditional-exponential} is reproduced on hardware | **yes at k=0**: TVD 0.030 and 0.049 over 10 000 shots on ibm_boston |
| The optimal assignment is the most frequent report | **yes**, on every device, amplified or not |
| Hardware selects a zone assignment at near-exact quality | **yes**, under IonQ's debiased aggregation: P(opt) 0.46 against 0.444 exactly |
| One amplification round raises acceptance on hardware | **no** — 0.125 against 0.138 at k=0 on ibm_boston |
| The full oracle circuit reproduces the accepted law on hardware | **no** — TVD 0.159 at best (ibm_boston, unequal demands) |

Two numbers must never be quoted as performance: IonQ's **sharpened
acceptance of 0.927**, which exceeds the noise-free 0.836 and is an artefact
of the aggregation, and the **sharpened TVD**, which does not measure the
accepted law. Use P(opt) and mean utility for sharpened output.

---

## Sec. V-B, replacing the quantum-processor paragraph

> Selected circuits are further executed on quantum processors: the
> superconducting devices ibm\_kingston and ibm\_boston, and the trapped-ion
> device IonQ Forte Enterprise, whose all-to-all connectivity requires no
> routing and executes the inter-cell circuit in $196$ two-qubit gates against
> $787$ after routing onto a heavy-hex lattice. Every reported assignment
> satisfies the original constraints on every device, since the acceptance
> test of Section \ref{sec:distribution} is applied to the decoded assignment
> rather than to the internal flags. Without amplification, ibm\_boston
> reproduces the accepted law of \eqref{eq:conditional-exponential}: over
> $10{,}000$ shots its distance to the exact law is $0.030$ for unit demands
> and $0.049$ for demands $(1,2,1,2)$, the optimal assignment is reported with
> probability $0.609$ against $0.632$ exactly, and an accepted sample carries
> $0.870\,J^\star$ against $0.878\,J^\star$. One amplification round adds
> $787$ two-qubit gates, and at present error rates the round costs more
> acceptance than it returns, so the accepted law is flattened; the selection
> the sampler is meant to make nonetheless survives, and under the debiased
> aggregation IonQ applies to trapped-ion jobs the optimum is reported with
> probability $0.46$ against $0.444$ exactly, an accepted sample carrying
> $0.806\,J^\star$. Job identifiers, device settings, transpiled gate counts
> and raw counts are provided in \cite{scqfrepo}.

Rests on: the six hardware rows of HW_RUNS.md, and the aggregation comparison
of job `01a0b994-…`.

**If a shorter paragraph is wanted**, the last sentence before the citation
can be dropped; the amplification sentence should not be, because the
repository publishes both the k=0 and k=1 runs.

---

## Table row for the noise table

The noise table of Sec. V-B reports the simulated conditions. One hardware row
per device, with the same two columns, extends it:

```latex
\midrule
ibm\_boston, no amplification$^b$ & 0.100 & 0.030 \\
ibm\_boston, one round$^b$        & 0.055 & 0.269 \\
IonQ Forte Enterprise$^{b,c}$     & 0.099 & 0.257 \\
```

with the footnotes

```latex
$^b$Measured on hardware; the simulated rows above use the unrouted circuit,
whereas a hardware run executes the routed one ($767$ two-qubit gates on
ibm\_boston, $196$ on IonQ).
$^c$Debiased, averaged aggregation. Under the plurality aggregation the same
job reports the optimum with probability $0.46$; see \cite{scqfrepo}.
```

The simulated and hardware rows are not the same circuit, which is why the
footnote is needed: comparing them directly would mix gate noise with routing.

---

## Abstract, Sec. I and the conclusion

These three passages already say the algorithm runs on hardware. Keep them,
and let them say what was measured:

* **Abstract.** "quantum hardware execution" stands. If the sentence lists
  what the results verify, keep "exact constraint satisfaction" (hardware
  supports it) and attach "reconstruction of the distribution over feasible
  assignments" to the simulation results, since on hardware it holds for the
  unamplified sampler only.
* **Sec. I.** "executed on a real NISQ processor" stands. Naming the devices
  costs four words and answers the first question a reviewer asks:
  "is executed on superconducting and trapped-ion processors".
* **Third contribution.** "executed on a real quantum processor to verify its
  operation under actual hardware conditions" stands as written.
* **Conclusion.** "Qiskit-based simulations and quantum-processor experiments"
  stands.

## Table I

"Gate-level QPU demo: $\bigcirc$" for the proposed method is now backed by
executions on three devices across two qubit technologies. No change.

---

## Answers to the questions a reviewer is likely to ask

**"Which circuit was executed?"** The inter-cell test circuit of Sec. V-B:
four UEs, two APs, $\mathcal{C}_z=\mathcal{G}_z\mathcal{P}_z$ with one
amplification round, plus the same circuit without amplification. It is not a
zone of the campus instance, and with two APs it decides both capacity
constraints from one counter and a two-sided comparison rather than the
general per-AP accumulation of Fig. \ref{fig:oracle-capacity}. HW_RUNS.md says
so in its first section.

**"Why is the amplified acceptance below the unamplified one?"** Because at
$767$–$787$ two-qubit gates present error rates remove more of the marked
branch than one round adds. The k=0 rows, which reach their exact values on
the same device, locate the loss in the oracle rather than in the sampler or
the readout.

**"Is the IonQ number produced by post-processing?"** Partly, and the
repository publishes both aggregations of the same job. The averaged
aggregation is the device's own law (acceptance 0.099, TVD 0.257); the
plurality aggregation is what IonQ applies to debiased jobs and is quoted only
through P(opt) and mean utility, never through acceptance.

**"Does the 2000-gate partition budget hold on these devices?"** No, and the
paper should not imply it. The measured decay for the executed circuit is
about $0.0042$ per two-qubit gate on ibm_kingston, against the $0.00118$ that
`model_cost.py` assumes; the budget is a fault-tolerant-era target, consistent
with the timing model of Sec. V-D. HW_RUNS.md records the discrepancy.

**"Why 500 shots on IonQ and 10 000 on IBM?"** IonQ bills by shot. At 500
shots the sampling floor on the accepted law is about 0.05, and the acceptance
estimate carries $\pm 0.03$; both are small against the effects being
reported.
