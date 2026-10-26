"""
The zone sampler on IonQ through Amazon Braket (Sec. V-B).
==========================================================
The Braket route to the same IonQ processors `ionq_oracle.py` reaches
directly: same circuits, same transpilation, same acceptance test
(`noise_oracle.law_from_counts`), same record format in `hw_runs/`.  It
exists because Braket bills the AWS account, not the IonQ project budget.

Braket prices IonQ per task and per shot (Forte: $0.30 + $0.08 per shot),
not per gate, so a shallow circuit costs the same as a deep one.  Debiasing
on Braket needs at least 2500 shots; without it any count is accepted.

Run:
    python braket_oracle.py --dry-run --only inter --demands 1,2,1,2 --cap 4 --k 0 --shots 500
    python braket_oracle.py --device Forte-1 --confirm-qpu --no-wait \\
        --only inter --demands 1,2,1,2 --cap 4 --k 0 --shots 500
    python braket_oracle.py --task <arn>          # collect a submitted task

AWS credentials: the usual boto3 chain (env AWS_ACCESS_KEY_ID /
AWS_SECRET_ACCESS_KEY, or ~/.aws/credentials), region us-east-1.
"""

import argparse
import glob
import json
import os
from datetime import datetime, timezone

from qiskit import transpile

from hw_oracle import OUT_DIR, case_key, n2q, report, save, select
from ionq_oracle import BASES, measured_map
from noise_oracle import law_from_counts

REGION = "us-east-1"
ARN = "arn:aws:braket:us-east-1::device/qpu/ionq/{}"
PRICE = {"task": 0.30, "shot": 0.08}          # Forte, USD
DEBIAS_MIN_SHOTS = 2500


def to_braket(tqc):
    """The transpiled circuit without its measurements, as a Braket circuit.

    Braket measures every qubit at the end, so the clbits are recovered from
    the full bitstrings with the circuit's own measure map, as for IonQ.
    """
    from qiskit_braket_provider.providers.adapter import to_braket as conv
    return conv(tqc.remove_final_measurements(inplace=False))


def marginal(counts, cmap):
    """Full-register Braket counts -> counts on the clbits, qiskit order.

    Braket bitstrings list qubit 0 first; qiskit keys list the last clbit
    first, which is the order `law_from_counts` reads.
    """
    out = {}
    for bits, n in counts.items():
        key = "".join(bits[q] for q in reversed(cmap))
        out[key] = out.get(key, 0) + n
    return out


def collect(record, task, meta, feas):
    res = task.result()
    counts = marginal(dict(res.measurement_counts), record["measure_map"])
    record["metrics"] = r = law_from_counts(counts, meta, feas)
    record["counts"] = counts
    record["task"] = {k: str(v) for k, v in task.metadata().items()
                      if k in ("status", "deviceArn", "shots", "createdAt",
                               "endedAt", "failureReason")}
    report("braket", r)
    save(record)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="Forte-1",
                    help="IonQ device name on Braket: Forte-1 | Forte-Enterprise-1")
    ap.add_argument("--shots", type=int)
    ap.add_argument("--k", type=int, default=1)
    ap.add_argument("--only", choices=("inter", "intra"))
    ap.add_argument("--basis", choices=tuple(BASES), default="rzz")
    ap.add_argument("--demands")
    ap.add_argument("--cap", type=int)
    ap.add_argument("--debias", action="store_true",
                    help=f"IonQ debiasing; Braket requires >= {DEBIAS_MIN_SHOTS} shots")
    ap.add_argument("--dry-run", action="store_true",
                    help="run on Braket's local simulator, free and offline")
    ap.add_argument("--no-wait", action="store_true")
    ap.add_argument("--task", help="collect a submitted task from its saved record")
    ap.add_argument("--confirm-qpu", action="store_true",
                    help="required to submit to a QPU (billed per task and shot)")
    a = ap.parse_args()

    if a.task:
        from braket.aws import AwsQuantumTask
        tid = a.task.rsplit("/", 1)[-1]
        path = glob.glob(os.path.join(OUT_DIR, f"*_{tid}.json"))[0]
        with open(path, encoding="utf-8") as f:
            record = json.load(f)
        (name, qc, meta, feas), = select(record["case"].split("-")[0],
                                         record["k"], record.get("demands"),
                                         record.get("cap"))
        print(f"{name}   k={record['k']}   task {tid}")
        collect(record, AwsQuantumTask(record["task_arn"]), meta, feas)
        print(f"  written to {os.path.relpath(path)}")
        return
    if a.shots is None or a.only is None:
        raise SystemExit("--shots and --only are required to submit")
    if a.debias and a.shots < DEBIAS_MIN_SHOTS:
        raise SystemExit(f"debiasing on Braket needs >= {DEBIAS_MIN_SHOTS} shots")
    if not a.dry_run and not a.confirm_qpu:
        cost = PRICE["task"] + PRICE["shot"] * a.shots
        raise SystemExit(f"QPU task costs ${cost:.2f}; add --confirm-qpu")

    demands = [int(d) for d in a.demands.split(",")] if a.demands else None
    for name, qc, meta, feas in select(a.only, a.k, demands, a.cap):
        tqc = transpile(qc, basis_gates=BASES[a.basis], optimization_level=3,
                        seed_transpiler=7)
        cmap = measured_map(tqc)
        circ = to_braket(tqc)
        where = "local simulator" if a.dry_run else f"ionq/{a.device}"
        print(f"{name}   k={a.k}   |F| = {len(feas)}   {qc.num_qubits} qubits   "
              f"{n2q(tqc)} two-qubit gates, depth {tqc.depth()}   on {where}   "
              f"{a.shots} shots", flush=True)

        if a.dry_run:
            from braket.devices import LocalSimulator
            res = LocalSimulator().run(circ, shots=a.shots).result()
            counts = marginal(dict(res.measurement_counts), cmap)
            report("braket-local", law_from_counts(counts, meta, feas))
            continue

        from braket.aws import AwsDevice
        from boto3 import Session
        from braket.aws import AwsSession
        device = AwsDevice(ARN.format(a.device),
                           aws_session=AwsSession(Session(region_name=REGION)))
        kw = {}
        if a.debias:
            from braket.error_mitigation import Debias
            kw["device_parameters"] = {"errorMitigation": Debias()}
        task = device.run(circ, shots=a.shots, **kw)
        tid = task.id.rsplit("/", 1)[-1]
        record = dict(
            case=name,
            case_key=f"{case_key(name, a.k, demands)}_{a.basis}_braket",
            demands=demands, cap=a.cap, k=a.k, shots=a.shots,
            backend=f"braket:ionq/{a.device}", simulator=False,
            job_id=tid, task_arn=task.id,
            submitted_utc=datetime.now(timezone.utc).isoformat(),
            logical_qubits=qc.num_qubits, two_qubit=n2q(tqc),
            depth=tqc.depth(), basis=a.basis, optimization_level=3,
            error_mitigation="debiasing" if a.debias else None,
            measure_map=cmap, lam=meta["lam"], n_feasible=len(feas),
            price_usd=PRICE["task"] + PRICE["shot"] * a.shots)
        path = save(record)                 # the task ARN is on disk before we wait
        print(f"  submitted task {tid}, saved {os.path.relpath(path)}", flush=True)
        if a.no_wait:
            print(f"  collect with: python braket_oracle.py --task {tid}")
            continue
        collect(record, task, meta, feas)


if __name__ == "__main__":
    main()
