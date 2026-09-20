"""
The zone sampler on IonQ (Sec. V-B).
====================================
The IonQ counterpart of `hw_oracle.py`: same circuits, same acceptance test
(`noise_oracle.law_from_counts`), same record format in `hw_runs/`.

IonQ traps are all-to-all connected, so there is no routing: the two-qubit
count sent to the device is the logical count (384 for the inter-cell
circuit at k=1, against 767 after routing onto a heavy-hex lattice).

IonQ bills by shot, so shot counts here are sized from the precision needed,
not set generously.  A QPU submission additionally requires `--confirm-qpu`;
the cloud simulator, with or without a device noise model, does not.

Run:
    python ionq_oracle.py --noise-model forte-1 --shots 500 \
        --only inter --demands 1,2,1,2 --cap 4 --k 1
    python ionq_oracle.py --backend qpu.forte-enterprise-1 --confirm-qpu         --no-wait ...                       # debiased; collect with --job <id>

API key: env IONQ_API_KEY, else ~/.ionq/api_key.
"""

import argparse
import glob
import json
import os
import pathlib
from datetime import datetime, timezone

from qiskit import transpile

from hw_oracle import OUT_DIR, case_key, n2q, report, save, select
from noise_oracle import law_from_counts

# Controlled phases of the QFT adder become one partial-angle ZZ each on IonQ,
# not two CNOTs; with optimization level 3 this halves the two-qubit count.
BASES = {"cx": ["cx", "rz", "sx", "x"],
         "rzz": ["rzz", "rz", "ry", "rx", "h", "x", "sx"]}


def api_key():
    key = os.getenv("IONQ_API_KEY")
    if key:
        return key.strip()
    return (pathlib.Path.home() / ".ionq" / "api_key").read_text().strip()


def job_details(job_id):
    """What IonQ recorded: status, compiled gate counts, timing, cost."""
    import requests
    d = requests.get(f"https://api.ionq.co/v0.4/jobs/{job_id}",
                     headers={"Authorization": f"apiKey {api_key()}"},
                     timeout=30).json()
    keep = ("status", "backend", "noise", "stats", "settings", "cost",
            "request", "response", "start", "end", "execution_duration_ms",
            "predicted_execution_duration_ms", "error_mitigation")
    return {k: d[k] for k in keep if k in d}


def measured_map(tqc):
    """clbit index -> qubit index, read off the circuit that was sent."""
    m = {}
    for inst in tqc.data:
        if inst.operation.name == "measure":
            m[tqc.find_bit(inst.clbits[0]).index] = tqc.find_bit(inst.qubits[0]).index
    return [m[c] for c in range(len(m))]


def raw_counts(job_id, shots, cmap):
    """Counts straight from IonQ's histogram, marginalized onto the clbits.

    The API returns the distribution over the whole register, keyed by the
    integer whose bit q is qubit q.  qiskit-ionq's own conversion resamples
    simulator counts client-side, so it is bypassed.  `sharpen=false` asks
    for the averaged debiased distribution rather than plurality voting.
    """
    import requests
    probs = requests.get(
        f"https://api.ionq.co/v0.4/jobs/{job_id}/results/probabilities",
        params={"sharpen": "false"},
        headers={"Authorization": f"apiKey {api_key()}"}, timeout=60).json()
    marg = {}
    for key, p in probs.items():
        x = int(key)
        bits = "".join(str((x >> q) & 1) for q in reversed(cmap))
        marg[bits] = marg.get(bits, 0.0) + p
    counts = {b: round(p * shots) for b, p in marg.items() if round(p * shots)}
    return counts, marg


def collect(record, job_id, meta, feas):
    counts, probs = raw_counts(job_id, record["shots"], record["measure_map"])
    record["metrics"] = r = law_from_counts(counts, meta, feas)
    record["sharpen"] = False
    record["counts"] = counts
    record["probabilities"] = probs
    record["job"] = job_details(job_id)
    report("ionq", r)
    save(record)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="simulator",
                    help="simulator | qpu.forte-enterprise-1 | qpu.forte-1")
    ap.add_argument("--noise-model", default="forte-1",
                    help="simulator only: ideal | forte-1 | aria-1 ...")
    ap.add_argument("--shots", type=int)
    ap.add_argument("--k", type=int, default=1)
    ap.add_argument("--only", choices=("inter", "intra"))
    ap.add_argument("--basis", choices=tuple(BASES), default="rzz")
    ap.add_argument("--demands")
    ap.add_argument("--cap", type=int)
    ap.add_argument("--no-wait", action="store_true")
    ap.add_argument("--job", help="collect a submitted job from its saved record")
    ap.add_argument("--confirm-qpu", action="store_true",
                    help="required to submit to a QPU (billed per shot)")
    a = ap.parse_args()

    if a.job:
        path = glob.glob(os.path.join(OUT_DIR, f"*_{a.job}.json"))[0]
        with open(path, encoding="utf-8") as f:
            record = json.load(f)
        (name, qc, meta, feas), = select(record["case"].split("-")[0],
                                         record["k"], record.get("demands"),
                                         record.get("cap"))
        if "measure_map" not in record:     # records saved before it was kept
            record["measure_map"] = measured_map(transpile(
                qc, basis_gates=BASES[record.get("basis", "cx")],
                optimization_level=record.get("optimization_level", 1),
                seed_transpiler=7))
        print(f"{name}   k={record['k']}   job {a.job}")
        collect(record, a.job, meta, feas)
        print(f"  written to {os.path.relpath(path)}")
        return
    if a.shots is None or a.only is None:
        raise SystemExit("--shots and --only are required to submit")

    qpu = a.backend.startswith("qpu.")
    if qpu and not a.confirm_qpu:
        raise SystemExit("QPU submission is billed per shot; add --confirm-qpu")

    from qiskit_ionq import IonQProvider
    provider = IonQProvider(api_key())
    backend = provider.get_backend(a.backend if qpu else "ionq_simulator")

    demands = [int(d) for d in a.demands.split(",")] if a.demands else None
    for name, qc, meta, feas in select(a.only, a.k, demands, a.cap):
        tqc = transpile(qc, basis_gates=BASES[a.basis], optimization_level=3,
                        seed_transpiler=7)
        where = a.backend if qpu else f"simulator ({a.noise_model} noise)"
        print(f"{name}   k={a.k}   |F| = {len(feas)}   {qc.num_qubits} qubits   "
              f"{n2q(tqc)} two-qubit gates, depth {tqc.depth()}   on {where}   "
              f"{a.shots} shots", flush=True)

        from qiskit_ionq import ErrorMitigation
        kw = dict(error_mitigation=ErrorMitigation.DEBIASING) if qpu else             dict(noise_model=a.noise_model)
        job = backend.run(tqc, shots=a.shots, **kw)
        record = dict(
            case=name, case_key=f"{case_key(name, a.k, demands)}_{a.basis}",
            demands=demands, cap=a.cap, k=a.k, shots=a.shots,
            backend=a.backend if qpu else f"ionq_simulator_{a.noise_model}",
            simulator=not qpu, job_id=job.job_id(),
            submitted_utc=datetime.now(timezone.utc).isoformat(),
            logical_qubits=qc.num_qubits, two_qubit=n2q(tqc),
            depth=tqc.depth(), basis=a.basis, optimization_level=3,
            error_mitigation="debiasing" if qpu else None,
            measure_map=measured_map(tqc),
            lam=meta["lam"], n_feasible=len(feas))
        path = save(record)                 # the job ID is on disk before we wait
        print(f"  submitted job {job.job_id()}, saved {os.path.relpath(path)}",
              flush=True)

        if a.no_wait:
            print(f"  collect with: python ionq_oracle.py --job {job.job_id()}")
            continue
        job.status()                        # block until done, then read raw
        job.result(sharpen=False)
        collect(record, job.job_id(), meta, feas)

if __name__ == "__main__":
    main()
