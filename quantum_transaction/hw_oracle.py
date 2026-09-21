"""
The zone sampler on a real device (Sec. V-B).
=============================================
`noise_oracle.py` asks what a noisy device costs by simulating one.  This
module runs the same two circuits on an IBM Quantum processor and scores the
returned counts with the same acceptance test (`noise_oracle.law_from_counts`),
so a hardware row sits beside the simulated rows without any change of metric.

Unlike the simulation, the hardware run is routed: the circuit is transpiled
onto the device's coupling map and native gates, and the routed two-qubit
count is recorded with the result.  Every submitted job is written to
`hw_runs/` with its job ID, backend, calibration timestamp, layout, options
and raw counts, so the numbers in the paper can be re-scored from the file
and the job can be looked up on the IBM Quantum dashboard.

Run:
    python hw_oracle.py --dry-run                 # transpile only, no QPU time
    python hw_oracle.py [--backend ibm_kingston] [--shots 10000] [--only inter]
    python hw_oracle.py --no-wait --backend ibm_boston   # submit, collect later
    python hw_oracle.py --only inter --demands 1,2,1,2 --cap 4 --k 1
    python hw_oracle.py --job <job_id>             # collect or re-score a job

Needs a saved account: QiskitRuntimeService.save_account(...).
"""

import argparse
import json
import os
from datetime import datetime, timezone

from qiskit.transpiler import generate_preset_pass_manager

import glob

from noise_oracle import cases, inter_case, law_from_counts

BACKEND = "ibm_kingston"
SHOTS = 10_000
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hw_runs")


def n2q(c):
    return sum(v for g, v in c.count_ops().items()
               if g in ("cx", "ecr", "cz", "rzz"))


def routed(qc, backend, seed=7):
    """The circuit as the device will execute it."""
    pm = generate_preset_pass_manager(optimization_level=3, backend=backend,
                                      seed_transpiler=seed)
    return pm.run(qc)


def calibration_time(backend):
    try:
        return backend.properties().last_update_date.isoformat()
    except Exception:
        return None


def job_metadata(job):
    """What IBM recorded about the run: device, queue and execution times.

    The instance CRN carries the IBM Cloud account ID, so only its region is
    kept -- this file is published with the repository.
    """
    m = job.metrics()
    backend = job.backend()
    crn = getattr(job, "instance", None) or ""
    return dict(
        backend=backend.name,
        simulator=bool(backend.configuration().simulator),
        status=str(job.status()),
        region=crn.split(":")[5] if crn.count(":") >= 5 else None,
        timestamps_utc=m.get("timestamps"),
        usage=m.get("usage"),
        versions=m.get("qiskit_version"))


def finish(record, job, meta, feas):
    """Score a finished job and attach its counts and metadata to `record`."""
    counts = job.result()[0].data.c.get_counts()
    record["metrics"] = r = law_from_counts(counts, meta, feas)
    record["job"] = job_metadata(job)
    record["counts"] = counts
    report("hardware", r)


def select(only, k, demands=None, cap=None):
    """The cases to run; `demands`/`cap` replace the inter-cell defaults."""
    todo = [c for c in cases(k) if not only or c[0].startswith(only)]
    if demands is not None:
        if only != "inter":
            raise SystemExit("--demands/--cap need --only inter")
        qc, meta, feas = inter_case(demands, cap, k)
        name = f"inter-cell, 4 UE x 2 AP, demands {tuple(demands)}, cap {cap}"
        todo = [(name, qc, meta, feas)]
    return todo


def case_key(name, k, demands):
    key = name.split(",")[0]
    if demands is not None:
        key += "_w" + "".join(str(d) for d in demands)
    return f"{key}_k{k}"


def save(record):
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(
        OUT_DIR, f"{record['backend']}_{record['case_key']}_{record['job_id']}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=1)
    return path


def report(name, r):
    print(f"  {name:<28} accept {r['accept']:.3f}   discard {r['discard']:.3f}"
          f"   TVD {r['tvd']:.3f}   infeas {r['infeas']:.4f}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default=BACKEND)
    ap.add_argument("--shots", type=int, default=SHOTS)
    ap.add_argument("--k", type=int, default=1)
    ap.add_argument("--only", choices=("inter", "intra"))
    ap.add_argument("--demands", help="inter-cell demands, e.g. 1,2,1,2")
    ap.add_argument("--cap", type=int, help="inter-cell AP capacity")
    ap.add_argument("--dry-run", action="store_true",
                    help="transpile and count gates, submit nothing")
    ap.add_argument("--job", help="collect or re-score a submitted job")
    ap.add_argument("--no-wait", action="store_true",
                    help="submit, save the job ID, and return without waiting")
    ap.add_argument("--no-dd", action="store_true",
                    help="disable dynamical decoupling, isolating idle decoherence")
    ap.add_argument("--seed", type=int, default=7,
                    help="transpiler seed; different seeds give different layouts")
    a = ap.parse_args()

    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

    service = QiskitRuntimeService()
    demands = [int(d) for d in a.demands.split(",")] if a.demands else None
    if (demands is None) != (a.cap is None):
        raise SystemExit("--demands and --cap go together")

    if a.job:
        found = glob.glob(os.path.join(OUT_DIR, f"*_{a.job}.json"))
        record = {}
        if found:                           # complete a record saved at submission
            with open(found[0], encoding="utf-8") as f:
                record = json.load(f)
            only = record["case"].split("-")[0]
            todo = select(only, record["k"], record.get("demands"),
                          record.get("cap"))
        else:
            todo = select(a.only, a.k, demands, a.cap)
        if len(todo) != 1:
            raise SystemExit("--job needs --only inter|intra (no saved record)")
        name, qc, meta, feas = todo[0]
        job = service.job(a.job)
        print(f"{name}   k={record.get('k', a.k)}   job {a.job} on "
              f"{job.backend().name}")
        finish(record, job, meta, feas)
        if found:
            save(record)
            print(f"  written to {os.path.relpath(found[0])}")
        return

    todo = select(a.only, a.k, demands, a.cap)
    backend = service.backend(a.backend)
    print("=" * 78)
    print(f"  Zone sampler on {backend.name} ({backend.num_qubits} qubits), "
          f"{a.shots:,} shots, k={a.k}")
    print(f"  calibration {calibration_time(backend)}")
    print("=" * 78, flush=True)

    for name, qc, meta, feas in todo:
        key = (case_key(name, a.k, demands) + ("_nodd" if a.no_dd else "")
               + (f"_s{a.seed}" if a.seed != 7 else ""))
        isa = routed(qc, backend, seed=a.seed)
        layout = isa.layout.final_index_layout() if isa.layout else None
        print()
        print(f"{name}   |F| = {len(feas)}   {qc.num_qubits} qubits   routed "
              f"{n2q(isa)} two-qubit gates, depth {isa.depth()}", flush=True)
        if a.dry_run:
            continue

        sampler = SamplerV2(mode=backend)
        sampler.options.dynamical_decoupling.enable = not a.no_dd
        record = dict(
            case=name, case_key=key, demands=demands, cap=a.cap,
            backend=backend.name, calibration=calibration_time(backend),
            submitted_utc=datetime.now(timezone.utc).isoformat(),
            shots=a.shots, k=a.k,
            logical_qubits=qc.num_qubits, routed_2q=n2q(isa),
            routed_depth=isa.depth(), layout=layout,
            options=dict(optimization_level=3,
                         dynamical_decoupling=not a.no_dd, seed_transpiler=a.seed),
            lam=meta["lam"], n_feasible=len(feas))
        job = sampler.run([isa], shots=a.shots)
        record["job_id"] = job.job_id()
        path = save(record)                 # the job ID is on disk before we wait
        print(f"  submitted job {job.job_id()}, saved {os.path.relpath(path)}",
              flush=True)
        if a.no_wait:
            print(f"  collect with: python hw_oracle.py --job {job.job_id()}")
            continue

        finish(record, job, meta, feas)
        save(record)
        print(f"  results written to {os.path.relpath(path)}", flush=True)

if __name__ == "__main__":
    main()
