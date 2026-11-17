#!/usr/bin/env python3
"""Fresh equivalent-target fits; one hard-deadline R child per controller."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--part", required=True, type=int, choices=[0, 1])
    args = parser.parse_args()
    repo = Path(os.environ["ALFAK2_REPO"]).resolve()
    root = repo / "validation/observation_design"
    source = Path(os.environ["ALFAK2_DIRECT_BATCH_SOURCE"])
    out = Path(os.environ["ALFAK2_DIRECT_BATCH_OUTPUT"])
    out.mkdir(exist_ok=True)
    spec = importlib.util.spec_from_file_location("control", root / "results/control_source_v1/bounded_run.py")
    control = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(control)
    validation = control.read_csv(root / "results/direct_validation_tests.csv")
    if not all(float(row["failed"]) == 0 and row["error"].lower() == "false" and
               float(row["warning"]) == 0 for row in validation):
        raise RuntimeError("Direct implementation equivalence validation not passed")
    smoke = json.loads((root / "results/direct_smoke_resources.json").read_text())
    if not smoke["passed_6GiB_300sec"] or smoke.get("numerical_status")!="passed_quadrature_refinement":
        raise RuntimeError("Direct smoke resource gate not passed")
    if (source / "direct_joint.R").read_bytes() != Path(smoke["implementation_source"]).read_bytes():
        raise RuntimeError("Batch implementation differs from validated implementation")
    manifest = control.read_csv(root / "results/dataset_manifest.csv")
    copied = out / "dataset_manifest.csv"
    if not copied.exists():
        shutil.copy2(root / "results/dataset_manifest.csv", copied)
    units = list(dict.fromkeys(row["unit_id"] for row in manifest))
    selected = {unit for index, unit in enumerate(units, 1) if index % 2 == args.part}
    order = ["control300", "depth600", "time600", "composition_anchor600", "composition600"]
    cases = sorted((row for row in manifest if row["unit_id"] in selected),
                   key=lambda row: (units.index(row["unit_id"]), order.index(row["design"])))
    destination = out / f"fits_{args.part}"
    destination.mkdir(exist_ok=True)
    if (destination / "status.csv").exists():
        raise RuntimeError("Direct part already owns results; do not overwrite attempts")
    statuses, metrics = [], []
    print(f"Direct controller PID {os.getpid()}, part {args.part}, {len(cases)} cases", flush=True)
    for row in cases:
        case_id = row["case_id"]
        environment = os.environ.copy()
        environment["ALFAK2_DIRECT_CASE"] = case_id
        environment.pop("ALFAK2_DESIGN_OUTPUT", None)
        (destination / "active_case.json").write_text(json.dumps({"case_id": case_id, "started_unix": time.time(),
            "parent_pid": os.getpid(), "budget_seconds": 300}, indent=2) + "\n")
        outcome = control.execute(["Rscript", str(source / "run_direct.R")], environment,
                                  str(destination / (case_id + ".log")), 300)
        case_out = out / "cases" / case_id
        if outcome["exit_code"] == 0 and (case_out / "status.csv").exists():
            status = control.read_csv(case_out / "status.csv")[0]
            metrics.extend(control.read_csv(case_out / "case_metrics.csv"))
        else:
            status = {key: row[key] for key in ["case_id", "unit_id", "landscape_id", "family", "design"]}
            label = "external_time_limit" if outcome["external_timeout"] else "child_process_error"
            status.update(status=label, raw_status=label, elapsed_seconds=outcome["child_seconds"],
                final_order="", conditional_nodes="", message=f"Child exit {outcome['exit_code']}; partial cache retained",
                result_path="", recovery="none")
        status.update(child_wall_seconds=outcome["child_seconds"], maximum_rss_bytes=outcome["maximum_rss_bytes"],
            resource_pass=not outcome["external_timeout"] and outcome["child_seconds"]<=300.3 and
            outcome["maximum_rss_bytes"]<=6*1024**3, computational_method="direct_joint_mode_hessian")
        statuses.append(status)
        control.write_csv(destination / "status.csv", statuses)
        control.write_csv(destination / "case_metrics.csv", metrics)
        if (case_out / "source_sha256.csv").exists() and not (destination / "source_sha256.csv").exists():
            shutil.copy2(case_out / "source_sha256.csv", destination / "source_sha256.csv")
        print(f"{len(statuses)}/30 {case_id} {status['status']} "
              f"fit={float(status['elapsed_seconds']):.2f}s peak={outcome['maximum_rss_bytes']/1024**3:.3f}GiB "
              f"resource_pass={status['resource_pass']}", flush=True)
    (destination / "active_case.json").unlink(missing_ok=True)
    (destination / "COMPLETE").write_text("All planned cases attempted; no failures removed.\n")


if __name__ == "__main__":
    main()
