#!/usr/bin/env python3
"""One sequential fitting child per controller; two parts only."""
import argparse
import importlib.util
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--part", type=int, choices=[0, 1], required=True)
    parser.add_argument("--phase", choices=["preflight", "main"], required=True)
    args = parser.parse_args()
    repo = Path(os.environ["ALFAK2_REPO"])
    root = repo / "validation/observation_design"
    source = Path(os.environ["ALFAK2_CONSTRAINED_SOURCE"])
    out = Path(os.environ["ALFAK2_CONSTRAINED_OUTPUT"])
    loader = importlib.util.spec_from_file_location("guard", root / "results/control_source_v1/bounded_run.py")
    guard = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(guard)
    cases = json.loads((out / "queue.json").read_text())[args.phase][args.part::2]
    destination = out / "workers" / f"{args.phase}_part_{args.part}"
    destination.mkdir(parents=True, exist_ok=False)
    statuses = []
    print("Owner", os.getpid(), args.phase, args.part, len(cases), "cases", flush=True)
    for case in cases:
        case_dir = out / "fits" / case
        if case_dir.exists():
            raise RuntimeError("Existing attempt must be preserved: " + case)
        env = os.environ.copy()
        env["ALFAK2_CONSTRAINED_CASE"] = case
        (destination / "active_case.json").write_text(json.dumps({"case": case, "owner": os.getpid()}))
        resource = guard.execute(["Rscript", str(source / "run_case.R")], env,
                                 str(destination / (case + ".log")), 300)
        path = case_dir / "status.csv"
        if resource["exit_code"] == 0 and path.exists():
            row = guard.read_csv(path)[0]
        else:
            row = dict(case_id=case, status="external_time_limit" if resource["external_timeout"] else "process_error",
                       validation_pass=False, elapsed_seconds=resource["child_seconds"], nodes="",
                       message=f"exit {resource['exit_code']}; no case or seed removed")
        row.update(child_wall_seconds=resource["child_seconds"], maximum_rss_bytes=resource["maximum_rss_bytes"],
                   resource_pass=not resource["external_timeout"] and resource["maximum_rss_bytes"] <= 6 * 1024**3
                   and resource["child_seconds"] <= 300.3)
        statuses.append(row)
        guard.write_csv(destination / "status.csv", statuses)
        print(case, row["status"], "wall", round(resource["child_seconds"], 2),
              "peakGiB", round(resource["maximum_rss_bytes"] / 1024**3, 3), flush=True)
    (destination / "active_case.json").unlink(missing_ok=True)
    (destination / "COMPLETE").write_text("All owned cases attempted; failures retained.\n")


if __name__ == "__main__":
    main()
