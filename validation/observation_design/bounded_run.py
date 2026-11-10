#!/usr/bin/env python3
"""Resume frozen v4 cases under an external, cumulative wall-clock deadline.

Only two instances (parts 0 and 1) may run. Each owns at most one R child.
Old results and unfinished node caches are preserved without rewriting them.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import signal
import sys
import time


def execute(command, environment, logfile, seconds):
    started = time.monotonic()
    child = os.fork()
    if child == 0:
        try:
            fd = os.open(logfile, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(fd, 1)
            os.dup2(fd, 2)
            os.close(fd)
            os.execvpe(command[0], command, environment)
        except BaseException:
            os._exit(127)
    limited = False
    while True:
        pid, status, usage = os.wait4(child, os.WNOHANG)
        if pid:
            break
        if time.monotonic() - started >= seconds:
            os.kill(child, signal.SIGKILL)
            _, status, usage = os.wait4(child, 0)
            limited = True
            break
        time.sleep(0.1)
    maximum_rss = usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    return {"exit_code": os.waitstatus_to_exitcode(status),
            "external_timeout": limited, "child_seconds": time.monotonic() - started,
            "maximum_rss_bytes": maximum_rss, "child_pid": child}


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--part", type=int, choices=[0, 1])
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "child.log")
            done = execute([sys.executable, "-c", "print('ok')"], os.environ.copy(), path, 3)
            cut = execute([sys.executable, "-c", "import time;time.sleep(3)"], os.environ.copy(), path, .2)
            assert done["exit_code"] == 0 and not done["external_timeout"]
            assert cut["external_timeout"] and cut["exit_code"] == -signal.SIGKILL and cut["child_seconds"] < 1
            print(json.dumps({"normal_exit": done, "hard_timeout": cut}), flush=True)
        return
    if args.part is None:
        parser.error("--part is required")
    repo = Path(os.environ.get("ALFAK2_REPO", Path(__file__).resolve().parents[2])).resolve()
    os.chdir(repo)
    root = repo / "validation/observation_design"
    out = root / "results"
    source = out / "computation_source_v4"
    manifest = read_csv(out / "dataset_manifest.csv")
    units = list(dict.fromkeys(row["unit_id"] for row in manifest))
    selected = {unit for index, unit in enumerate(units, 1) if index % 2 == args.part}
    design_order = ["control300", "depth600", "time600", "composition_anchor600", "composition600"]
    cases = sorted((row for row in manifest if row["unit_id"] in selected),
                   key=lambda row: (units.index(row["unit_id"]), design_order.index(row["design"])))
    destination = out / f"fits_bounded_{args.part}"
    destination.mkdir(exist_ok=True)
    if (destination / "status.csv").exists():
        raise RuntimeError("Bounded part already exists; inspect ownership before resuming it")
    old = out / f"fits_{args.part}"
    statuses = read_csv(old / "status.csv")
    metrics = read_csv(old / "case_metrics.csv")
    for row in statuses:
        row["original_elapsed_seconds"] = row["elapsed_seconds"]
        row["maximum_rss_bytes"] = ""
        row["resource_pass"] = float(row["elapsed_seconds"]) <= 300
        row["result_path"] = str(old / (row["case_id"] + "_result.rds"))
        row["execution"] = "original_v4"
    completed = {row["case_id"] for row in statuses}
    transition = json.loads((out / "external_deadline_transition.json").read_text())
    partial = next(row for row in transition if row["part"] == args.part)
    remaining = [row for row in cases if row["case_id"] not in completed]
    write_csv(destination / "status.csv", statuses)
    write_csv(destination / "case_metrics.csv", metrics)
    print(f"Bounded owner PID {os.getpid()}, part {args.part}, {len(remaining)} remaining cases", flush=True)
    for index, row in enumerate(remaining):
        case_id = row["case_id"]
        prior = partial["seconds_since_last_status"] if index == 0 else 0.0
        budget = max(0.0, 300.0 - prior)
        child_part = f"bounded_child_{args.part}_{case_id}"
        child_out = out / ("fits_" + child_part)
        environment = os.environ.copy()
        environment.update(ALFAK2_DESIGN_SOURCE=str(source), ALFAK2_DESIGN_PART=child_part,
                           ALFAK2_DESIGN_ONLY=case_id, ALFAK2_REPO=str(repo))
        active = dict(case_id=case_id, prior_seconds=prior, remaining_budget_seconds=budget,
                      parent_pid=os.getpid(), started_unix=time.time())
        (destination / "active_case.json").write_text(json.dumps(active, indent=2) + "\n")
        result = execute(["Rscript", str(source / "run.R")], environment,
                         str(destination / (case_id + ".log")), budget)
        status_path = child_out / "status.csv"
        if result["exit_code"] == 0 and status_path.exists():
            status = read_csv(status_path)[0]
            status["elapsed_seconds"] = float(status["elapsed_seconds"]) + prior
            metrics.extend(read_csv(child_out / "case_metrics.csv"))
            status["result_path"] = str(child_out / (case_id + "_result.rds"))
        else:
            status = {key: row[key] for key in ["case_id", "unit_id", "landscape_id", "family", "design"]}
            status.update(status="external_time_limit" if result["external_timeout"] else "child_process_error",
                          elapsed_seconds=prior + result["child_seconds"], final_order="",
                          conditional_nodes="", message="Hard child deadline; partial caches retained" if result["external_timeout"] else
                          f"Child exit code {result['exit_code']}; inspect case log", result_path="")
        status.update(original_elapsed_seconds=prior, child_wall_seconds=result["child_seconds"],
                      maximum_rss_bytes=result["maximum_rss_bytes"], execution="external_deadline_v1",
                      resource_pass=not result["external_timeout"] and
                      prior + result["child_seconds"] <= 300.3 and result["maximum_rss_bytes"] <= 6 * 1024**3)
        statuses.append(status)
        write_csv(destination / "status.csv", statuses)
        write_csv(destination / "case_metrics.csv", metrics)
        print(f"{len(statuses)}/30 {case_id} {status['status']} "
              f"total={float(status['elapsed_seconds']):.2f}s peak={result['maximum_rss_bytes']/1024**3:.3f}GiB "
              f"resource_pass={status['resource_pass']}", flush=True)
    (destination / "active_case.json").unlink(missing_ok=True)
    (destination / "COMPLETE").write_text("All planned cases attempted; failures retained.\n")


if __name__ == "__main__":
    main()
