#!/usr/bin/env python3
"""Preflight gate, then exactly two fitting controllers; successful preflights reused."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    src = Path(os.environ["ALFAK2_CONSTRAINED_SOURCE"])
    out = Path(os.environ["ALFAK2_CONSTRAINED_OUTPUT"])
    if not (out / "noise/COMPLETE").exists():
        raise RuntimeError("Finish process-noise computation before budgeted fit timing")
    for phase in ["preflight", "main"]:
        children = []
        streams = []
        for part in range(2):
            stream = (out / f"{phase}_controller_{part}.log").open("w")
            streams.append(stream)
            child = subprocess.Popen([sys.executable, str(src / "controller.py"), "--phase", phase, "--part", str(part)],
                                     stdout=stream, stderr=subprocess.STDOUT)
            children.append(child)
        print(phase, "owners", [p.pid for p in children], flush=True)
        while any(p.poll() is None for p in children):
            time.sleep(1)
        for stream in streams:
            stream.close()
        if any(p.returncode != 0 for p in children):
            raise RuntimeError("Controller failed; preserve outputs and inspect logs")
        statuses = []
        for part in range(2):
            p = out / "workers" / f"{phase}_part_{part}" / "status.csv"
            if p.exists():
                with p.open() as stream:
                    statuses.extend(csv.DictReader(stream))
        if phase == "preflight" and not all(r["validation_pass"].lower() == "true" and
                                           r["resource_pass"].lower() == "true" for r in statuses):
            (out / "PREFLIGHT_FAILED").write_text(json.dumps(statuses, indent=2))
            raise RuntimeError("Preflight gate failed; full batch not started")
        print(phase, "complete", len(statuses), flush=True)
    (out / "BATCH_ATTEMPT_COMPLETE").write_text("All planned main cases attempted; inspect statuses.\n")


if __name__ == "__main__":
    main()
