#!/usr/bin/env python3
"""Exactly one R fitting child per controller, with external resource accounting."""
import argparse
import importlib.util
import json
import os
from pathlib import Path


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--mode",choices=["reference","primary"],required=True)
    parser.add_argument("--part",type=int,choices=[0,1],required=True)
    parser.add_argument("--cases",nargs="*")
    args=parser.parse_args()
    repo=Path(os.environ["ALFAK2_REPO"]);root=repo/"validation/observation_design"
    src=Path(os.environ["ALFAK2_ADAPTIVE_SOURCE"]);out=Path(os.environ["ALFAK2_ADAPTIVE_OUTPUT"])
    out.mkdir(exist_ok=True)
    spec=importlib.util.spec_from_file_location("guard",root/"results/control_source_v1/bounded_run.py")
    guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
    if args.cases:
        cases=args.cases
    elif args.mode=="reference":
        cases=[f"n{n}_{profile}_beta{beta}_counts" for n in [4,8] for profile in ["balanced","sparse"] for beta in ["0","0p005"]]
    else:
        old=guard.read_csv(root/"results/analysis_input/all_fit_status.csv")
        cases=[row["case_id"] for row in old if row["status"]!="passed_quadrature_refinement"]
    cases=[case for i,case in enumerate(cases) if i%2==args.part]
    destination=out/f"{args.mode}_part_{args.part}";destination.mkdir(exist_ok=True)
    if (destination/"COMPLETE").exists():raise RuntimeError("Completed controller cannot be overwritten")
    statuses=[]
    print("Adaptive controller",os.getpid(),args.mode,"cases",cases,flush=True)
    for case in cases:
        if (out/args.mode/case/"status.csv").exists():raise RuntimeError("Existing case attempt must be retained")
        env=os.environ.copy();env.update(ALFAK2_ADAPTIVE_CASE=case,ALFAK2_ADAPTIVE_MODE=args.mode)
        (destination/"active_case.json").write_text(json.dumps({"case":case,"owner_pid":os.getpid()}))
        result=guard.execute(["Rscript",str(src/"run_case.R")],env,str(destination/(case+".log")),300)
        p=out/args.mode/case/"status.csv"
        if result["exit_code"]==0 and p.exists():row=guard.read_csv(p)[0]
        else:row=dict(case_id=case,mode=args.mode,status="external_time_limit" if result["external_timeout"] else "process_error",
            validation_pass=False,elapsed_seconds=result["child_seconds"],nodes="",message=f"exit {result['exit_code']}")
        row.update(child_wall_seconds=result["child_seconds"],maximum_rss_bytes=result["maximum_rss_bytes"],
            resource_pass=not result["external_timeout"] and result["maximum_rss_bytes"]<=6*1024**3 and result["child_seconds"]<=300.3)
        statuses.append(row);guard.write_csv(destination/"status.csv",statuses)
        print(case,row["status"],"validation",row["validation_pass"],"wall",round(result["child_seconds"],2),
            "peakGiB",round(result["maximum_rss_bytes"]/1024**3,3),flush=True)
    (destination/"active_case.json").unlink(missing_ok=True)
    (destination/"COMPLETE").write_text("All owned cases attempted; no failures removed.\n")


if __name__=="__main__":main()
