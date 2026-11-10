#!/usr/bin/env python3
import argparse
import importlib.util
import os
from pathlib import Path

parser=argparse.ArgumentParser();parser.add_argument("--index",type=int,choices=[1,2],required=True)
args=parser.parse_args();repo=Path(os.environ["ALFAK2_REPO"]);root=repo/"validation/observation_design"
src=Path(os.environ["ALFAK2_IS_SOURCE"]);out=Path(os.environ["ALFAK2_IS_OUTPUT"])
spec=importlib.util.spec_from_file_location("guard",root/"results/control_source_v1/bounded_run.py")
guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
marker=out/f"reference_{args.index}_resource.csv"
if marker.exists():raise RuntimeError("Existing reference attempt cannot be overwritten")
env=os.environ.copy();env["ALFAK2_IS_INDEX"]=str(args.index)
print("High-precision reference owner",os.getpid(),args.index,flush=True)
r=guard.execute(["Rscript",str(src/"extend_reference.R")],env,str(out/f"reference_{args.index}.log"),float("inf"))
r['memory_pass']=r['maximum_rss_bytes']<=6*1024**3
guard.write_csv(marker,[r]);print(r,flush=True)
