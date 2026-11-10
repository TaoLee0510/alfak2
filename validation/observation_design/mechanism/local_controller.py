#!/usr/bin/env python3
import argparse
import importlib.util
import os
from pathlib import Path


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--index",type=int,choices=[1,2],required=True)
    args=parser.parse_args();repo=Path(os.environ["ALFAK2_REPO"]);root=repo/"validation/observation_design"
    out=Path(os.environ["ALFAK2_MECHANISM_OUTPUT"]);source=Path(os.environ["ALFAK2_MECHANISM_LOCAL_SOURCE"])
    dest=out/"localization";dest.mkdir(exist_ok=True)
    marker=dest/f"controller_{args.index}_resource.csv"
    if marker.exists():raise RuntimeError("Existing localization attempt must be retained")
    loader=importlib.util.spec_from_file_location("guard",root/"results/control_source_v1/bounded_run.py")
    guard=importlib.util.module_from_spec(loader);loader.loader.exec_module(guard)
    env=os.environ.copy();env["ALFAK2_MECHANISM_INDEX"]=str(args.index)
    print("Localization owner",os.getpid(),"index",args.index,flush=True)
    result=guard.execute(["Rscript",str(source/"localize.R")],env,str(dest/f"localize_{args.index}.log"),float("inf"))
    result["case_index"]=args.index;result["memory_pass"]=result["maximum_rss_bytes"]<=6*1024**3
    guard.write_csv(marker,[result]);print(result,flush=True)


if __name__=="__main__":main()
