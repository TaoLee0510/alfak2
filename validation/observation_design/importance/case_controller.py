#!/usr/bin/env python3
import argparse, importlib.util, os, json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--part',type=int,choices=[0,1],required=True);a=p.parse_args()
repo=Path(os.environ['ALFAK2_REPO']);root=repo/'validation/observation_design';src=Path(os.environ['ALFAK2_IS_SOURCE'])
spec=importlib.util.spec_from_file_location('guard',root/'results/control_source_v1/bounded_run.py')
guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
rows=guard.read_csv(Path(a.manifest));selected=[r for i,r in enumerate(rows) if i%2==a.part]
for r in selected:
 out=Path(r['output']);out.mkdir(parents=True,exist_ok=False)
 env=os.environ.copy();env.update(ALFAK2_IS_CASE=r['file'],ALFAK2_IS_CONTEXT=r['context'],ALFAK2_IS_CASE_OUTPUT=str(out),ALFAK2_IS_SEED=r['seed'])
 (out/'active.json').write_text(json.dumps({'owner_pid':os.getpid(),'case_id':r['case_id']}))
 z=guard.execute(['Rscript',str(src/'run_case.R')],env,str(out/'run.log'),300)
 z.update(case_id=r['case_id'],resource_pass=not z['external_timeout'] and z['child_seconds']<=300.3 and z['maximum_rss_bytes']<=6*1024**3)
 guard.write_csv(out/'resource.csv',[z]);(out/'active.json').unlink()
 print(z,flush=True)
