#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('manifest');a=p.parse_args()
rows=list(csv.DictReader(open(a.manifest)));done=[];bad=[];active=[]
for r in rows:
 out=Path(r['output']);resource=out/'resource.csv';status=out/'status.csv'
 if resource.exists():
  rr=next(csv.DictReader(resource.open()));ss=next(csv.DictReader(status.open())) if status.exists() else {}
  valid=rr['resource_pass'].lower()=='true' and ss.get('status')=='passed_joint_importance'
  done.append((r['case_id'],rr,ss));
  if not valid:bad.append(r['case_id'])
 elif (out/'active.json').exists():active.append(r['case_id'])
print(json.dumps(dict(total=len(rows),completed=len(done),passed=len(done)-len(bad),failed=bad,active=active,
 max_seconds=max([float(r['child_seconds']) for _,r,_ in done],default=0),
 max_rss_gib=max([float(r['maximum_rss_bytes'])/1024**3 for _,r,_ in done],default=0)),ensure_ascii=False))
