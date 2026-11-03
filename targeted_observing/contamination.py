"""Bright-end supplemented risk report at actual planned visits, never a clean-field certificate."""
import argparse
import json
from pathlib import Path
from .catalog import BASE
from .bright_stars import assess

def main():
    p=argparse.ArgumentParser();p.add_argument('--campaign',default=str(BASE/'data/campaign_20261103.json'))
    p.add_argument('--plans',nargs='*',default=[]);p.add_argument('--output',default=str(BASE/'products/pro_review/bright_stars.json'));a=p.parse_args()
    campaign=json.loads(Path(a.campaign).read_text());visits=dict(campaign.get('planned_visits',{}))
    for filename in a.plans:
        plan=json.loads(Path(filename).read_text());visits[plan['date']]=[r for b in plan['blocks'] if b['mode']=='trojan' for r in b['rows']]
    if not visits:raise ValueError('Supply --plans or a campaign containing planned_visits; do not substitute daytime hourly extrema')
    result=assess(campaign,visits);out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'risk_rows':len(result['risks']),'hip_rows':result['hip_rows'],'output':str(out)}))
if __name__=='__main__':main()
