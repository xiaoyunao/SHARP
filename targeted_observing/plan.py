"""CLI for local plan generation; intentionally has no publish or execute option."""
import argparse
import json
from pathlib import Path
from .catalog import BASE
from .geometry import setup_iers
from .scheduler import make_plan
from .validate import validate_plan
from .products import write_products,plot_plan,summarize


def main():
    p=argparse.ArgumentParser();p.add_argument('--date',required=True);p.add_argument('--config',default=str(BASE/'defaults.json'));p.add_argument('--iers',default=str(BASE/'cache/iers.all'));p.add_argument('--campaign');p.add_argument('--previous');p.add_argument('--shift-minutes',type=float,default=0);p.add_argument('--output',required=True)
    a=p.parse_args();c=json.loads(Path(a.config).read_text());iers=setup_iers(a.iers)
    campaign=json.loads(Path(a.campaign).read_text()) if a.campaign else None
    previous=json.loads(Path(a.previous).read_text())['blocks'] if a.previous else None
    plan=make_plan(c,a.date,campaign,previous,a.shift_minutes);plan['iers']=iers
    if campaign:
        from .bright_stars import assess
        plan['bright_star_risks']=assess(campaign,{a.date:[r for b in plan['blocks'] if b['mode']=='trojan' for r in b['rows']]})
    qa=validate_plan(plan,campaign);write_products(plan,qa,a.output);plot_plan(plan,campaign,a.output)
    print(json.dumps({'date':a.date,'validation':qa,'skips':plan['skips'],'fields':summarize(plan)},indent=2))
    if not qa['passed']:raise SystemExit(1)

if __name__=='__main__':main()
