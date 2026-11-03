"""Recheck the entire frozen campaign night by night."""
import argparse
import json
from datetime import date,timedelta
from pathlib import Path
from .catalog import BASE,write_csv
from .geometry import setup_iers,Sky,unix
from .scheduler import trojan


def main():
    p=argparse.ArgumentParser();p.add_argument('--campaign',default=str(BASE/'data/campaign_20261103.json'));p.add_argument('--config',default=str(BASE/'defaults.json'));p.add_argument('--iers',default=str(BASE/'cache/iers.all'));p.add_argument('--output',default=str(BASE/'products/campaign_20261103/visibility_14days.csv'));a=p.parse_args()
    c=json.loads(Path(a.config).read_text());setup_iers(a.iers);s=Sky(c);camp=json.loads(Path(a.campaign).read_text());rows=[]
    d=date.fromisoformat(camp['start_date']);stop=date.fromisoformat(camp['stop_date'])
    while d<stop:
        night=s.night(d.isoformat());b=trojan(s,night['dark'],camp['fields'],camp['tracks'])
        rows.append({'date':d.isoformat(),'eligible_block':bool(b),'start_local':b['not_before_local'] if b else '',
            'duration_minutes':(unix(b['hard_deadline_utc'])-unix(b['not_before_utc']))/60 if b else 0,'exposures':len(b['rows']) if b else 0})
        d+=timedelta(days=1)
    Path(a.output).parent.mkdir(parents=True,exist_ok=True);write_csv(a.output,rows)
    print(json.dumps({'nights':len(rows),'eligible_nights':sum(r['eligible_block'] for r in rows)},indent=2))

if __name__=='__main__':main()
