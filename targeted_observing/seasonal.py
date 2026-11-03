"""Representative seasonal dry runs plus clearly labelled approximate population scope."""
import json
from pathlib import Path
import numpy as np
from .catalog import BASE,approximate,write_csv
from .geometry import setup_iers,Sky,local_time
from .scheduler import make_plan
from .validate import validate_plan
from .products import write_products,plot_plan


def main():
    c=json.loads((BASE/'defaults.json').read_text());eop=setup_iers(BASE/'cache/iers.all')
    cat=json.loads((BASE/'cache/trojans.json').read_text());out=BASE/'products/seasonal';out.mkdir(parents=True,exist_ok=True)
    rows=[];s=Sky(c)
    for d in ['2026-12-21','2027-03-20','2027-06-21','2027-09-23']:
        plan=make_plan(c,d);plan['iers']=eop;qa=validate_plan(plan);write_products(plan,qa,out/d);plot_plan(plan,None,out/d)
        night=s.night(d);r=approximate(cat,d+'T16:00:00');r=[x for x in r if x['V_proxy']<=18.5]
        m=s.metrics([x['ra_deg'] for x in r],[x['dec_deg'] for x in r],np.arange(*night['dark'],300))
        for cloud in ['L4','L5']:
            idx=np.array([x['cloud_proxy']==cloud for x in r]);mask=m['alt']>32.5
            consecutive=[]
            for flags in mask:
                run=best=0
                for flag in flags:
                    run=run+1 if flag else 0;best=max(best,run)
                consecutive.append(max(0,(best-1)*300))
            n=int(np.sum(idx&(np.array(consecutive)>=7200)))
            rows.append({'date':d,'cloud_proxy':cloud,'approx_V18p5_geometric_2h_count':n,
                'evening_start_local':local_time(night['evening'][0]),'morning_end_local':local_time(night['morning'][1]),
                'twilight_blocks':','.join(b['block_id'] for b in plan['blocks']),
                'validation_passed':qa['passed'],
                'limits':'two-body/H-G ranking only; no lunar/crowding cuts for population count; not selectable target count; no future frozen Trojan field'})
        print(d,qa,flush=True)
    write_csv(out/'seasonal_summary.csv',rows)
    if not all(r['validation_passed'] for r in rows):raise SystemExit(1)

if __name__=='__main__':main()
