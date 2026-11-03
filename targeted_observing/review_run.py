"""Reproduce the three-night, three-campaign review using cached public inputs."""
import json
from pathlib import Path
from .catalog import BASE
from .geometry import Sky,setup_iers
from .scheduler import make_plan,trojan
from .validate import validate_plan
from .products import write_products,plot_plan
from .bright_stars import assess

def main():
    config=json.loads((BASE/'defaults.json').read_text());setup_iers(BASE/'cache/iers.all')
    base=BASE/'products/pro_review';comparison=json.loads((base/'campaigns/comparison.json').read_text());results=[]
    for item in comparison['candidates']:
        camp=json.loads(Path(item['campaign_file']).read_text());previous=None
        for date in ['2026-11-03','2026-11-04','2026-11-06']:
            plan=make_plan(config,date,camp,previous)
            plan['bright_star_risks']=assess(camp,{date:[r for b in plan['blocks'] if b['mode']=='trojan' for r in b['rows']]})
            qa=validate_plan(plan,camp);out=base/'plans'/item['label']/date
            write_products(plan,qa,out);plot_plan(plan,camp,out)
            if not qa['passed']:raise ValueError(qa)
            previous=plan['blocks'];result={'candidate':item['label'],'night':date,'exposures':qa['exposures'],
                'anchors':{b['block_id']:b.get('retained_anchor_count') for b in plan['blocks']}}
            print(result,flush=True);results.append(result)
        # Audit each of the 14 reserved nights independently, including targets leaving the field.
        sky=Sky(config)
        for date in camp['planned_visits']:
            b=trojan(sky,sky.night(date)['dark'],camp['fields'],camp['tracks'])
            qa=validate_plan({'config':config,'date':date,'blocks':[b]},camp)
            if not qa['passed']:raise ValueError((item['label'],date,qa))
    (base/'verification.json').write_text(json.dumps({'plans':results,'campaign_nights_audited':42,'passed':True},indent=2)+'\n')
if __name__=='__main__':main()
