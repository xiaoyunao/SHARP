"""Read source-preserving VizieR tables; never translate fitted sigma into LCDB U."""
import csv
import hashlib
import json
from pathlib import Path
import math
from .catalog import BASE,download

TESS_URL='https://vizier.idia.ac.za/viz-bin/asu-tsv?-source=J%2FA%2BA%2F693%2FA66%2Ftess&-out.all=&-out.max=unlimited'
HIP_URL='https://vizier.idia.ac.za/viz-bin/asu-tsv?-source=I%2F239%2Fhip_main&Vmag=%3C6&-out=HIP,RAICRS,DEICRS,Vmag,pmRA,pmDE&-out.max=unlimited'

def table(path):
    lines=[s for s in Path(path).read_text().splitlines() if s and not s.startswith('#')]
    if len(lines)<4:raise ValueError('Empty or invalid VizieR table')
    names=lines[0].split('\t')
    return [{k:v.strip() for k,v in zip(names,row)} for row in csv.reader(lines[3:],delimiter='\t') if len(row)==len(names)]

def tess_records(path):
    rows=table(path);records={};zero=0
    if len(rows)!=9912:raise ValueError('TESS snapshot must contain 9912 rows; possible truncation/schema change')
    for raw in rows:
        p,lo,hi,w=[float(raw[k]) for k in ['Per','b_Per','B_Per','e_Per']]
        if not all(math.isfinite(x) for x in [p,lo,hi,w]) or not lo<=p<=hi or not math.isclose(hi-lo,w,abs_tol=3e-8):raise ValueError('Invalid TESS period bounds')
        key=str(int(raw['MPC'])) if raw['MPC'] else raw['Name']
        if key in records:raise ValueError('Duplicate TESS identifier')
        zero+=w==0
        records[key]={'period_h':p,'lower_h':lo,'upper_h':hi,'interval_width_h':w,
            'zero_width':w==0,'other_periods_h':[float(x) for x in raw['PerOthers'].split(';') if x],
            'raw':raw,'period_kind':'synodic','error_kind':'full interval width, not +/- sigma',
            'source':'J/A+A/693/A66/tess','source_url':TESS_URL,'query_status':'matched',
            'dependence':'Reanalysis of TSSYS-DR1; not an independent observing data set'}
    return records,{'rows':len(rows),'numbered':sum(bool(r['MPC']) for r in rows),'zero_width':zero,'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'url':TESS_URL}

def main():
    cache=BASE/'cache'
    for url,name in [(TESS_URL,'tess2025.tsv'),(HIP_URL,'hip_bright_pm.tsv')]:download(url,cache/name)
    records,meta=tess_records(cache/'tess2025.tsv')
    cat=json.loads((cache/'trojans.json').read_text());idx=cat['fields'].index('pdes');ids={r[idx] for r in cat['data']}
    candidates=list(csv.DictReader((BASE/'data/target_audit_20261103.csv').open()))
    selected={r['pdes']:records.get(r['pdes'],{'query_status':'not_in_this_table'}) for r in candidates}
    meta['tjn_matches']=len(ids&records.keys());meta['candidate_matches']=sum(k in records for k in selected)
    meta['hip_sha256']=hashlib.sha256((cache/'hip_bright_pm.tsv').read_bytes()).hexdigest();meta['hip_url']=HIP_URL
    (BASE/'data/review_sources_20261103.json').write_text(json.dumps({'sources':meta,'TESS2025_candidates':selected},indent=2)+'\n')
    print(meta)

if __name__=='__main__':main()
