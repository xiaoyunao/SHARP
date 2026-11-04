"""Evidence-preserving period cross-match. Absence never means unknown period."""
import csv
import json
from pathlib import Path


def load_periods(cache):
    cache=Path(cache);out={}
    def entry(k):return out.setdefault(str(k),{'LCDB':None,'Li2025':None,'models2026':None,'DAMIT':None,'TESS':[], 'TESS_scope':'DR1 full table; later releases not exhaustive'})
    for line in (cache/'lcdb_summary.txt').read_text().splitlines()[5:]:
        if not line[:7].strip().isdigit():continue
        period=line[145:156].strip();quality=line[187:190].strip()
        try:p=float(period)
        except ValueError:p=None
        entry(str(int(line[:7])))['LCDB']={'period_h':p,'U':quality,'release':'2023-10-01'}
    for row in csv.DictReader((cache/'li2025.csv').open()):
        entry(row['Designation'])['Li2025']={'period_h':float(row['Best_Period [hours]']),
                                             'quality':row['Reliability_Flag'],'comment':row['Comments']}
    for name in ['spin_results','spin_results_partial']:
        for line in (cache/name).read_text().splitlines()[1:]:
            if not line[:6].strip().isdigit():continue
            k=str(int(line[:6]));p=float(line[6:20])
            entry(k)['models2026']={'period_h':p,'solution':name}
    ids={}
    for row in csv.DictReader((cache/'damit_asteroids.csv').open(encoding='utf-8-sig')):
        if row['number']:ids[row['id']]=str(int(row['number']))
    for row in csv.DictReader((cache/'damit_asteroid_models.csv').open(encoding='utf-8-sig')):
        k=ids.get(row['asteroid_id'])
        if k:
            entry(k)['DAMIT']=entry(k)['DAMIT'] or []
            entry(k)['DAMIT'].append({'model_id':row['id'],'period_h':float(row['period']) if row['period'] else None})
    for line in (cache/'tssys_release.merge').read_text().splitlines():
        cols=line.split()
        if not cols or not cols[0].isdigit():continue
        entry(str(int(cols[0])))['TESS'].append({'period_h':float(cols[3]),'type':cols[2],
            'amplitude_mag':float(cols[4]),'points':int(cols[5]),'sector':int(cols[6]),'release':'TSSYS-DR1'})
    from .review_data import tess_records
    path=cache/'tess2025.tsv'
    if path.exists():
        tess,_=tess_records(path)
        for k,r in tess.items():entry(k)['TESS2025']=r
        for ev in out.values():ev.setdefault('TESS2025',{'query_status':'not_in_this_table'})
    else:
        for ev in out.values():ev['TESS2025']={'query_status':'source_unavailable'}
    reviews=json.loads((Path(__file__).parent/'data/manual_science_review.json').read_text())
    for k,review in reviews.items():entry(k)['manual_review']=review
    return out


def assess(evidence):
    """Candidate review flags do not prove unresolved aliases or source independence."""
    if not evidence:return {'role':'literature_unchecked','weight':0.1,'periods':[],'review_reasons':['missing_evidence']}
    periods=[]
    for source in ['LCDB','Li2025','models2026']:
        r=evidence.get(source)
        if r and r.get('period_h'):periods.append({'source':source,**r})
    for source in ['DAMIT','TESS']:
        for r in evidence.get(source) or []:
            if r.get('period_h'):periods.append({'source':source,**r})
    tess=evidence.get('TESS2025',{})
    if tess.get('query_status')=='matched':periods.append({'source':'TESS2025',**tess})
    values=[r['period_h'] for r in periods];reasons=[]
    if values and max(values)/min(values)>1.1:reasons.append('numerical_period_difference_requires_quality_and_literature_review')
    if tess.get('other_periods_h'):reasons.append('TESS_alternative_periods')
    if tess.get('zero_width'):reasons.append('TESS_zero_width_not_infinite_precision')
    if tess.get('interval_width_h',0)>0.1*tess.get('period_h',float('inf')):reasons.append('TESS_broad_interval')
    lc=evidence.get('LCDB') or {};model=evidence.get('models2026') or (evidence.get('DAMIT') or [None])[0]
    if lc.get('U','').startswith('3') and lc.get('period_h'):role='control'
    elif model:role='known_other'  # A newer model may have resolved an older discrepancy.
    elif values:role='unknown_or_weak'  # Remains a review candidate, not an automatic high-value conflict.
    else:role='literature_unchecked'
    weights={'unknown_or_weak':1.,'alias_review':1.,'slow_candidate':.7,'control':.2,'known_other':.1,'literature_unchecked':.1}
    manual=evidence.get('manual_review')
    if manual:
        role=manual['role']
        if role not in weights:raise ValueError('Unsupported manually reviewed science role')
        reasons.append('manual_review: '+manual['reason'])
    if any(p>100 for p in values):reasons.append('long_period_evidence; dense_single_night_science_requires_justification')
    return {'role':role,'weight':weights[role],'periods':periods,'review_reasons':reasons,'manual_review':manual,
            'classification':'conservative automatic triage; alias/slow science role requires explicit review'}

def evaluate(evidence):
    a=assess(evidence)
    return a['role'],a['weight']
