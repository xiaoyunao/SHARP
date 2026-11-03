"""Evidence-preserving period cross-match. Absence never means unknown period."""
import csv
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
    return out


def evaluate(evidence):
    if evidence is None:return 'unverified',1.0
    lc=evidence['LCDB'];li=evidence['Li2025'];model=evidence['models2026']
    periods=[x['period_h'] for x in [lc,li,model] if x and x['period_h']]
    periods += [x['period_h'] for x in evidence.get('TESS',[]) if x['period_h']]
    if len(periods)>1 and max(periods)/min(periods)>1.1:
        # A low-quality published estimate is a candidate discrepancy, not a demonstrated conflict.
        return 'period_discrepancy_candidate',2.0
    if lc and lc['U'].startswith('3') and lc['period_h']:
        return 'control_candidate',1.0
    if model or evidence['DAMIT']:return 'published_model',.6
    if li or (lc and lc['period_h']):return 'weak_or_published_period',1.5
    return 'unverified',1.0
