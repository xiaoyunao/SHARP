"""Joint A/B choice with explicit old fields and local spherical repair candidates."""
import numpy as np
from datetime import datetime,timedelta
from zoneinfo import ZoneInfo
from .geometry import coord,offset

def plan(sky,window,side,previous=None):
    from .scheduler import twilight_candidates,field,visits,check_rows,block_record,closed_neighbor_tour
    c=sky.config;tc=c['twilight'];a,b=window
    names=['A','B'] if side=='EV' else ['B','A'];slots=[(a,(a+b)/2),((a+b)/2,b)]
    cycles=tc['visits']
    if cycles==10:names=['A'];slots=[(a,b)]
    pools=[];audit=[];future=[]
    if not previous:
        night_date=datetime.fromtimestamp(a,ZoneInfo('Asia/Shanghai')).date()
        if side=='AM':night_date-=timedelta(days=1)
        for days in c['twilight'].get('anchor_lookahead_days',[1,3]):
            n=sky.night((night_date+timedelta(days=days)).isoformat())
            future.append(n['evening' if side=='EV' else 'morning'])
    for label,(start,end) in zip(names,slots):
        from .geometry import unix
        name=f'{side}_{label}';old=next((x for x in previous or [] if x['block_id']==name),None)
        oldids={f['field_id'] for f in old['fields']} if old else set()
        pool=[];raw=[]
        if not any(max(start,unix(r['start_utc']))<min(end,unix(r['end_utc'])) for r in c['reservations']):
            if old:
                raw.append((old['fields'],'previous_group'))
                rr=visits(old['fields'],start,end,tc['exposure_s'],cycles,c,name,'twilight')
                for i,f in enumerate(old['fields']):
                    if check_rows(sky,[f],[r for r in rr if r['field_id']==f['field_id']],'twilight',True):continue
                    for size in [.15,.3,.45]:
                        for dx,dy in [(size,size),(size,0),(0,size),(-size,size),(size,-size),(-size,0),(0,-size),(-size,-size)]:
                            p=offset(coord(f['ra_deg'],f['dec_deg']),dx,dy)
                            fs=list(old['fields']);fs[i]=field(p.ra.deg,p.dec.deg,version=c['version'])
                            raw.append((fs,'local_repair_of_'+f['field_id']))
            raw.extend((fs,'new_grid') for _,_,_,fs in twilight_candidates(sky,start,end,side,c['version'])[:8])
        for fs,reason in raw:
            if not closed_neighbor_tour(fs,c):continue
            rr=visits(fs,start,end,tc['exposure_s'],cycles,c,name,'twilight')
            if len(rr)!=4*cycles:continue
            q=check_rows(sky,fs,rr,'twilight',True)
            if not q:continue
            ext=any(v['elongation_max_deg']>tc['elongation_deg'][1] for v in q.values())
            block=block_record(name,fs,rr,start,end,q,ext)
            block['retained_anchor_count']=sum(f['field_id'] in oldids for f in fs)
            block['selection_reason']=reason
            block['replaced_field_ids']=sorted(oldids-{f['field_id'] for f in fs})
            block['replacement_reason']='previous group not feasible at scheduled visits or joint A/B incompatibility' if block['replaced_field_ids'] else None
            future_passes=0
            for fa,fb in future:
                index=names.index(label);lo=fa+index*(fb-fa)/len(names);hi=lo+(fb-fa)/len(names)
                future_rows=visits(fs,lo,hi,tc['exposure_s'],cycles,c,name,'twilight')
                if len(future_rows)==4*cycles and check_rows(sky,fs,future_rows,'twilight',True):future_passes+=1
            block['anchor_lookahead_valid_nights']=future_passes
            pool.append(block)
        pool.sort(key=lambda x:(-x['retained_anchor_count'],-x['anchor_lookahead_valid_nights'],x['extended_search'],-min(q['min_center_alt_deg'] for q in x['quality'].values())))
        pools.append(pool[:12])
        audit.append({'block_id':name,'generated_candidates':len(raw),'feasible_candidates':len(pool),
                      'alternatives':[{'field_ids':[f['field_id'] for f in x['fields']],'retained':x['retained_anchor_count'],'reason':x['selection_reason']} for x in pool[:12]]})
    if len(pools)==1:
        chosen=pools[0][:1]
        for block in chosen:block['joint_selection_audit']=audit
        return chosen
    choices=[]
    for x in pools[0]+[None]:
        for y in pools[1]+[None]:
            selected=[z for z in [x,y] if z]
            if x and y:
                f,g=x['fields'],y['fields']
                sep=coord([z['ra_deg'] for z in f],[z['dec_deg'] for z in f])[:,None].separation(coord([z['ra_deg'] for z in g],[z['dec_deg'] for z in g])[None,:]).deg
                if sep.min()<2.6:continue
                el=sky.metrics([z['ra_deg'] for z in f+g],[z['dec_deg'] for z in f+g],[(slots[1][0]+slots[1][1])/2])['elong'].ravel()
                if (side=='EV' and el[:4].mean()>el[4:].mean()) or (side=='AM' and el[:4].mean()<el[4:].mean()):continue
            choices.append(((len(selected),sum(z['retained_anchor_count'] for z in selected),sum(z['anchor_lookahead_valid_nights'] for z in selected),-sum(z['extended_search'] for z in selected)),selected))
    chosen=max(choices,key=lambda z:z[0])[1]
    for block in chosen:block['joint_selection_audit']=audit
    return chosen
