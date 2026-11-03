"""Independent product audit: chronology, budgets, constraints and target retention."""
from collections import Counter
import numpy as np
from .geometry import Sky,unix,coord,footprint,in_field
from .catalog import interpolate
from .scheduler import slew


def validate_plan(plan,campaign=None):
    c=plan['config'];sky=Sky(c);errors=[];nrows=0;minimum=90;prior=-np.inf
    expected_night=sky.night(plan['date'])
    for block in sorted(plan['blocks'],key=lambda b:unix(b['not_before_utc'])):
        name=block['block_id'];rows=block['rows'];fs={f['field_id']:f for f in block['fields']}
        start,end=unix(block['not_before_utc']),unix(block['hard_deadline_utc'])
        if start<prior-.002:errors.append(name+': overlapping blocks')
        prior=end
        window=expected_night['dark'] if name=='TJ' else expected_night['evening' if name.startswith('EV') else 'morning']
        if start<window[0]-.02 or end>window[1]+.02:errors.append(name+': outside solar window')
        if name=='TJ' and (c['trojan']['block_s']==0 or end-start>c['trojan']['block_s']+.01):errors.append('TJ: configured budget exceeded or disabled')
        if name=='TJ' and (end-start>7200.01 or end-start<3599.99):errors.append('TJ: block duration')
        for r in c['reservations']:
            if max(start,unix(r['start_utc']))<min(end,unix(r['end_utc'])):errors.append(name+': reservation overlap')
        if name!='TJ':
            route=block['fields']
            for i,f in enumerate(route):
                g=route[i-1]
                if coord(f['ra_deg'],f['dec_deg']).separation(coord(g['ra_deg'],g['dec_deg'])).deg>c['adjacent_max_deg']:
                    errors.append(name+': closed tour exceeds neighbor threshold')
        counts=Counter(r['field_id'] for r in rows)
        if name!='TJ' and (len(counts)!=4 or set(counts.values())!={5}):errors.append(name+': incomplete five-visit group')
        if name=='TJ' and len(set(counts.values()))!=1:errors.append('TJ: incomplete cycle')
        if name=='TJ' and end-start<5400-.01 and len(fs)!=1:errors.append('TJ: short block must use primary field')
        for fid in counts:
            seq=[r['visit'] for r in rows if r['field_id']==fid]
            if seq!=list(range(1,len(seq)+1)):errors.append(name+': duplicate or unordered visit')
        previous=None;ready=start
        for r in rows:
            nrows+=1
            f=fs[r['field_id']];a,b=unix(r['start_utc']),unix(r['end_utc'])
            if r['ra_deg']!=f['ra_deg'] or r['dec_deg']!=f['dec_deg'] or r['field_version']!=f['version']:
                errors.append(name+': exposure/field definition mismatch')
            if abs(unix(r['release_utc'])-b-c['readout_s'])>.003:errors.append(name+': readout release mismatch')
            needed=slew(previous,f,c)+(c['filter_change_s'] if previous is None else 0)
            if abs(unix(r['overhead_start_utc'])-ready)>.003:errors.append(name+': overhead chronology')
            if a<ready+needed-.003:errors.append(name+': missing transition budget')
            if abs((a+b)/2-unix(r['mid_utc']))>.002:errors.append(name+': incorrect midpoint')
            if abs(b-a-r['exposure_s'])>.002:errors.append(name+': exposure duration')
            ready=b+c['readout_s']
            if ready>end+.002:errors.append(name+': deadline overrun')
            if name=='TJ' and ready>end-c['trojan']['handover_s']+.002:errors.append(name+': missing handover')
            if r['filter']!=c['filter']:errors.append(name+': filter mismatch')
            previous=f
        # Vectorized dense samples of every exposure, independent of planner's whole-span grid.
        all_ra=[];all_dec=[];all_t=[]
        for r in rows:
            a,b=unix(r['start_utc']),unix(r['end_utc'])
            ts=np.linspace(a,b,int(np.ceil((b-a)/5))+1)
            all_t.extend(ts);all_ra.extend([r['ra_deg']]*len(ts));all_dec.extend([r['dec_deg']]*len(ts))
        frame,sun,moon,fraction=sky.environment(np.array(all_t))
        center=coord(all_ra,all_dec).transform_to(frame)
        # Circumscribed cap encloses the complete 3x3 footprint, regardless of PA.
        radius=2.122;guard=.022
        lower=center.alt.deg-radius-guard
        minimum=min(minimum,float(lower.min()))
        if np.any(lower<=c['min_alt_deg']):errors.append(name+': footprint below altitude limit')
        if np.any(center.alt.deg<c['center_alt_deg']+guard):errors.append(name+': center margin')
        limits=np.where(fraction>=c['bright_moon_fraction'],c['bright_moon_sep_deg'],c['moon_sep_deg'])
        if np.any((center.separation(moon).deg-radius<limits+guard)&(moon.alt.deg> -guard)):errors.append(name+': moon clearance')
        sunlim=c['trojan']['sun_alt_deg'] if name=='TJ' else c['twilight']['sun_alt_deg']
        if np.any(sun.alt.deg>sunlim+.002):errors.append(name+': solar altitude')
        if name!='TJ':
            el=center.separation(sun).deg
            hi=c['twilight']['extended_elongation_deg'][1] if block['extended_search'] else c['twilight']['elongation_deg'][1]
            if el.min()<c['twilight']['elongation_deg'][0] or el.max()>hi:errors.append(name+': solar elongation')
        if name=='TJ':
            if not campaign:errors.append('TJ: campaign required for target validation')
            else:
                if not campaign['start_date']<=plan['date']<campaign['stop_date']:errors.append('TJ: expired field version')
                for f in fs.values():
                    for k in f['core_targets']:
                        track=campaign['tracks'][k]
                        if not np.all(in_field(track['RA'],track['DEC'],f['ra_deg'],f['dec_deg'],c['fov_deg'],c['edge_margin_deg'])):errors.append(k+': fails full campaign retention')
                        ts=[unix(r['mid_utc']) for r in rows if r['field_id']==f['field_id']]
                        ra,dec=interpolate(track,ts)
                        if not np.all(in_field(ra,dec,f['ra_deg'],f['dec_deg'],c['fov_deg'],c['edge_margin_deg'])):errors.append(k+': leaves actual exposure footprint')
    return {'passed':not errors,'errors':errors,'exposures':nrows,'blocks':len(plan['blocks']),
            'min_footprint_alt_bound_deg':minimum if nrows else None,
            'sample_step_s':5,'whole_footprint_bound':'center altitude minus circumscribed cap minus time-grid guard',
            'on_sky_validation':False}
