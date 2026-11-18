"""Simulation-only acceleration; preserve the production route equations.

Numba is optional outside this annual experiment. Production never imports it.
Final checks can interpolate the exact Astropy selection grid; direct checks on
representative dates and interpolation error certificates accompany the result.
"""
import math
import numpy as np
from numba import njit


@njit(cache=True)
def factor(sun,alt,elong,phase,moonsep,moonalt,dark,enhancement):
    x=min(10.,max(0.,-sun-5));mu=11.40+1.567*x-.064*x*x
    if sun < -15:mu=20.67+(dark-20.67)*min(1.,max(0.,(-sun-15)/3))
    excess=max(0.,10**(.4*(dark-mu))-1.)
    if sun<=-18:excess=0.
    return (1+excess*(1+enhancement*math.exp(-elong/35))+(8*phase*math.exp(-moonsep/40) if moonalt>0 else 0))/math.sin(math.radians(max(alt,5.)))


@njit(cache=True)
def dose(sun,alt,elong,phase,moonsep,moonalt,c):
    if sun<=c[9]:return c[1]
    anchor=factor(-10.,25.,40.,0.,180.,-10.,c[6],c[7])
    value=c[5]*anchor/factor(sun,alt,elong,phase,moonsep,moonalt,c[6],c[7])
    return min(c[1],max(c[0],math.floor(value/c[2])*c[2]))


@njit(cache=True)
def route(times,alt,elong,moonsep,sun,phase,moonalt,ha,positions,radii,ids,bounds,c,cycles):
    out=np.empty((300,15));n=0;t=bounds[0];last=-1;visit=0;done=False
    while not done and (cycles<0 or visit<cycles):
        visit+=1
        for i in ids:
            if last<0:turn=c[10]
            elif last==i:turn=0.
            else:
                dra=abs((positions[i,0]-positions[last,0]+180)%360-180)
                ddec=abs(positions[i,1]-positions[last,1]);turn=9.11+max(11.95*math.sqrt(dra),12.30*math.sqrt(ddec))
            change=c[4] if last<0 else 0.;begin=t+turn+change
            if begin+c[0]+c[3]+c[11]>bounds[1]:
                if cycles>=0:return out[:0]
                done=True;break
            early=max(bounds[0],begin-c[12]-c[13]*n);late=min(begin+c[11]+c[1],bounds[1]);ex=c[1]
            for q in (early,begin,late):
                if q<times[0]-.001 or q>times[-1]+.001:return out[:0]
                ex=min(ex,dose(np.interp(q,times,sun),np.interp(q,times,alt[i]),np.interp(q,times,elong[i]),
                               np.interp(q,times,phase),np.interp(q,times,moonsep[i]),np.interp(q,times,moonalt),c))
            finish=begin+ex;end=finish+c[11]
            if finish+c[3]+c[11]>bounds[1]:
                if cycles>=0:return out[:0]
                done=True;break
            if early<bounds[0]-.001 or end>bounds[1]+.001:return out[:0]
            lo=np.searchsorted(times,early,side='right');hi=np.searchsorted(times,end,side='left')
            for k in range(lo-2,hi):
                q=early if k==lo-2 else end if k==lo-1 else times[k]
                a=np.interp(q,times,alt[i]);m=np.interp(q,times,moonsep[i]);up=np.interp(q,times,moonalt)
                ph=np.interp(q,times,phase);el=np.interp(q,times,elong[i]);s=np.interp(q,times,sun)
                if not (a-radii[i]>.15+c[8] and (up<-.15 or m-radii[i]>(40 if ph>=.75 else 30)+.15)
                        and el<=c[14] and s<=c[15]+.003 and abs(np.interp(q,times,ha[i]))<c[19]):return out[:0]
            mid=begin+ex/2;s=np.interp(mid,times,sun);a=np.interp(mid,times,alt[i]);el=np.interp(mid,times,elong[i])
            sk=factor(s,a,el,np.interp(mid,times,phase),np.interp(mid,times,moonsep[i]),np.interp(mid,times,moonalt),c[6],c[7])
            out[n]=np.array([i,visit,ex,t,begin,mid,finish,finish+c[3],turn,c[3],change,s,a,el,sk]);n+=1
            if n>=300:return out[:0]
            t=finish+c[3];last=i
    for i in ids:
        count=0;previous=0.
        for j in range(n):
            if int(out[j,0])!=i:continue
            if count and not c[17]<=out[j,5]-previous<=c[18]:return out[:0]
            previous=out[j,5];count+=1
        if count<c[16]:return out[:0]
    return out[:n]


def simulate_route(g,ids,cycles,c):
    if not hasattr(g,'_annual_numeric'):
        g._annual_numeric=(np.asarray([[f['ra_deg'],f['dec_deg']] for f in g.fields]),np.asarray([f['radius_deg'] for f in g.fields]))
    positions,radii=g._annual_numeric
    values=np.array([c[k] for k in ('min_exposure_s','max_exposure_s','exposure_quantum_s','readout_s','filter_change_s',
        'anchor_exposure_s','dark_sky_R_mag_arcsec2','near_sun_enhancement','min_footprint_alt_deg','dark_sun_deg',
        'initial_slew_s','late_tolerance_s','initial_early_tolerance_s','early_per_frame_s','max_elongation_deg',
        'sun_start_deg','minimum_visits','cadence_min_s','cadence_max_s')]+[c.get('max_tracking_hour_angle_h',4.95)],dtype=float)
    m=g.m;rows=route(g.times,m['alt'],m['elong'],m['moon_sep'],m['sun_alt'],m['moon_fraction'],m['moon_alt'],g.ha,positions,radii,
                    np.asarray(ids,dtype=np.int64),np.asarray(g.bounds),values,-1 if cycles is None else cycles)
    if not len(rows):return None
    result=[]
    for row in rows:
        i=int(row[0]);f=g.fields[i]
        result.append(dict(field_id=f['field_id'],ra_deg=f['ra_deg'],dec_deg=f['dec_deg'],field_index=i,visit=int(row[1]),
            filter='r',exposure_s=int(row[2]),overhead_start_unix=row[3],start_unix=row[4],mid_unix=row[5],end_unix=row[6],
            release_unix=row[7],slew_s=row[8],readout_s=row[9],filter_change_s=row[10],sun_alt_deg=row[11],
            alt_deg=row[12],elongation_deg=row[13],sky_factor=row[14]))
    return result


class GridSky:
    def __init__(self,sky,shared,yso):
        self.location=sky.location;self.shared=shared;self.yso=yso
        self.lookup={}
        for g in (shared,yso):
            for i,f in enumerate(g.fields):self.lookup[(f['ra_deg'],f['dec_deg'])]=(g,i)

    def metrics(self,ra,dec,seconds):
        q=np.atleast_1d(seconds)
        if len(q)==1 and len(np.atleast_1d(ra))==1:
            exact=getattr(self,'scalar_exact',{}).get((float(np.atleast_1d(ra)[0]),float(np.atleast_1d(dec)[0]),float(q[0])))
            if exact is not None:return exact
        out={key:[] for key in ('alt','elong','moon_sep')}
        for a,d in zip(np.atleast_1d(ra),np.atleast_1d(dec)):
            g,i=self.lookup[(float(a),float(d))]
            if q.min()<g.times[0] or q.max()>g.times[-1]:raise ValueError('Grid interpolation outside coverage')
            for key in out:out[key].append(np.interp(q,g.times,g.metrics[key][i]))
        out={key:np.asarray(v) for key,v in out.items()}
        for key in ('sun_alt','moon_alt','moon_fraction'):out[key]=np.interp(q,self.shared.times,self.shared.metrics[key])
        return out
