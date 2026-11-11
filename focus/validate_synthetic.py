"""Independent analytic Gaussian/Moffat/donut and contaminated-source checks."""
import json
import numpy as np
from measure_focus import measure_tile
rng=np.random.default_rng(7)
y,x=np.indices((800,800));results=[]
for kind,width in [('gaussian',2.),('gaussian',4.),('moffat',3.),('ring',5.)]:
    data=rng.normal(500,4,x.shape)
    for cy in [100.2,300.4,500.1,700.3]:
        for cx in [100.3,300.2,500.4,700.1]:
            r2=(x-cx)**2+(y-cy)**2
            if kind=='gaussian':psf=np.exp(-r2/(2*width**2));expected=2*width*np.sqrt(2*np.log(2))
            elif kind=='moffat':psf=(1+r2/width**2)**-3.5;expected=2*width*np.sqrt(2**(1/2.5)-1)
            else:psf=np.exp(-(np.sqrt(r2)-width)**2/(2*.8**2));expected=None
            data+=2000*psf
    rows,_,stats=measure_tile(data.astype('float32'),0,0,0)
    measured=float(np.median([r['d50'] for r in rows]))
    assert len(rows)==16,(kind,stats)
    if expected:assert abs(measured/expected-1)<.06,(kind,expected,measured)
    if kind=='ring':assert np.mean([r['donut'] for r in rows])>.9
    results.append(dict(kind=kind,width=width,expected_d50=expected,measured_d50=measured,n=len(rows)))
assert results[1]['measured_d50']>results[0]['measured_d50']*1.9
# Saturated source and a one-pixel cosmic ray must not enter accepted catalogs.
data=np.random.default_rng(8).normal(500,4,(800,800)).astype('float32')
data+=np.minimum(100000*np.exp(-((x-300.2)**2+(y-300.3)**2)/8),51000) # Real detector clips below header SATURATE.
data[500,500]=50000
rows,_,stats=measure_tile(data,0,0,0)
assert len(rows)==0,stats
print(json.dumps(dict(status='passed',analytic_cases=results,contamination=stats),indent=2))
