"""Shared weather ensembles and explicitly hypothetical r-band performance.

Latent weather is generated in advance for reproducibility. A scheduler may read
only the current/past row; future weather must not influence target selection.
"""
import json
import numpy as np
from scipy.special import ndtr
from .catalog import BASE


def ensemble(nights, realizations, seed):
    prior=json.loads((BASE/'data/weather_prior.json').read_text())
    monthly=np.mean(list(prior['observable_percent_by_year'].values()),axis=0)/100
    p=prior['scenario_parameters_not_measured'];rng=np.random.default_rng(seed)
    z=rng.normal(size=(realizations,len(nights)))
    rho=p['latent_ar1_rho']
    for j in range(1,len(nights)):z[:,j]=rho*z[:,j-1]+np.sqrt(1-rho*rho)*z[:,j]
    clear=ndtr(z)<np.array([monthly[int(day[5:7])-1] for day in nights])[None,:]
    shape=(realizations,len(nights),3)  # EV, Trojan, AM: one common night state.
    partial=rng.random(shape)<p['partial_cloud_probability_given_observable']
    return {'observable':np.repeat(clear[:,:,None],3,axis=2),
            'reserved':rng.random(shape)<p['external_block_loss_probability'],
            'cloud_start_fraction':rng.uniform(.05,.65,shape),
            'cloud_duration_fraction':partial*rng.uniform(.15,.5,shape),
            'delay_s':rng.uniform(0,p['delay_max_s'],shape),
            'extra_readout_s':rng.exponential(p['extra_readout_mean_s'],shape),
            'extra_slew_s':rng.exponential(p['extra_slew_mean_s'],shape),
            'seeing_arcsec':rng.lognormal(np.log(p['seeing_median_arcsec']),.2,shape),
            'point_seed':rng.integers(0,2**31-1,size=shape,dtype=np.int64),
            'zeropoint_standard_normal':rng.normal(size=(realizations,len(nights)))}


def random_error(magnitude,airmass,moon_fraction,moon_separation,seeing,exposure_s,scenario):
    """Scenario calibration at r=18, X=1, dark sky, 2 arcsec, 120s.

    Photon/sky approximation includes extinction and lunar-background penalty;
    coefficients are sensitivity assumptions, not a fitted Xinglong sky model.
    """
    mag=np.asarray(magnitude);x=np.maximum(1,np.asarray(airmass))
    moon=np.maximum(0,np.asarray(moon_fraction))*np.exp(-np.asarray(moon_separation)/40)
    sky=1+8*moon
    flux=10**(-.4*(mag-18+.15*(x-1)))
    variance=.5/flux+.5*sky*(np.asarray(seeing)/2)**2/flux**2
    return np.sqrt(.005**2+scenario['sigma_r18_120s']**2*variance*120/np.asarray(exposure_s))


def lightcurve_quality(times,night_indices,errors,zeropoints,valid,thresholds):
    """Count distinct nights; never convert sampling compliance into a period claim."""
    t=np.asarray(times);night=np.asarray(night_indices);err=np.asarray(errors)
    clean=np.asarray(valid,dtype=bool)&np.isfinite(err)&np.isfinite(t)
    good=clean&(err<=.05)&(np.abs(np.asarray(zeropoints))<=thresholds['max_zeropoint_mag'])
    good_nights=[]
    for n in np.unique(night[good]):
        ts=t[good&(night==n)]
        if len(ts)>=thresholds['points_per_night'] and np.ptp(ts)>=thresholds['night_span_s']:
            good_nights.append(int(n))
    baseline=(max(good_nights)-min(good_nights)) if good_nights else 0
    median=float(np.median(err[clean])) if clean.any() else None
    fraction=float(np.mean(err[clean]<=.05)) if clean.any() else 0
    passed=(int(good.sum())>=thresholds['points'] and len(good_nights)>=thresholds['nights']
            and baseline>=thresholds['baseline_days'] and median is not None
            and median<=thresholds['median_random_error_mag'] and fraction>=thresholds['fraction_below_005'])
    return {'quality_lightcurve':bool(passed),'good_points':int(good.sum()),'clean_points':int(clean.sum()),
            'good_nights':len(good_nights),'baseline_days':baseline,'median_random_error_mag':median,
            'fraction_below_005':fraction,'period_reliable':None}
