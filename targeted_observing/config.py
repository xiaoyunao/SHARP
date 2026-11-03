"""Fail explicitly on unsupported v0 geometry instead of silently ignoring options."""
import math


def validate_config(c):
    if c['fov_deg']!=3 or c['twilight']['fields_per_group']!=4:
        raise ValueError('v0 implements 3-degree footprints and four-field twilight groups only')
    if c['filter'] not in {'g','r','i','z'}:raise ValueError('Unsupported science filter')
    for k in ['readout_s','adjacent_slew_s','far_slew_s','filter_change_s','edge_margin_deg']:
        if not math.isfinite(c[k]) or c[k]<0:raise ValueError('Invalid nonnegative parameter '+k)
    if not 0<c['trojan'].get('min_core_retention_fraction',.8)<=1:raise ValueError('Invalid retention fraction')
    budget=c['trojan']['block_s']
    if not math.isfinite(budget) or (budget!=0 and not 3600<=budget<=7200):raise ValueError('Trojan block must be 0 or 60--120 minutes')
    if not math.isfinite(c['adjacent_max_deg']) or c['adjacent_max_deg']<=0:raise ValueError('Invalid neighbor threshold')
    if c['twilight']['window_s']!=3600 or c['twilight']['visits']!=5:raise ValueError('v0 twilight requires one hour and five visits')
    if c['min_alt_deg']<30 or c['center_alt_deg']<32.5:raise ValueError('Cannot relax v0 altitude safety limits')
    if c['trojan']['sun_alt_deg']!=-18 or c['twilight']['sun_alt_deg']!=-12:raise ValueError('v0 solar crossing levels are fixed')
    if not 0<=c['edge_margin_deg']<1.5:raise ValueError('Invalid field margin')
    for section in ['twilight','trojan']:
        if not math.isfinite(c[section]['exposure_s']) or c[section]['exposure_s']<=0:raise ValueError('Exposure must be positive')
    if not math.isfinite(c['trojan']['handover_s']) or c['trojan']['handover_s']<0 or (budget and c['trojan']['handover_s']>=budget):raise ValueError('Invalid handover')
