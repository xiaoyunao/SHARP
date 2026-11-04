"""Dry-run execution gate. Returns decisions only; never sends telescope commands."""
import copy
from .geometry import Sky,unix,isot
from .scheduler import check_rows


def gate(plan,block_id,row_index,ready_utc,campaign=None):
    """Call after real slew/filter readiness; early readiness waits, late readiness revalidates.

    The caller must track completed row indices, actual telescope state and readout completion.
    This module supplies a timing/geometry contract, not an observatory driver.
    """
    block=next(b for b in plan['blocks'] if b['block_id']==block_id)
    if not 0<=row_index<len(block['rows']):raise ValueError('Invalid row index')
    row=block['rows'][row_index];now=unix(ready_utc);scheduled=unix(row['start_utc'])
    if now<scheduled:return {'action':'wait','until_utc':row['start_utc']}
    delay=now-scheduled
    remaining=copy.deepcopy(block['rows'][row_index:])
    for r in remaining:
        for key in ['overhead_start_utc','start_utc','mid_utc','end_utc','release_utc']:
            r[key]=isot(unix(r[key])+delay)
    cutoff=unix(block['hard_deadline_utc'])
    if block_id=='TJ':cutoff-=plan['config']['trojan']['handover_s']
    if unix(remaining[-1]['release_utc'])>cutoff:
        return {'action':'abort_remaining_block','reason':'cannot finish remaining sequence before deadline'}
    present={r['field_id'] for r in remaining};fields=[f for f in block['fields'] if f['field_id'] in present]
    if block_id=='TJ' and campaign is None:raise ValueError('Campaign required for Trojan execution gate')
    q=check_rows(Sky(plan['config']),fields,remaining,block['mode'],block['extended_search'],
                 campaign['tracks'] if campaign and block_id=='TJ' else None)
    if q is None:return {'action':'abort_remaining_block','reason':'actual timing violates geometry or target retention'}
    for f in fields:
        selected=[r for r in remaining if r['field_id']==f['field_id']]
        sampling=q[f['field_id']].get('target_sampling',{})
        for i,r in enumerate(selected):
            r['target_states']={v['object_id']:{'in_field':v['nights']['night']['in_field'][i],
                'photometric_usable':None} for v in sampling.values()}
    return {'action':'dry_run_exposure_eligible','exposure':remaining[0],'delay_s':delay,
            'target_sampling_at_actual_times':{fid:v.get('target_sampling',{}) for fid,v in q.items()},
            'hardware_command_sent':False,'remaining_sequence_revalidated':True}
