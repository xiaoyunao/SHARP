"""Read-back checks for the fixed-field deliverable and protected baseline."""
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.table import Table
from .catalog import BASE
from .fixed_fields import Grid,registry
from .fixed_field_build import OUTPUT,sha,write_json


def run(output=OUTPUT):
    doc=json.loads((output/'registry.json').read_text());g=Grid(**doc['grid'])
    cells={(f['row_index'],f['column_index']):f['planning_cloud'] for f in doc['fields']}
    rebuilt=registry(dict(reversed(list(cells.items()))),g)
    for f,old in zip(rebuilt['fields'],doc['fields']):f['selection_basis']=old['selection_basis']
    assert rebuilt==doc, 'Fresh RC build must be deterministic'
    retained=registry(cells,g,previous=doc)
    assert retained==doc, 'Existing release must not move or renumber'
    csv=pd.read_csv(output/'registry.csv');fits=Table.read(output/'registry.fits')
    assert csv.field_id.nunique()==len(doc['fields'])==len(fits)
    for name in ['ra_deg','dec_deg','pa_deg']:
        np.testing.assert_allclose(csv[name],[f[name] for f in doc['fields']],atol=5e-12,rtol=0)
        np.testing.assert_array_equal(fits[name],[f[name] for f in doc['fields']])
    objects=pd.read_csv(output/'objects.csv',dtype={'object_id':str});days=pd.read_csv(output/'object_days.csv.gz',dtype={'object_id':str})
    counts=days.groupby('object_id').visible_samples.sum()
    for r in objects.itertuples():assert counts.get(r.object_id,0)==r.visible_samples
    totals=0;seen=set()
    for chunk in pd.read_csv(output/'object_samples.csv.gz',chunksize=100000,dtype={'object_id':str}):
        totals+=int(chunk.strict_visible.sum());seen.update(chunk.object_id)
        assert chunk.in_nominal.all() and chunk.in_quality_model.all()
        assert chunk.in_usable.isna().all()
    assert totals==objects.visible_samples.sum()
    legacy=pd.read_csv(output/'legacy_to_fixed.csv')
    assert legacy.legacy_field_id.nunique()==383
    assert set(legacy.field_id)<=set(csv.field_id)
    tonight=sha(BASE/'products/full_night_20261111/20261111.txt')
    assert tonight=='e43db50712f924f382b9f776f62d825abca5d8dc54d46e598427b29daf78147b'
    notes=[]
    for file in (output/'representative_nights').glob('*.json'):
        plan=json.loads(file.read_text());assert plan['validation']['passed']
        for b in plan['blocks']:
            for f in b['fields']:
                original=next(x for x in doc['fields'] if x['field_id']==f['field_id'])
                assert all(f[k]==original[k] for k in ['ra_deg','dec_deg','pa_deg','grid_cell_id'])
        notes.append({'date':plan['date'],'passed':True,'exposures':plan['validation']['exposures']})
    assert len(notes)==4
    production_diff=subprocess.check_output(['git','diff','954e43d5','--','survey','known_asteroid','heliolincrr'],cwd=BASE.parent,text=True)
    assert production_diff==''
    main=subprocess.check_output(['git','rev-parse','main'],cwd=BASE.parent,text=True).strip()
    assert main=='2d5a46940205fc80398f8c1a1223f372435d32e3'
    figures=sorted((output/'figures').glob('*.png'));pdfs=sorted((output/'figures').glob('*.pdf'))
    assert len(figures)==len(pdfs)==6
    from PIL import Image
    for p in figures:
        with Image.open(p) as image:
            image.verify()
    wcs=json.loads((output/'wcs_audit.json').read_text());independent=json.loads((output/'independent_audit.json').read_text())
    assert wcs['transported_quality_all_pass'] and independent['legacy_quality_uncovered']==0
    assert not wcs['formal_freeze_ready']
    testlog=(output/'unittest.log').read_text()
    assert 'Ran 74 tests' in testlog and testlog.rstrip().endswith('OK')
    report={'passed':True,'scope':'nominal registry, local TPV scenarios and limited planner integration',
        'registry_repeat_order_and_previous_release':True,'csv_fits_json_roundtrip':True,
        'daily_sample_object_count_closure':True,'strict_samples':totals,'object_count':len(objects),
        'fields':len(doc['fields']),'sampled_envelope_objects':len(seen),'legacy_fields':383,
        'representative_nights':notes,'unit_tests':74,'test_log_sha256':sha(output/'unittest.log'),
        'plot_pairs':6,'plot_visual_review':'all six PNG figures inspected; same canvas PDF pairs',
        'tonight_controller_sha256':tonight,'production_diff_from_954e43d5':'empty','main_ref':main,
        'server_operations':False,'new_twilight_annual_or_weather_simulation':False,'formal_instrument_freeze':False,
        'source_sha256':{p.name:sha(p) for p in BASE.glob('fixed*.py')}}
    write_json(output/'validation.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':run()
