"""Publish a reviewable local artifact and small versioned registry snapshot."""
import argparse
import json
from pathlib import Path
import shutil
import zipfile
from .catalog import BASE
from .fixed_field_build import OUTPUT,sha,write_json


def deliver(output=OUTPUT):
    validation=json.loads((output/'validation.json').read_text())
    if not validation['passed']:raise ValueError('Validate deliverable before packaging')
    if validation['source_sha256']!={p.name:sha(p) for p in BASE.glob('fixed*.py')}:
        raise ValueError('Source changed since validation')
    # Exclude the large derived sampling cache; source ephemerides stay in the existing cache.
    members=sorted(p for p in output.rglob('*') if p.is_file() and p.suffix not in ('.zip','.npz') and p.name!='MANIFEST_SHA256.json')
    manifest={str(p.relative_to(output)):sha(p) for p in members}
    write_json(output/'MANIFEST_SHA256.json',manifest)
    archive=output/'SHARP_TJN_fixed_fields_RC1_20261111.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in members+[output/'MANIFEST_SHA256.json']:z.write(p,arcname='SHARP_TJN_fixed_fields_RC1_20261111/'+str(p.relative_to(output)))
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError('ZIP CRC failed')
        import hashlib
        for name,h in manifest.items():
            if hashlib.sha256(z.read('SHARP_TJN_fixed_fields_RC1_20261111/'+name)).hexdigest()!=h:raise ValueError('Archive hash mismatch')
    data=BASE/'data/tjn_fixed_20261111';data.mkdir(exist_ok=True)
    shutil.copy2(output/'registry.json',data/'registry.json')
    docs=BASE/'docs/tjn_fixed_20261111';docs.mkdir(exist_ok=True)
    selected=['REPORT.md','summary.json','independent_audit.json','wcs_audit.json','row_geometry_audit.csv','wcs_samples.csv','wcs_row_scenarios.csv',
        'representative_nights.csv','special_targets.csv','objects.csv','seed_vs_final_objects.csv','field_occupancy.csv','overlap_occupancy.csv',
        'TJ_test_retention.csv','TJ_test_provenance.json','provenance.json','validation.json','legacy_to_fixed.csv','neighbors.csv','TJL4.csv','TJL5.csv']
    for name in selected:shutil.copy2(output/name,docs/name)
    shutil.copytree(output/'figures',docs/'figures',dirs_exist_ok=True)
    write_json(docs/'delivery.json',{'archive_local':str(archive),'archive_sha256':sha(archive),'archive_bytes':archive.stat().st_size,
        'archive_members':len(members)+1,'registry_sha256':sha(data/'registry.json'),
        'source_commit':'see git commit containing this delivery; baseline 7747cc16',
        'snapshot_files_sha256':{str(p.relative_to(docs)):sha(p) for p in docs.rglob('*') if p.is_file() and p.name!='delivery.json'}})
    print(json.dumps({'archive':str(archive),'size_mb':archive.stat().st_size/1e6,'members':len(members)+1},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=OUTPUT);args=parser.parse_args();deliver(args.output)
