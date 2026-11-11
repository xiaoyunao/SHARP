"""Audit saved TPV WCS without mistaking stamps for the whole detector mask."""
import copy
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from .catalog import BASE
from .geometry import coord
from .fixed_fields import Grid,cell_boundary
from .fixed_field_build import OUTPUT,sha,write_json

SOURCES=[
 'outputs/exhibition_20261019/revision/hapke_native_280_metadata.json',
 'outputs/exhibition_20261019/data/unknown_example_00001NM.json',
 'outputs/exhibition_20261019/expanded/figure_data/candidate_00001L5_native_280_metadata.json',
 'outputs/exhibition_20261019/expanded/figure_data/candidate_00001mX_native_280_metadata.json']


def extract():
    saved=[]
    for name in SOURCES:
        path=BASE.parent/name;doc=json.loads(path.read_text())
        for frame in doc['frames']:
            header=copy.deepcopy(frame.get('crop_wcs_header',frame.get('wcs_header')))
            header['CRPIX1']+=frame['crop_origin_xy0'][0];header['CRPIX2']+=frame['crop_origin_xy0'][1]
            row=frame.get('row',frame)
            saved.append({'source_local':name,'source_sha256':sha(path),
                'source_image':frame.get('image_path',frame.get('l1_path')),
                'full_wcs_header':header,'detector_shape':[9232,9216],
                'shape_evidence':'raw_header NAXIS in 00001NM; same instrument shape assumed for remaining snapshots',
                'catalog_origin':1,
                'catalog_pixel_xy0':[row['X_Win']-1,row['Y_Win']-1],
                'catalog_ra_dec':[row['RA_Win'],row['DEC_Win']]})
    return saved


def run(output=OUTPUT):
    saved_path=BASE/'data/tjn_wcs_evidence_20261111.json'
    samples=json.loads(saved_path.read_text()) if saved_path.exists() else extract()
    write_json(saved_path,samples)
    doc=json.loads((output/'registry.json').read_text());grid=Grid(**doc['grid']);js=sorted({f['row_index'] for f in doc['fields']})
    diagnostics=[];coverage=[];boundaries={}
    for i,s in enumerate(samples):
        w=WCS(s['full_wcs_header']);ny,nx=s['detector_shape'];cx,cy=(nx-1)/2,(ny-1)/2
        c=coord(*w.all_pix2world(cx,cy,0));x=coord(*w.all_pix2world(cx+10,cy,0));y=coord(*w.all_pix2world(cx,cy+10,0))
        p=w.all_world2pix([s['catalog_ra_dec']],0)[0];res=np.linalg.norm(p-s['catalog_pixel_xy0'])
        # Existing metadata uses its explicit origin=0 conversion, independently verify it.
        q=np.linspace(0,1,257);xx=np.r_[q*(nx-1),np.full(257,nx-1),(1-q)*(nx-1),np.zeros(257)]
        yy=np.r_[np.zeros(257),q*(ny-1),np.full(257,ny-1),(1-q)*(ny-1)]
        rr,dd=w.all_pix2world(xx,yy,0);px,py=w.all_world2pix(rr,dd,0)
        scale=proj_plane_pixel_scales(w)*3600
        diagnostics.append({'sample':i,'source_image':s['source_image'],'scale_x_arcsec':scale[0],'scale_y_arcsec':scale[1],
            'x_pa_deg':float(c.position_angle(x).deg)%360,'y_pa_deg':float(c.position_angle(y).deg)%360,
            'parity':int(np.sign(np.linalg.det(w.pixel_scale_matrix))),
            'catalog_reconstruction_error_pix':float(res),
            'dense_boundary_roundtrip_error_pix':float(np.max(np.hypot(px-xx,py-yy))),
            'boundary_samples':len(xx),'observed_dec_deg':float(c.dec.deg)})
        boundaries[str(i)]={'ra_deg':rr.tolist(),'dec_deg':dd.tolist(),'status':'measured distortion, assumed full science rectangle; mask unknown'}
        for j in js:
            # Repointing this measured optical model is an engineering scenario, not a high-Dec measurement.
            h=copy.deepcopy(s['full_wcs_header']);h['CRVAL1']=0.;h['CRVAL2']=j*grid.sy;h['LATPOLE']=j*grid.sy
            transported=WCS(h);a,d=cell_boundary(j,grid)
            xp,yp=transported.all_world2pix(a,d,0)
            edge=np.minimum.reduce([xp,yp,nx-1-xp,ny-1-yp])
            margin_pix=(grid.margin+grid.guard)*3600/1.15
            coverage.append({'sample':i,'row_index':j,'dec_deg':j*grid.sy,'min_edge_pix':float(edge.min()),
                'required_quality_plus_guard_pix':margin_pix,'nominal_pass':bool(edge.min()>=0),
                'quality_scenario_pass':bool(edge.min()>=margin_pix),
                'scope':'transported measured TPV; not observed field orientation or template mask'})
    pd.DataFrame(diagnostics).to_csv(output/'wcs_samples.csv',index=False)
    pd.DataFrame(coverage).to_csv(output/'wcs_row_scenarios.csv',index=False)
    write_json(output/'measured_wcs_boundaries.json',boundaries)
    result={'samples':len(samples),'complete_TPV':True,'origin_zero_reconstruction_max_pix':max(r['catalog_reconstruction_error_pix'] for r in diagnostics),
        'dense_boundary_roundtrip_max_pix':max(r['dense_boundary_roundtrip_error_pix'] for r in diagnostics),
        'transported_quality_all_pass':all(r['quality_scenario_pass'] for r in coverage),
        'minimum_transported_edge_pix':min(r['min_edge_pix'] for r in coverage),
        'actual_science_mask':None,'template_common_region':None,'pointing_repeatability':None,
        'hardware_PA_zero_certified':False,'formal_freeze_ready':False,
        'missing_inputs':['full detector bad-pixel/vignette mask','template common crop','repeated-pointing and PA range at TJN low/high declination'],
        'evidence_sha256':sha(saved_path)}
    write_json(output/'wcs_audit.json',result)
    if result['origin_zero_reconstruction_max_pix']>1e-3 or not result['transported_quality_all_pass']:raise ValueError(result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
