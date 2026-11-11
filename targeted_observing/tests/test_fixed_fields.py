import copy
import unittest
import numpy as np
from astropy.wcs import WCS
from targeted_observing.fixed_fields import Grid,registry,FixedFields,tan_xy,from_tan,edge_distance,cell_boundary


class FixedFieldContracts(unittest.TestCase):
    def test_wrap_high_dec_and_inverse(self):
        g=Grid();j,k=g.assign(np.array([359.99,.01,180]),np.array([0,0,60]))
        self.assertEqual(k[0],k[1]);self.assertEqual(j[2],24)
        ra,dec=from_tan(359.9,60,np.array([-.5,0,.5]),np.array([.3,0,-.3]))
        x,y=tan_xy(ra,dec,359.9,60)
        np.testing.assert_allclose(x,[-.5,0,.5],atol=1e-12);np.testing.assert_allclose(y,[.3,0,-.3],atol=1e-12)
        w=WCS(naxis=2);w.wcs.crval=[359.9,60];w.wcs.crpix=[1,1];w.wcs.cdelt=[1,1];w.wcs.ctype=['RA---TAN','DEC--TAN']
        np.testing.assert_allclose(w.all_world2pix(np.c_[ra,dec],0),np.c_[x,y],atol=1e-10)

    def test_coverage_and_original_pitch_regression(self):
        g=Grid()
        for j in range(-12,25):
            a,d=cell_boundary(j,g,samples=301)
            self.assertGreater(edge_distance(a,d,0,j*g.sy,g).min(),g.margin+g.guard)
        old=Grid(sx=2.67264,sy=2.67496)
        a,d=cell_boundary(20,old)
        self.assertLess(edge_distance(a,d,0,20*old.sy,old).min(),old.margin)

    def test_registry_append_reorder_retirement(self):
        cells={(1,3):'L4',(2,4):'L5',(3,5):'L4'}
        first=registry(cells);self.assertEqual(first,registry(dict(reversed(list(cells.items())))))
        first['fields'][0]['status']='RETIRED';before=copy.deepcopy(first)
        with self.assertRaises(ValueError):registry({(-1,2):'L4'},previous=first)
        second=registry({(-1,2):'L4',(1,3):'L5'},previous=first,release='later')
        self.assertEqual(first,before);self.assertEqual(second['fields'][:3],before['fields'])
        self.assertEqual(second['fields'][3]['field_id'],'TJL4_003')
        with self.assertRaises(ValueError):registry(cells,Grid(sx=2.4),previous=first)

    def test_overlap_unique_and_template_identity(self):
        doc=registry({(0,0):'L4',(0,1):'L5'})
        r=FixedFields(doc);hits=r.locate([1.25],[0],quality=True)[0]
        self.assertEqual(len(hits),2)
        fid=hits[0]
        self.assertNotEqual(r.template_key(fid,'r','a'),r.template_key(fid,'g','a'))
        with self.assertRaises(ValueError):r.template_key(fid,'W','a')
        doc['fields'][0]['status']='RETIRED'
        self.assertNotIn(fid,r.locate([1.25],[0],quality=True)[0])

    def test_dynamic_query_moves_objects_not_fields(self):
        from astropy.time import Time
        g=Grid();r=FixedFields(registry({(0,0):'L4',(0,1):'L4'}));before=copy.deepcopy(r.document)
        t0,t1=Time(['2026-11-01','2026-12-15']).jd
        tracks={'demo':{'datetime_jd':[t0,t1],'RA':[0.,2.5],'DEC':[0.,0.],'V':[18.,18.]}}
        first=r.at(tracks,Time(t0,format='jd').unix)
        later=r.at(tracks,Time(t1,format='jd').unix)
        self.assertEqual(first[0]['field_id'],'TJL4_001')
        self.assertEqual(later[0]['field_id'],'TJL4_002')
        self.assertEqual(r.document,before)
        self.assertEqual(first,r.at(tracks,Time(t0,format='jd').unix))

    def test_saved_tpv_origin_and_edges(self):
        from pathlib import Path
        import json
        path=Path(__file__).resolve().parents[1]/'data/tjn_wcs_evidence_20261111.json'
        for row in json.loads(path.read_text()):
            w=WCS(row['full_wcs_header'])
            xy=w.all_world2pix([row['catalog_ra_dec']],0)[0]
            np.testing.assert_allclose(xy,row['catalog_pixel_xy0'],atol=1e-6)
            points=np.array([[0,0],[9215,0],[9215,9231],[0,9231],[4608,4616]])
            np.testing.assert_allclose(w.all_world2pix(w.all_pix2world(points,0),0),points,atol=1e-6)

    def test_ids_beyond_999(self):
        cells={(j,k):'L4' for j in range(8) for k in range(140) if k<Grid().row(j)[1]}
        doc=registry(cells)
        self.assertIn('TJL4_1000',[r['field_id'] for r in doc['fields']])

if __name__=='__main__':unittest.main()
