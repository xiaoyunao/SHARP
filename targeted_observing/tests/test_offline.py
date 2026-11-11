import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from astropy.table import Table
from astropy.time import Time
from targeted_observing.unknown_pipeline.exposure_timing import header_mid_mjd
from targeted_observing.unknown_pipeline.make_tracklet_linreproj import pair_two_exposures, image_center_and_radius
from targeted_observing.offline_twilight import (
    load_config, load_grid, slew, route_order, exposure, write_plan)


class OfflineContracts(unittest.TestCase):
    def setUp(self):self.c=load_config()

    def test_variable_midpoints_and_midnight(self):
        for seconds in (5,10,30):
            h={'EXPSTA':'2026-11-11T23:59:58','EXPOSURE':seconds}
            expected=Time(h['EXPSTA']).mjd+seconds/2/86400
            self.assertAlmostEqual(header_mid_mjd(h),expected,places=10)
        h={'EXPSTA':'2026-11-11T23:59:58','EXPEND':'2026-11-12T00:00:08','EXPOSURE':30}
        self.assertAlmostEqual(header_mid_mjd(h),Time('2026-11-12T00:00:03').mjd,places=10)

    def test_missing_and_invalid_time_do_not_assume_30s(self):
        self.assertIsNone(header_mid_mjd({'EXPSTA':'2026-11-11T00:00:00'}))
        self.assertIsNone(header_mid_mjd({'EXPSTA':'2026-11-11T00:00:10','EXPEND':'2026-11-11T00:00:05','EXPOSURE':30}))
        self.assertAlmostEqual(header_mid_mjd({'MJD':61000,'EXPTIME':10}),61000+5/86400,places=10)

    def test_original_mp_grid_not_equatorial_replacement(self):
        fields,meta=load_grid()
        self.assertEqual(len(fields),1855)
        self.assertEqual(fields[0]['field_id'],'MP0001')
        self.assertAlmostEqual(fields[0]['ra_deg'],7.590927324098)
        self.assertEqual(meta['original_fits_sha256'],'13f394877ce26e50e5cee5053a4ad845f412be4520421e46931b3e157d9e0225')

    def test_slew_wrap_and_zero(self):
        a={'field_id':'a','ra_deg':359,'dec_deg':20}
        b={'field_id':'b','ra_deg':1,'dec_deg':20}
        self.assertAlmostEqual(slew(a,b),9.11+11.95*np.sqrt(2))
        self.assertEqual(slew(a,a),0)

    def test_route_never_worse_than_same_set_ra_order(self):
        f=[{'field_id':str(i),'ra_deg':r,'dec_deg':d} for i,(r,d) in enumerate(
            [(359,1),(1,4),(2,0),(358,3),(3,2),(357,-1)])]
        path,baseline,best=route_order(list(range(6)),f)
        self.assertEqual(set(path),set(range(6)))
        self.assertLessEqual(best,baseline+1e-8)
        self.assertAlmostEqual(best,sum(slew(f[path[i-1]],f[path[i]]) for i in range(6)))

    def test_theory_is_monotone_and_dark_is_30(self):
        values=[exposure(h,30,40,0,180,-10,self.c) for h in (-10,-11,-12,-13,-14,-18,-20)]
        self.assertEqual(values,sorted(values))
        self.assertEqual(values[0],5)
        self.assertEqual(values[-2:],[30,30])

    def test_high_speed_actual_pairer_and_original_lower_bound(self):
        start=61000;end=start+8/1440
        one=Table({'RA':[10.],'DEC':[0.],'objID':[1],'MAG_PSF':[18.]})
        for speed in (2.,30.,150.,700.,800.):
            two=Table({'RA':[10+speed*(end-start)*24/3600],'DEC':[0.],
                       'objID':[2],'MAG_PSF':[18.]})
            result=pair_two_exposures(one,two,3,750,mjd1=start,mjd2=end)
            self.assertEqual(result is not None,3<=speed<=750)
            if result is not None:self.assertAlmostEqual(float(result['v_arcsec_hr'][0]),speed,places=6)
            old=pair_two_exposures(one,two,3,63,mjd1=start,mjd2=end)
            self.assertEqual(old is not None,3<=speed<=63)

    def test_timing_primary_hdu_with_extension_wcs(self):
        from astropy.io import fits
        primary=fits.PrimaryHDU()
        primary.header['EXPSTA']='2026-11-11T12:00:00'
        primary.header['EXPEND']='2026-11-11T12:00:05'
        extension=fits.BinTableHDU()
        for k,v in {'CRVAL1':20.,'CRVAL2':10.,'CRPIX1':1.,'CRPIX2':1.,
                    'CTYPE1':'RA---TAN','CTYPE2':'DEC--TAN',
                    'CD1_1':-.0003,'CD1_2':0.,'CD2_1':0.,'CD2_2':.0003}.items():
            extension.header[k]=v
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'test.fits';fits.HDUList([primary,extension]).writeto(path)
            *_,mid=image_center_and_radius(path)
            self.assertAlmostEqual(mid,Time('2026-11-11T12:00:02.5').mjd,places=10)

    def test_new_processing_entry_uses_configured_speed_only(self):
        from targeted_observing.offline_processing import command
        cmd=command('20261111','/tmp/input','/tmp/output',self.c,2)
        self.assertEqual(cmd[cmd.index('--vmin')+1],'3')
        self.assertEqual(cmd[cmd.index('--vmax')+1],'750')
        self.assertNotIn('submit',str(cmd))
        self.assertIn('targeted_observing/unknown_pipeline/',cmd[1])
        self.assertNotIn('/heliolincrr/',cmd[1])

    def test_static_export_has_real_R_and_no_waits(self):
        row={'field_id':'MP0001','ra_deg':7.590927324098,'dec_deg':-16.93414391,
             'exposure_s':5,'start_unix':1800000000,'mid_unix':1800000002.5,'end_unix':1800000005}
        block={'side':'EV','rows':[row],'window_unix':[1800000000,1800003600],'reason':'test'}
        plan={'date':'2026-11-11','config':self.c,'blocks':[block]}
        with tempfile.TemporaryDirectory() as tmp:
            write_plan(plan,tmp)
            tokens=(Path(tmp)/'EV.txt').read_text().strip().split('\t')
            self.assertEqual(len(tokens),8)
            self.assertEqual(tokens[0],'MP_0001')
            self.assertEqual(tokens[4],'R');self.assertEqual(tokens[6],'5')
            self.assertNotIn('WAIT',(Path(tmp)/'EV.txt').read_text())


if __name__=='__main__':unittest.main()
