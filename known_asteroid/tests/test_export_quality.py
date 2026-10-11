"""Known-report quality and local schema gates; never send network requests."""
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from astropy.table import Table, MaskedColumn

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import export_ades as exporter


class KnownExportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.input = self.root / 'matched.fits'
        self.out = self.root / 'report.psv'
        env = patch.dict(os.environ, {'MPC_ADES_VALIDATOR_PYTHON':os.environ.get('MPC_ADES_VALIDATOR_PYTHON',sys.executable),
                                     'MPC_SUBMISSION_HOLD_FILE':str(self.root / 'HOLD')})
        env.start();self.addCleanup(env.stop)
        self.args = exporter.build_argparser().parse_args(['--fits',str(self.input),'--out',str(self.out)])

    def table(self, n=4):
        return Table({'number':[12345.] * n, 'name':['Test'] * n,
                      'epoch':61318 + np.arange(n) / 144.,
                      'RA_Win':120 + np.arange(n) * .001,'DEC_Win':[10.] * n,
                      'RAErr_Win':[1e-5] * n,'DECErr_Win':[1e-5] * n,
                      'Mag_Aper4':[18.] * n,'MagErr_Aper4':[.05] * n,
                      'Flag':[0] * n,'Flux_Aper4':[1000.] * n,'FluxErr_Aper4':[50.] * n,
                      'source_file':[f'frame{i}.fits' for i in range(n)],'objID':np.arange(n)})

    def export(self, table):
        table.write(self.input, overwrite=True)
        return exporter.export_psv(self.args)

    def audit(self):
        return json.loads(Path(str(self.out)+'.audit.json').read_text())

    def test_valid_report_preserves_uncertainties(self):
        self.assertEqual(self.export(self.table(3)),3)
        self.assertIn('|0.0360|0.0360|',self.out.read_text())
        self.assertEqual(self.audit()['validation']['status'],'local_validated')
        self.assertEqual(exporter.validate_payload(self.out.read_text())['observations'],3)

    def test_incident_roundoff_and_bad_flag_are_excluded_without_floor(self):
        table=self.table();table['RAErr_Win'][-1]=9.11565933847669e-9
        table['DECErr_Win'][-1]=8.486015445896555e-9
        table['Flag'][-1]=4;table['Mag_Aper4'][-1]=99;table['Flux_Aper4'][-1]=-1
        self.assertEqual(self.export(table),3)
        self.assertNotIn('|0.0000|',self.out.read_text())
        with Path(str(self.out)+'.excluded.csv').open() as f:rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),1)
        self.assertIn('nonzero_or_missing_flag',rows[0]['reason'])
        self.assertIn('uncertainty_rounds_to_zero',rows[0]['reason'])
        self.assertIn('invalid_photometry',rows[0]['reason'])

    def test_flag_zero_does_not_bypass_error_and_photometry_gates(self):
        cases=[('RAErr_Win',1e-9),('RAErr_Win',0),('DECErr_Win',-1e-5),
               ('RAErr_Win',float('nan')),('MagErr_Aper4',.0001),
               ('Flux_Aper4',0),('FluxErr_Aper4',-1),('Mag_Aper4',99)]
        for col,value in cases:
            with self.subTest(col=col,value=value):
                table=self.table();table[col][-1]=value
                self.assertEqual(self.export(table),3)

    def test_masked_uncertainty_not_treated_as_valid_storage_value(self):
        table=self.table();table['RAErr_Win']=MaskedColumn(table['RAErr_Win'],mask=[0,0,0,1])
        self.assertEqual(self.export(table),3)

    def test_minimum_observations_after_quality_cuts(self):
        table=self.table(3);table['Flag'][-1]=4
        self.assertEqual(self.export(table),0)
        self.assertEqual(self.audit()['objects_below_minimum_after_filtering'],1)

    def test_minimum_observations_after_dedup(self):
        table=self.table(3);table['epoch'][2]=table['epoch'][1]
        self.assertEqual(self.export(table),0)
        self.assertEqual(self.audit()['duplicate_rows'],1)

    def test_quality_filter_precedes_duplicate_ranking(self):
        table=self.table(4);table['epoch'][3]=table['epoch'][2]
        table['Flag'][3]=4;table['RAErr_Win'][3]=table['DECErr_Win'][3]=1e-9
        self.assertEqual(self.export(table),3)
        self.assertNotIn(f'{table["RA_Win"][3]:.10f}',self.out.read_text())

    def test_minimum_observations_after_history_filter(self):
        table=self.table(3)
        self.args.dedup_history_root=str(self.root / 'history')
        epoch=exporter.mjd_to_obs_time_z(np.array(table['epoch']))[0]
        with patch.object(exporter,'load_submitted_history_keys',return_value={('12345',epoch)}):
            self.assertEqual(self.export(table),0)
        self.assertEqual(self.audit()['history_removed'],1)

    def test_bad_schema_does_not_replace_existing_report(self):
        self.out.write_text('original evidence')
        self.args.obs_code='INVALID'
        with self.assertRaisesRegex(ValueError,'Local ADES validation failed'):
            self.export(self.table(3))
        self.assertEqual(self.out.read_text(),'original evidence')

    def test_schema_rejects_serialized_zero_negative_and_bad_coordinates(self):
        self.export(self.table(3));text=self.out.read_text()
        for bad in [text.replace('|0.0360|','|0.0000|'),text.replace('|0.0360|','|-0.0010|'),
                    text.replace('120.0000000000','400.0000000000'),text.replace('|0.050|','|0.000|')]:
            with self.assertRaisesRegex(ValueError,'schema validation failed'):
                exporter.validate_payload(bad)

    def test_validate_cli_is_local_and_hold_still_blocks_formal_transport(self):
        self.table(3).write(self.input)
        argv=['export_ades.py','--fits',str(self.input),'--out',str(self.out),'--validate']
        with patch.object(sys,'argv',argv),patch.object(exporter,'run_curl') as send:
            exporter.main();send.assert_not_called()
        (self.root/'HOLD').touch()
        with patch.object(exporter.subprocess,'run') as run:
            with self.assertRaisesRegex(RuntimeError,'submission paused'):
                exporter.run_curl(exporter.SUBMIT_URL,self.args)
            run.assert_not_called()

    def test_missing_quality_column_blocks_export(self):
        table=self.table(3);table.remove_column('Flag')
        with self.assertRaisesRegex(KeyError,'Flag'):
            self.export(table)

    def test_formal_transport_uses_validated_snapshot_when_file_changes(self):
        self.out.write_text('frozen PSV')
        def validate(payload):
            self.assertEqual(payload,'frozen PSV')
            self.out.write_text('changed after validation')
            return {'status':'local_validated'}
        reply=SimpleNamespace(returncode=0,stdout='Submission ID example',stderr='')
        with patch.object(exporter,'validate_payload',side_effect=validate),patch.object(exporter.subprocess,'run',return_value=reply) as run:
            exporter.run_curl(exporter.SUBMIT_URL,self.args)
        self.assertEqual(run.call_args.kwargs['input'],'frozen PSV')
        self.assertIn('source=<-',run.call_args.args[0])

    def test_invalid_snapshot_blocks_formal_transport(self):
        self.out.write_text('invalid PSV')
        with patch.object(exporter,'validate_payload',side_effect=ValueError('invalid schema')),patch.object(exporter.subprocess,'run') as run:
            with self.assertRaisesRegex(ValueError,'invalid schema'):
                exporter.run_curl(exporter.SUBMIT_URL,self.args)
            run.assert_not_called()


if __name__=='__main__':unittest.main()
