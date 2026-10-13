import csv
import json
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from coupling.batch import discover,preflight,run_batch,load_manifest
from coupling.power import STIMULI


def fixture(root,name='session'):
    root=root/name;root.mkdir()
    fs=1000;n=15000
    (root/'info.rhd').write_bytes(struct.pack('<Ihhf',0xc6912702,1,0,fs))
    np.arange(n,dtype='<i4').tofile(root/'time.dat')
    rng=np.random.default_rng(13)
    for c in range(10):
        rng.integers(-100,100,n,dtype=np.int16).astype('<i2').tofile(root/f'amp-C-{c:03d}.dat')
    for din in range(8):
        x=np.zeros(n,dtype='<u2')
        for k in range(3 if din<6 else 4): x[3000+din*10+k*2500:3005+din*10+k*2500]=1
        x.tofile(root/f'board-DIN-{din:02d}.dat')
    info=dict(file_type='one file per channel',electrode_layout={'gc':[[0,1,2,3,4]],'pc':[[5,6,7,8,9]]},
        taste_params=dict(dig_in_nums=list(range(6)),tastes=list(STIMULI),trial_count=[3]*6),
        dig_ins=dict(nums=list(range(8)),trial_counts=[3]*6+[4,4]),recording_params=dict(sampling_rate=fs),laser_params=dict(dig_in_nums=[]))
    (root/'session.info').write_text(json.dumps(info));(root/'session.info.bak').write_text('{}')
    (root/'session.params').write_text(json.dumps(dict(sampling_rate=fs,bandpass_lower_cutoff=300,bandpass_upper_cutoff=3000)))
    with (root/'session_electrode_layout.csv').open('w') as h:
        w=csv.writer(h);w.writerow(['filename','electrode_ind','CAR_group','original_CAR_group'])
        for c in range(10):w.writerow([f'amp-C-{c:03d}.dat',c,'gc' if c<5 else 'pc','gc' if c<5 else 'pc'])
    return dict(path=str(root),session_id=name,animal_id='TG22')


class BatchTests(unittest.TestCase):
    def test_preflight_mapping_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            s=fixture(Path(tmp));r=preflight(s,{})
            self.assertEqual(r['errors'],[])
            self.assertEqual(r['facts']['channels_by_region'],dict(GC=5,PC=5))
            self.assertEqual(r['config']['channels'][5]['region'],'PC')
            self.assertEqual(r['config']['stimulus_by_din']['4'],'oIA')
            self.assertEqual(r['facts']['ignored_spike_bandpass_hz'],[300,3000])
            self.assertTrue(preflight(s,{'misspelled_setting':1})['errors'])

    def test_discovery_backup_ignored_and_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);animal=p/'TG22';animal.mkdir();fixture(animal)
            sessions=discover(p,p/'manifest.json')
            self.assertEqual(len(sessions),1);self.assertEqual(sessions[0]['animal_id'],'TG22')
            with self.assertRaises(FileExistsError):discover(p,p/'manifest.json')

    def test_missing_files_and_rate_conflict(self):
        with tempfile.TemporaryDirectory() as tmp:
            s=fixture(Path(tmp));p=Path(s['path'])
            (p/'amp-C-000.dat').unlink();(p/'amp-C-005.dat').unlink()
            (p/'session.params').write_text('{"sampling_rate":30000}')
            r=preflight(s,{})
            self.assertIsNone(r['config'])
            self.assertTrue(any('amp-C-000.dat' in x and 'amp-C-005.dat' in x for x in r['errors']))
            self.assertTrue(any('disagrees with RHD' in x for x in r['errors']))

    def test_mapping_disagreement_and_opto_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            s=fixture(Path(tmp));p=Path(s['path'])/'session.info';v=json.loads(p.read_text())
            v['electrode_layout']={'gc':[[0,1,3]],'pc':[[2,4,5]]};v['laser_params']['dig_in_nums']=[7];p.write_text(json.dumps(v))
            r=preflight(s,{})
            self.assertTrue(any('region disagreement' in x for x in r['errors']))
            self.assertTrue(any('Optogenetic' in x for x in r['errors']))

    def test_timestamp_gap_and_trial_count_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            s=fixture(Path(tmp));p=Path(s['path'])
            t=np.arange(15000,dtype='<i4');t[500:]+=1;t.tofile(p/'time.dat')
            v=json.loads((p/'session.info').read_text());v['taste_params']['trial_count'][0]=30
            (p/'session.info').write_text(json.dumps(v));r=preflight(s,{})
            self.assertTrue(any('Timestamp gaps' in x for x in r['errors']))
            self.assertTrue(any('actual trial count' in x for x in r['errors']))

    def test_check_is_read_only_and_isolates_sessions(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);a=fixture(p,'a');b=fixture(p,'b')
            (Path(b['path'])/'info.rhd').unlink()
            manifest=p/'batch.json';manifest.write_text(json.dumps(dict(version=1,sessions=[a,b])))
            before={str(f):(f.stat().st_size,f.stat().st_mtime_ns) for s in [a,b] for f in Path(s['path']).iterdir()}
            rows=run_batch(manifest,p/'report')
            self.assertEqual([r['status'] for r in rows],['ready','blocked'])
            for f,stamp in before.items():self.assertEqual((Path(f).stat().st_size,Path(f).stat().st_mtime_ns),stamp)
            with self.assertRaises(FileExistsError):run_batch(manifest,p/'report')
            with self.assertRaises(ValueError):run_batch(manifest,Path(a['path'])/'out')

    def test_duplicate_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);s=fixture(p);m=p/'m.json';m.write_text(json.dumps(dict(version=1,sessions=[s,s])))
            with self.assertRaises(ValueError):load_manifest(m)

    def test_processing_failure_does_not_stop_next_session(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);sessions=[fixture(p,'a'),fixture(p,'b')];m=p/'m.json'
            m.write_text(json.dumps(dict(version=1,sessions=sessions)))
            def fake_prepare(source,config,out):
                if Path(source).name=='a': raise ValueError('intentional test failure')
                out.mkdir();return out/'cleaned_lfp.npz'
            def fake_analyze(source,out,**kw):
                out.mkdir();(out/'retention.csv').write_text('accepted_trials\n3\n')
                (out/'granger_diagnostics.csv').write_text('passed\nTrue\nFalse\n')
            with patch('coupling.batch.prepare',side_effect=fake_prepare),patch('coupling.lfp_coupling.analyze',side_effect=fake_analyze):
                rows=run_batch(m,p/'out',execute=True)
            self.assertEqual([r['status'] for r in rows],['failed','completed'])
            self.assertEqual(rows[1]['granger_fits_passed'],1)

if __name__=='__main__':unittest.main()
