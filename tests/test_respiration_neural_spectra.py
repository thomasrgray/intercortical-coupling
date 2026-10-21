"""Synthetic integration tests for respiration–LFP–sorted-spike spectral workflow."""
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

ENTRY = Path(__file__).resolve().parents[1] / "respiration_sniffing" / "03_neural_spectra.py"
spec = importlib.util.spec_from_file_location("neural_spectra", ENTRY)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def make_fixture(root):
    session=root/"TG22_Feb2"
    session.mkdir()
    labels=["BA","S","oIA","IA","oBA","BK"]*2
    dins={"BA":0,"S":1,"oIA":2,"IA":3,"oBA":4,"BK":5}
    times=np.arange(54000)/1000
    raw=(32768+400*np.sin(2*np.pi*5*times)).astype("<u2")
    adc=session/"board-ADC-00.dat"
    raw.tofile(adc)
    trials=session/"trial_info_frame.csv"
    counts={}
    keys=[]
    with trials.open("w",newline="") as handle:
        w=csv.DictWriter(handle,fieldnames=["start_taste","end_taste","taste",
                                            "dig_in_num_taste","abs_trial_num"])
        w.writeheader()
        for i,stim in enumerate(labels):
            din=dins[stim]
            k=counts.get(din,0)
            keys.append(f"DIN{din}:trial{k}")
            counts[din]=k+1
            start=2500+i*4000
            w.writerow(dict(start_taste=start,end_taste=start+500,
                            taste=stim,dig_in_num_taste=din,abs_trial_num=i))
    t=np.arange(-2000,5000)/1000
    rng=np.random.default_rng(12)
    lfp=np.stack([np.stack([25*np.sin(2*np.pi*(5+ch)*t)
                          +rng.normal(0,2,len(t)) for _ in labels])
                  for ch in range(2)]).astype("float32")
    provenance=dict(source=str(session),reference="none (synthetic)",
                    filtering=dict(fs_hz=1000),qc_method="synthetic accepted",
                    event_alignment="TTL first-high sample")
    archive=root/"cleaned_lfp.npz"
    np.savez_compressed(archive,lfp_uv=lfp,time_s=t,stimulus=np.array(labels),
                        event_id=np.arange(len(labels)),electrode_id=np.array([10,40]),
                        region=np.array(["GC","PC"]),good_channel=np.ones(2,bool),
                        good_trial=np.ones((2,len(labels)),bool),
                        session_id="TG22_Feb2",animal_id="TG22",
                        provenance_json=json.dumps(provenance))
    spike_t=t
    spikes=rng.binomial(1, .03, size=(len(labels),2,len(spike_t))).astype("uint8")
    spike_file=session/"sorted_spikes.npz"
    np.savez(spike_file,spikes=spikes,time_s=spike_t,
             unit_id=np.array([0,1]),region=np.array(["GC","PC"]),
             good_unit=np.array([True,True]),trial_id=np.asarray(keys),
             stimulus=np.asarray(labels),stratum=np.array(["laser:off"]*len(labels)),
             trial_order=np.arange(len(labels)),good_trial=np.ones(len(labels),bool))
    config=root/"spike_config.json"
    config.write_text(json.dumps({"defaults":{}}))
    return adc,trials,archive,spike_file,config


class NeuralSpectraTests(unittest.TestCase):
    def test_spectral_recovery(self):
        t=np.arange(3000)/1000
        f,p=module.welch_epoch(np.sin(2*np.pi*5*t))
        band=(f>=1)&(f<=12)
        self.assertEqual(float(f[band][np.argmax(p[band])]),5.)
        norm=module.normalized_psd(f,p)
        self.assertAlmostEqual(float(np.trapz(norm[band],f[band])),1.,places=4)

    def test_trial_id_validation(self):
        trials=[dict(trial="0",din=4),dict(trial="1",din=0),dict(trial="2",din=4)]
        self.assertEqual(module.validate_trial_table(trials),
                         ["DIN4:trial0","DIN0:trial0","DIN4:trial1"])
        with self.assertRaises(ValueError):
            module.validate_trial_table(trials[::-1])

    def test_end_to_end_with_and_without_sorted_spikes(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            adc,trials,lfp,spikes,config=make_fixture(root)
            rows=module.run(adc,trials,lfp,root/"with_spikes",
                            spikes=spikes,spike_config=config,raw_fs=1000)
            self.assertEqual(len(rows),10)
            self.assertEqual({r["kind"] for r in rows},
                             {"respiration","lfp","firing_rate"})
            self.assertGreater((root/"with_spikes"/"respiration_neural_spectra.pdf").stat().st_size,1000)
            self.assertGreater((root/"with_spikes"/"spectra.csv").stat().st_size,500)
            second=module.run(adc,trials,lfp,root/"no_spikes",raw_fs=1000)
            self.assertEqual({r["kind"] for r in second},{"respiration","lfp"})
            with self.assertRaises(FileExistsError):
                module.run(adc,trials,lfp,root/"no_spikes",raw_fs=1000)

    def test_mismatched_lfp_session_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            adc,trials,lfp,_,_=make_fixture(root)
            with np.load(lfp,allow_pickle=False) as d:
                payload={key:d[key] for key in d.files}
            labels=payload["stimulus"].copy()
            labels[0]="S"
            payload["stimulus"]=labels
            broken=root/"bad_lfp.npz"
            np.savez_compressed(broken,**payload)
            with self.assertRaisesRegex(ValueError,"stimulus order"):
                module.run(adc,trials,broken,root/"bad_results",raw_fs=1000)


if __name__=="__main__":
    unittest.main()
