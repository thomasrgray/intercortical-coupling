import json
from types import SimpleNamespace
import h5py
import numpy as np
import pytest
from coupling.phase_coherence import load_input
from coupling.julia_lfp import remove_lfp_artifacts as clean


def make_session(root):
    rng = np.random.default_rng(0)
    (root/'a.info').write_text(json.dumps({'file_type': 'one file per channel',
        'electrode_layout': {'GC': [[0, 1]], 'PC': [[2, 3]]},
        'taste_params': {'tastes': ['waterIA', 'IA'], 'dig_in_nums': [0, 1]}}))
    (root/'a.params').write_text(json.dumps({'spike_array_durations': [2000, 5000]}))
    folder = root/'cleaned_LFP_artifacts'
    folder.mkdir()
    rows = []
    with h5py.File(root/'a.h5', 'w') as h:
        h.create_dataset('Parsed_LFP_channels', data=np.arange(4))
        for din in range(2):
            data = rng.normal(size=(4, 4, 7000))
            h.create_dataset(f'Parsed_LFP/dig_in_{din}_LFPs', data=data)
            original = data.transpose(1, 0, 2)
            kept = original[[0, 2, 3]]
            np.save(folder/f'taste_{din:02d}_cleaned.npy', kept)
            rows.append(dict(taste_index=din, retained_trial_indices=[0, 2, 3],
                             input_shape=list(original.shape), output_shape=list(kept.shape)))
    report = dict(approved_lfp_positions_by_region={'GC': [0, 1], 'PC': [2, 3]},
                  excluded_taste_indices=[], tastes=rows)
    (folder/'artifact_removal_report.json').write_text(json.dumps(report))
    return folder, report


def test_julia_adapter_keeps_original_ids_qc_and_routes(tmp_path):
    folder, report = make_session(tmp_path)
    spec = {'input_backend': 'julia', 'preprocessing': 'synthetic unreferenced 1 kHz epochs'}
    dat = load_input(tmp_path, spec)
    assert dat['counts'] == {'IA': 3, 'oIA': 3}
    assert dat['good'].tolist() == [True, False, True, True]*2
    assert dat['event_id'][2] == 'dig_in_0_LFPs:trial2'
    assert dat['time_s'][0] == -2
    # Reject changed trial identity, even when the shapes and counts still match.
    report['tastes'][0]['retained_trial_indices'] = [1, 2, 3]
    (folder/'artifact_removal_report.json').write_text(json.dumps(report))
    with pytest.raises(ValueError, match='no longer match'):
        load_input(tmp_path, spec)


def test_missing_trial_ids_requires_reclean(tmp_path):
    folder, report = make_session(tmp_path)
    del report['tastes'][0]['retained_trial_indices']
    (folder/'artifact_removal_report.json').write_text(json.dumps(report))
    with pytest.raises(ValueError, match='Original trial IDs'):
        load_input(tmp_path, {'input_backend':'julia', 'preprocessing':'synthetic'})


def test_all_rejected_trials_stay_rejected(tmp_path, monkeypatch):
    dat = SimpleNamespace(get_lfps=lambda: None, lfp_array_list=[np.zeros((3, 2, 10))])
    monkeypatch.setattr(clean.ephys_data, 'ephys_data', lambda _: dat)
    monkeypatch.setattr(clean, 'get_taste_names_in_lfp_order', lambda _: ['BA'])
    monkeypatch.setattr(clean, 'validate_dataset', lambda *a, **kw: {'good_lfp_positions_by_region': {'GC':[0], 'PC':[1]}})
    monkeypatch.setattr(clean.lfp_processing, 'return_good_lfp_trial_inds', lambda *a, **kw: np.zeros(3, bool))
    clean.remove_dataset_artifacts(str(tmp_path), 3)
    report = json.loads((tmp_path/'cleaned_LFP_artifacts/artifact_removal_report.json').read_text())
    assert report['tastes'][0]['retained_trial_indices'] == []
    assert report['tastes'][0]['retained_trials'] == 0
