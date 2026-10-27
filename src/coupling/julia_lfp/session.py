"""Load Julia-cleaned trials for phase coherence without intermediate archives."""
import json
from pathlib import Path
import numpy as np
import tables
from .lfp_common import ephys_data, get_lfp_epoch_timing, get_taste_names_in_lfp_order, CLEANED_LFP_FOLDER
from .lfp_session_labels import canonicalize_stimulus_names


def load_phase_session(session, spec):
    session = Path(session)
    declaration = spec.get('preprocessing', '')
    if not declaration or 'REPLACE_ME' in str(declaration):
        raise ValueError('Declare the actual filtering, reference and extraction in preprocessing')
    dat = ephys_data.ephys_data(str(session))
    dat.get_lfps()
    dat.get_lfp_channels()
    timing = get_lfp_epoch_timing(dat)
    names = canonicalize_stimulus_names(get_taste_names_in_lfp_order(dat))
    folder = session/CLEANED_LFP_FOLDER
    report_path = folder/'artifact_removal_report.json'
    report = json.loads(report_path.read_text())
    reports = {row['taste_index']: row for row in report['tastes']}
    if len(reports) != len(dat.lfp_array_list):
        raise ValueError('Artifact report does not match session; rerun cleaning')
    groups = report.get('approved_lfp_positions_by_region')
    if groups is None:
        raise ValueError('Artifact report lacks channel QC identities; rerun coupling-lfp-clean')
    regions = np.full(len(dat.parsed_lfp_channels), 'excluded', dtype='U8')
    for label, positions in groups.items():
        canonical = label.upper()
        if canonical not in {'GC', 'PC'}:
            continue
        for position in positions:
            if not 0 <= position < len(regions) or regions[position] != 'excluded':
                raise ValueError('Invalid/overlapping QC channel positions')
            regions[position] = canonical
    # Verify the reported regional identities against current electrode metadata.
    dat.get_region_electrodes()
    from ..phase_coherence import flatten_ids
    for label in ('GC', 'PC'):
        ids = {electrode for key, groups_ in dat.region_electrode_dict.items()
               if key.upper() == label for electrode in flatten_ids(groups_)}
        if any(int(dat.parsed_lfp_channels[i]) not in ids for i in np.flatnonzero(regions == label)):
            raise ValueError('Region metadata changed since cleaning; rerun QC')
    with tables.open_file(dat.hdf5_path, 'r') as h:
        nodes = [n._v_name for n in h.list_nodes('/Parsed_LFP')
                 if n._v_name.endswith('_LFPs') or n._v_name.startswith('DIN')]
    arrays, masks, labels, event_ids = [], [], [], []
    for index, original in enumerate(dat.lfp_array_list):
        if index in report.get('excluded_taste_indices', []):
            continue
        row = reports[index]
        if 'retained_trial_indices' not in row:
            raise ValueError('Original trial IDs are missing; rerun coupling-lfp-clean')
        retained = np.asarray(row['retained_trial_indices'], dtype=int)
        if retained.ndim != 1 or len(np.unique(retained)) != len(retained) or np.any(retained < 0) or np.any(retained >= len(original)):
            raise ValueError('Invalid retained trial indices')
        cleaned = np.load(folder/f'taste_{index:02d}_cleaned.npy', allow_pickle=False)
        if list(original.shape) != row['input_shape'] or list(cleaned.shape) != row['output_shape']:
            raise ValueError('Stale cleaned trial shapes; rerun cleaning')
        if cleaned.shape != original[retained].shape or not np.array_equal(cleaned, original[retained], equal_nan=True):
            raise ValueError('Cleaned trials no longer match original trial IDs; rerun cleaning')
        arrays.append(original.transpose(1, 0, 2))
        mask = np.zeros(len(original), dtype=bool)
        mask[retained] = True
        masks.append(np.broadcast_to(mask, (len(regions), len(mask))))
        labels.extend([names[index]] * len(original))
        event_ids.extend([f'{nodes[index]}:trial{i}' for i in range(len(original))])
    if not arrays:
        raise ValueError('No selected Julia LFP conditions')
    x = np.concatenate(arrays, axis=1)
    return dict(x=x, time_s=np.arange(x.shape[-1])/timing['fs']-timing['pre_stim_s'],
                electrode_id=np.asarray(dat.parsed_lfp_channels), region=regions,
                stimulus=np.asarray(labels), event_id=np.asarray(event_ids),
                good_channel=np.isin(regions, ['GC', 'PC']), good_trial=np.concatenate(masks, axis=1),
                provenance={'declared_preprocessing': declaration, 'input_backend': 'julia',
                            'event_alignment': spec.get('event_alignment', 'declared parsed-trial event'),
                            'artifact_report': str(report_path), 'artifact_report_contents': report,
                            'timing': timing, 'trial_identity_check': 'cleaned values equal original indexed trials'})
