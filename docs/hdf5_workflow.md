# HDF5-only analysis archives

All analysis archives previously stored as NPZ now use **real HDF5**, through
`coupling.array_store`. This includes prepared LFPs, power, phase coherence,
respiration, spike/transition null distributions, Julia statistical contrasts,
figure source values, and cross-session summaries. All consumers, demonstrations,
configuration examples and tests use the new format. This is a breaking storage
change: regenerate old outputs; renaming an old file does not convert its contents.

Existing single-array `.npy` caches and the Julia STFT `.pkl` export are retained.
They are not NPZ files. No analysis code loads or writes NPZ. The archive tool
still recognizes old NPZ files so it can move previous runs out of the way.

## What Abu's repository does

Verified against `abuzarmahmood/Mahmood_2023_BLA_GC` at
`9c54acb57e227dc9d13597750165b3ebc7be68e9`:

- [Updated phase aggregation](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/blob/9c54acb57e227dc9d13597750165b3ebc7be68e9/lfp_phase_coherence/updated/lfp_coherence_aggregate.py)
  reads session directories, calls `ephys_data(...).get_stft(..., dat_type=['phase'])`,
  and writes bootstrap/shuffle arrays into the session HDF5 through PyTables.
- [Original setup](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/blob/9c54acb57e227dc9d13597750165b3ebc7be68e9/lfp_phase_coherence/original/lfp_coherence_setup.py)
  uses a separate `AM_LFP_analyses.h5` file.

Neither phase workflow requires NPZ. Our derived analysis files are kept separate
from recording HDF5, and their schema marker distinguishes them from recordings.
Changing storage does not itself change the coherence statistic, STFT, bootstrap,
normalization, or animal-level inference. This is not a claim that all methods
in the repositories are identical.

## Existing workflows

Commands retain their arguments, with `.h5` replacing `.npz` in archive paths:
`cleaned_lfp.h5`, `session_power.h5`, `pressure_1khz.h5`,
`phase_coherence.h5`, `spectra.h5`, and transition-posterior `GC.h5`/`PC.h5`.
The new Julia power commands remain session-based. The older power interfaces
remain named `coupling-power-legacy` and `coupling-compare-power-legacy` but now
also read/write HDF5. Phase and spike commands can still read recording HDF5
files directly. No conversion of original recordings is necessary.

Historical local configs contain output paths from earlier runs. Merely replacing
an extension is insufficient: point them to newly generated output directories,
or use the Julia session backend below, which needs no prepared-LFP archive.
`configs/phase_coherence.h5.json` replaces the older archive-input example name.

## Julia-cleaned data → phase coherence

The new `input_backend: "julia"` path consumes the **same cleaned trials and
approved channels** used by the Julia power workflow. It reads the recording
HDF5, `cleaned_LFP_artifacts/*.npy`, and its artifact-removal report directly.
Original per-DIN trial indices are retained, so rejection does not renumber events.
It checks that cleaned values equal their indexed original trials and that regional
electrode assignments still agree. Old reports without trial IDs must be regenerated.

QC changes are explicit: all-rejected trials remain rejected; no approved channels
is an error. Power plots skip empty conditions; phase coherence reports an error
when a selected condition lacks jointly accepted trials or lacks the required two
approved channels per region for its intra-region control.

```bash
cd "$HOME/intercortical-coupling"
python -m pip install -e '.[lfp,phase,spikes]'
cp configs/phase_coherence.julia.example.json configs/phase_coherence.julia.local.json
```

Edit `defaults.preprocessing` in that local configuration to describe the actual
filtering, referencing and extraction used for the existing `/Parsed_LFP` arrays.
Do not infer those settings from the filename. Override per-session values in
`sessions` when acquisition/extraction differs. The example deliberately refuses
to run with the placeholder. Labels and regions come from `.info`, and epoch
onset follows the Julia timing helper; do not add old `input` overrides.

Preview and archive previous outputs, then run:

```bash
python scripts/archive_results.py --data-list configs/fin_inter_list.txt
python scripts/archive_results.py --data-list configs/fin_inter_list.txt --apply

coupling-lfp-pipeline \
  --data-list configs/fin_inter_list.txt \
  --output-dir "results/julia_hdf5_$(date +%Y%m%d_%H%M%S)" \
  --stages validate clean plot average phase \
  --phase-config configs/phase_coherence.julia.local.json \
  --phase-bootstrap 500 --seed 0 --jobs 1
```

Use `--stages validate clean plot average compare phase` to include the Julia
animal-level power comparisons too. An insufficient-animal comparison returns
nonzero and stops later stages; run it separately if descriptive and phase outputs
should complete regardless. Phase output appears in the run's `phase/` directory.
For a preflight after cleaning:

```bash
coupling-phase-figure configs/fin_inter_list.txt \
  configs/phase_coherence.julia.local.json results/phase_check --check
```

The runner retains its default stages unless explicitly supplied; phase is opt-in.

## Storage contract and safety

Derived stores have root attributes `schema="intercortical-coupling.arrays"` and
`schema_version=1`, with named root datasets. Use `save_arrays(path, **arrays)` and
`with load_arrays(path) as arrays:`. Numeric/complex/boolean arrays and Unicode
labels retain their shapes and dtypes. JSON metadata stays a string dataset.
Objects and pickles are rejected; writing is atomic and cannot overwrite an
unrecognized recording HDF5. HDF5 strings can be inspected with ordinary h5py or
PyTables. The archive script permits schema-marked derived HDF5 while continuing
to reject recording or unrecognized HDF5 files inside selected output folders.

The storage migration does not unify every preprocessing algorithm. Respiration
realignment still needs continuous LFP and ADC caches; transition correlations
still need independently fitted transition posteriors. Those stages now use
HDF5 archives, but their scientific prerequisites remain unchanged. Original
recordings and pre-existing results are not modified automatically.
