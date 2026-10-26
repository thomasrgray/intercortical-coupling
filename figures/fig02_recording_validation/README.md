# Figure 2 — simultaneous spikes, alpha LFP, and respiration

Implemented: one descriptive figure per stimulus, from the **same accepted events**:

- **A:** sorted single-unit GC and PC rasters for the first accepted example trial.
- **B:** GC and PC 8–12 Hz LFP from that same trial, one accepted electrode per region.
- **C:** filtered nasal pressure for the same trial, with a 200 ms scale bar.
- **D:** firing rate per unit, averaged over units and matched trials, in 50 ms bins.

Histological placement and waveform/drift-quality panels remain separate inputs;
these plots do not establish anatomical placement or sorting quality by themselves.
Panel letters are provisional. No synthetic output is study data.

## Run on a prepared session

Use `prepared` and `respiration_audit` from the **same** run of the
[respiration_sniffing workflow](../../respiration_sniffing/README.md). That workflow
creates these caches in steps 2 and 3. The standalone sniffing step's CSV summaries
are not pressure waveforms; Figure 2 reads the audited continuous pressure cache.

```bash
cd "$HOME/intercortical-coupling"
python -m pip install -e '.[spikes]'

# Replace this with the completed run directory, not the raw session directory.
RUN="$PWD/results/YOUR_COMPLETED_RESPIRATION_LFP_RUN"
python -m coupling.figure2 \
  "$RUN/prepared" "$RUN/respiration_audit" \
  "$RUN/fig02_recording_validation" --check
python -m coupling.figure2 \
  "$RUN/prepared" "$RUN/respiration_audit" \
  "$RUN/fig02_recording_validation"
```

`coupling-figure2` is the equivalent installed command. `--check` validates and
reports retained trial counts without writing plots. Existing output directories
are never overwritten. If multiple HDF5 or `.info` files exist, pass their exact
session-local filenames with `--h5 session.h5 --info session.info`.

## Signal sources and matching

| Signal | Source | Interpretation |
| --- | --- | --- |
| Sorted spikes | Session `*.h5`/`*.hdf5`, `/sorted_units/unitNNN/times` | Absolute raw acquisition sample indices |
| Unit assignment | HDF5 `/unit_descriptor` and session `*.info` → `electrode_layout.gc`/`.pc` | Descriptor row matches numeric unit ID; retain `single_unit` flags |
| LFP | `prepared/cleaned_lfp.npz` | Existing filtering, within-region reference and channel/trial QC |
| Pressure | `respiration_audit/pressure_1khz.npz` | Existing 0.5–20 Hz filtered counts |
| Event clock | `prepared/events.csv`, audited respiration events and provenance | DIN identity, onset and rounded 1 kHz epoch center |

Verified against upstream [blech_make_arrays.py](https://github.com/katzlabbrandeis/blech_clust/blob/master/blech_make_arrays.py)
and [ephys_data.py](https://github.com/katzlabbrandeis/blech_clust/blob/master/utils/ephys_data/ephys_data.py).
`blech_clust` also stores `/spike_trains/dig_in_N/spike_array` as trial × unit ×
millisecond. Figure 2 deliberately uses the original sorted times: prebuilt arrays
can have different epoch lengths and omitted cutoff trials. No positional trial
join is needed. Absolute times are binned into 1 ms **counts**, preserving multiple
spikes per bin; zero is the same rounded center used for LFP and pressure.

The source session is inferred from prepared provenance. ADC/LFP source, clock,
event identity, cache provenance and source-file freshness use the existing
spectra validator. Ambiguous files, inconsistent unit IDs, overlapping regions,
or no eligible single units in either region are errors. Unit IDs are descriptor
row indices, not electrode IDs. HDF5/region file hashes are recorded.

For each matched odor, use the same accepted LFP electrodes across routes.
Select each region's electrode by maximizing the smaller of the two route-specific
valid trial counts, then their sum (ties: archive order). Other stimuli use
maximum valid trial coverage. Intersect GC/PC/pressure accepted trials per stimulus.
Plot the first shared event. Selection is based on QC coverage, never alpha
amplitude or apparent coupling. Export channel choices and counts. Regional CAR
waveforms are **not averaged**, because that can cancel the signal.

Alpha filtering is fourth-order Butterworth 8–12 Hz, applied forward/backward to
full −2 to +5 s epochs; display only −0.5 to +2.5 s to reduce boundary artifacts.
This offline filter is acausal: do not use it to infer response onset or direction.
Pressure polarity is retained in provenance, not shown in titles; unknown polarity does not label inhalation.
DIN zero is a command timestamp, not measured odor arrival.

Firing rates are descriptive equal-unit averages on the same retained trials.
Laser conditions are not separated in this recording-validation view; do not use
these pooled summaries as optogenetic contrasts. Single-unit flags alone do not
establish waveform isolation or stability. No animal-level inference is performed.
An electrode may have valid LFP without a sorted unit, and the illustrated LFP
electrode need not carry the displayed units.

## Outputs

`01_figure2.png/pdf`, etc. follow the stimulus order in the prepared archive.
`retention.csv` maps each stimulus to accepted counts and selected electrodes;
stimuli with no shared trials are recorded but have no figure.
`units.csv` lists retained unit IDs and electrode/region assignments.
`pressure_trial_qc.csv` preserves pressure exclusions.
`plotted_data.npz` contains the time axes, event IDs, pressure, alpha traces,
spike counts, and trial-level population rates. `provenance.json` records source
hashes, inherited processing, selection/filtering rules and limitations.

Synthetic tests verify exact spike timing/counts, region/identity failures, alpha
passband rejection, cache mismatch rejection and complete PNG/PDF/data export.
Biological validation still requires running on the local recordings and inspecting
the figure, sorting QC, assignments and retention tables.

### Region names in older session metadata

The default loader accepts `gc`/`pc`, `GC`/`PC`, and mixed-case variants.
Numeric CAR groups are not interpreted as anatomical regions. If `.info` lacks
anatomical keys, add `--unit-regions-from-prepared` to both check and generation
commands to use the electrode-ID/region mapping already used for LFP processing.
This requires HDF5 `unit_descriptor.electrode_number` to use the same electrode
numbering as the prepared archive. No channel-index or unit-index substitution
is performed. All mapped electrodes are used for unit assignment, including
electrodes excluded by LFP QC; unit eligibility still uses `single_unit`. The
chosen mapping is recorded in provenance. Existing prepared/respiration caches
can be reused without rerunning preprocessing.

### Orthonasal–retronasal comparison

`BA_ortho_vs_retro.png/pdf` compares oBA (left) with BA (right);
`IA_ortho_vs_retro.png/pdf` compares oIA with IA. Columns share time and row-specific
y-axis scales, unit ordering and LFP electrodes. Each column shows its own first
accepted trial for rasters/traces, followed by its matched-trial mean firing rate.
These are descriptive examples, not simultaneous trials across routes or evidence
of a significant route difference. Missing conditions are explicitly labeled.
The existing per-stimulus figures also put rasters first.


## Preferred version restored

The active code is the blue/orange rasters-first version from commit `e27e8c2`,
without background overlays or envelopes. Both per-stimulus and ortho/retro
comparison figures remain available. Later experiments are preserved in
[the Figure 2 archive](../../archive/figure2/README.md), along with instructions
for archiving generated images and associated results on the workstation.
