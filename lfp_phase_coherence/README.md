# LFP phase coherence and coupling analyses

This workflow directory follows the analysis-oriented organization of
[Mahmood_2023_BLA_GC](https://github.com/abuzarmahmood/Mahmood_2023_BLA_GC/tree/9c54acb57e227dc9d13597750165b3ebc7be68e9/lfp_phase_coherence).
Shared numerical implementations remain in `src/coupling` so there is only one
implementation to test. [Methods and profiles](methods.md) describe the existing
repository/paper differences. Power has its own [workflow](../lfp_power/README.md).

## Run from the repository root

```bash
python -m pip install -e '.[phase]'
RUN_STAMP=$(date +%Y%m%d_%H%M%S)
coupling-phase-figure configs/fin_inter_list.phase_ready.txt \
  configs/phase_coherence.routes.local.json "results/phase_check_${RUN_STAMP}" --check && \
coupling-phase-figure configs/fin_inter_list.phase_ready.txt \
  configs/phase_coherence.routes.local.json "results/phase_routes_${RUN_STAMP}" \
  --profile repository --bootstrap 500 --seed 0
```

Update the session list and local config to point at your actual inputs first.
Check mode validates data/mapping; it does not certify recording quality or guarantee
that RAM/disk will suffice. Sessions run sequentially and failures are recorded in
`status.json`; a failed session does not stop the remaining sessions. Use a new output
folder when rerunning. Existing progress reporting is retained. Power context is retained only for selected/control
channels, and previous-session tensors are released before loading the next recording.

## Export completed sessions to Figure 6

The following reads only sessions marked completed; it does not rerun any analyses:

```bash
python - "results/phase_routes_${RUN_STAMP}" <<'PY'
import json
import sys
from pathlib import Path
from coupling.figure6 import export
run = Path(sys.argv[1]).resolve()
for row in json.loads((run / 'status.json').read_text()):
    if row['status'] == 'completed':
        export(Path('figures/fig06_lfp_coupling') / ('exports_' + run.name) / row['output_key'],
               phase_session=row['output'])
PY
```

Or export one completed session folder (the folder containing `conditions.json`):

```bash
coupling-export-figure6 --phase-session results/RUN/SESSION_KEY \
  --output figures/fig06_lfp_coupling/exports_SESSION
```

`python lfp_phase_coherence/03_export_figure6.py` accepts the same arguments.
An existing destination is never overwritten. Invalid inputs or failed rendering leave
no final export folder. Numerical values and SHA-256 source hashes accompany each export.
Missing matched pairs are reported in the manifest; no available pair is an error.

## Interpretation and reference adaptation

- The `repository` profile retains the existing Mahmood-style representative-channel
  selection, phase-vector coherence, replacement-shuffle controls and band definitions.
- Comparisons are oBA versus BA and oIA versus IA. The combined route estimate gives
  the two within-odor coherence estimates equal weight; it never pools opposite odor
  phase lags before estimating coherence. S and BK remain separate controls in the analysis.
- Clean exports contain time–frequency coherence, signed ortho-minus-retro maps, and
  band traces. Both route maps share a 0–1 scale; difference maps use −1 to 1.
- Band shading is pointwise 95% bootstrap percentiles, explicitly different from the
  ±3 bootstrap-SD shading of the full diagnostic plots. Neither is a corrected route
  significance test. No significance stars are added.
- These are single-session descriptive panels. Animals, not channels or repeated
  sessions, must define replication for the final paper's inferential statistics.
- `coupling-lfp-figure` / `02_spectral_coupling.py` retains magnitude-squared coherence
  and diagnostic-gated Granger as a separate workflow. Phase coherence is not MSC.
  This initial clean exporter covers power and phase coherence, not Granger panels.

## Tests

```bash
python -m pip install -e '.[phase,spikes]'
python -m unittest discover -s tests -v
```

Tests cover numerical recovery, mapping/QC, synthetic HDF5 input, condition weighting,
export shape/grid validation, path safety, failed-render cleanup and PDF/SVG/PNG writing.
No test can guarantee crash-free execution on arbitrary recordings or available memory.
