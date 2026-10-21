# LFP power spectrogram generation

Run from the repository root after `python -m pip install -e '.[phase]'`.
The numbered scripts are thin entry points into the tested `src/coupling` package.
Existing console commands, Python imports and configuration paths are unchanged.

| Stage | Command | Output |
| --- | --- | --- |
| Raw Intan → QC + cleaned trials | `coupling-prepare --help` / `python lfp_power/01_prepare.py --help` | Cleaned LFP NPZ, rejection masks, preprocessing provenance |
| Cleaned trials → power | `coupling-power CLEANED.npz results/power_RUN` | `session_power.npz`, QC counts, diagnostic plots |
| Compare completed power sessions | `coupling-compare-power --help` | Session comparisons; see [details](power_comparison.md) |
| Paper panels | `coupling-export-figure6 --power results/power_RUN/session_power.npz --output figures/fig06_lfp_coupling/exports_RUN` | PDF, SVG, PNG, source values and manifest |

Power is squared STFT magnitude, normalized per accepted trial and frequency against
−1.5 to −0.5 s, then averaged with equal trial and channel weights. This is **not**
Mahmood's original median-amplitude/baseline-ratio overview. Figure labels retain the
actual normalization. No filtering, reference or trial rejection is changed for appearance.
Power panels use the event label from preprocessing provenance rather than assuming TTL.

`python lfp_power/demo_power.py` creates explicitly synthetic input under `results/`.
The old `scripts/demo_power.py` remains a compatibility entry point.

Command to run coupling-summarize-LFP: output is to figures/fig06_lfp_coupling/summary_YYYYMMDD_HHMMSS, where YYYYMMDD_HHMMSS is the current date and time.

coupling-summarize-lfp \
   --phase-run "$LFP_RUN/phase" \
   --power-root "$LFP_RUN" \
   --output "figures/fig06_lfp_coupling/summary_$(date +%Y%m%d_%H%M%S)"