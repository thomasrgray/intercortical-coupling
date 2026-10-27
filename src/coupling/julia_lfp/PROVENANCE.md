# Source provenance

Ported from https://github.com/thomasrgray/julia_blech_clust_testing at
`1a4a4eb6346b00226e907228f4a2c3f658223907` (MIT; LICENSE included).

Includes all five LFP workflow scripts, shared helpers, stimulus/animal labels,
and the matching ephys_data, lfp_processing and BAKS dependencies. Source plots,
QC thresholds, normalization, reductions, and statistical contrasts are retained.
Integration changes: package-relative imports, repository defaults, comment-aware
session lists, lazy imports for unrelated GUI/spike ANOVA methods, read-only
loading of existing Parsed_LFP/channel metadata, an error instead of an interactive
prompt for ambiguous HDF5 files, and bounded STFT workers (COUPLING_LFP_JOBS,
default 1). Worker failures now stop instead of continuing with partial arrays.
The orchestrator and archive utility are new. Upstream tests are included.

Existing Parsed_LFP is reused. If absent, the upstream extractor writes derived
LFP nodes into the session HDF5; this is not a wholly read-only workflow.
