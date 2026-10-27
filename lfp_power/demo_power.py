"""Generate clearly synthetic, already-cleaned demo input, never a study result."""

from coupling.array_store import save_arrays
from pathlib import Path
import json
import numpy as np
from coupling.power import STIMULI
out=Path("results/demo_input.h5")
if out.exists():
    raise FileExistsError(out)
out.parent.mkdir(exist_ok=True)
rng=np.random.default_rng(20261011)
t=np.arange(-2000,2500)/1000
labels=np.repeat(STIMULI,3)
x=rng.normal(size=(2,len(labels),len(t)))
for c in range(2):
    x[c] += (t>0)*np.sin(2*np.pi*(12+8*c)*t)*1.5
save_arrays(out,lfp_uv=x,time_s=t,stimulus=labels,electrode_id=[12,59],
 region=np.array(["GC","PC"]),good_channel=np.ones(2,bool),good_trial=np.ones((2,len(labels)),bool),
 session_id="SYNTHETIC_DEMO",animal_id="SYNTHETIC",
 provenance_json=json.dumps(dict(source="synthetic seeded Gaussian + sinusoid",reference="none; synthetic",
 filtering="none; synthetic",qc_method="synthetic finite data; all accepted",event_alignment="synthetic TTL zero")))
print(out)
