"""Compatibility entry point; canonical demo lives in lfp_phase_coherence/demo_lfp_coupling.py."""
from pathlib import Path
import runpy

if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parents[1] / "lfp_phase_coherence/demo_lfp_coupling.py"), run_name="__main__")
