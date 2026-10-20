"""Compatibility entry point; canonical demo lives in lfp_power/demo_power.py."""
from pathlib import Path
import runpy

if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parents[1] / "lfp_power/demo_power.py"), run_name="__main__")
