"""Smoke tests for the publication-style sniff figure and 200-ms bar."""
import tempfile
import unittest
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# 02_figure is a numbered workflow entry point, load from filesystem.
import importlib.util

MODULE = Path(__file__).resolve().parents[1] / "respiration_sniffing" / "02_figure.py"
spec = importlib.util.spec_from_file_location("sniff_figure", MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class FigureTests(unittest.TestCase):
    def test_bar_exactly_200ms(self):
        fig, ax = plt.subplots()
        ax.plot([-.5, 2.5], [0, 1])
        ax.set_xlim(-.5, 2.5)
        left, right = mod.scale_bar(ax)
        self.assertAlmostEqual(right-left, .2)
        plt.close(fig)

    def test_four_panel_figure(self):
        t = np.arange(-500, 2500)/1000
        rng = np.random.default_rng(42)
        synthetic = {
            "oBA":[np.sin(2*np.pi*5*t) + .05*rng.standard_normal(len(t)) for _ in range(3)],
            "BA":[np.sin(2*np.pi*4*t+.3) + .05*rng.standard_normal(len(t)) for _ in range(3)]
        }
        with tempfile.TemporaryDirectory() as root:
            dest = Path(root)/"reference_style.pdf"
            mod.figure(synthetic,t,dest)
            self.assertTrue(dest.exists())
            self.assertGreater(dest.stat().st_size,1000)
