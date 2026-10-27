
from coupling.array_store import load_arrays, save_arrays
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from coupling.lfp_summary import animal_average, read_power, read_phase, summarize


class LFPSummaryTests(unittest.TestCase):
    def test_equal_animal_weighting(self):
        values = np.array([[0.], [.4], [.8]])
        animals, per, grand = animal_average(values, ["A","A","B"])
        np.testing.assert_allclose(per[:,0], [.2,.8])
        np.testing.assert_allclose(grand, [.5])
        with self.assertRaises(ValueError):
            animal_average(values, ["A","A",""])

    def fixture(self, root, n=2):
        run = root/"phase"; run.mkdir()
        power = root/"power"; power.mkdir()
        statuses = []
        t, f = np.linspace(-1,1,8), np.array([4.,8.,12.])
        for i in range(n):
            sid = f"session_{i}"; animal = f"animal_{i}"
            folder = run/sid; folder.mkdir()
            statuses.append(dict(status="completed", session=f"/raw/{sid}",
                                 output_key=sid, animal_id=animal))
            (folder/"provenance.json").write_text(json.dumps(dict(
                animal_id=animal, profile="repository", settings={"window_samples":500})))
            conditions = {}
            for label, value in [("ortho", .6+i*.1), ("retro", .3+i*.1)]:
                target = folder/label; target.mkdir()
                save_arrays(target/"phase_coherence.h5", time_s=t, frequency_hz=f,
                         band=["theta"], band_limits_hz=[[3,8]],
                         alignment_label="Time from stimulus DIN onset (s)",
                         selected_electrode_id=[1,33],
                         bootstrap_band=np.full((3,1,len(t)), value))
                conditions[label] = dict(status="completed", directory=label)
            (folder/"conditions.json").write_text(json.dumps(conditions))
            p = power/sid; p.mkdir()
            x = np.empty((2,4,len(f),len(t)))
            # Noncanonical order; code must resolve labels.
            for region in range(2):
                for stim in range(4):
                    x[region,stim] = region*10+stim+i
            save_arrays(p/"session_power.h5", time_s=t, frequency_hz=f,
                     region=["PC","GC"], stimulus=["IA","oIA","BA","oBA"], z_power=x,
                     animal_id=animal, session_id=sid, settings_json=json.dumps({"baseline_s":[-1.5,-.5]}))
            (p/"provenance.json").write_text(json.dumps({
                "preprocessing":{"event_alignment":"stimulus DIN onset"}}))
        statuses.append(dict(status="failed", session="/raw/failed", error="QC failed"))
        (run/"status.json").write_text(json.dumps(statuses))
        return run, power

    def test_power_label_mapping_equal_odors(self):
        with tempfile.TemporaryDirectory() as tmp:
            run, root = self.fixture(Path(tmp), n=1)
            records = read_power(root, "route", [])
            # GC ortho = average GC oIA=11 and GC oBA=13.
            np.testing.assert_allclose(records[0]["values"][0,0], 12.)
            np.testing.assert_allclose(records[0]["values"][1,0], 11.)

    def test_end_to_end_saved_figures_and_exclusions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); run, power = self.fixture(root)
            out = summarize(root/"out", run, power)
            for stem in ("all_sessions_power","animal_mean_power",
                         "all_sessions_phase","animal_mean_phase"):
                for ext in ("png","pdf","svg"):
                    self.assertGreater((out/f"{stem}.{ext}").stat().st_size, 100)
            with load_arrays(out/"phase_summary_values.h5") as d:
                np.testing.assert_allclose(d["equal_animal_mean"][0], .65)
                np.testing.assert_allclose(d["equal_animal_mean"][1], .35)
            manifest = json.loads((out/"summary_manifest.json").read_text())
            self.assertEqual(len(manifest["excluded"]), 1)
            self.assertEqual(len(manifest["phase_sessions"]), 2)
            with self.assertRaises(FileExistsError):
                summarize(out, run, power)

    def test_grid_mismatch_fails_before_publication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); run, power = self.fixture(root)
            path = run/"session_1/retro/phase_coherence.h5"
            with load_arrays(path) as d:
                values = dict(d)
            values["time_s"] = values["time_s"]+.01
            save_arrays(path, **values)
            with self.assertRaisesRegex(ValueError, "time_s"):
                summarize(root/"out", run, power)
            self.assertFalse((root/"out").exists())

    def test_duplicate_power_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); run, power = self.fixture(root)
            path = power/"session_1/session_power.h5"
            with load_arrays(path) as d:
                values = dict(d)
            values.update(animal_id="animal_0", session_id="session_0")
            save_arrays(path, **values)
            with self.assertRaisesRegex(ValueError, "Duplicate power"):
                read_power(power, "route", [])

    def test_render_failure_is_atomic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); run, power = self.fixture(root, n=1)
            with patch("coupling.lfp_summary.save", side_effect=OSError("disk full")):
                with self.assertRaises(OSError):
                    summarize(root/"out", run, power)
            self.assertFalse((root/"out").exists())
            self.assertFalse(list(root.glob(".lfp-summary-*")))

    def test_missing_condition_recorded_not_replaced(self):
        with tempfile.TemporaryDirectory() as tmp:
            run, power = self.fixture(Path(tmp), n=1)
            path = run/"session_0/conditions.json"
            conditions = json.loads(path.read_text())
            conditions["retro"]["status"] = "skipped"
            path.write_text(json.dumps(conditions))
            excluded = []
            self.assertEqual(read_phase(run, "route", [], excluded), [])
            self.assertEqual(len(excluded), 2)


if __name__ == "__main__":
    unittest.main()
