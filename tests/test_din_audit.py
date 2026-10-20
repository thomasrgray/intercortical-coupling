import hashlib
import json
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from coupling.din_audit import read_pulses, audit, classify


class DinAuditTests(unittest.TestCase):
    def write(self, path, x):
        np.asarray(x, dtype="<u2").tofile(path)
        return path

    def test_chunk_edges_and_widths(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp)/"din.dat", [0,0,1,1,1,1,0,0,1,0])
            on, off, good, n = read_pulses(path, chunk=3)
            np.testing.assert_array_equal(on, [2,8])
            np.testing.assert_array_equal(off, [6,9])
            self.assertTrue(good.all())
            self.assertEqual(n, 10)

    def test_censored_boundaries_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp)/"din.dat", [1,1,0,1,1,0,1])
            on, off, good, n = read_pulses(path, chunk=2)
            np.testing.assert_array_equal(on, [0,3,6])
            np.testing.assert_array_equal(off, [2,5,7])
            np.testing.assert_array_equal(good, [False,True,False])
            self.write(path, [1,1,1])
            self.assertFalse(read_pulses(path, chunk=1)[2].any())

    def test_idle_and_invalid(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp)/"din.dat", [0,0,0])
            self.assertEqual(len(read_pulses(path)[0]), 0)
            self.write(path, [0,2,0])
            with self.assertRaises(ValueError):
                read_pulses(path)
            path.write_bytes(b"\x00")
            with self.assertRaises(ValueError):
                read_pulses(path)
            path.write_bytes(b"")
            with self.assertRaises(ValueError):
                read_pulses(path)

    def test_heuristic_does_not_assign_odor(self):
        self.assertEqual(classify(np.array([2000.]), "oBA", 100), ("long", "consistent"))
        self.assertEqual(classify(np.array([15.]), "BA", 100), ("short", "consistent"))
        self.assertEqual(classify(np.array([15.]), "oBA", 100), ("short", "REVIEW"))
        self.assertEqual(classify(np.array([2000.]), "unknown", 100), ("long", "unassigned"))

    def fixture(self, root, name):
        session = root/name
        session.mkdir()
        (session/"info.rhd").write_bytes(struct.pack("<Ihhf", 0xc6912702, 1, 0, 1000.))
        (session/"session.info").write_text(json.dumps({
            "taste_params": {"dig_in_nums": [0,2], "tastes": ["BA","oBA"]}}))
        x = np.zeros(2500, dtype="<u2")
        x[10:25] = 1
        self.write(session/"board-DIN-00.dat", x)
        x[:] = 0
        x[100:2100] = 1
        x[-10:] = 1
        self.write(session/"board-DIN-02.dat", x)
        return session

    def test_end_to_end_comparison_and_no_raw_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = self.fixture(root, "a")
            b = self.fixture(root, "b")
            self.write(b/"board-DIN-07.dat", [0,2,0])
            listing = root/"sessions.txt"; listing.write_text("# sessions\na\nb\n")
            before = {p: hashlib.sha256(p.read_bytes()).hexdigest()
                      for folder in (a,b) for p in folder.iterdir()}
            rows = audit(listing, root/"output")
            r = next(r for r in rows if r.get("din") == 2 and r["session"] == "a")
            self.assertEqual(r["median_ms"], 2000.)
            self.assertEqual(r["complete_pulses"], 1)
            self.assertEqual(r["partial_pulses"], 1)
            self.assertTrue(any(r["status"] == "error" for r in rows))
            for name in ("pulse_widths.csv", "din_duration_summary.csv",
                         "duration_comparison_01.png", "duration_comparison_01.pdf", "audit.json"):
                self.assertTrue((root/"output"/name).is_file(), name)
            for p, digest in before.items():
                self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), digest)
            with self.assertRaises(FileExistsError):
                audit(listing, root/"output")
            with self.assertRaises(ValueError):
                audit(listing, a/"output")

    def test_bad_metadata_still_measures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); a = self.fixture(root, "a")
            (a/"session.info").write_text("{}")
            listing = root/"list.txt"; listing.write_text("a\n")
            with patch("coupling.din_audit.plot_results"):
                rows = audit(listing, root/"out")
            self.assertEqual(rows[0]["median_ms"], 15.)
            self.assertEqual(rows[0]["stimulus"], "")
            self.assertEqual(rows[0]["status"], "warning")

    def test_duplicate_and_missing_din(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); a = self.fixture(root, "a")
            (a/"board-DIN-02.dat").unlink()
            self.write(a/"board-DIN-0.dat", [0,1,0])
            listing = root/"list.txt"; listing.write_text("a\n")
            with patch("coupling.din_audit.plot_results"):
                rows = audit(listing, root/"out")
            self.assertTrue(all(r["status"] == "error" for r in rows))
            self.assertTrue(any("Declared DIN" in r["notes"] for r in rows))


if __name__ == "__main__":
    unittest.main()
