"""Semantic rejection controls; full-domain arithmetic is run separately."""
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from checking.exact_audit import mono3 as q
from checking.exact_audit.runner import validate_c13_conjugators


class ConjugatorRejectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.valid = json.loads((ROOT / "src/checking/exact_audit/data/n3-exceptional-conjugators.json").read_text())
        first = cls.valid[0]
        C = q.mul(q.diag(first["E_phase_zeta"]), q.perm(first["alpha"]))
        cls.H = q.conj(C, (tuple(first["pi"]), tuple(first["phase_zeta"])))
        # These unit controls use the supplied rays as a fixture. The accepted
        # full-domain run separately regenerates and binds the exceptional set.
        validate_c13_conjugators(cls.valid, cls.valid, cls.H)

    def test_missing_record(self):
        with self.assertRaisesRegex(ValueError, "count"):
            validate_c13_conjugators(self.valid[:-1], self.valid, self.H)

    def test_duplicate_record(self):
        bad = deepcopy(self.valid)
        bad[-1] = deepcopy(bad[0])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            validate_c13_conjugators(bad, self.valid, self.H)

    def test_corrupted_exact_phase(self):
        bad = deepcopy(self.valid)
        bad[0]["E_phase_zeta"][0] = (bad[0]["E_phase_zeta"][0] + 1) % 8
        with self.assertRaisesRegex(ValueError, "Clifford|equality"):
            validate_c13_conjugators(bad, self.valid, self.H)

    def test_wrong_dimension(self):
        bad = deepcopy(self.valid)
        bad[0]["alpha"] = bad[0]["alpha"][:-1]
        with self.assertRaisesRegex(ValueError, "dimension"):
            validate_c13_conjugators(bad, self.valid, self.H)

    def test_non_integer_phase(self):
        bad = deepcopy(self.valid)
        bad[0]["E_phase_zeta"][0] = True
        with self.assertRaisesRegex(ValueError, "integer"):
            validate_c13_conjugators(bad, self.valid, self.H)

    def test_invalid_phase_range(self):
        bad = deepcopy(self.valid)
        bad[0]["phase_zeta"][0] = 8
        with self.assertRaisesRegex(ValueError, "range"):
            validate_c13_conjugators(bad, self.valid, self.H)

    def test_nonaffine_conjugator(self):
        bad = deepcopy(self.valid)
        bad[0]["alpha"][0], bad[0]["alpha"][1] = bad[0]["alpha"][1], bad[0]["alpha"][0]
        with self.assertRaisesRegex(ValueError, "permutation"):
            validate_c13_conjugators(bad, self.valid, self.H)

    def test_disabled_assertions_rejected(self):
        cp = subprocess.run([sys.executable, "-B", "-O", "-c",
                             "import sys;sys.path.insert(0,sys.argv[1]);import checking.exact_audit.runner",
                             str(ROOT / "src")], capture_output=True, text=True)
        self.assertNotEqual(cp.returncode, 0)
        self.assertIn("Run without -O", cp.stderr)


if __name__ == "__main__":
    unittest.main()
