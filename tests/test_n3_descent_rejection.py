"""Reject omitted C11/C12 domains and disabled proof-critical assertions."""
from copy import deepcopy
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from checking.n3_descent.runner import run, validate_c11, validate_c12


class DescentRejectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c11 = run("C11")["result"]
        cls.c12 = run("C12")["result"]

    def test_incomplete_ray_domain(self):
        bad = deepcopy(self.c11)
        bad['non_Pauli_rays'] -= 1
        with self.assertRaisesRegex(ValueError, "domain incomplete"):
            validate_c11(bad)

    def test_omitted_phase_family(self):
        bad = deepcopy(self.c12)
        bad['families'] = bad['families'][:-1]
        with self.assertRaisesRegex(ValueError, "family incomplete"):
            validate_c12(bad)

    def test_omitted_target_rows(self):
        bad = deepcopy(self.c12)
        bad['target_rows_checked'] -= 1
        with self.assertRaisesRegex(ValueError, "domain incomplete"):
            validate_c12(bad)

    def test_false_lattice_conclusion(self):
        bad = deepcopy(self.c12)
        bad['all_integer_lattice_inclusions'] = False
        with self.assertRaisesRegex(ValueError, "implication failed"):
            validate_c12(bad)

    def test_disabled_assertions_rejected(self):
        cp = subprocess.run([sys.executable, "-B", "-O", "-c",
                             "import sys;sys.path.insert(0,sys.argv[1]);import checking.n3_descent.runner",
                             str(ROOT / "src")], capture_output=True, text=True)
        self.assertNotEqual(cp.returncode, 0)
        self.assertIn("Run without -O", cp.stderr)


if __name__ == "__main__":
    unittest.main()
