"""Mathematical primitive controls and domain/conclusion rejection tests for (C19)."""
from copy import deepcopy
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from checking.root_words import binary as f
from checking.root_words.runner import run, validate


class RootWordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c19 = run('C18')['result']

    def test_binary_characteristic_polynomial(self):
        # All 512 binary 3x3 matrices satisfy their independently computed Leibniz characteristic
        # polynomial (Cayley-Hamilton), which covers the shared exact primitives.
        for packed in range(1 << 9):
            M = tuple((packed >> (3 * j)) & 7 for j in range(3))
            polynomial = f.characteristic_polynomial(M)
            value, current = (0, 0, 0), f.identity(3)
            for k in range(4):
                if polynomial >> k & 1:
                    value = tuple(a ^ b for a, b in zip(value, current))
                current = f.multiply(M, current)
            self.assertEqual(value, (0, 0, 0))

    def test_irreducibility_controls(self):
        self.assertTrue(f.polynomial_irreducible(0b1001001))      # z^6+z^3+1
        self.assertFalse(f.polynomial_irreducible(0b1000001))     # z^6+1 = (z^3+1)^2
        self.assertTrue(f.polynomial_irreducible(0b111))          # z^2+z+1

    def test_order_bound_is_enforced(self):
        with self.assertRaisesRegex(ValueError, 'finite bound'):
            f.order(((2, 1, 4)), 1)

    def test_changed_root_domain(self):
        bad = deepcopy(self.c19)
        bad['case']['edges'][4] = [3, 5]
        with self.assertRaisesRegex(ValueError, 'root domain'):
            validate(bad)

    def test_changed_order(self):
        bad = deepcopy(self.c19)
        bad['case']['order_A'] = 4
        with self.assertRaisesRegex(ValueError, 'order of A'):
            validate(bad)

    def test_changed_characteristic_polynomial(self):
        bad = deepcopy(self.c19)
        bad['case']['characteristic_polynomial_bits'] ^= 2
        with self.assertRaisesRegex(ValueError, 'characteristic polynomial'):
            validate(bad)

    def test_unknown_key(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported'):
            run('C19')

    def test_optimized_mode_cannot_bypass_failures(self):
        cp = subprocess.run([sys.executable, '-B', '-O', '-c',
                             "import sys;sys.path.insert(0,sys.argv[1]);from checking.root_words.runner import validate;"
                             "validate({'case':{'name':'D4','dimension':4,'polar_rank':4,'edges':[]}})",
                             str(ROOT / 'src')], capture_output=True, text=True)
        self.assertNotEqual(cp.returncode, 0)
        self.assertIn('root domain', cp.stderr)


if __name__ == '__main__':
    unittest.main()
