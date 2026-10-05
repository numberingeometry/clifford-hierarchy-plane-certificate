from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/checking'))
import copy
import unittest

import c14


class ConventionAndRejectionTests(unittest.TestCase):
    def test_hstar_input_phases_and_branch_order(self):
        permutation, phases = c14.hstar()
        self.assertEqual((permutation[1], phases[1]), (3, 1))
        self.assertEqual((permutation[3], phases[3]), (1, 3))
        self.assertEqual((permutation[5], phases[5]), (5, 2))
        self.assertEqual((permutation[6], phases[6]), (7, 0))

    def test_pauli_composition_is_right_to_left_with_phase(self):
        X = tuple(j ^ 4 for j in range(8)), (0,)*8
        Z = tuple(range(8)), tuple(2*((j>>2)&1) for j in range(8))
        XZ, ZX = c14.compose(X,Z), c14.compose(Z,X)
        self.assertEqual(XZ[0], ZX[0])
        self.assertTrue(all((a-b)%4 == 2 for a,b in zip(XZ[1],ZX[1])))

    def test_reject_noncommuting_clifford_and_nonhermitian_ray(self):
        Zy = tuple(range(8)), tuple(2*((j>>1)&1) for j in range(8))
        self.assertFalse(c14.commutes(Zy,c14.hstar()))
        Sy = tuple(range(8)), tuple((j>>1)&1 for j in range(8))
        self.assertIsNone(c14.hermitian_representative(Sy))

    def test_non_real_scalar_hermitian_lift(self):
        imaginary_identity = tuple(range(8)), (1,)*8
        self.assertEqual(c14.hermitian_representative(imaginary_identity), (tuple(range(8)), (0,)*8))

    def test_results_reject_omission_and_wrong_family(self):
        cresult = c14.run()
        bad = copy.deepcopy(cresult)
        bad["complete_domain_size"] -= 1
        with self.assertRaises(c14.CheckFailure):
            c14.validate_summary(bad)
        bad = copy.deepcopy(cresult)
        bad["claimed_family_equals_hermitian_rays"] = False
        with self.assertRaises(c14.CheckFailure):
            c14.validate_summary(bad)


if __name__ == "__main__":
    unittest.main()
