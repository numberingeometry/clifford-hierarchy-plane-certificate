from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/checking'))
import copy
import unittest
import c16


class C16Tests(unittest.TestCase):
    def test_y_phase_is_hermitian_and_involutory(self):
        Yx = c16.pauli((4<<3)|4)
        self.assertEqual((Yx[0][0],Yx[1][0]),(4,1))
        self.assertEqual((Yx[0][4],Yx[1][4]),(0,3))
        self.assertEqual(c16.mono_inverse(Yx),Yx)
        self.assertEqual(c16.mono_compose(Yx,Yx),(c16.IDENTITY,(0,)*8))

    def test_projector_with_y_generator_avoids_unphased_pauli_error(self):
        basis = (36,2,1)
        projector = c16.projector_numerator(basis,(0,0,0))
        # |+y,0,0><+y,0,0| has imaginary off-diagonal entries.
        self.assertEqual(projector[4],(0,-4))
        self.assertEqual(projector[32],(0,4))

    def test_support_rejects_nonclifford_phase(self):
        CS = c16.IDENTITY,tuple(((j>>1)&1)*(j&1) for j in range(8))
        self.assertFalse(c16.is_clifford(CS))
        self.assertTrue(c16.is_third_level(CS))

    def test_complete_run_and_reject_omitted_frame_or_control(self):
        result = c16.run()
        bad = copy.deepcopy(result)
        bad["frame_affine_pairs_checked"] -= 1344
        with self.assertRaises(c16.CheckFailure):
            c16.validate_summary(bad)
        bad = copy.deepcopy(result)
        bad["counterexample_product_flag_indices"] = []
        with self.assertRaises(c16.CheckFailure):
            c16.validate_summary(bad)
        bad = copy.deepcopy(result)
        bad["counterexample_common_flag_indices"] = [0]
        with self.assertRaises(c16.CheckFailure):
            c16.validate_summary(bad)


if __name__ == "__main__":
    unittest.main()
