from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/checking'))
import copy
import unittest
import c15


class C15Tests(unittest.TestCase):
    def test_coordinate_order_and_right_composition(self):
        self.assertEqual(c15.TAU0,(0,1,2,3,4,5,7,6))
        flip_x = tuple(j^4 for j in range(8))
        self.assertNotEqual(c15.compose(c15.TAU0,flip_x),c15.compose(flip_x,c15.TAU0))
        self.assertEqual(c15.compose(c15.TAU0,c15.TAU0),c15.IDENTITY)

    def test_affine_plane_translation_uses_linear_part(self):
        alpha = tuple(j^7 for j in range(8))
        config = 2,frozenset((4,5,6,7))
        self.assertEqual(c15.action(alpha,config),(2,frozenset((0,1,2,3))))

    def test_complete_census_rejects_changed_orbits_and_omitted_domain(self):
        result = c15.run()
        bad = copy.deepcopy(result)
        bad["representative_orbit_sizes"]["D_z_eq_1_u_y"] -= 1
        with self.assertRaises(c15.CheckFailure):
            c15.validate_summary(bad)
        bad = copy.deepcopy(result)
        bad["complete_permutation_domain"] -= 1
        with self.assertRaises(c15.CheckFailure):
            c15.validate_summary(bad)


if __name__ == "__main__":
    unittest.main()
