"""Controls for (C18): the hierarchy arithmetic of U_E6 against an independent dense reference, and
rejection of altered caches, failure chains and incomplete reports. No full hierarchy run."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
import dense_reference as ref  # noqa: E402
from checking import e2_hierarchy as subject  # noqa: E402
from checking import e6_level  # noqa: E402

SELECTED = (0, 1, 2, 3, 4, 16, 31, 32, 33, 64, 65, 98, 255, 257, 511, 682, 1023)


def hierarchy_fixture():
    item = {"in_C5": True, "in_C4": False, "exact": True,
            "lower_failure_chain_independently_verified": True,
            "lower_failure_chain": {"level": 4, "pauli_label": 1, "next": {"level": 3}},
            "complete_root_Pauli_count": 1024, "full_Pauli_domain_at_levels_at_least_four": 1024,
            "generator_count_at_levels_two_and_three": 10,
            "root_words": [[1, 24, 2, 48, 6], [8, 24]], "control_block_order": ["I", "A", "B", "AB"],
            "unitarity_checked_exactly": True, "root_lifts_square_to_declared_Paulis": True,
            "upper_coverage": {"every_distinct_upper_C4_node_checked_all_1024_Paulis": True,
                               "every_upper_C3_and_C2_node_checked_all_ten_generators": True,
                               "distinct_memoized_operators_by_level": {"1": 100, "2": 40, "3": 20, "4": 1024, "5": 1},
                               "computed_conjugations_by_level": {"2": 400, "3": 200, "4": 1048576, "5": 1024}}}
    return {"status": "HIERARCHY_PASS", "exact": True, "results": {"E6": item}}


def block_fixture():
    return {"status": "EXACT_BLOCK_AUDIT_PASS", "exact": True, "blocks_checked": 3,
            "all_registered_blocks_covered": True, "unitarity_blocks_checked": 3,
            "signed_Pauli_images_checked": 192, "Pauli_action_domain_per_block": 64}


class DenseReferenceControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gate = subject.gates()["E6"]
        cls.dense = ref.gate("E6")

    def test_root_lifts_square_to_paulis_and_are_unitary(self):
        for symbols in sorted({word for side in ref.WORDS["E6"] for word in side}):
            with self.subTest(symbols=symbols):
                root = ref.root(symbols)
                self.assertEqual(ref.product(root, root), ref.pauli_word(symbols))
                self.assertEqual(ref.product(root, ref.adjoint(root)), ref.identity(8))

    def test_gate_is_unitary_and_matches_the_independent_construction(self):
        self.assertEqual(ref.product(self.dense, ref.adjoint(self.dense)), ref.identity(32))
        self.assertEqual(subject.to_dense(self.gate), self.dense)

    def test_conjugations_match_dense_arithmetic_including_control_shifts(self):
        for label in SELECTED:
            with self.subTest(label=label):
                self.assertEqual(subject.to_dense(subject.conjugate_pauli(self.gate, label)),
                                 ref.conjugate(self.dense, label))
        shifted, shifted_dense = subject.conjugate_pauli(self.gate, 1), ref.conjugate(self.dense, 1)
        for label in (2, 33, 65, 98, 257, 682, 1023):
            with self.subTest(shifted_label=label):
                self.assertEqual(subject.to_dense(subject.conjugate_pauli(shifted, label)),
                                 ref.conjugate(shifted_dense, label))

    def test_all_1024_paulis_bind_to_the_tensor_convention(self):
        unit = subject.conjugate_pauli(self.gate, 0)
        self.assertEqual(subject.to_dense(unit), ref.identity(32))
        for label in range(1024):
            with self.subTest(label=label):
                image = subject.conjugate_pauli(unit, label)
                self.assertEqual(subject.to_dense(image), ref.pauli(label))
                self.assertTrue(subject.is_pauli(image))

    def test_common_scalar_key_preserves_relative_block_phase(self):
        shift, blocks = self.gate
        expected = subject.canonical(self.gate)
        for exponent, scalar in enumerate(((1, 0), (0, 1), (-1, 0), (0, -1))):
            altered = (shift, tuple((ident, (phase+exponent) % 4) for ident, phase in blocks))
            self.assertEqual(subject.canonical(altered), expected)
            self.assertEqual(subject.to_dense(altered), ref.global_phase(self.dense, scalar))
        altered = (shift, (blocks[0], (blocks[1][0], (blocks[1][1]+1) % 4), *blocks[2:]))
        self.assertNotEqual(subject.canonical(altered), expected)
        self.assertNotEqual(subject.to_dense(altered), self.dense)


class CacheAndChainRejection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gate = subject.gates()["E6"]

    def test_corrupted_cached_action_rejected_by_block_audit(self):
        actions = list(subject._ACTIONS)
        identity_action = list(actions[0])
        label, sign = identity_action[3]
        identity_action[3] = (label, (sign+2) % 4)
        actions[0] = tuple(identity_action)
        with patch.object(subject, "_ACTIONS", actions):
            with self.assertRaisesRegex(ValueError, "stored Clifford action does not match"):
                subject.audit_exact_blocks()

    def test_boolean_phase_and_cached_sign_rejected(self):
        shift, blocks = self.gate
        with self.assertRaisesRegex(ValueError, "invalid coherent Clifford reference"):
            subject.canonical((shift, ((blocks[0][0], False), *blocks[1:])))
        actions = list(subject._ACTIONS)
        identity_action = list(actions[0])
        identity_action[3] = (identity_action[3][0], False)
        actions[0] = tuple(identity_action)
        with patch.object(subject, "_ACTIONS", actions):
            with self.assertRaisesRegex(ValueError, "malformed Pauli/sign fields"):
                subject.audit_exact_blocks()

    def test_level_four_failure_chain_rejects_alterations(self):
        chain = subject.Engine().failure(self.gate, 4)
        self.assertTrue(subject.verify_failure(self.gate, chain, 4))
        changed = deepcopy(chain)
        changed["level"] = 3
        with self.assertRaisesRegex(ValueError, "level changed or skipped"):
            subject.verify_failure(self.gate, changed, 4)
        changed = deepcopy(chain)
        del changed["next"]["next"]
        with self.assertRaises(ValueError):
            subject.verify_failure(self.gate, changed, 4)
        changed = deepcopy(chain)
        changed["pauli_label"] = True
        with self.assertRaisesRegex(ValueError, "malformed failure-chain Pauli label"):
            subject.verify_failure(self.gate, changed, 4)
        changed = deepcopy(chain)
        changed["next"]["next"]["next"]["non_pauli_dense"]["gaussian_integer_pairs"][0][0] += 1
        with self.assertRaisesRegex(ValueError, "leaf operator changed"):
            subject.verify_failure(self.gate, changed, 4)

    def test_chain_ending_in_a_pauli_rejected(self):
        identity = subject.pauli_operator(0)
        leaf = {"level": 1, "non_pauli_dense": {
            "denominator_power": 0, "gaussian_integer_pairs": [list(v) for v in ref.identity(32)[1]]}}
        for level in (2, 3, 4):
            leaf = {"level": level, "pauli_label": 0, "next": leaf}
        with self.assertRaisesRegex(ValueError, "leaf is Pauli"):
            subject.verify_failure(identity, leaf, 4)


class ReportRejection(unittest.TestCase):
    def test_complete_reports_pass(self):
        e6_level.validate_hierarchy(hierarchy_fixture())
        e6_level.validate_block_audit(block_fixture())

    def test_level_or_domain_alterations_rejected(self):
        for change, message in (
                (lambda r: r["results"]["E6"].update(in_C5=False), "level-five conclusion"),
                (lambda r: r["results"].update(A7=deepcopy(r["results"]["E6"])), "not the gate U_E6"),
                (lambda r: r["results"]["E6"]["root_words"][1].append(2), "root words"),
                (lambda r: r["results"]["E6"].update(complete_root_Pauli_count=1023), "root Pauli domain"),
                (lambda r: r["results"]["E6"]["upper_coverage"]["computed_conjugations_by_level"].update({"4": 1048575}),
                 "level-four coverage"),
                (lambda r: r["results"]["E6"]["upper_coverage"]["computed_conjugations_by_level"].update({"3": 199}),
                 "lower group domain")):
            report = hierarchy_fixture()
            change(report)
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                e6_level.validate_hierarchy(report)

    def test_partial_block_audit_rejected(self):
        report = block_fixture()
        report["signed_Pauli_images_checked"] -= 1
        with self.assertRaisesRegex(ValueError, "signed Pauli cache audit"):
            e6_level.validate_block_audit(report)
        report = block_fixture()
        report["all_registered_blocks_covered"] = False
        with self.assertRaisesRegex(ValueError, "block audit incomplete"):
            e6_level.validate_block_audit(report)


if __name__ == "__main__":
    unittest.main()
