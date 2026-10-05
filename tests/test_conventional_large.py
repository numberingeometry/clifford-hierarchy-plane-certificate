from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
"""Fail-closed controls for C2 gating and fresh C6 outputs; no census rerun."""
from copy import deepcopy
import importlib.util
import os
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

try:
    from checking import conventional_large as checker
except ImportError:
    import conventional_large as checker


class GateControls(unittest.TestCase):
    def test_explicit_zig_executable_and_cc_get_cc_subcommand(self):
        for executable in ("zig", "zig.exe", str(Path("example-toolchain") / "zig.exe")):
            self.assertEqual(checker.choose_compiler(executable), [executable, "cc"])
            with patch.dict(os.environ, {"CC": executable}):
                self.assertEqual(checker.choose_compiler(None), [executable, "cc"])
        argv = ["custom-compiler", "--tool-prefix"]
        self.assertEqual(checker.choose_compiler(argv), argv)
        self.assertEqual(checker.choose_compiler("clang"), ["clang"])

    def test_unflagged_c2_cannot_start_or_claim_full_acceptance(self):
        with patch.object(checker, "isolated_package", side_effect=AssertionError("execution started")):
            record = checker.run("C2", mode="full", allow_expensive=False)
        self.assertEqual(record["status"], "NOT_RUN")
        self.assertFalse(record["accepted"])
        self.assertFalse(record["fresh_mathematical_acceptance"])
        self.assertGreater(record["observed_seconds"], 29000)


class FreshC2Bindings(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="paper-i-c2-binding-test-")
        self.addCleanup(temporary.cleanup)
        self.package = Path(temporary.name) / "package"
        self.run = Path(temporary.name) / "run"
        selected = {
            "src/ge4_check.c": b"selected exact kernel\n",
            "src/fast_helpers_v3.h": b"selected helpers\n",
            "data/roots-n4-all.txt": b"selected complete roots\n",
            "tools/recheck_locks.py": b"selected exact replay tool\n",
        }
        for relative, payload in selected.items():
            path = self.package / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            target = self.run / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
        checker.write_json(self.package / "manifest.json", {"files": {
            relative: checker.sha(self.package / relative) for relative in selected}})
        self.verifier = SimpleNamespace(LOCK_INPUTS={
            "src/ge4_check.c", "src/fast_helpers_v3.h", "data/roots-n4-all.txt"})
        self.report = {
            "source_sha256": {path.relative_to(self.package).as_posix(): checker.sha(path)
                              for path in (self.package / "src").rglob("*") if path.is_file()},
            "package_inputs": {relative: checker.sha(self.package / relative) for relative in selected},
            "workspace_inputs": {path.relative_to(self.run).as_posix(): checker.sha(path)
                                 for path in self.run.rglob("*") if path.is_file()},
        }

    def validate(self, report=None):
        return checker.validate_fresh_lock_input_bindings(
            self.package, self.run, self.verifier, self.report if report is None else report)

    def test_original_selected_package_bindings_pass(self):
        self.assertEqual(self.validate(), {
            relative: checker.sha(self.package / relative) for relative in self.verifier.LOCK_INPUTS})

    def test_self_consistent_modified_fresh_source_rejected(self):
        relative = "src/ge4_check.c"
        (self.run / relative).write_bytes(b"different arithmetic\n")
        # A kernel and fresh reports replaced consistently after execution
        # must still be rejected against the package.
        for field in ("source_sha256", "package_inputs", "workspace_inputs"):
            self.report[field][relative] = checker.sha(self.run / relative)
        with self.assertRaisesRegex(ValueError, "fresh lock inputs differ"):
            self.validate()

    def test_modified_actual_tool_rejected_even_with_original_snapshot(self):
        (self.run / "tools/recheck_locks.py").write_bytes(b"different acceptance tool\n")
        with self.assertRaisesRegex(ValueError, "fresh workspace inputs changed"):
            self.validate()

    def test_modified_tool_snapshot_rejected(self):
        relative = "tools/recheck_locks.py"
        (self.run / relative).write_bytes(b"different acceptance tool\n")
        self.report["workspace_inputs"][relative] = checker.sha(self.run / relative)
        with self.assertRaisesRegex(ValueError, "census workspace_inputs binding"):
            self.validate()

    def test_missing_extra_or_changed_report_bindings_rejected(self):
        for field in ("source_sha256", "package_inputs", "workspace_inputs"):
            for operation in ("missing", "extra", "changed"):
                report = deepcopy(self.report)
                relative = next(iter(report[field]))
                if operation == "missing":
                    del report[field][relative]
                elif operation == "extra":
                    report[field]["unknown-input"] = "0" * 64
                else:
                    report[field][relative] = "0" * 64
                with self.subTest(field=field, operation=operation), self.assertRaisesRegex(
                        ValueError, "census " + field + " binding"):
                    self.validate(report)

    def test_changed_selected_package_rejected_against_manifest(self):
        relative = "tools/recheck_locks.py"
        replacement = b"different selected tool\n"
        (self.package / relative).write_bytes(replacement)
        (self.run / relative).write_bytes(replacement)
        for field in ("package_inputs", "workspace_inputs"):
            self.report[field][relative] = checker.sha(self.package / relative)
        with self.assertRaisesRegex(ValueError, "package input changed"):
            self.validate()


class FreshC6Controls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture = os.environ.get("C6_TEST_WORKSPACE")
        if not fixture:
            raise unittest.SkipTest("Set C6_TEST_WORKSPACE to a fresh C6 execution workspace")
        cls.work = Path(fixture)
        cls.result = checker.load_json(cls.work / "symplectic-results.json")
        cls.common = checker.load_json(cls.work / "conventional/evidence/compact/common-input.json")
        cls.third = checker.load_json(cls.work / "conventional/evidence/compact/third-level.json")
        cls.bindings = checker.load_json(cls.work / "physical-bindings.json")
        spec = importlib.util.spec_from_file_location("c6_test_exact_kernel",
                    cls.work / "conventional/src/compact/compact_core.py")
        cls.q = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.q)

    def test_wrong_family_allocation_rejected_despite_correct_global_totals(self):
        record = deepcopy(self.result)
        # Transfer one qualifying common disposition between families. All
        # global totals and per-record partition identities remain unchanged.
        first = next(row for row in record["cases"][:360] if row["common"] < 1451520)
        second = next(row for row in record["cases"][360:] if row["common"] > 0)
        first["common"] += 1
        first["unipotence_pass"] += 1
        second["common"] -= 1
        second["unipotence_pass"] -= 1
        with self.assertRaisesRegex(ValueError, "family acceptance allocation"):
            checker.validate_c6_result(record)

    def test_omitted_duplicate_or_filtered_record_rejected(self):
        record = deepcopy(self.result)
        record["cases"].pop()
        with self.assertRaises(ValueError):
            checker.validate_c6_result(record)
        record = deepcopy(self.result)
        record["cases"][1]["id"] = record["cases"][0]["id"]
        with self.assertRaises(ValueError):
            checker.validate_c6_result(record)
        record = deepcopy(self.result)
        record["cases"][0]["stationarity_excluded"] = 1
        with self.assertRaises(ValueError):
            checker.validate_c6_result(record)

    def test_malformed_counts_rejected(self):
        for wrong in (-1, True, 1.0, "1"):
            record = deepcopy(self.result)
            record["cases"][0]["bad"] = wrong
            with self.subTest(wrong=wrong), self.assertRaises(ValueError):
                checker.validate_c6_result(record)

    def test_noninteger_header_or_zero_filter_count_rejected(self):
        record = deepcopy(self.result)
        record["symplectic_matrices"] = 1451520.0
        with self.assertRaises(ValueError):
            checker.validate_c6_result(record)
        record = deepcopy(self.result)
        record["cases"][0]["balanced"] = False
        with self.assertRaises(ValueError):
            checker.validate_c6_result(record)

    def test_corrupted_binary_or_json_binding_rejected(self):
        cases = self.common["cases"] + self.third["cases"]
        payload = (self.work / "merged-input.bin").read_bytes()
        for malformed in (payload[:-1], payload + b"x"):
            with self.assertRaises(ValueError):
                checker.parse_binary(malformed, cases, self.q)
        malformed = bytearray(payload)
        malformed[1493 + 20] = 0  # first identity test has a zero first column
        with self.assertRaises(ValueError):
            checker.parse_binary(malformed, cases, self.q)
        changed = deepcopy(cases)
        changed[0]["tests"][0][0] = 2
        with self.assertRaises(ValueError):
            checker.parse_binary(payload, changed, self.q)

    def test_corrupted_physical_endpoint_rejected_by_exact_reconstruction(self):
        index = next(i for i, row in enumerate(self.common["cases"]) if row["returns"])
        record = deepcopy(self.common["cases"][index])
        binding = self.bindings["common_input"][index]
        a = tuple(tuple(values) for values in binding["A"])
        b = tuple(tuple(values) for values in binding["B"])
        checker.validate_common_physical_record(record, self.q, a, b)
        record["returns"][0][0] ^= 1
        with self.assertRaisesRegex(ValueError, "physical common-block endpoints"):
            checker.validate_common_physical_record(record, self.q, a, b)


if __name__ == "__main__":
    unittest.main()
