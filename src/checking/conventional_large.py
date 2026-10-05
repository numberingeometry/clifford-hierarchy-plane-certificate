"""Facts (C7) and (C15): exact recomputation with the compiled C kernels of the four-qubit package.

(C7) recomputes all 1,010 gates of the census and must be requested with --allow-expensive; it is never
sampled. (C15) regenerates the pair domains, binds every record to the gates it describes, and enumerates
Sp(6,2). Keys: C2 = (C7), C6 = (C15).
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import runpy
import shutil
import struct
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
CENSUS_SHA256 = "0f57793b16c11879417f3a9816b672c276abf62f451217f1b23b882229bb142a"


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            need(key not in result, "duplicate JSON key: " + key)
            result[key] = value
        return result
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def owner_path(owner_root):
    return (Path(owner_root) if owner_root is not None else
            Path(__file__).resolve().parents[2] / "four-qubit")


def availability(code, owner_root=None):
    need(code in ("C2", "C6"), "supported large conventional facts are C2 and C6")
    owner = owner_path(owner_root)
    needed = ["manifest.json", "scripts/verify.py", "src/ge4_check.c"] if code == "C2" else [
        "manifest.json", "src/compact/build_common_input.py", "src/compact/build_third_level.py",
        "src/compact/check_symplectic.c", "src/compact/compact_core.py", "src/completion_core.py"]
    return {"status": "NOT_RUN", "accepted": False, "route": "recompute", "fact": code, "exact": True,
            "availability": "available" if all((owner / rel).is_file() for rel in needed) else "unavailable",
            "fresh_mathematical_acceptance": False,
            "expected_domain": {"nondirect_roots": 1010, "charts": 2295,
                                "point_permutations": 163216, "all_phase_families": True}
                               if code == "C2" else
                               {"common_input_records": 360, "third_level_records": 2076,
                                "symplectic_matrices": 1451520},
            "observed_seconds": 29509.2396364 if code == "C2" else 38.163903,
            "resource_note": "Observed times on one worker, not a runtime promise. The C6 time is the group "
                             "test alone, without the construction of the domains.",
            "full_execution_requires_flag": code == "C2",
            "command": "python -m checking.conventional_large " + code +
                       (" --allow-expensive --workers 1" if code == "C2" else ""),
            "analytic_domain_reduction": "paper premise; not established by finite execution"}


def choose_compiler(compiler):
    def executable_argv(value):
        executable = str(value)
        return [executable, "cc"] if Path(executable).name.lower() in ("zig", "zig.exe") else [executable]
    if compiler:
        return [str(item) for item in compiler] if isinstance(compiler, (tuple, list)) else executable_argv(compiler)
    if os.environ.get("CC"):
        return executable_argv(os.environ["CC"])
    for name in ("clang", "gcc", "zig"):
        if shutil.which(name):
            return [name, "cc"] if name == "zig" else [name]
    if importlib.util.find_spec("ziglang"):
        return [sys.executable, "-m", "ziglang", "cc"]
    raise RuntimeError("Install Clang, GCC or ziglang, set CC, or supply compiler explicitly")


def isolated_package(owner, work_directory):
    need(not sys.flags.optimize, "Run without -O/-OO: preserved exact sources use assertions")
    need(sys.byteorder == "little", "the binary kernel requires a little-endian host")
    manifest = load_json(owner / "manifest.json")
    files = manifest.get("files")
    need(isinstance(files, dict) and files, "missing package manifest")
    if work_directory is None:
        work = Path(tempfile.mkdtemp(prefix="paper-i-conventional-"))
    else:
        work = Path(work_directory)
        work.mkdir(parents=True, exist_ok=False)
    package = work / "conventional"
    package.mkdir()
    for relative, digest in files.items():
        need(isinstance(relative, str) and "\\" not in relative and ":" not in relative
             and not Path(relative).is_absolute() and ".." not in Path(relative).parts,
             "unsafe manifest path")
        source = owner / relative
        expected = digest.get("sha256") if isinstance(digest, dict) else digest
        need(sha(source) == expected, "package input changed: " + relative)
        target = package / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    shutil.copyfile(owner / "manifest.json", package / "manifest.json")
    (package / "census").mkdir()
    shutil.copyfile(owner / "census/replay.py", package / "census/replay.py")
    need(sha(package / "census/replay.py") == CENSUS_SHA256, "census runner changed")
    spec = importlib.util.spec_from_file_location("public_conventional_integrity", package / "scripts/verify.py")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    verifier.validate_integrity(package)
    (work / "logs").mkdir()
    return work, package, verifier


def execute(work, label, command, cwd=None, env=None):
    started = time.perf_counter()
    with (work / "logs" / (label + ".txt")).open("w", encoding="utf-8") as stream:
        environment = dict(os.environ, ZIG_GLOBAL_CACHE_DIR=str(work / "tool-cache"),
                           ZIG_LOCAL_CACHE_DIR=str(work / "local-tool-cache"), PYTHONDONTWRITEBYTECODE="1")
        if env is not None:
            environment.update(env)
        process = subprocess.run(command, cwd=cwd or work, env=environment, stdout=stream,
                                 stderr=subprocess.STDOUT, text=True)
    need(process.returncode == 0, "execution failed: " + label)
    return {"check": label, "exit": process.returncode, "seconds": time.perf_counter() - started}


def compile_abi(work, cc, steps):
    source = work / "abi.c"
    source.write_text("#include <stdint.h>\n#include <limits.h>\n"
                      "int main(void){unsigned x=1;return CHAR_BIT!=8 || sizeof(unsigned)!=4 || "
                      "sizeof(uint64_t)!=8 || *(unsigned char*)&x!=1;}\n", encoding="ascii")
    executable = work / ("abi.exe" if os.name == "nt" else "abi")
    steps.append(execute(work, "compiler-version", cc + ["--version"]))
    steps.append(execute(work, "compile-abi", cc + [str(source), "-o", str(executable)]))
    steps.append(execute(work, "binary-abi", [str(executable)]))


def validate_fresh_lock_input_bindings(package, run, verifier, report):
    """Bind a C2 census to the package, including its copied tools."""
    inputs = {relative: sha(run / relative) for relative in verifier.LOCK_INPUTS}
    need(inputs == {relative: sha(package / relative) for relative in verifier.LOCK_INPUTS},
         "fresh lock inputs differ from the package")
    source = {path.relative_to(package).as_posix(): sha(path)
              for path in (package / "src").rglob("*")
              if path.is_file() and "__pycache__" not in path.parts}
    need(report.get("source_sha256") == source, "census source_sha256 binding")
    files = load_json(package / "manifest.json")["files"]
    package_inputs = {relative: sha(package / relative) for relative in files}
    for relative, digest in files.items():
        expected = digest.get("sha256") if isinstance(digest, dict) else digest
        need(package_inputs[relative] == expected, "package input changed: " + relative)
    need(report.get("package_inputs") == package_inputs, "census package_inputs binding")
    # The census runner snapshots exactly these initial workspace files.
    # Extra fresh outputs (binaries, logs and evidence) are checked separately.
    workspace = {}
    for directory in ("src", "data", "tools"):
        for path in (package / directory).rglob("*"):
            if (path.is_file() and "__pycache__" not in path.parts
                    and not path.name.endswith(".pyc")):
                workspace[path.relative_to(package).as_posix()] = sha(path)
    need(report.get("workspace_inputs") == workspace, "census workspace_inputs binding")
    need({relative: sha(run / relative) for relative in workspace} == workspace,
         "fresh workspace inputs changed after snapshot")
    return inputs


def validate_fresh_locks(package, run, verifier):
    report = load_json(run / "run.json")
    need(report.get("status") == "PASS" and report.get("mode") == "locks"
         and report.get("fresh_all_roots") is True
         and report.get("reused_canonical_records") is False, "not a fresh full lock census")
    need(report.get("runner_sha256") == CENSUS_SHA256
         and report.get("manifest_sha256") == sha(package / "manifest.json"), "census source binding")
    steps = report.get("steps", [])
    expected_steps = ["compiler-version", "compile-abi", "binary-abi", "compile-lock", "coverage",
                      "selftest", "testfamily", "fasttests", "testtiers", "all-root-locks"]
    need([step.get("check") for step in steps] == expected_steps
         and all(step.get("exit") == 0 for step in steps), "incomplete exact lock execution")
    result = load_json(run / "evidence/recheck-all.json")
    need(result.get("status") == "PASS" and result.get("scope") == "all"
         and result.get("roots") == sorted(verifier.ROOTS)
         and result.get("new_roots_completed") == 1010
         and result.get("cached_roots_reused") == 0 and result.get("failures") == [],
         "incomplete, failed or cached root census")
    inputs = validate_fresh_lock_input_bindings(package, run, verifier, report)
    binary = run / "src" / ("ge4_check.exe" if os.name == "nt" else "ge4_check")
    need(result.get("inputs") == inputs and result.get("binary_sha256") == sha(binary),
         "root census inputs changed")
    lock_rows = (run / "data/lock-n4.txt").read_text(encoding="utf-8").splitlines()
    for root in sorted(verifier.ROOTS):
        base = run / "evidence/recheck" / ("root-%04d" % root)
        record = load_json(base.with_suffix(".json"))
        need(record.get("root") == root and record.get("status") == "PASS"
             and record.get("inputs") == inputs and record.get("binary_sha256") == sha(binary)
             and record.get("result_sha256") == sha(base.with_suffix(".txt")), "root binding failure")
        fields = lock_rows[root].split()
        verifier._root_log(base.with_suffix(".txt").read_text(encoding="utf-8"), root,
                           fields[-1], int(fields[2]), int(fields[4]))
    return report


def c2(work, package, verifier, cc, workers, mode):
    need(1 <= workers <= 64, "workers must be between one and sixty-four")
    if mode == "full":
        # The census runner is loaded unchanged; this driver supplies only the
        # compiler command and the task arguments.
        spec = importlib.util.spec_from_file_location("c2_census_runner", package / "census/replay.py")
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        runner.compiler = lambda: list(cc)
        previous_argv = list(sys.argv)
        previous_env = dict(os.environ)
        sys.argv = [str(package / "census/replay.py"), "locks", "--workers", str(workers), "--directory", "fresh-c2"]
        os.environ.update(PYTHONUTF8="1", PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1",
                          ZIG_GLOBAL_CACHE_DIR=str(work / "tool-cache"), ZIG_LOCAL_CACHE_DIR=str(work / "local-tool-cache"))
        try:
            with (work / "logs/complete-root-census.txt").open("w", encoding="utf-8") as log:
                with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                    runner.main()
        finally:
            sys.argv[:] = previous_argv
            os.environ.clear()
            os.environ.update(previous_env)
        report = validate_fresh_locks(package, package / "runs/fresh-c2", verifier)
        return {"status": "PASS", "accepted": True, "route": "recompute", "fact": "C2", "exact": True,
                "fresh_mathematical_acceptance": True, "fresh_roots": 1010,
                "cached_roots_reused": 0, "workers": workers,
                "seconds": report["elapsed_seconds"], "all_phase_families": True,
                "analytic_domain_reduction": "paper premise; not established by finite execution"}
    steps = []
    compile_abi(work, cc, steps)
    executable = work / ("ge4_check.exe" if os.name == "nt" else "ge4_check")
    steps.append(execute(work, "compile-lock", cc + ["-O2", "-std=c99", "-DFAST_TESTS",
                         "-o", str(executable), str(package / "src/ge4_check.c")]))
    steps.append(execute(work, "selftest", [str(executable), "selftest"], cwd=package / "data"))
    steps.append(execute(work, "testfamily", [str(executable), "testfamily"], cwd=package / "data"))
    need("selftest: lagrangians 2295 planes 10795 cliffords ok" in
         (work / "logs/selftest.txt").read_text(encoding="utf-8"), "missing mathematical selftest verdict")
    need("FAMILY TEST PASSED: repaired invariant correct on all 8 cases, legacy accumulator wrong where expected"
         in (work / "logs/testfamily.txt").read_text(encoding="utf-8"), "universal-family regression failed")
    lock_rows = (package / "data/lock-n4.txt").read_text(encoding="utf-8").splitlines()
    for root in (18, 400):
        output = work / ("sample-root-%04d.txt" % root)
        steps.append(execute(work, "sample-root-" + str(root), [str(executable), "lock", "roots-n4-all.txt",
                             "pointdomain-n4.bin", str(output), str(root), str(root),
                             "fingerprints2-n4.txt", "rigid3", "inside", "guided"], cwd=package / "data"))
        fields = lock_rows[root].split()
        verifier._root_log(output.read_text(encoding="utf-8"), root, fields[-1], int(fields[2]), int(fields[4]))
    return {"status": "INTEGRATION_PASS", "accepted": False, "route": "recompute", "fact": "C2", "exact": True,
            "fresh_mathematical_acceptance": False, "sample_roots": [18, 400],
            "pending_roots": 1008, "all_phase_families": "only the two sampled roots",
            "steps": steps, "kernel_binary_sha256": sha(executable)}


def normalized(value):
    return json.loads(json.dumps(value))


def parse_binary(payload, expected_cases, q):
    position = 0

    def take(fmt):
        nonlocal position
        size = struct.calcsize(fmt)
        need(position + size <= len(payload), "truncated binary")
        value = struct.unpack_from(fmt, payload, position)
        position += size
        return value

    def symplectic(columns):
        need(len(columns) == 6 and all(0 < value < 64 for value in columns)
             and all(q.om(columns[i], columns[j]) == q.om(1 << i, 1 << j)
                     for i in range(6) for j in range(6)), "non-symplectic binary test")

    need(take("<I")[0] == 135, "complete Lagrangian header count")
    for mask, basis in zip(q.LMASK, q.BASES):
        encoded_mask, *encoded_basis = take("<Q3B")
        need(encoded_mask == mask and tuple(encoded_basis) == basis, "Lagrangian header semantics")
    need(take("<I")[0] == len(expected_cases), "binary case coverage")
    tests = 0
    for index, expected in enumerate(expected_cases):
        nt, ns, ne, nr, nb = take("<5I")
        need(0 < nt <= 4096 and ns <= 64 and ne <= 135 and nr <= 4096 and nb <= 4096,
             "binary record sizes exceed the format")
        record = {"tests": [list(take("<6B")) for _ in range(nt)],
                  "station": [], "edges": [], "returns": [], "balanced": []}
        need(record["tests"][0] == [1, 2, 4, 8, 16, 32], "missing identity test")
        for matrix in record["tests"]:
            symplectic(matrix)
        tests += nt
        for _ in range(ns):
            matrix = list(take("<6B"))
            symplectic(matrix)
            count = take("<I")[0]
            need(count <= 135, "stationary index count")
            indices = list(take("<" + "I" * count))
            need(all(index < 135 for index in indices) and len(indices) == len(set(indices)),
                 "stationary indices")
            record["station"].append([matrix, indices])
        for _ in range(ne):
            mask, *basis = take("<Q3B")
            need(mask in q.LMASK and len(q.span(basis)) == 8
                 and all(0 < value < 64 for value in basis)
                 and all(not q.om(a, b) for a in basis for b in basis), "binary physical edge")
            record["edges"].append([mask, basis])
        for name, count in (("returns", nr), ("balanced", nb)):
            for _ in range(count):
                pair = list(take("<2B"))
                need(all(0 < value < 64 for value in pair), "binary return label")
                record[name].append(pair)
        need(all(record[key] == expected[key] for key in record),
             "binary differs from fresh physical constructor: " + str(index))
    need(position == len(payload), "trailing binary data")
    return tests


def validate_common_physical_record(record, q, a, b):
    mixed = [q.mul(q.mul(b, pauli), q.inv(a)) for pauli in q.PAULIS]
    decompositions = [q.decompose(operator) for operator in mixed]
    need(all(item is not None and q.is_c3(operator)
             for item, operator in zip(decompositions, mixed)), "non-third-level physical mixed block")
    g0i = q.inverse(decompositions[0][0])
    tests, station = set(), set()
    for operator, (g, core) in zip(mixed, decompositions):
        need(q.symp(q.mul(operator, q.inv(core))) == g, "exact physical decomposition binding")
        gg = q.prod(g0i, g)
        tests.update(q.prod(gg, e) for e in q.STAR_E[core])
        station.add((gg, tuple(i for i, j in q.STAR_CHARTS[core])))
    ordered = sorted(tests)
    ordered.remove(q.I)
    ordered.insert(0, q.I)
    edges, returns, balanced = q.witnesses(a, b)
    physical = {"tests": ordered, "station": sorted(station),
                "edges": [(q.LMASK[i], tuple(q.app(g0i, v) for v in q.BASES[j])) for i, j in edges],
                "returns": sorted(set((i, q.app(g0i, j)) for i, j in returns)),
                "balanced": sorted(set((i, q.app(g0i, j)) for i, j in balanced))}
    need(all(normalized(value) == record[key] for key, value in physical.items()),
         "physical common-block endpoints differ from binary obligation")
    return {"name": record["name"], "A": normalized(a), "B": normalized(b),
            "normalization_symplectic": list(g0i), "mixed_blocks_checked": 64,
            "edges": record["edges"], "returns": record["returns"]}


def physical_bindings(common, third, common_namespace, third_namespace):
    q = common_namespace["q"]
    inputs = {name: (a, mark) for name, a, mark in common_namespace["alist"]}
    common_bindings = []
    for record in common["cases"]:
        name, parameters = record["name"].rsplit("_", 1)
        need(name in inputs and len(parameters) == 5, "common physical name binding")
        a, mark = inputs[name]
        perm, x, y, z, d = map(int, parameters)
        need(perm in range(8) and x in range(2) and y in range(2) and z in range(2) and d in range(4),
             "common physical parameter domain")
        b = (q.POINTS[perm], q.mono(0, x, y, z, d, mark)[1])
        common_bindings.append(validate_common_physical_record(record, q, a, b))
    # Export and check one explicit symplectic conjugator for every generated
    # correction-image group; lifting it to a Clifford gate (up to a scalar) is the
    # standard correspondence between Sp(6,2) and the Clifford group used in the paper.
    q = third_namespace["q"]
    frame_witnesses = {}
    for name, core in third_namespace["bases"]:
        corrections = q.m.group(core)
        for lift in third_namespace["lifts"] if name != "I" else [q.I]:
            li = q.inverse(lift)
            for linear, inverse in zip(third_namespace["linears"], third_namespace["lininv"]) if name != "I" else [(q.I, q.I)]:
                frame = q.prod(lift, linear)
                fi = q.prod(inverse, li)
                group = tuple(sorted(q.prod(q.prod(frame, e), fi) for e in corrections))
                frame_witnesses.setdefault(group, {"core": name, "frame": list(frame)})
    need(len(frame_witnesses) == 1396, "physical conjugate-group witness coverage")
    definitions = third["conjugate_groups"]
    for definition in definitions:
        group = tuple(tuple(matrix) for matrix in definition["group"])
        need(group in frame_witnesses and frame_witnesses[group]["core"] == definition["kind"],
             "conjugate-group origin binding")
        need(list(third_namespace["invariant_charts"](group)) == definition["charts"]
             and sorted(third_namespace["fixed"](group)) == definition["fixed"],
             "physical conjugate-group charts or fixed axes changed")
    cores = dict(third_namespace["alist"])
    third_bindings = []
    for record in third["cases"]:
        name, index, kind = record["name"].split("_", 2)
        a = cores[name]
        definition = definitions[int(index)]
        need(definition["kind"] == kind, "third-level physical core binding")
        ea = q.m.group(a)
        xa = third_namespace["invariant_charts"](tuple(ea))
        ra = third_namespace["fixed"](ea)
        ed, xd = definition["group"], definition["charts"]
        need(all(set(third_namespace["lag_action"](e)[i] for i in xd) == set(xd) for e in ea),
             "third-level chart invariance hypothesis")
        tests = sorted({q.prod(e, f) for e in ed for f in ea})
        tests.remove(q.I)
        tests.insert(0, q.I)
        physical = {"tests": tests, "station": [(f, xd) for f in ea],
                    "edges": [(q.LMASK[i], q.BASES[i]) for i in set(xa) & set(xd)],
                    "returns": [(r, r) for r in ra & set(definition["fixed"])], "balanced": []}
        need(all(normalized(value) == record[key] for key, value in physical.items()),
             "physical third-level endpoints differ from binary obligation")
        rr = q.rows(a)
        need(all(q.one_label(rr[r]) == r for r in ra), "third-level returned axis operator meaning")
        third_bindings.append({"name": record["name"], "A": normalized(a),
                               "other_group_origin": frame_witnesses[tuple(tuple(e) for e in ed)],
                               "edges": record["edges"], "returns": record["returns"]})
    return q, {"common_input": common_bindings, "third_level": third_bindings,
               "conjugate_group_witnesses": [dict(frame_witnesses[tuple(tuple(e) for e in definition["group"])],
                                                    group=definition["group"]) for definition in definitions]}


def validate_c6_result(result):
    need(all(type(result.get(field)) is int for field in ("symplectic_matrices", "unipotent", "bad")),
         "binary census count types")
    need(result.get("symplectic_matrices") == 1451520 and result.get("unipotent") == 262144
         and result.get("bad") == 0 and len(result.get("cases", [])) == 2436, "complete binary group domain")
    totals = {name: 0 for name in ("unipotence_pass", "common", "returned", "bad")}
    for index, record in enumerate(result["cases"]):
        need(type(record.get("id")) is int
             and all(type(record.get(field)) is int for field in
                     ("unipotence_pass", "common", "returned", "balanced", "stationarity_excluded", "bad")),
             "binary record count types")
        need(record.get("id") == index and record.get("bad") == 0
             and record.get("balanced") == record.get("stationarity_excluded") == 0
             and record.get("first") == [0] * 6
             and record.get("common", -1) + record.get("returned", -1) == record.get("unipotence_pass"),
             "binary implication failure or undeclared filter at record " + str(index))
        for name in totals:
            need(type(record[name]) is int and 0 <= record[name] <= 1451520, "binary count type/domain")
            totals[name] += record[name]
    need(totals == {"unipotence_pass": 18903552, "common": 12979712, "returned": 5923840, "bad": 0},
         "complete pair-domain aggregate changed")
    families = {}
    for family, begin, end, expected in (
            ("common-input", 0, 360, {"unipotence_pass": 572416, "common": 565248, "returned": 7168, "bad": 0}),
            ("third-level", 360, 2436, {"unipotence_pass": 18331136, "common": 12414464, "returned": 5916672, "bad": 0})):
        actual = {name: sum(row[name] for row in result["cases"][begin:end]) for name in totals}
        need(actual == expected, "physical family acceptance allocation changed: " + family)
        families[family] = dict(actual, records=end-begin)
    return totals, families



def c6(work, package, cc):
    started = time.perf_counter()
    steps = []
    compile_abi(work, cc, steps)
    evidence = package / "evidence/compact"
    evidence.mkdir(parents=True)
    original_path = list(sys.path)
    old_module = sys.modules.pop("compact_core", None)
    namespaces = []
    sys.path.insert(0, str(package / "src/compact"))
    try:
        for name in ("build_common_input", "build_third_level"):
            tick = time.perf_counter()
            with (work / "logs" / (name + ".txt")).open("w", encoding="utf-8") as log:
                with contextlib.redirect_stdout(log):
                    namespaces.append(runpy.run_path(str(package / "src/compact" / (name + ".py"))))
            steps.append({"check": name, "exit": 0, "seconds": time.perf_counter() - tick})
        common = load_json(evidence / "common-input.json")
        third = load_json(evidence / "third-level.json")
        need(len(common["cases"]) == 360 and len(common["excluded"]) == 4760
             and len(common["cases"]) + len(common["excluded"]) == 5120
             and common["decompositions"] == 13228, "complete common-input physical domain")
        need(len(third["cases"]) == 2076 and len(third["conjugate_groups"]) == 1396
             and third["excluded"] == 11884 and third["excluded"] + len(third["cases"]) == 13960,
             "complete third-level physical domain")
        for data in (common, third):
            need(len({row["name"] for row in data["cases"]}) == len(data["cases"]), "duplicate physical case")
        q, bindings = physical_bindings(common, third, namespaces[0], namespaces[1])
        write_json(work / "physical-bindings.json", bindings)
        a, b = [(evidence / (family + ".bin")).read_bytes() for family in ("common-input", "third-level")]
        parse_binary(a, common["cases"], q)
        parse_binary(b, third["cases"], q)
        need(a[:1489] == b[:1489], "physical Lagrangian headers differ")
        cases = common["cases"] + third["cases"]
        merged = a[:1489] + struct.pack("<I", 2436) + a[1493:] + b[1493:]
        need(len(merged) == 889083 and parse_binary(merged, cases, q) == 50059,
             "merged constructor binding or test domain changed")
    finally:
        sys.path[:] = original_path
        sys.modules.pop("compact_core", None)
        if old_module is not None:
            sys.modules["compact_core"] = old_module
    binary_input = work / "merged-input.bin"
    binary_input.write_bytes(merged)
    executable = work / ("check_symplectic.exe" if os.name == "nt" else "check_symplectic")
    steps.append(execute(work, "compile-symplectic", cc + ["-O3", "-std=c99",
                         str(package / "src/compact/check_symplectic.c"), "-o", str(executable)]))
    for label, malformed in (("truncated", merged[:-1]), ("trailing", merged + b"x")):
        path = work / (label + ".bin")
        path.write_bytes(malformed)
        process = subprocess.run([str(executable), str(path)], capture_output=True, text=True)
        need(process.returncode != 0 and label in process.stderr, "binary parser accepted " + label)
    steps.append(execute(work, "complete-symplectic-census", [str(executable), str(binary_input)]))
    result = load_json(work / "logs/complete-symplectic-census.txt")
    totals, families = validate_c6_result(result)
    write_json(work / "symplectic-results.json", result)
    return {"status": "PASS", "accepted": True, "route": "recompute", "fact": "C6", "exact": True,
            "fresh_mathematical_acceptance": True, "fresh_constructor_records": [360, 2076],
            "physical_binding_records": [len(bindings["common_input"]), len(bindings["third_level"])],
            "conjugate_group_witnesses": len(bindings["conjugate_group_witnesses"]),
            "symplectic_matrices": 1451520, "unipotent": 262144,
            "records": 2436, "tests": 50059, "totals": totals, "per_family": families, "steps": steps,
            "binary_negative_controls": ["truncated", "trailing"],
            "seconds": time.perf_counter() - started, "independent_kernel": False,
            "physical_binding_sha256": sha(work / "physical-bindings.json"),
            "merged_input_sha256": sha(binary_input), "kernel_binary_sha256": sha(executable),
            "analytic_domain_reduction": "paper premise; not established by finite execution"}


def run(code, *, owner_root=None, work_directory=None, compiler=None, workers=1,
        allow_expensive=False, mode="full"):
    info = availability(code, owner_root)
    need(mode in ("availability", "integration", "full"), "unknown execution mode")
    if mode == "availability" or (code == "C2" and mode == "full" and not allow_expensive):
        return info
    need(info["availability"] == "available", "package sources unavailable")
    need(code != "C6" or mode == "full", "C6 integration must use its complete bounded exact domain")
    owner = owner_path(owner_root)
    cc = choose_compiler(compiler)
    work, package, verifier = isolated_package(owner, work_directory)
    record = c2(work, package, verifier, cc, workers, mode) if code == "C2" else c6(work, package, cc)
    verifier.validate_integrity(package)
    record.update({"workspace_name": work.name, "package_manifest_sha256": sha(package / "manifest.json"),
                   "checker_sha256": sha(__file__), "python": sys.version.split()[0],
                   "compiler_invocation_kind": "installed executable" if len(cc) == 1 else "installed argv",
                   "artifact_integrity": "PASS", "formal_verification": False})
    write_json(work / "acceptance.json", record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("code", choices=("C2", "C6"))
    parser.add_argument("--mode", choices=("availability", "integration", "full"), default="full")
    parser.add_argument("--owner-root", type=Path)
    parser.add_argument("--work-directory", type=Path)
    parser.add_argument("--compiler")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--allow-expensive", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args.code, owner_root=args.owner_root, work_directory=args.work_directory,
                         compiler=args.compiler, workers=args.workers, mode=args.mode,
                         allow_expensive=args.allow_expensive), indent=2))


if __name__ == "__main__":
    main()
