#!/usr/bin/env python3
"""Integrity of the four-qubit package and the format of a census root log.

Integrity is relative to manifest.json: every listed file must be present with its SHA-256 hash, and no
unlisted file may appear in the protected directories. It is not a signature and not a recomputation.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys

sys.dont_write_bytecode = True
PROTECTED = ("src", "data", "tools", "scripts")
SOURCES = {
    "src/check_coverage.py", "src/completion_core.py", "src/ge4_check.c",
    "src/fast_helpers_v3.h", "src/fast_tests.h", "tools/recheck_locks.py",
    *{"src/compact/" + x for x in (
        "build_common_input.py", "build_third_level.py",
        "check_block_cover.py", "check_core_levels.py",
        "check_symplectic.c", "compact_core.py")},
}
DATA = {"data/" + x for x in (
    "roots-n4-all.txt", "lock-n4.txt", "pointdomain-n4.bin",
    "fingerprints2-n4.txt", "atlas-witnesses.json", "phase-forests/manifest.json")}
DATA |= {f"data/phase-forests/class-{i:02d}-even.bin" for i in range(36)}
DATA |= {f"data/phase-forests/class-{i:02d}-odd.bin" for i in range(29, 36)}
LOCK_INPUTS = {"src/ge4_check.c", "src/fast_helpers_v3.h", "src/fast_tests.h",
               "data/roots-n4-all.txt", "data/lock-n4.txt", "data/pointdomain-n4.bin",
               "data/fingerprints2-n4.txt"}
LARGE = {0, 2, 90, 92, 669, 674, 677, 682, 744, 751}
ROOTS = set(range(1020)) - LARGE
REQUIRED = SOURCES | DATA | {"scripts/verify.py", "scripts/audit_block_edges.py"}


class ValidationError(ValueError):
    pass


def need(condition, message):
    if not condition:
        raise ValidationError(message)


def _guard():
    need(not sys.flags.optimize, "optimized Python (-O/-OO) is not supported")


def _path(root, relative):
    need(isinstance(relative, str) and "\\" not in relative, "nonportable path")
    p = PurePosixPath(relative)
    need(not p.is_absolute() and p.parts and not any(x in ("..", ".") for x in p.parts)
         and ":" not in relative, "unsafe path: " + relative)
    result = root.joinpath(*p.parts)
    need(result.resolve().is_relative_to(root.resolve()), "path escapes package: " + relative)
    need(not result.is_symlink(), "symlink not accepted: " + relative)
    return result


def _sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _pairs(pairs):
    d = {}
    for k, v in pairs:
        need(k not in d, "duplicate JSON key: " + k)
        d[k] = v
    return d


def _json(path):
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)


def _equal(value, expected, name):
    need(type(value) is type(expected) and value == expected, "unexpected " + name)


def validate_integrity(root):
    """Require every listed source and data file with its hash, and reject extra protected files."""
    _guard()
    root = Path(root)
    manifest = _json(root / "manifest.json")
    files = manifest.get("files")
    need(isinstance(files, dict), "manifest.files must be a path-to-hash mapping")
    need(REQUIRED <= set(files), "manifest omits required source or data files: "
         + ", ".join(sorted(REQUIRED - set(files))[:8]))
    for name, entry in files.items():
        digest = entry.get("sha256") if isinstance(entry, dict) else entry
        need(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest),
             "invalid manifest hash: " + name)
        file = _path(root, name)
        need(file.is_file(), "missing protected file: " + name)
        need(_sha(file) == digest, "hash mismatch: " + name)
    actual = set()
    for directory in PROTECTED:
        for file in (root / directory).rglob("*"):
            if file.is_file() and "__pycache__" not in file.parts:
                actual.add(file.relative_to(root).as_posix())
    need(actual <= set(files), "unexpected protected files: "
         + ", ".join(sorted(actual - set(files))[:8]))
    return {"checked_files": len(files), "kind": "manifest integrity, not recomputation"}


def _root_log(text, rid, plane, level, sector):
    """Check one census root log against the plane, level and sector listed in data/lock-n4.txt."""
    lines = text.splitlines()
    need(len(lines) == 1 and lines[0].startswith("ROOT "), "partial/extra root log")
    need(not any(s in text for s in ("INCOMPLETE", "FIBER_TOO_LARGE", "FAIL")), "failed root log")
    tokens = lines[0].split()
    for key, value in {"ROOT": str(rid), "level": str(level), "sector": str(sector),
                       "survivors": "1", "skel_checked": "374580720",
                       "incomplete_families": "0", "transport_failures": "0",
                       "ray_mismatch": "0", "kernel": "publication-lowtier-v1",
                       "scope": "full", "guided": "1", "charts": "0:2294"}.items():
        need(tokens.count(key) == 1, "missing/duplicate log field " + key)
        i = tokens.index(key)
        need(i + 1 < len(tokens) and tokens[i + 1] == value, "root log mismatch: " + key)
    _equal(re.findall(r"\(\d+,\d+,\d+\)", text), [plane], "proposed plane")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        print(json.dumps({"status": "VALIDATED", "integrity": validate_integrity(args.root)}, indent=2))
    except (ValidationError, OSError, ValueError, KeyError, TypeError) as exc:
        print(json.dumps({"status": "REJECTED", "error": str(exc)}), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
