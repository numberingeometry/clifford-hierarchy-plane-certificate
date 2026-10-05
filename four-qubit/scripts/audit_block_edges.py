"""Replay every block-cover edge with its exact scalar and transported mark.

The workspace must be a disposable certificate copy with src/ and evidence/.
The underlying checker rewrites evidence/compact/block-cover.json there.
This script also writes a new audit record inside that workspace; never point
it at the package itself.
"""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy
import sys
import time

if not __debug__:
    raise SystemExit("Run without Python -O: assertions are verification checks")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path,
                        help="disposable certificate root containing src/ and evidence/")
    parser.add_argument("--output", default="evidence/compact/block-edge-audit.json",
                        help="report path relative to, and contained in, the workspace")
    args = parser.parse_args()
    workspace = args.workspace.resolve(strict=True)
    if not workspace.is_dir():
        parser.error("workspace is not a directory")
    if workspace == Path(__file__).resolve().parents[1]:
        parser.error("refusing to run in this script's canonical package tree")
    if Path(args.output).is_absolute():
        parser.error("output must be a path relative to the workspace")
    output = (workspace / args.output).resolve()
    try:
        output.relative_to(workspace)
    except ValueError:
        parser.error("output escapes the workspace")
    checker = workspace / "src/compact/check_block_cover.py"
    required = [checker, workspace / "src/compact/compact_core.py",
                workspace / "src/completion_core.py",
                workspace / "src/compact/check_core_levels.py",
                workspace / "evidence/compact/core-levels.json"]
    for path in required:
        if not path.is_file():
            parser.error("missing required input: " + path.relative_to(workspace).as_posix())
    bindings = {path.relative_to(workspace).as_posix(): sha(path) for path in required}
    start = time.perf_counter()
    saved_argv = sys.argv
    try:
        sys.argv = [str(checker)]
        with contextlib.redirect_stdout(io.StringIO()):
            data = runpy.run_path(str(checker))
    finally:
        sys.argv = saved_argv
    q = data["q"]
    params, states, moves = data["PARAMS"], data["STATE"], data["MOVES"]
    mat, multiply, unmat = data["mat_of_mono"], data["matmul"], data["mono_of_mat2"]

    def scale(unitary, exponent):
        return unitary[0], tuple((e + exponent) % 8 for e in unitary[1])

    count = odd = changed = 0
    for index, a in enumerate(states):
        p = params[index][5]
        assert all(a[1][n] % 2 == (n & p).bit_count() % 2 for n in range(8))
        for move in moves:
            if move["kind"] == "monomial":
                ci = (q.inv(data["sp_mono"](p)) if move["name"] == "root_Sp"
                      else move["Ci"])
                image = q.mul(q.mul(move["Co"], a), ci)
            else:
                image = unmat(multiply(multiply(move["Co2"], mat(a)), move["Ci2"]))
                if image is None:
                    continue
            canonical = data["canon"](image)
            if canonical is None:
                continue
            target, alpha, s = canonical
            b, new_p = states[target], params[target][5]
            residue = q.mul(q.inv(b), q.mul(alpha, image))
            assert residue in data["ABSORB"]
            diagonal = scale(residue, -s)
            assert diagonal[1][0] == 0 and all(e % 2 == 0 for e in diagonal[1])
            assert q.symp(diagonal) is not None
            # L = alpha Co; R = Ci diagonal^-1, with s retained exactly.
            assert q.mul(q.mul(alpha, image), q.inv(diagonal)) == scale(b, s)
            if move["kind"] == "monomial":
                right = q.mul(ci, q.inv(diagonal))
                transported = q.conj(q.inv(right), q.pauli(0, p))
            else:
                # Co2 = sqrt(2) Ci^dagger for each target-conjugation move.
                transported = unmat(multiply(
                    multiply(mat(diagonal), multiply(
                        multiply(move["Co2"], mat(q.pauli(0, p))), move["Ci2"])),
                    mat(q.inv(diagonal))))
            assert transported == scale(q.pauli(0, new_p), 4 * (s % 2)), (
                params[index], move["name"], params[target], s)
            count += 1
            odd += s % 2
            changed += p != new_p
    assert len(states) == 512 and len(data["ORBS"]) == len(data["CORE_INDEX"]) == 30
    assert (count, odd, changed) == (52320, 768, 36096)
    assert all(sha(workspace / name) == digest for name, digest in bindings.items())
    record = {
        "status": "PASS",
        "scope": "all applicable block-cover edges; exact scalar and mark interface",
        "forms": len(states),
        "applicable_edges_including_self": count,
        "odd_scalar_edges": odd,
        "mark_coordinate_changes": changed,
        "orbit_count": len(data["ORBS"]),
        "listed_cores": len(data["CORE_INDEX"]),
        "exact_identity": "L A R = zeta8^s A_next",
        "mark_identity": "R^dagger Z^p R = (-1)^s Z^p_next",
        "inputs": bindings,
        "audit_source_sha256": sha(Path(__file__)),
        "seconds": round(time.perf_counter() - start, 3),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
