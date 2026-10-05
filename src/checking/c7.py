"""Fact (C10): the case rank E = 2 in the proof of Lemma app:return-balanced. Exact recomputation.

Three parts: the cycle type of Z_1 h (and of the equivalent 16-point word), the 512 two-target phase
tables, and the 64 parameter choices with their unipotence conditions. The monomial and chart arithmetic
is that of the module compact_core.py of the four-qubit package, loaded unchanged and checked
against its hash.
Internal key: C7.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True


KERNEL_HASHES = {
    "src/compact/compact_core.py": "b2d52370fdab007d4efe9f711460f6af40ab59d8edc44b0c472b4f3878978856",
    "src/completion_core.py": "75cf2068e4328a04089832a7f256a1a331702d2156b265155548643aa8b6aab7",
}
EXPECTED_SURVIVORS = [
    {"parameters": [0, 0, 0, 0, 0, 0], "common_charts": [15]},
    {"parameters": [0, 0, 0, 1, 0, 1], "common_charts": [20]},
]
EXPECTED_CYCLES = [[0, 1], [2, 3], [4, 9, 8, 5, 12, 13], [6, 11, 10, 7, 14, 15]]


def need(condition, message):
    if not condition:
        raise ValueError(message)


def load_kernel(owner_root):
    need(not sys.flags.optimize, "Run without -O/-OO: the unchanged kernel uses assertions")
    for relative, expected in KERNEL_HASHES.items():
        path = owner_root / relative
        need(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
             "selected exact kernel changed: " + relative)
    path = owner_root / "src/compact/compact_core.py"
    spec = importlib.util.spec_from_file_location("paper_i_c7_compact_core", path)
    need(spec is not None and spec.loader is not None, "cannot load selected exact kernel")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(owner_root=None):
    started = time.perf_counter()
    if owner_root is None:
        owner_root = Path(__file__).resolve().parents[2] / "four-qubit"
    q = load_kernel(Path(owner_root))
    need(len(q.LAGS) == 135 and len(set(q.LAGS)) == 135, "three-qubit Lagrangian domain")
    need(q.LAGS[15] == q.span((8, 2, 4)), "chart 15 must be <Z_c,X_y,X_z>")
    need(q.LAGS[20] == q.span((8, 2, 36)), "chart 20 must be <Z_c,X_y,Y_z>")

    def unipotent(a):
        b = tuple(a[j] ^ (1 << j) for j in range(6))
        for _ in range(3):
            b = tuple(q.app(b, v) for v in b)
        return not any(b)

    def hadamard(i):
        cols = list(q.I)
        cols[i], cols[i + 3] = cols[i + 3], cols[i]
        return tuple(cols)

    def cz(i, j):
        cols = list(q.I)
        cols[i] ^= 1 << (j + 3)
        cols[j] ^= 1 << (i + 3)
        return tuple(cols)

    def cx(i, j):
        cols = list(q.I)
        cols[i] ^= 1 << j
        cols[j + 3] ^= 1 << (i + 3)
        return tuple(cols)

    hh = q.prod(hadamard(0), hadamard(1))
    C = q.prod(q.prod(q.prod(hh, cz(0, 1)), hh), cx(1, 2))
    need(C == (1, 6, 4, 10, 17, 49), "physical C symplectic action")

    # The obstruction as a word on the sixteen basis states: pp X_s pp^-1 X_c X_s.
    pp = []
    for atom in range(16):
        c, s, y, z = ((atom >> j) & 1 for j in range(4))
        zz = z ^ (s & y)
        yy = y ^ (c & zz)
        pp.append(c + 2 * s + 4 * yy + 8 * zz)
    need(sorted(pp) == list(range(16)), "physical routing word is not a permutation")
    inverse = [pp.index(atom) for atom in range(16)]
    obstruction = [pp[inverse[atom ^ 3] ^ 2] for atom in range(16)]
    need(sorted(obstruction) == list(range(16)), "obstruction is not a permutation")
    cycles, seen = [], set()
    for atom in range(16):
        if atom in seen:
            continue
        cycle, next_atom = [], atom
        while next_atom not in seen:
            seen.add(next_atom)
            cycle.append(next_atom)
            next_atom = obstruction[next_atom]
        need(next_atom == atom, "cycle does not close at its first atom")
        cycles.append(cycle)
    need(cycles == EXPECTED_CYCLES, "six-cycle obstruction changed")
    need(sorted(map(len, cycles)) == [2, 2, 6, 6] and seen == set(range(16)),
         "complete six-cycle coverage")

    # The same obstruction in the form stated in the paper. In the basis |+-> (x) |s,y,z>,
    # h = W X_s W^dagger = diag(M, tau M tau) with M of permutation sigma(s,y,z) = (s+1, y, y+z)
    # and tau(s,y,z) = (s, y+z, z); Z_1 exchanges the two blocks. Points are 8*b + (s + 2y + 4z).
    def point(s, y, z):
        return s % 2 + 2 * (y % 2) + 4 * (z % 2)

    triples = [((n & 1), (n >> 1) & 1, (n >> 2) & 1) for n in range(8)]
    sigma = [point(s + 1, y, y + z) for s, y, z in triples]
    tau = [point(s, y + z, z) for s, y, z in triples]
    blocks = (sigma, [tau[sigma[tau[n]]] for n in range(8)])
    z1h = [8 * (1 - (p >> 3)) + blocks[p >> 3][p & 7] for p in range(16)]
    need(sorted(z1h) == list(range(16)), "Z_1 h is not a permutation")

    def cycle_type(permutation):
        lengths, done = [], set()
        for start in range(len(permutation)):
            if start not in done:
                length, current = 0, start
                while current not in done:
                    done.add(current)
                    current = permutation[current]
                    length += 1
                lengths.append(length)
        return sorted(lengths)

    need(cycle_type(z1h) == [2, 2, 6, 6], "Z_1 h does not have cycle type 2^2 6^2")
    square = [tau[sigma[tau[sigma[n]]]] for n in range(8)]
    need(square == [point(s, y + z, y) for s, y, z in triples] and cycle_type(square) == [1, 1, 3, 3],
         "(Z_1 h)^2 on the first block is not (s,y,z) -> (s,y+z,y) of cycle type 1^2 3^2")

    # All normalized eighth-root two-bit phase tables, exact exponents modulo eight.
    phase_cases = []
    for f1, f2, f3 in product(range(8), repeat=3):
        phases = (0, f1, f2, f3)
        diagonal = (tuple(range(8)), tuple(phases[atom >> 1] for atom in range(8)))
        marked = (diagonal[0], tuple((diagonal[1][atom] + 2 * (atom & 1)) % 8
                                    for atom in range(8)))
        stationary = [i for i, j in q.charts(marked)
                      if all(q.app(C, v) in q.LAGS[i] for v in q.BASES[j])]
        clifford = q.symp(diagonal) is not None
        need(not stationary or clifford, "stationary phase table is not Clifford: " + str(phases))
        phase_cases.append({"phase": list(phases), "stationary_charts": stationary,
                            "clifford": clifford})
    need(len(phase_cases) == 512 and {tuple(row["phase"]) for row in phase_cases}
         == {(0,) + values for values in product(range(8), repeat=3)},
         "complete 512-table domain")
    phase_survivors = sum(bool(row["stationary_charts"]) for row in phase_cases)
    clifford_tables = sum(row["clifford"] for row in phase_cases)
    need(phase_survivors == 16 and clifford_tables == 32, "balanced table exact counts")

    # Every binary parameter choice, both required unipotence hypotheses.
    survivors, parameter_cases = [], []
    for a, b, c, d, e, f in product((0, 1), repeat=6):
        parameters = [a, b, c, d, e, f]
        a0 = (tuple(range(8)), tuple(((atom & 1) + a * ((atom >> 1) & 1)
              + b * ((atom >> 2) & 1) + 2 * c * ((atom >> 1) & 1) * ((atom >> 2) & 1)) % 8
              for atom in range(8)))
        mark = (tuple(range(8)), tuple((2 * d * ((atom >> 1) & 1)
                + 2 * e * ((atom >> 2) & 1)
                + 4 * f * ((atom >> 1) & 1) * ((atom >> 2) & 1)) % 8
                for atom in range(8)))
        mark_action = q.symp(mark)
        need(mark_action is not None, "binary phase correction is not Clifford")
        F = q.prod(C, mark_action)
        Fi = q.inverse(F)
        corrections = q.m.group(a0)
        need(corrections and len(corrections) == len(set(corrections)), "correction-image domain")
        tests = [(unipotent(q.prod(F, g)),
                  unipotent(q.prod(q.prod(q.prod(F, g), Fi), g))) for g in corrections]
        accepted = all(first and second for first, second in tests)
        common = []
        if accepted:
            common = [i for i, j in q.charts(a0) if i == j
                      and all(q.app(F, v) in q.LAGS[i] for v in q.BASES[i])]
            need(common, "surviving parameter choice has no common chamber")
            survivors.append({"parameters": parameters, "common_charts": common})
        parameter_cases.append({"parameters": parameters, "hypothesis_pass": accepted,
                                "correction_matrices": len(corrections), "common_charts": common})
    need(len(parameter_cases) == 64 and {tuple(row["parameters"]) for row in parameter_cases}
         == set(product((0, 1), repeat=6)), "complete 64-parameter domain")
    need(survivors == EXPECTED_SURVIVORS, "exact balanced parameter survivors or chambers changed")
    # The restored balanced rotations centralize both displayed common chambers.
    for sign in (1, -1):
        rotation = (tuple(range(8)), tuple(sign * (atom & 1) % 8 for atom in range(8)))
        for chart in (15, 20):
            for label in q.LAGS[chart]:
                pauli = q.pauli(label & 7, label >> 3)
                need(q.conj(rotation, pauli) == pauli, "restored rotation fails to centralize chamber")

    return {"status": "PASS", "route": "recompute", "exact": True,
            "fact": "balanced-rank-two-branch", "kernel_sha256": dict(KERNEL_HASHES),
            "six_cycle_atoms": 16, "six_cycle_obstruction": cycles,
            "z1h_permutation": z1h, "z1h_cycle_type": cycle_type(z1h),
            "obstruction_permutation": obstruction, "balanced_phase_domain": 512,
            "balanced_phase_survivors": phase_survivors, "clifford_phase_tables": clifford_tables,
            "balanced_phase_checks": phase_cases, "balanced_parameter_cases": 64,
            "parameter_checks": parameter_cases, "balanced_survivors": survivors,
            "restored_rotation_signs_checked": 2, "seconds": time.perf_counter() - started,
            "arithmetic": "exact permutations, modulo-eight monomials and binary matrices",
            "independent_kernel": False,
            "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--owner-root", type=Path, help="local selected conventional source root")
    parser.add_argument("--output", type=Path, help="write a new execution record")
    args = parser.parse_args()
    record = run(args.owner_root)
    encoded = json.dumps(record, indent=2) + "\n"
    if args.output is not None:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
