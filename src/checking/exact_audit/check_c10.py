"""Fact (C14): reduction of B_0 to 256 gates, checked exactly.

Claims:
 1. All 10,752 affine or Toffoli-type permutations of F_2^3 lie in a left affine coset of the
    eight listed permutations (the identity and one Toffoli-type permutation per radical direction).
 2. For the identity and six of the seven radical directions, the 32 listed phases
    i^(a'xy+b'xz+c'yz+d'xyz) represent all 32 classes of Z/4-valued phases modulo left diagonal
    Cliffords; for the direction (1,1,1) they represent only 16.
 3. For the direction (1,1,1), for every phase class and every one of the twenty non-third-level
    cores A, some B_0 P_v A^dagger lies outside C_3(3).
 4. Consistency: the listed domain reproduces the 360 common-input pairs of (C15)(a).
The eight permutations and twenty cores are derived here from the parameter definitions in
four-qubit/src/compact/{compact_core,build_common_input}.py. Key: C10.
"""
import itertools, time
from .mono3 import *


if not __debug__:
    raise RuntimeError("Run without -O: exact audit assertions must be enabled.")

def recompute():
    t0 = time.perf_counter()
    out = {"fact": "C10", "manuscript": "327e9e7"}

    # Listed permutations: same deterministic rule as compact_core.POINTS.
    POINTS = [tuple(range(8))]
    RADICAL = [0]
    for k in range(1, 8):
        g = next(g for g in LINEARS if g[k] == 4)          # linear g with g(k) = e_z
        POINTS.append(pcomp(pinv(g), pcomp(TOFFOLI, g)))    # g^-1 o Toffoli o g, radical k
        RADICAL.append(k)
    assert len(set(POINTS)) == 8

    # 1. Affine or Toffoli-type permutations.
    AFF = set(AFFINE)
    toff_type = {pcomp(a, pcomp(TOFFOLI, b)) for a in AFFINE for b in AFFINE}
    S = AFF | toff_type
    coset_of = {}
    for s in S:
        reps = [i for i, l in enumerate(POINTS) if pcomp(s, pinv(l)) in AFF]
        coset_of[s] = reps
    out["claim1"] = {"affine": len(AFF), "toffoli_type": len(toff_type), "total": len(S),
                     "every_one_in_a_left_affine_coset": all(len(r) >= 1 for r in coset_of.values()),
                     "cosets_disjoint": all(len(r) == 1 for r in coset_of.values()),
                     "holds": len(S) == 10752 and all(len(r) >= 1 for r in coset_of.values())}

    # 2. Phase classes: g ~ g + e o l, e = c + sum l_j x_j + 2 q (diagonal Clifford phases, Z/4).
    def diag_cliff_phases():
        fs = set()
        for c, l0, l1, l2, q0, q1, q2 in itertools.product(range(4), range(4), range(4), range(4), (0, 1), (0, 1), (0, 1)):
            f = []
            for n in range(8):
                x, y, z = bit(n, 0), bit(n, 1), bit(n, 2)
                f.append((c + l0 * x + l1 * y + l2 * z + 2 * (q0 * x * y + q1 * x * z + q2 * y * z)) % 4)
            fs.add(tuple(f))
        return fs

    E = diag_cliff_phases()
    assert len(E) == 2048

    def class_map(l):
        """Canonical class representative for every g in (Z/4)^8 modulo {e o l}."""
        El = [tuple(e[l[x]] for x in range(8)) for e in E]
        rep = {}
        for g in itertools.product(range(4), repeat=8):
            if g in rep:
                continue
            for e in El:
                rep[tuple((g[x] + e[x]) % 4 for x in range(8))] = g
        return rep

    LISTED = []
    for a, b, c, d in itertools.product((0, 1), (0, 1), (0, 1), range(4)):
        g = []
        for n in range(8):
            x, y, z = bit(n, 0), bit(n, 1), bit(n, 2)
            g.append((a * x * y + b * x * z + c * y * z + d * x * y * z) % 4)
        LISTED.append(tuple(g))

    reps_by_dir = {}
    classes_hit = {}
    for i, l in enumerate(POINTS):
        rep = class_map(l)
        allreps = sorted(set(rep.values()))
        reps_by_dir[i] = allreps
        classes_hit[RADICAL[i]] = {"classes": len(allreps), "listed_classes_hit": len({rep[g] for g in LISTED})}
    full = [k for k, v in classes_hit.items() if v["listed_classes_hit"] == 32]
    short = {k: v["listed_classes_hit"] for k, v in classes_hit.items() if v["listed_classes_hit"] != 32}
    out["claim2"] = {"by_radical_direction(0=identity; k=x+2y+4z)": classes_hit,
                     "complete_directions": full, "incomplete_directions": short,
                     "holds": sorted(full) == [0, 1, 2, 3, 4, 5, 6] and short == {7: 16}}

    # Cores (build_common_input.alist), with their marks p_A.
    CORES = []
    for sign in (1, 3):
        CORES += [((0, 0, 0, 0, sign, 4), 4), ((1, 0, 0, 0, sign, 1), 1), ((1, 0, 0, 0, sign, 4), 4)]
    for eta in (0, 1):
        CORES += [((1, 0, 0, cross, 2 * eta, mark), mark) for mark, cross in ((1, 1), (2, 1), (4, 0), (4, 1))]
    for sign in (1, 3):
        CORES += [((t, 0, 0, 0, sign, 0), 0) for t in (0, 1)]
    CORES += [((1, 0, 0, 1, 2 * eta, 0), 0) for eta in (0, 1)]
    assert len(CORES) == 20
    core_gates = [(mono(*params), mark) for params, mark in CORES]
    out["cores_all_non_third_level"] = all(not is_c3(A) for A, _ in core_gates)

    def all_derivs_c3(B0, A):
        Ai = inv(A)
        return all(is_c3(mul(mul(B0, PAULIS[v]), Ai)) for v in range(64))

    def gate(l, g, mark):
        T = mono(0, 0, 0, 0, 0, mark)
        return mul(perm(l), mul(diag([2 * x for x in g]), T))

    # 3. Direction (1,1,1): every class, every core, some derivative outside C_3.
    i111 = RADICAL.index(7)
    bad = []
    for A, mark in core_gates:
        for g in reps_by_dir[i111]:
            if all_derivs_c3(gate(POINTS[i111], g, mark), A):
                bad.append([CORES[core_gates.index((A, mark))][0], g])
    out["claim3"] = {"direction": [1, 1, 1], "classes": len(reps_by_dir[i111]), "cores": 20,
                     "counterexamples": bad, "holds": not bad}

    # 4. Consistency with the released 360 qualifying pairs.
    qual = 0
    for A, mark in core_gates:
        for l in POINTS:
            for g in LISTED:
                if all_derivs_c3(gate(l, g, mark), A):
                    qual += 1
    out["claim4"] = {"domain": 20 * 8 * 32, "qualifying_pairs": qual, "holds": qual == 360}

    out["holds"] = all(out[k]["holds"] for k in ("claim1", "claim2", "claim3", "claim4")) and out["cores_all_non_third_level"]
    out["seconds"] = round(time.perf_counter() - t0, 1)
    return out
