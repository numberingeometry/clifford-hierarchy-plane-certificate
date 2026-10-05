"""Fact (C2): the n=3 ray table and the 256 conjugators.

Rays M|x> = zeta8^lam i^f(x) |pi x>, pi a standard-flag involution, f + f o pi = -lam (mod 4),
tr M = 0, modulo M ~ -M.  Columns: #Fix, #pi, #rays, #in C3, #Clifford-conjugate to [H_*].
For every non-C3 ray we search an affine-diagonal Clifford conjugator (affine permutation
alpha, diagonal Clifford E) with (E alpha) M (E alpha)^dagger = +-H_* exactly.
Also checks (BA)^4 = -I and (BA)^2 not in mu4 P(2), i.e. H_* in C4 \\ C3 via the controlled-Clifford lemma,
and directly H_* in C4 \\ C3.
"""
import itertools, time
from .mono3 import *


if not __debug__:
    raise RuntimeError("Run without -O: exact audit assertions must be enabled.")

def recompute():
    t0 = time.perf_counter()

    def idx(x, y, z):
        return x | (y << 1) | (z << 2)

    # Standard flag permutations pi(x,y,z)=(x+a, y+b0+b1 x, z+c0+c1 x+c2 y+c3 xy).
    flag = []
    for a, b0, b1, c0, c1, c2, c3 in itertools.product((0, 1), repeat=7):
        p = []
        for n in range(8):
            x, y, z = bit(n, 0), bit(n, 1), bit(n, 2)
            p.append(idx(x ^ a, y ^ b0 ^ (b1 & x), z ^ c0 ^ (c1 & x) ^ (c2 & y) ^ (c3 & x & y)))
        flag.append(tuple(p))
    invols = sorted({p for p in flag if all(p[p[n]] == n for n in range(8))})

    # H_* = diag_x(A, B) on (y, z); x = qubit 0 = control.  Build as a monomial.
    def Hstar():
        p, f = [0] * 8, [0] * 8
        for n in range(8):
            x, y, z = bit(n, 0), bit(n, 1), bit(n, 2)
            if x == 0:   # A = P+(Z_z) X_y + P-(Z_z) Y_y ; Y|0>=i|1>, Y|1>=-i|0>
                if z == 0:
                    p[n], f[n] = idx(0, y ^ 1, z), 0
                else:
                    p[n], f[n] = idx(0, y ^ 1, z), (2 if y == 0 else 6)
            else:        # B = P+(Z_y) Z_z + P-(Z_y) X_z
                if y == 0:
                    p[n], f[n] = n, (0 if z == 0 else 4)
                else:
                    p[n], f[n] = idx(1, y, z ^ 1), 0
        return tuple(p), tuple(f)

    H = Hstar()
    assert mul(H, H) == ID, "H_* must be an involution"

    # Diagonal Clifford phase functions in zeta units: 2(c + l.x) + 4 q.
    DC = set()
    for c, l0, l1, l2, q0, q1, q2 in itertools.product(range(4), range(4), range(4), range(4), (0, 1), (0, 1), (0, 1)):
        DC.add(tuple((2 * (c + l0 * bit(n, 0) + l1 * bit(n, 1) + l2 * bit(n, 2))
                      + 4 * (q0 * bit(n, 0) * bit(n, 1) + q1 * bit(n, 0) * bit(n, 2) + q2 * bit(n, 1) * bit(n, 2))) % 8
                     for n in range(8)))
    DCG = [diag(e) for e in DC]

    def conjugator(M):
        """Affine-diagonal Clifford C = E alpha with C M C^-1 = +-H_*, or None."""
        targets = {H, (H[0], tuple((v + 4) % 8 for v in H[1]))}
        for al in AFFINE:
            Ma = conj(perm(al), M)
            if Ma[0] != H[0]:
                continue
            for E in DCG:
                if conj(E, Ma) in targets:
                    return {"alpha": list(al), "E_phase_zeta": list(E[1])}
        return None

    rows = {}
    conjugators = []
    for p in invols:
        fix = sum(p[n] == n for n in range(8))
        r = rows.setdefault(fix, {"pi": 0, "rays": 0, "in_C3": 0, "conj_to_Hstar": 0})
        r["pi"] += 1
        seen = set()
        for lam in (0, 1):
            for f4 in itertools.product(range(4), repeat=8):
                if any((f4[n] + f4[p[n]] + lam) % 4 for n in range(8)):
                    continue
                ph = tuple((lam + 2 * f4[n]) % 8 for n in range(8))
                # trace = sum over fixed points of zeta^ph ; exact zero test
                if any(zeta_vec([ph[n] for n in range(8) if p[n] == n])):
                    continue
                key = min(ph, tuple((v + 4) % 8 for v in ph))
                if key in seen:
                    continue
                seen.add(key)
                M = (p, ph)
                assert mul(M, M) == ID
                r["rays"] += 1
                if is_c3(M):
                    r["in_C3"] += 1
                else:
                    cj = conjugator(M)
                    if cj is not None:
                        r["conj_to_Hstar"] += 1
                        conjugators.append({"pi": list(p), "phase_zeta": list(ph), **cj})

    expected = {8: [1, 35, 35, 0], 6: [4, 160, 160, 0], 4: [10, 480, 480, 0], 2: [12, 768, 512, 256], 0: [17, 4352, 4352, 0]}
    table_ok = all([rows[k]["pi"], rows[k]["rays"], rows[k]["in_C3"], rows[k]["conj_to_Hstar"]] == v for k, v in expected.items())
    noncl = sum(r["rays"] - r["in_C3"] for r in rows.values())

    # H_* level: in C4 \ C3.
    Hlevel = {"in_C3": is_c3(H), "in_C4": all(is_c3(conj(H, g)) for g in GENS)}

    out = {"lemma": "gauss-lem:exceptional", "manuscript": "327e9e7",
           "flag_involutions": len(invols), "table": {str(k): rows[k] for k in sorted(rows, reverse=True)},
           "non_C3_rays": noncl, "conjugators_found": len(conjugators),
           "table_matches_manuscript": table_ok, "Hstar": Hlevel,
           "holds": table_ok and len(conjugators) == noncl == 256 and Hlevel == {"in_C3": False, "in_C4": True},
           "seconds": round(time.perf_counter() - t0, 1)}
    return out, conjugators, H
