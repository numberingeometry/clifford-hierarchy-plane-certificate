"""Fact (C3): Cliffords commuting with H_*. Exact census, without tolerances.

Basis order is |xyz>, x the most significant bit. Operators act on columns, and products compose
right to left. A monomial operator is a pair (permutation, phase), phase[j] the exponent of i on input
column j. The domain is the 1,344 affine permutations times the 512 diagonal Clifford phases, modulo
scalars. The reduction from arbitrary commuting Cliffords to this domain is part of the proof in the
paper. Internal key: C14.
"""

from itertools import product
import json
import time


class CheckFailure(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise CheckFailure(message)


def bits(j):
    return ((j >> 2) & 1, (j >> 1) & 1, j & 1)


def rank_binary(columns):
    basis = {}
    for value in columns:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def affine_permutations():
    out = []
    for columns in product(range(1, 8), repeat=3):
        if rank_binary(columns) != 3:
            continue
        linear = tuple(
            (columns[0] if x else 0) ^ (columns[1] if y else 0) ^ (columns[2] if z else 0)
            for x, y, z in map(bits, range(8))
        )
        out.extend(tuple(value ^ translation for value in linear) for translation in range(8))
    require(len(out) == len(set(out)) == 1344, "affine domain is incomplete or has duplicates")
    return out


def diagonal_cliffords():
    out = []
    for a in product(range(4), repeat=3):
        for q in product(range(2), repeat=3):
            out.append(tuple(
                (a[0]*x + a[1]*y + a[2]*z + 2*(q[0]*x*y + q[1]*x*z + q[2]*y*z)) % 4
                for x, y, z in map(bits, range(8))
            ))
    require(len(out) == len(set(out)) == 512, "diagonal-Clifford domain is incomplete or has duplicates")
    require(all(phases[0] == 0 for phases in out), "scalar normalization changed")
    return out


def compose(left, right):
    p, f = left
    q, g = right
    return tuple(p[q[j]] for j in range(8)), tuple((g[j] + f[q[j]]) % 4 for j in range(8))


def adjoint(operator):
    p, f = operator
    inverse = [0] * 8
    for j, value in enumerate(p):
        inverse[value] = j
    return tuple(inverse), tuple(-f[inverse[j]] % 4 for j in range(8))


def hstar():
    p, f = [], []
    for x, y, z in map(bits, range(8)):
        if x == 0:
            p.append((y ^ 1)*2 + z)
            f.append(z*(1 + 2*y) % 4)
        else:
            p.append(4 + y*2 + (z ^ y))
            f.append(2*(1-y)*z % 4)
    result = tuple(p), tuple(f)
    identity = tuple(range(8)), (0,)*8
    require(compose(result, result) == identity, "H* must square to identity")
    require(adjoint(result) == result, "H* must be Hermitian")
    return result


def commutes(left, right):
    return compose(left, right) == compose(right, left)


def hermitian_representative(operator):
    p, f = operator
    if any(p[p[j]] != j for j in range(8)):
        return None
    sums = {(f[j] + f[p[j]]) % 4 for j in range(8)}
    if len(sums) != 1:
        return None
    scalar_square = next(iter(sums))
    # Multiply by zeta_8^(-scalar_square), then store phases as mu_8 exponents.
    phases8 = tuple((2*f[j] - scalar_square) % 8 for j in range(8))
    require(all((phases8[j] + phases8[p[j]]) % 8 == 0 for j in range(8)), "Hermitian phase lift failed")
    return p, phases8


def claimed_family():
    return {
        abc: (tuple(range(8)), tuple(4*(a*x+b*x*y+c*(1+x)*z) % 8 for x, y, z in map(bits, range(8))))
        for abc in product(range(2), repeat=3)
        for a, b, c in (abc,)
    }


def walsh_numerators(phases8):
    require(all(e in (0, 4) for e in phases8), "claimed Hermitian diagonal has a non-real phase")
    signs = [1 if e == 0 else -1 for e in phases8]
    return {
        z: sum(signs[j] * (-1 if (j & z).bit_count() % 2 else 1) for j in range(8))
        for z in range(8)
    }


def validate_summary(summary):
    expected = {
        "affine_permutations": 1344,
        "diagonal_cliffords_mod_scalar": 512,
        "complete_domain_size": 688128,
        "commuting_cliffords_mod_scalar": 64,
        "hermitian_rays": 8,
        "traceless_hermitian_rays": 3,
    }
    for key, value in expected.items():
        require(summary.get(key) == value, f"C14 {key}: expected {value}, received {summary.get(key)!r}")
    require(summary.get("claimed_family_equals_hermitian_rays") is True, "C14 claimed-family comparison failed")
    require(summary.get("every_claimed_representative_commutes") is True, "C14 reverse inclusion failed")


def run():
    started = time.perf_counter()
    H = hstar()
    affines, diagonals = affine_permutations(), diagonal_cliffords()
    commuting, hermitian = [], []
    checked = 0
    for p in affines:
        # A commuting monomial operator must first commute at permutation level.
        permutation_commutes = all(p[H[0][j]] == H[0][p[j]] for j in range(8))
        for f in diagonals:
            checked += 1
            if not permutation_commutes:
                continue
            if all((H[1][j]+f[H[0][j]]-f[j]-H[1][p[j]]) % 4 == 0 for j in range(8)):
                operator = p, f
                require(commutes(operator, H), "optimized commutation disagrees with direct composition")
                commuting.append(operator)
                representative = hermitian_representative(operator)
                if representative is not None:
                    hermitian.append(representative)
    claimed = claimed_family()
    require(len(set(hermitian)) == len(hermitian), "Hermitian representative census has duplicate rays")
    same_family = set(hermitian) == set(claimed.values())
    reverse = all(commutes((p, tuple(e//2 for e in phases8)), H) for p, phases8 in claimed.values())
    traces = []
    representatives = []
    for abc, (p, phases8) in sorted(claimed.items()):
        coeffs = walsh_numerators(phases8)
        if coeffs[0] == 0:
            traces.append(abc)
        representatives.append({
            "parameters_a_b_c": list(abc),
            "permutation": list(p),
            "phases_mu8": list(phases8),
            "pauli_diagonal_coefficients_numerator_over_8": {str(z): value for z, value in coeffs.items() if value},
        })
    require(traces == [(0,1,1),(1,0,0),(1,1,1)], "traceless parameter list changed")
    expected_supports = {(0,1,1): {1,2,5,6}, (1,0,0): {4}, (1,1,1): {1,2,5,6}}
    for abc in traces:
        support = {z for z, value in walsh_numerators(claimed[abc][1]).items() if value}
        require(support == expected_supports[abc], f"Pauli support changed for {abc}")
        require(rank_binary(support) == (1 if abc == (1,0,0) else 3), f"Pauli support span changed for {abc}")
    result = {
        "fact": "cliffords-commuting-with-hstar",
        "route": "recompute",
        "status": "accepted",
        "scope": "complete monomial-Clifford commutant domain; analytic reduction to this domain remains a paper premise",
        "affine_permutations": len(affines),
        "diagonal_cliffords_mod_scalar": len(diagonals),
        "complete_domain_size": checked,
        "commuting_cliffords_mod_scalar": len(commuting),
        "hermitian_rays": len(hermitian),
        "traceless_hermitian_rays": len(traces),
        "claimed_family_equals_hermitian_rays": same_family,
        "every_claimed_representative_commutes": reverse,
        "hermitian_representatives_mod_sign": representatives,
        "seconds": time.perf_counter()-started,
    }
    validate_summary(result)
    return result


if __name__ == "__main__":
    print(json.dumps(run(), sort_keys=True))
