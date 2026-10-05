"""Exact finite checks for the conductor-descent lemmas of the paper.

Run ``python n3_descent_check.py`` with Python 3.10 or newer.
No third-party packages are used.
Pauli label v=a+8b denotes the Hermitian Pauli i^(a.b) X^a Z^b.
There are 135 Pauli Lagrangians. Each has seven nonzero quadratic
polar forms and four balanced linear parts, giving 28 non-Pauli
involution rays, each with a unique containing Lagrangian.
The script enumerates the Clifford orbit of a CCZ-conjugated X,
constructs stabilizer conjugators for the second generator and the
two types of third-generator subgroups, and checks the integer
row-lattice implications for the phase equations. It also checks
the cubic Boolean-difference rank and the nonaffine route used in
the proof.
"""

from itertools import combinations

if not __debug__:
    raise RuntimeError("Run without -O: exact descent assertions must be enabled.")


def bits(x):
    return x.bit_count()


def symp(v, w):
    a, b = v & 7, v >> 3
    c, d = w & 7, w >> 3
    return (bits(a & d) + bits(b & c)) & 1


def pauli_phase(v, w):
    a, b = v & 7, v >> 3
    c, d = w & 7, w >> 3
    exponent = bits(a & b) + bits(c & d) - bits((a ^ c) & (b ^ d)) + 2 * bits(b & c)
    return (1, 1j, -1, -1j)[exponent & 3]


def canonical(terms):
    terms = tuple(sorted((v, int(c.real if isinstance(c, complex) else c))
                         for v, c in terms if c))
    return terms if terms[0][1] > 0 else tuple((v, -c) for v, c in terms)


def span(basis):
    return tuple(sorted({sum_bits(basis, mask)
                         for mask in range(1 << len(basis))}))


def sum_bits(basis, mask):
    v = 0
    for i, b in enumerate(basis):
        if mask >> i & 1:
            v ^= b
    return v


def lagrangians():
    found = {}
    for basis in combinations(range(1, 64), 3):
        if any(symp(v, w) for v, w in combinations(basis, 2)):
            continue
        space = span(basis)
        if len(space) == 8:
            found.setdefault(space, basis)
    assert len(found) == 135, len(found)
    return found


def quadratic(x, polar, linear):
    out = bits(x & linear) & 1
    out ^= ((x & 1) and (x & 2) and (polar & 1)) != 0
    out ^= ((x & 1) and (x & 4) and (polar & 2)) != 0
    out ^= ((x & 2) and (x & 4) and (polar & 4)) != 0
    return int(out)


def candidate(basis, polar, linear):
    q = [quadratic(x, polar, linear) for x in range(8)]
    if sum((-1) ** z for z in q) != 0:
        return None
    out = []
    for mask in range(8):
        coeff = sum((-1) ** (q[x] ^ (bits(mask & x) & 1)) for x in range(8))
        if not coeff:
            continue
        assert abs(coeff) == 4
        label, phase = 0, 1
        for i in range(3):
            if mask >> i & 1:
                phase *= pauli_phase(label, basis[i])
                label ^= basis[i]
        assert phase in (1, -1)
        out.append((label, coeff // 4 * int(phase)))
    assert len(out) == 4
    return canonical(out)


def all_candidates():
    found = set()
    for basis in lagrangians().values():
        for polar in range(1, 8):
            for linear in range(8):
                row = candidate(basis, polar, linear)
                if row:
                    found.add(row)
    assert len(found) == 3780, len(found)
    print('non-Pauli involution rays:', len(found))
    found.update(((v, 2),) for v in range(1, 64))
    assert len(found) == 3843
    print('Pauli rays: 63')
    return found


def multiply(g, h):
    out = {}
    for v, c in g:
        for w, d in h:
            label = v ^ w
            out[label] = out.get(label, 0) + c * d * pauli_phase(v, w)
    return tuple(sorted((v, z / 2) for v, z in out.items() if z))


def ray(g):
    return canonical(g)


def nonisotropic(*gs):
    labels = [v for g in gs for v, _ in g]
    return any(symp(v, w) for v, w in combinations(labels, 2))


def commuting_candidate(g, h, candidates):
    gh = multiply(g, h)
    if gh != multiply(h, g):
        return False
    assert all(z.imag == 0 and z.real.is_integer() for _, z in gh)
    return ray(gh) in candidates


def clifford_map(images):
    """Signed Pauli-label action from the images of X0,X1,X2,Z0,Z1,Z2."""
    out = []
    roots = (1, 1j, -1, -1j)
    for v in range(64):
        a, b = v & 7, v >> 3
        label = 0
        phase = roots[bits(a & b) & 3]
        for i in range(6):
            if v >> i & 1:
                phase *= pauli_phase(label, images[i])
                label ^= images[i]
        assert phase in (1, -1)
        out.append((label, int(phase.real)))
    return tuple(out)


def clifford_generators():
    images = tuple(1 << i for i in range(6))
    out = []
    for j in range(3):
        h = list(images)
        h[j], h[j + 3] = h[j + 3], h[j]
        out.append(clifford_map(h))
        s = list(images)
        s[j] |= 1 << (j + 3)
        out.append(clifford_map(s))
    for control in range(3):
        for target in range(3):
            if control != target:
                cnot = list(images)
                cnot[control] |= 1 << target
                cnot[target + 3] |= 1 << (control + 3)
                out.append(clifford_map(cnot))
    return out


def apply_map(clifford, candidate):
    return canonical((clifford[v][0], c * clifford[v][1])
                     for v, c in candidate)


def compose(left, right):
    return tuple((left[w][0], sign * left[w][1]) for w, sign in right)


def inverse(clifford):
    out = [None] * 64
    for v, (w, sign) in enumerate(clifford):
        out[w] = (v, sign)
    assert all(entry is not None for entry in out)
    return tuple(out)


def orbit_with_conjugators(start, generators):
    identity = clifford_map(tuple(1 << i for i in range(6)))
    conjugators = {start: identity}
    queue = [start]
    for candidate in queue:
        for gate in generators:
            image = apply_map(gate, candidate)
            if image not in conjugators:
                conjugators[image] = compose(gate, conjugators[candidate])
                queue.append(image)
    return conjugators


def stabilizer_generators(start, orbit, generators):
    inverses = {candidate: inverse(gate) for candidate, gate in orbit.items()}
    stabilizers = set()
    for candidate, conjugator in orbit.items():
        for gate in generators:
            image = apply_map(gate, candidate)
            stabilizer = compose(inverses[image], compose(gate, conjugator))
            assert apply_map(stabilizer, start) == start
            stabilizers.add(stabilizer)
    return stabilizers


def third_subgroup(third, pair):
    return frozenset([ray(third)] + [ray(multiply(third, h)) for h in pair])


def subgroup_orbit(start, generators):
    found = {start}
    queue = [start]
    for subgroup in queue:
        for gate in generators:
            image = frozenset(apply_map(gate, h) for h in subgroup)
            if image not in found:
                found.add(image)
                queue.append(image)
    return found


def twice_matrix(candidate):
    roots = (1, 1j, -1, -1j)
    matrix = [[0] * 8 for _ in range(8)]
    for v, coefficient in candidate:
        a, b = v & 7, v >> 3
        phase = roots[bits(a & b) & 3]
        for x in range(8):
            sign = -1 if bits(b & x) & 1 else 1
            matrix[x ^ a][x] += coefficient * phase * sign
    return matrix


def ccz_conjugate_pauli_label(candidate):
    matrix = twice_matrix(candidate)
    conjugate = [[matrix[y][x] * (-1 if x == 7 else 1) *
                  (-1 if y == 7 else 1) for x in range(8)]
                 for y in range(8)]
    for v in range(1, 64):
        pauli = twice_matrix(((v, 2),))
        if conjugate == pauli or conjugate == [[-entry for entry in row]
                                              for row in pauli]:
            return v
    raise AssertionError('CCZ conjugate is not Pauli')


def check_c11():
    cubic_check()
    candidates = all_candidates()
    # CCZ X_j CCZ = X_j CZ on the other two qubits.
    A = []
    for j in range(3):
        others = [i for i in range(3) if i != j]
        terms = []
        for mask in range(4):
            z = sum((1 << others[i]) for i in range(2) if mask >> i & 1)
            terms.append(((1 << j) | (z << 3), -1 if mask == 3 else 1))
        A.append(canonical(terms))
        assert A[-1] in candidates
    generators = clifford_generators()
    first_orbit = orbit_with_conjugators(A[0], generators)
    assert set(first_orbit) == {g for g in candidates if len(g) == 4}
    print('non-Pauli Clifford orbit:', len(first_orbit))
    first_stabilizers = stabilizer_generators(A[0], first_orbit, generators)
    print('first-generator stabilizer generators:', len(first_stabilizers))
    Z2 = ((1 << (3 + 2), 2),)
    assert Z2 in candidates
    second = [g for g in candidates if g != A[0] and
              commuting_candidate(A[0], g, candidates) and
              nonisotropic(A[0], g)]
    print('second after A0:', len(second))
    assert A[1] in second
    second_orbit = orbit_with_conjugators(A[1], first_stabilizers)
    assert set(second_orbit) == set(second)
    assert all(apply_map(conjugator, A[0]) == A[0]
               for conjugator in second_orbit.values())
    print('second-generator stabilizer orbit:', len(second_orbit))
    subgroup = {ray(A[0]), ray(A[1]), ray(multiply(A[0], A[1]))}
    third = [g for g in candidates if g not in subgroup and
             all(commuting_candidate(g, h, candidates) for h in subgroup) and
             nonisotropic(A[0], A[1], g)]
    print('third after A0,A1:', len(third))
    print('third Pauli:', sum(len(g) == 1 for g in third))
    print('third non-Pauli:', sum(len(g) == 4 for g in third))
    types = {}
    for g in third:
        members = [ray(g), ray(multiply(g, A[0])), ray(multiply(g, A[1])),
                   ray(multiply(g, multiply(A[0], A[1])))]
        num_pauli = sum(len(h) == 1 for h in members)
        types[num_pauli] = types.get(num_pauli, 0) + 1
    print('third by Pauli kernel:', types)
    assert len(second) == 48
    assert types == {0: 8, 1: 4}
    assert Z2 in third
    assert A[2] in third
    third_subgroups = {third_subgroup(g, subgroup) for g in third}
    assert len(third_subgroups) == 3
    pair_stabilizers = stabilizer_generators(A[1], second_orbit,
                                             first_stabilizers)
    assert all(apply_map(gate, A[0]) == A[0] for gate in pair_stabilizers)
    pauli_groups = {group for group in third_subgroups
                    if any(len(g) == 1 for g in group)}
    nonpauli_groups = third_subgroups - pauli_groups
    assert subgroup_orbit(third_subgroup(Z2, subgroup), pair_stabilizers) == pauli_groups
    assert subgroup_orbit(third_subgroup(A[2], subgroup), pair_stabilizers) == nonpauli_groups
    assert {ccz_conjugate_pauli_label(g) for group in pauli_groups for g in group} == {
        a | (4 << 3) for a in range(4)}
    assert {ccz_conjugate_pauli_label(g) for group in nonpauli_groups for g in group} == {
        (a | 4) | (b << 3) for a in range(4) for b in (0, 4)}
    print('third-subgroup stabilizer orbits:', len(pauli_groups), len(nonpauli_groups))
    route_check()
    return {
        'lagrangians': 135, 'non_Pauli_rays': sum(len(g) == 4 for g in candidates),
        'Pauli_rays': sum(len(g) == 1 for g in candidates), 'total_candidate_rays': len(candidates),
        'Clifford_generators': len(generators), 'first_orbit': len(first_orbit),
        'first_stabilizer_generators': len(first_stabilizers), 'second_candidates': len(second),
        'second_orbit': len(second_orbit), 'third_candidates': len(third),
        'third_Pauli': sum(len(g) == 1 for g in third),
        'third_non_Pauli': sum(len(g) == 4 for g in third), 'third_by_Pauli_kernel': types,
        'third_subgroups': len(third_subgroups), 'pauli_subgroup_orbit': len(pauli_groups),
        'nonpauli_subgroup_orbit': len(nonpauli_groups),
        'cubic_difference_matrices': 1024, 'cubic_difference_rank': 7,
        'nonaffine_relative_action_order': 3,
        'signed_action_orbit_equalities_checked': True,
        'CCZ_normal_form_label_sets_checked': True,
    }


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            j = row.bit_length() - 1
            if j in pivots:
                row ^= pivots[j]
            else:
                pivots[j] = row
                break
    return len(pivots)


def cubic_check():
    def q(x, lower):
        x0, x1, x2 = x & 1, (x >> 1) & 1, (x >> 2) & 1
        monomials = (1, x0, x1, x2, x0 * x1, x0 * x2, x1 * x2)
        return (x0 * x1 * x2) ^ (sum((lower >> i & 1) * monomials[i]
                                        for i in range(7)) & 1)
    for lower in range(128):
        for p in range(8):
            rows = []
            for a in range(8):
                row = sum((q(x ^ a, lower) ^ q(x, lower) ^
                           (dot(p, a) & dot(p, x))) << x for x in range(8))
                assert bits(row) % 2 == 0
                rows.append(row)
            assert binary_rank(rows) == 7
    print('cubic difference matrices: 1024; each has even rows and rank 7')


def route_check():
    J = [[0, 1], [1, 0]]
    S = [[0, 0, *J[0]], [0, 0, *J[1]],
         [*J[0], 1, 0], [*J[1], 0, 1]]
    def product(A, B):
        return [[sum(A[i][k] * B[k][j] for k in range(4)) & 1
                 for j in range(4)] for i in range(4)]
    S2 = product(S, S)
    I = [[int(i == j) for j in range(4)] for i in range(4)]
    assert all((S2[i][j] ^ S[i][j] ^ I[i][j]) == 0
               for i in range(4) for j in range(4))
    assert S != I and product(S2, S) == I
    print('nonaffine route: relative symplectic action has order 3')


def rank_f2(T):
    A = [row[:] for row in T]
    rank = 0
    for j in range(3):
        pivot = next((i for i in range(rank, 3) if A[i][j]), None)
        if pivot is None:
            continue
        A[rank], A[pivot] = A[pivot], A[rank]
        for i in range(3):
            if i != rank and A[i][j]:
                A[i] = [x ^ y for x, y in zip(A[i], A[rank])]
        rank += 1
    return rank


def symmetric_matrices():
    for mask in range(64):
        a, b, c, d, e, f = [(mask >> i) & 1 for i in range(6)]
        yield [[a, b, c], [b, d, e], [c, e, f]]


def apply_T(T, a):
    return sum((sum(T[i][j] * ((a >> j) & 1) for j in range(3)) & 1) << i
               for i in range(3))


def dot(a, b):
    return bits(a & b) & 1


def row(a, k, x):
    out = [0] * 8
    for i, c in ((x ^ a ^ k, 1), (x ^ a, -1), (x ^ k, -1), (x, 1)):
        out[i] += c
    return out


def row_lattice_basis(rows):
    A = [r[:] for r in rows]
    pivots = []
    rank = 0
    for j in range(8):
        pivot = next((i for i in range(rank, len(A)) if A[i][j]), None)
        if pivot is None:
            continue
        A[rank], A[pivot] = A[pivot], A[rank]
        for i in range(rank + 1, len(A)):
            while A[i][j]:
                quotient = A[rank][j] // A[i][j]
                A[rank] = [x - quotient * y for x, y in zip(A[rank], A[i])]
                A[rank], A[i] = A[i], A[rank]
        if A[rank][j] < 0:
            A[rank] = [-x for x in A[rank]]
        assert all(A[i][j] == 0 for i in range(rank + 1, len(A)))
        pivots.append(j)
        rank += 1
    assert all(not any(r) for r in A[rank:])
    return A[:rank], pivots


def lattice_contains(basis, pivots, target):
    v = target[:]
    for r, j in zip(basis, pivots):
        if v[j] % r[j]:
            return False
        q = v[j] // r[j]
        v = [x - q * y for x, y in zip(v, r)]
    return not any(v)


def phase_rows(T, model):
    out = []
    for a in range(8):
        c = apply_T(T, a)
        ks = [k for k in range(8) if dot(k, c) == 0] if model == 'N' else [
            ((c & 1) << 1) | ((c & 2) >> 1)]
        for k in ks:
            for x in range(8):
                out.append(row(a, k, x))
    return out


def phase_check():
    nonsingular = [T for T in symmetric_matrices() if rank_f2(T) == 3]
    assert len(nonsingular) == 28
    rank_two = []
    for p in range(4):
        p0, p1 = p & 1, (p >> 1) & 1
        b = 1 ^ (p0 & p1)
        rank_two.append([[p0, b, 0], [b, p1, 0], [0, 0, 0]])
    assert len(rank_two) == 4 and all(rank_f2(T) == 2 for T in rank_two)
    families = [('N', nonsingular, (0, 1, 2)),
                ('P', rank_two, (0, 1)),
                ('P', nonsingular, (0, 1))]
    targets_checked = 0
    equation_rows = 0
    family_records = []
    for model, family, js in families:
        family_records.append({'model': model, 'systems': len(family), 'directions': list(js)})
        for T in family:
            rows = phase_rows(T, model)
            equation_rows += len(rows)
            basis, pivots = row_lattice_basis(rows)
            for j in js:
                for x in range(8):
                    target = [0] * 8
                    target[x ^ (1 << j)] += 4
                    target[x] -= 4
                    assert lattice_contains(basis, pivots, target)
                    targets_checked += 1
            if model == 'P' and T == [[0, 1, 0], [1, 0, 0], [0, 0, 0]]:
                for j in js:
                    for x in range(8):
                        target = [0] * 8
                        target[x ^ (1 << j)] += 2
                        target[x] -= 2
                        assert lattice_contains(basis, pivots, target)
                        targets_checked += 1
    print('phase families: 28, 4, 28; all lattice inclusions verified')
    return {'families': family_records, 'phase_systems': sum(len(family) for _, family, _ in families),
            'equation_rows_generated': equation_rows, 'target_rows_checked': targets_checked,
            'all_integer_lattice_inclusions': True}


if __name__ == '__main__':
    raise SystemExit('Use the maintained route-A runner to select C11 or C12.')
