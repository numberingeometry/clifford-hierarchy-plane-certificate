"""Exact three-qubit monomial algebra over Z[zeta8].

A monomial gate is u = (p, f): u|x> = zeta8^f[x] |p[x]>, x = x0 + 2 x1 + 4 x2,
with qubit 0 (called x) the low bit, then y, then z.  Phases f are integers mod 8.
Pauli labels v = a + 8 b mean P_v = i^(a.b) X^a Z^b (Z acts first).
Every test below is exact: no floating point is used.
"""
import itertools

N = 8
ID = (tuple(range(N)), (0,) * N)


def bit(n, j):
    return (n >> j) & 1


def mul(u, w):
    """u*w (w acts first)."""
    pu, fu = u
    pw, fw = w
    return tuple(pu[pw[x]] for x in range(N)), tuple((fw[x] + fu[pw[x]]) % 8 for x in range(N))


def inv(u):
    p, f = u
    q = [0] * N
    g = [0] * N
    for x in range(N):
        q[p[x]] = x
        g[p[x]] = (-f[x]) % 8
    return tuple(q), tuple(g)


def conj(u, w):
    return mul(mul(u, w), inv(u))


def pauli(v):
    a, b = v & 7, v >> 3
    return (tuple(x ^ a for x in range(N)),
            tuple((2 * bin(a & b).count("1") + 4 * bin(b & x).count("1")) % 8 for x in range(N)))


PAULIS = [pauli(v) for v in range(64)]
GENS = [pauli(v) for v in (1, 2, 4, 8, 16, 32)]


def diag(f):
    return tuple(range(N)), tuple(x % 8 for x in f)


def perm(p):
    return tuple(p), (0,) * N


def pauli_label(u):
    """Label v if u is a scalar multiple of P_v, else None (projective equality)."""
    p, f = u
    a = p[0]
    if any(p[x] != x ^ a for x in range(N)):
        return None
    b = 0
    for j in range(3):
        d = (f[1 << j] - f[0]) % 8
        if d not in (0, 4):
            return None
        b |= (d // 4) << j
    if any((f[x] - f[0]) % 8 != 4 * (bin(b & x).count("1") % 2) for x in range(N)):
        return None
    return a + 8 * b


def is_pauli(u):
    return pauli_label(u) is not None


def is_clifford(u):
    return all(is_pauli(conj(u, g)) for g in GENS)


def is_c3(u):
    return all(is_clifford(conj(u, g)) for g in GENS)


def zeta_vec(exps_signs):
    """Sum of +-zeta^e as a vector in the Q-basis 1,zeta,zeta^2,zeta^3 (zeta^4=-1)."""
    v = [0, 0, 0, 0]
    for e in exps_signs:
        e %= 8
        v[e % 4] += 1 if e < 4 else -1
    return tuple(v)


def pauli_support(u):
    """Labels v with tr(P_v^dagger u) != 0, computed exactly."""
    p, f = u
    out = []
    for v in range(64):
        a, b = v & 7, v >> 3
        terms = [f[x] + 4 * (bin(b & x).count("1") % 2) for x in range(N) if p[x] == x ^ a]
        if any(zeta_vec(terms)):
            out.append(v)
    return out


def affine(p):
    return all(p[n] == p[0] ^ ((p[1] ^ p[0]) if n & 1 else 0) ^ ((p[2] ^ p[0]) if n & 2 else 0)
               ^ ((p[4] ^ p[0]) if n & 4 else 0) for n in range(N))


ALL_PERMS = list(itertools.permutations(range(N)))
AFFINE = [p for p in ALL_PERMS if affine(p)]
LINEARS = [p for p in AFFINE if p[0] == 0]
TOFFOLI = tuple(n ^ (4 * bit(n, 0) * bit(n, 1)) for n in range(N))  # (x,y,z)->(x,y,z+xy)


def pcomp(p, q):
    """p o q."""
    return tuple(p[q[x]] for x in range(N))


def pinv(p):
    q = [0] * N
    for x in range(N):
        q[p[x]] = x
    return tuple(q)


def mono(t, a, b, c, d, p):
    """pi^t i^(a xy + b xz + c yz + d xyz) T_p as in the certificate's parametrization."""
    pp = TOFFOLI if t else tuple(range(N))
    ph = []
    for n in range(N):
        x, y, z = bit(n, 0), bit(n, 1), bit(n, 2)
        ph.append((2 * (a * x * y + b * x * z + c * y * z + d * x * y * z) + (bin(n & p).count("1") % 2)) % 8)
    return pp, tuple(ph)


def span_dim(vs):
    lead = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h not in lead:
                lead[h] = v
                break
            v ^= lead[h]
    return len(lead)


def is_subspace(s):
    s = set(s)
    return 0 in s and all((u ^ w) in s for u in s for w in s)
