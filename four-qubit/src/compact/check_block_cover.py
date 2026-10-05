"""Exact orbit check of the block cover used for (C12).

Every rooted block form

    A = pi_0 i^(a xy + b xz + c yz + d xyz) T_p ,   a,b,c in F_2, d in Z/4, pi_0 in {1, pi: z -> z+xy}, p in F_2^3

is carried to one of thirty listed canonical cores.  This script verifies it by an exact orbit computation on
the 2*2^3*4*8 = 512 rooted forms of the parametrization `compact_core.mono(t,a,b,c,d,p)`, under the moves

  * the three translations of the inputs
  * conjugation by S_z
  * linear changes of the controls -- arbitrary for pi_0 = 1, radical-line preserving for pi_0 = pi --
    including the four target shears z -> z + lambda(x,y)
  * root-orientation changes by S_p
  * the two target-Clifford identities H_z pi H_z = Q and R_z (pi Q) R_z^dagger = CS_xy^{-1} Q, R_z = (H S^dagger)_z

Every move is applied as an exact conjugation/multiplication of the actual monomial: the monomial moves through
`completion_core`'s exact (permutation, phase) arithmetic, the two identity moves through an exact Z[zeta_8] matrix
layer (H_z and R_z are not monomial).  The action on (t,a,b,c,d,p) is COMPUTED by re-reading the image back through
the parametrization, never typed in; the left affine output correction is DERIVED from the image permutation and
asserted affine; where a move leaves the parametrized family the image is rejected (the move does not apply there)
and that fact is recorded.

Quotients (the omitted-phase argument, in computational form): the absorbed group is the input-side diagonal
Clifford group times the eighth-root scalars, exponents s + 2*lambda.n + 4*mu(n) (4096 monomials, each asserted
Clifford).  Marks therefore survive only as linear characters p in F_2^3, the cubic coefficient survives in full as
d in Z/4, and route/orientation scalars drop out.  The 512 parametrized forms are asserted to be a transversal of
that group, so the reduced parameters of an image are unique.

Verified: every orbit meets the thirty listed cores; the orbit -> cores map; all 14 marked and 6 unmarked
non-third-level cores and all 10 third-level cores are reached; every state fixed by all of its applicable moves is
a listed core; `is_c3` and dim Gamma_A (and the parity of d) are constant on orbits; dim Gamma_A of the twenty
marked cores reproduces rows 3/4/5/6 of the cover table; `is_c3` of all 512 forms agrees with the bound
record evidence/compact/core-levels.json.

Run without -O from the package root:  python src/compact/check_block_cover.py
                                       python src/compact/check_block_cover.py --negative-control
Writes evidence/compact/block-cover.json (status PASS only if every assertion holds).
"""
import sys
if not __debug__:
    raise RuntimeError("Run without Python -O: assertions are certificate checks")
from pathlib import Path
from itertools import product
import argparse, hashlib, json, time

t0 = time.perf_counter()
SRC = Path(__file__).resolve().parent
ROOT = SRC.parents[1]                      # package root
sys.path.insert(0, str(SRC))
import compact_core as q                   # exact small monomials; q.m is completion_core
EVID = ROOT / 'evidence/compact'


class CoverFailure(AssertionError):
    """A theorem-bearing check of the cover failed (raised so the negative control can catch it)."""


def require(cond, tag):
    if not cond:
        raise CoverFailure(tag)


# --------------------------------------------------------------------------------------------------
# exact Z[zeta_8] and 8x8 matrices -- only needed for H_z and R_z, which are not monomial
# --------------------------------------------------------------------------------------------------
ZERO = (0, 0, 0, 0)


def zeta(k):
    k %= 8
    s, k = (1, k) if k < 4 else (-1, k - 4)
    return tuple(s if i == k else 0 for i in range(4))


def zadd(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


def zmul(a, b):
    c = [0] * 7
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    c[i + j] += x * y
    return (c[0] - c[4], c[1] - c[5], c[2] - c[6], c[3])


def mat_of_mono(u):
    p, f = u
    M = [[ZERO] * 8 for _ in range(8)]
    for x in range(8):
        M[p[x]][x] = zeta(f[x])
    return tuple(tuple(r) for r in M)


def matmul(A, B):
    out = [[ZERO] * 8 for _ in range(8)]
    for i in range(8):
        Ai, Oi = A[i], out[i]
        for k in range(8):
            a = Ai[k]
            if a == ZERO:
                continue
            Bk = B[k]
            for j in range(8):
                b = Bk[j]
                if b != ZERO:
                    Oi[j] = zadd(Oi[j], zmul(a, b))
    return tuple(tuple(r) for r in out)


SCALED = {tuple(2 * c for c in zeta(k)): k for k in range(8)}


def mono_of_mat2(M):
    """Read a matrix known to be 2*(monomial) back as (perm, phase); None if it is not of that shape."""
    p, f = [0] * 8, [0] * 8
    for x in range(8):
        nz = [(i, M[i][x]) for i in range(8) if M[i][x] != ZERO]
        if len(nz) != 1 or nz[0][1] not in SCALED:
            return None
        p[x], f[x] = nz[0][0], SCALED[nz[0][1]]
    if sorted(p) != list(range(8)):
        return None
    return (tuple(p), tuple(f))


def lift(single, bit):
    """Embed a 2x2 matrix single[out_bit][in_bit] on the given qubit of three."""
    M = [[ZERO] * 8 for _ in range(8)]
    for n in range(8):
        for r in range(2):
            M[(n & ~(1 << bit)) | (r << bit)][n] = single[r][(n >> bit) & 1]
    return tuple(tuple(r) for r in M)


ONE, MI, PI_ = zeta(0), zeta(6), zeta(2)                       # 1, -i, i
SQRT2_HZ = lift([[ONE, ONE], [ONE, zeta(4)]], 2)               # sqrt2 * H on qubit z
SQRT2_RZ = lift([[ONE, MI], [ONE, PI_]], 2)                    # sqrt2 * (H S^dagger) on qubit z
SQRT2_RZD = lift([[ONE, ONE], [PI_, MI]], 2)                   # sqrt2 * R^dagger

# --------------------------------------------------------------------------------------------------
# the 512 rooted block forms, and the absorbed group
# --------------------------------------------------------------------------------------------------
IDP = tuple(range(8))
IDM = (IDP, (0,) * 8)
PIM = q.mono(1, 0, 0, 0, 0, 0)                                 # the Toffoli permutation, zero phase
PIP = PIM[0]
QM = q.mono(0, 0, 0, 0, 2, 0)                                  # Q = CCZ
PIQ = q.mono(1, 0, 0, 0, 2, 0)                                 # pi Q
CSINVQ = q.mono(0, 3, 0, 0, 2, 0)                              # CS_xy^{-1} Q

PARAMS = [tp for tp in product(range(2), range(2), range(2), range(2), range(4), range(8))]
STATE = [q.mono(*tp) for tp in PARAMS]
INDEX = {s: i for i, s in enumerate(STATE)}
assert len(PARAMS) == 512 and len(INDEX) == 512

BITS = [tuple((n >> j) & 1 for j in range(3)) for n in range(8)]


def diag(expo):
    return (IDP, tuple(e % 8 for e in expo))


ABSORB = set()
for lam in product(range(4), repeat=3):
    for mu in product(range(2), repeat=3):
        for s in range(8):
            ABSORB.add(diag([s + 2 * (lam[0] * x + lam[1] * y + lam[2] * z)
                             + 4 * (mu[0] * x * y + mu[1] * x * z + mu[2] * y * z)
                             for x, y, z in BITS]))
assert len(ABSORB) == 4096                                     # 4^3 lambda * 2^3 mu * 8 scalars, all distinct
for dmon in ABSORB:
    assert q.symp(dmon) is not None                            # every absorbed element really is a Clifford
# it is a group: generated by the S_j, the CZ_jk and the scalar, closed under those generators and under inverses
ABSORB_GENS = ([diag([2 * b[j] for b in BITS]) for j in range(3)]
               + [diag([4 * b[j] * b[k] for b in BITS]) for j, k in ((0, 1), (0, 2), (1, 2))]
               + [diag([1] * 8)])
for g in ABSORB_GENS:
    assert g in ABSORB
for dmon in ABSORB:
    assert q.inv(dmon) in ABSORB
    for g in ABSORB_GENS:
        assert q.mul(g, dmon) in ABSORB
# and it is ALL of the input-side diagonal Clifford group times mu_8: of the 8^7 diagonal monomials normalized to
# phase zero at n = 0, exactly the 512 normalized members of ABSORB are Clifford.  Being a group inside the
# Clifford group, ABSORB makes "Clifford" constant on each of the 8^7/512 = 4096 cosets of its normalized part,
# so testing one representative of each settles all 2,097,152 of them: exactly the trivial coset is Clifford.
COSETS = 0
for c1 in range(2):
    for c2 in range(2):
        for c4 in range(2):
            for c12 in range(4):
                for c13 in range(4):
                    for c23 in range(4):
                        for c123 in range(8):
                            rep = diag([c1 * x + c2 * y + c4 * z + c12 * x * y + c13 * x * z + c23 * y * z
                                        + c123 * x * y * z for x, y, z in BITS])
                            assert (q.symp(rep) is not None) == (not (c1 or c2 or c4 or c12 or c13 or c23 or c123))
                            COSETS += 1
assert COSETS == 4096


# every affine permutation matrix -- the only left corrections the moves are allowed to derive -- is a Clifford
AFFINE_SET = {(pm, (0,) * 8) for pm in (tuple(g[0][n] ^ s for n in range(8)) for g in q.LINEARS for s in range(8))}
assert len(AFFINE_SET) == 1344
for am in AFFINE_SET:
    assert q.affine(am[0]) and q.symp(am) is not None


def _inv_perm(p):
    out = [0] * 8
    for x in range(8):
        out[p[x]] = x
    return tuple(out)


def fit_perm(rho):
    """The derived affine output correction: the unique alpha with alpha . rho in {1, pi}, if affine."""
    ir = _inv_perm(rho)
    for t, target in ((0, IDP), (1, PIP)):
        alpha = tuple(target[ir[k]] for k in range(8))
        if q.affine(alpha):
            assert (alpha, (0,) * 8) in AFFINE_SET       # hence a Clifford, by the assertion above
            return t, (alpha, (0,) * 8)
    return None


def canon(v):
    """Read an exact monomial back through the parametrization; None if it leaves the family.

    Returns (state index, derived affine output correction, absorbed scalar).  The parameters are computed from
    the image (mark from the linear part mod two, quadratic/cubic from `compact_core.signature`) and then VERIFIED
    by exact membership of the residue in the absorbed diagonal-Clifford group; the scalar is that residue's
    eighth-root part, odd exactly when the image is an odd power of zeta_8 times a parametrized form.
    """
    fp = fit_perm(v[0])
    if fp is None:
        return None
    t, alpha = fp
    w = q.mul(alpha, v)                                        # permutation part is now exactly pi_t
    G = w[1]
    p = 0
    for j in range(3):
        if (G[1 << j] - G[0]) % 2:
            p |= 1 << j
    res = q.mul(q.inv(q.mono(t, 0, 0, 0, 0, p)), w)            # strip the mark; the rest must be even
    assert res[0] == IDP
    A, B, C, D = q.signature(res[1])                           # xy, xz, yz mod 4 and xyz mod 8
    if (A | B | C | D) & 1:
        return None
    cand = q.mono(t, (A // 2) % 2, (B // 2) % 2, (C // 2) % 2, (D // 2) % 4, p)
    absorbed = q.mul(q.inv(cand), w)
    if absorbed not in ABSORB:
        return None
    return INDEX[cand], alpha, absorbed[1][0]


# the family is a transversal of the absorbed group: distinct parameters are never absorbed into each other
INV_STATE = [q.inv(u) for u in STATE]
for i in range(512):
    ui = INV_STATE[i]
    for j in range(i + 1, 512):
        assert q.mul(ui, STATE[j]) not in ABSORB
    assert canon(STATE[i]) == (i, IDM, 0)

# --------------------------------------------------------------------------------------------------
# the generating moves
# --------------------------------------------------------------------------------------------------
SZ = diag([2 * z for _, _, z in BITS])                         # S_z = i^z
SZD = q.inv(SZ)

MOVES = []


def add(name, kind, cite, **kw):
    mv = {'name': name, 'kind': kind, 'citation': cite}
    mv.update(kw)
    MOVES.append(mv)
    return mv


for j, ax in enumerate('xyz'):
    add('translate_' + ax, 'monomial', 'input translation',
        Co=IDM, Ci=q.pauli(1 << j, 0),
        note='input translation X_%s with the derived affine output correction' % ax)
add('conj_Sz', 'monomial', 'conjugation by S_z', Co=SZ, Ci=SZD,
    note='conjugation by S_z')
add('conj_Sz_dagger', 'monomial', 'conjugation by S_z^dagger', Co=SZD, Ci=SZ,
    note='inverse of the previous move')
add('root_Sp', 'monomial', 'root-orientation change', Co=IDM, Ci=None,
    state_dep='right multiplication by S_p^{-1}, i.e. T_p -> T_p^{-1} (root-orientation change)')

LIN = list(q.LINEARS)                                          # the 168 linear coordinate changes
assert len(LIN) == 168
FIX_Z = [g for g in LIN if g[0][4] == 4]                       # those preserving the radical line <e_z>
SHEAR = [g for g in FIX_Z if g[0][1] in (1, 5) and g[0][2] in (2, 6)]
CTRL = [g for g in FIX_Z if g[0][1] in (1, 2, 3) and g[0][2] in (1, 2, 3)]
assert len(FIX_Z) == 24 and len(SHEAR) == 4 and len(CTRL) == 6
for k, g in enumerate(LIN):
    kind = 'shear' if g in SHEAR else ('control' if g in CTRL else ('radical' if g in FIX_Z else 'general'))
    add('linear_%03d' % k, 'monomial', 'linear change of the inputs',
        Co=g, Ci=q.inv(g), linear=kind, perm=list(g[0]))

add('identity_Hz', 'matrix', 'target-Clifford identity',
    Co2=SQRT2_HZ, Ci2=SQRT2_HZ, note='conjugation by H_z; H_z pi H_z = Q')
add('identity_Rz', 'matrix', 'target-Clifford identity',
    Co2=SQRT2_RZ, Ci2=SQRT2_RZD, note='conjugation by R_z = (H S^dagger)_z; R_z (pi Q) R_z^dagger = CS_xy^{-1} Q')
add('identity_Rz_dagger', 'matrix', 'target-Clifford identity',
    Co2=SQRT2_RZD, Ci2=SQRT2_RZ, note='inverse of the previous move')

IDENTITY_MOVES = ('identity_Hz', 'identity_Rz', 'identity_Rz_dagger')


def sp_mono(p):
    """S_p = i^(parity of p & n): the Clifford square root of the mark's square."""
    return diag([2 * (bin(n & p).count('1') % 2) for n in range(8)])


# The closed form of each named move's action on (t,a,b,c,d,p).  These are NOT the definition of the move -- the
# move is the exact two-sided Clifford change above -- they are read off the computed action and asserted below, so
# the record can state the action in the proof's own language and be wrong only if the assertion fires.
ACTION = {
    'translate_x': ('(t, a, b, c+d, -d, p)', lambda s: (s[0], s[1], s[2], (s[3] + s[4]) % 2, (-s[4]) % 4, s[5])),
    'translate_y': ('(t, a, b+d, c, -d, p)', lambda s: (s[0], s[1], (s[2] + s[4]) % 2, s[3], (-s[4]) % 4, s[5])),
    'translate_z': ('(t, a+d, b, c, -d, p)', lambda s: (s[0], (s[1] + s[4]) % 2, s[2], s[3], (-s[4]) % 4, s[5])),
    'conj_Sz': ('identity on the affine forms; (1, a+1, b, c, d+2, p) on the Toffoli forms',
                lambda s: s if s[0] == 0 else (1, (s[1] + 1) % 2, s[2], s[3], (s[4] + 2) % 4, s[5])),
    'conj_Sz_dagger': ('identity on the affine forms; (1, a+1, b, c, d-2, p) on the Toffoli forms -- the same map, '
                       'because d is taken modulo four',
                       lambda s: s if s[0] == 0 else (1, (s[1] + 1) % 2, s[2], s[3], (s[4] - 2) % 4, s[5])),
    'root_Sp': ('the identity on the parametrization: S_p is itself an absorbed diagonal Clifford, so the mark '
                'survives only as the linear character p', lambda s: s),
}
IDENTITY_DOMAIN = 'b = c = 0 (zero cross functional), d even, and the mark tangent or absent (p_z = 0)'


def apply_move(mv, i):
    u = STATE[i]
    if mv['kind'] == 'monomial':
        Ci = q.inv(sp_mono(PARAMS[i][5])) if mv['name'] == 'root_Sp' else mv['Ci']
        r = canon(q.mul(q.mul(mv['Co'], u), Ci))
    else:
        w = mono_of_mat2(matmul(matmul(mv['Co2'], mat_of_mono(u)), mv['Ci2']))
        r = None if w is None else canon(w)
    if r is not None and mv.get('affine_exit_only') and PARAMS[i][0] == 1 and PARAMS[r[0]][0] == 1:
        return None          # control: use the identities only to leave the Toffoli family
    return r


# every Clifford a move multiplies by is asserted to be one: the monomial conjugators through `symp`, the two
# non-monomial target Cliffords by conjugating all 64 Paulis in the exact matrix layer and reading a Pauli back
GENERATORS_CHECKED = 0
for mv in MOVES:
    if mv['kind'] == 'monomial':
        for c in (mv['Co'], mv['Ci']):
            if c is not None:
                assert q.symp(c) is not None, mv['name']
                GENERATORS_CHECKED += 1
for p in range(8):
    assert q.symp(sp_mono(p)) is not None
    GENERATORS_CHECKED += 1
for nm, (Cs, Cds) in (('H_z', (SQRT2_HZ, SQRT2_HZ)), ('R_z', (SQRT2_RZ, SQRT2_RZD)),
                      ('R_z^dagger', (SQRT2_RZD, SQRT2_RZ))):
    for v in range(64):
        w = mono_of_mat2(matmul(matmul(Cs, mat_of_mono(q.PAULIS[v])), Cds))
        assert w is not None and q.m.label(w) is not None, (nm, v)
    GENERATORS_CHECKED += 1


# --------------------------------------------------------------------------------------------------
# structural cross-checks of the two layers and of the two identities
# --------------------------------------------------------------------------------------------------
for u in (STATE[0], PIM, QM, PIQ, STATE[357], STATE[199]):
    for w in (PIQ, SZ, STATE[411]):
        assert mat_of_mono(q.mul(u, w)) == matmul(mat_of_mono(u), mat_of_mono(w))     # the two layers agree


def scale2(M):
    return tuple(tuple(tuple(2 * c for c in e) for e in row) for row in M)


assert matmul(matmul(SQRT2_HZ, mat_of_mono(PIM)), SQRT2_HZ) == scale2(mat_of_mono(QM))        # H_z pi H_z = Q
assert matmul(matmul(SQRT2_RZ, mat_of_mono(PIQ)), SQRT2_RZD) == scale2(mat_of_mono(CSINVQ))   # R_z(pi Q)R_z^+
assert mono_of_mat2(scale2(mat_of_mono(PIQ))) == PIQ                                          # round trip

# --------------------------------------------------------------------------------------------------
# invariants of the 512 states
# --------------------------------------------------------------------------------------------------
GAMMA, ISC3 = [], []
for u in STATE:
    S = [v for v in range(64) if q.symp(q.conj(u, q.PAULIS[v])) is not None]
    st = set(S)
    assert 0 in st and all((v ^ w) in st for v in S for w in S)        # Gamma_A is a subspace
    n = len(st)
    assert n and not (n & (n - 1))
    GAMMA.append(n.bit_length() - 1)
    ISC3.append(q.is_c3(u))
assert all((g == 6) == c for g, c in zip(GAMMA, ISC3))                 # is_c3 <=> dim Gamma_A = 6

levels = json.loads((EVID / 'core-levels.json').read_text(encoding='utf-8'))
assert levels['status'] == 'PASS' and levels['rooted_forms'] == 512
rec_c3 = {tuple(f['parameters']): f['third_level'] for f in levels['forms']}
assert len(rec_c3) == 512
for i, tp in enumerate(PARAMS):
    assert rec_c3[tp] == ISC3[i], tp                                   # agrees with the bound record

CORES = [(r['core'], tuple(r['parameters']), r['marked'], r['third_level']) for r in levels['listed_cores']]
assert len(CORES) == 30
assert sum(1 for _, _, mk, c3 in CORES if not c3 and mk) == 14
assert sum(1 for _, _, mk, c3 in CORES if not c3 and not mk) == 6
assert sum(1 for _, _, _, c3 in CORES if c3) == 10
CORE_INDEX = {INDEX[q.mono(*prm)]: nm for nm, prm, _, _ in CORES}
assert len(CORE_INDEX) == 30
# rows 3/4/5/6 of the cover table (dim Gamma_A) for the twenty marked cores
TABLE_ROW = {}
for nm, prm, mk, _ in CORES:
    if not mk:
        continue
    t, a, b, c, d, p = prm
    if d % 2:
        TABLE_ROW[nm] = 3 if (t == 0 or p == 1) else 4
    else:
        TABLE_ROW[nm] = 5 if t else 6
assert len(TABLE_ROW) == 20 and sorted(TABLE_ROW.values()) == [3] * 4 + [4] * 2 + [5] * 8 + [6] * 6


# --------------------------------------------------------------------------------------------------
# orbits
# --------------------------------------------------------------------------------------------------
def orbits(moves, core_index):
    parent = list(range(512))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)

    stats = {mv['name']: {'applies': 0, 'moves_state': 0, 'applies_t0': 0, 'applies_t1': 0,
                          'corrections': set(), 'odd_scalar': 0} for mv in moves}
    edges = set()
    for i in range(512):
        for mv in moves:
            r = apply_move(mv, i)
            if r is None:
                continue
            j, alpha, scalar = r
            s = stats[mv['name']]
            s['applies'] += 1
            s['applies_t%d' % PARAMS[i][0]] += 1
            s['corrections'].add(alpha[0])
            s['odd_scalar'] += scalar % 2
            if j != i:
                s['moves_state'] += 1
                edges.add((i, j))
                union(i, j)
    for i, j in edges:                                                 # the relation is symmetric
        require((j, i) in edges, 'move relation not symmetric')
    comp = {}
    for i in range(512):
        comp.setdefault(find(i), []).append(i)
    orbs = [sorted(v) for v in comp.values()]
    orbs.sort(key=lambda o: (-len(o), o[0]))
    # theorem-bearing checks
    coreless = [(k, o) for k, o in enumerate(orbs) if not any(i in core_index for i in o)]
    if coreless:
        exc = CoverFailure('orbit %d meets no listed core (%d coreless class(es), %d states in all)'
                           % (coreless[0][0], len(coreless), sum(len(o) for _, o in coreless)))
        exc.coreless = [{'states': len(o), 'parameters': [list(PARAMS[i]) for i in o]} for _, o in coreless]
        raise exc
    for k, o in enumerate(orbs):
        hit = [core_index[i] for i in o if i in core_index]
        require(hit, 'orbit %d meets no listed core' % k)
        require(len({ISC3[i] for i in o}) == 1, 'is_c3 not constant on orbit %d' % k)
        require(len({GAMMA[i] for i in o}) == 1, 'dim Gamma_A not constant on orbit %d' % k)
        require(len({PARAMS[i][4] % 2 for i in o}) == 1, 'parity of d not constant on orbit %d' % k)
    outgoing = {i for i, _ in edges}
    fixed = [i for i in range(512) if i not in outgoing]        # no applicable move changes the state
    for i in fixed:
        require(i in core_index, 'terminal state %s is not a listed core' % (PARAMS[i],))
    return orbs, stats, edges, fixed


def report(orbs, core_index):
    out = []
    for k, o in enumerate(orbs):
        hit = sorted(core_index[i] for i in o if i in core_index)
        out.append({'orbit': k, 'states': len(o), 'dim_gamma': GAMMA[o[0]], 'third_level': ISC3[o[0]],
                    'cubic_parity': PARAMS[o[0]][4] % 2, 'permutations': sorted({PARAMS[i][0] for i in o}),
                    'marks': sorted({PARAMS[i][5] for i in o}), 'cores': hit,
                    'representative_parameters': list(PARAMS[min(i for i in o if i in core_index)])})
    return out


# the computed action of each named move, checked against the closed form the record states
for mv in MOVES:
    dom = [(i, apply_move(mv, i)) for i in range(512)]
    dom = [(i, r[0]) for i, r in dom if r is not None]
    if mv['name'] in ACTION:
        text, f = ACTION[mv['name']]
        assert len(dom) == 512, mv['name']
        for i, j in dom:
            assert f(PARAMS[i]) == PARAMS[j], (mv['name'], PARAMS[i], PARAMS[j])
    elif mv['name'].startswith('linear_'):
        # the move is u -> Lambda u Lambda^{-1} with Lambda the printed `permutation` L, so the phase polynomial
        # is substituted by n -> L^{-1} n and the mark by p -> L^{-T} p.  BOTH are asserted here, on every state.
        lam_inv = mv['Ci'][0]                                  # the permutation of Lambda^{-1}
        cols = [lam_inv[1], lam_inv[2], lam_inv[4]]            # the columns L^{-1} e_k
        for i, j in dom:
            u = STATE[i]
            v = q.mul(q.mul(mv['Co'], u), mv['Ci'])
            assert all(v[1][n] == u[1][lam_inv[n]] for n in range(8)), mv['name']      # n -> L^{-1} n
            pk = PARAMS[i][5]
            pred = sum((bin(cols[k] & pk).count('1') % 2) << k for k in range(3))      # p -> L^{-T} p
            assert PARAMS[j][5] == pred, (mv['name'], PARAMS[i], PARAMS[j], pred)
    elif mv['name'] in IDENTITY_MOVES:
        want = {i for i, s in enumerate(PARAMS) if s[2] == 0 and s[3] == 0 and s[4] % 2 == 0 and s[5] < 4}
        assert {i for i, _ in dom} == want and len(want) == 32, mv['name']
        # the four stationary Toffoli forms leave the Toffoli family under these moves
        for i, j in dom:
            assert PARAMS[i][2] == PARAMS[j][2] == 0 and PARAMS[i][3] == PARAMS[j][3] == 0
            assert PARAMS[i][5] == PARAMS[j][5] and PARAMS[j][4] % 2 == 0

# the linear changes act on the mark alone when the phase polynomial is zero: transitively on the seven nonzero
# marks of an affine form, with the tangent/transverse split on a Toffoli form
MARK_ORBITS = {}
for t in (0, 1):
    for p0 in range(1, 8):
        reach = set()
        for mv in MOVES:
            if not mv['name'].startswith('linear_'):
                continue
            r = apply_move(mv, INDEX[q.mono(t, 0, 0, 0, 0, p0)])
            if r is not None:
                s = PARAMS[r[0]]
                assert s[:5] == (t, 0, 0, 0, 0), s
                reach.add(s[5])
        MARK_ORBITS[(t, p0)] = sorted(reach)
require(all(MARK_ORBITS[(0, p)] == list(range(1, 8)) for p in range(1, 8)),
        'the linear changes are not transitive on the nonzero marks of an affine form')
require(all(MARK_ORBITS[(1, p)] == [1, 2, 3] for p in (1, 2, 3)) and
        all(MARK_ORBITS[(1, p)] == [4, 5, 6, 7] for p in (4, 5, 6, 7)),
        'the radical-preserving changes do not split the nonzero marks into tangent and transverse')

ORBS, STATS, EDGES, FIXED = orbits(MOVES, CORE_INDEX)
ORBIT_OF = {i: k for k, o in enumerate(ORBS) for i in o}
TABLE = report(ORBS, CORE_INDEX)
require(sum(len(o) for o in ORBS) == 512, 'orbits do not partition the 512 forms')
reached = sorted(n for row in TABLE for n in row['cores'])
require(len(reached) == 30 and len(set(reached)) == 30, 'the thirty listed cores are not each reached once')
CORE_BY_NAME = {nm: (prm, mk, c3) for nm, prm, mk, c3 in CORES}
for tag, sel in (('14 marked non-third-level', lambda v: v[1] and not v[2]),
                 ('6 unmarked non-third-level', lambda v: not v[1] and not v[2]),
                 ('10 third-level', lambda v: v[2])):
    want = sorted(nm for nm, v in CORE_BY_NAME.items() if sel(v))
    require(all(nm in reached for nm in want), 'the %s cores are not all reached' % tag)
# cores sharing an orbit must share the cover-table row
for row in TABLE:
    rows = {TABLE_ROW[n] for n in row['cores'] if n in TABLE_ROW}
    require(len(rows) <= 1, 'orbit %d mixes rows of the cover table' % row['orbit'])
for nm, prm, mk, _ in CORES:
    if mk:
        require(GAMMA[INDEX[q.mono(*prm)]] == TABLE_ROW[nm],
                'dim Gamma_A of %s does not match the cover table' % nm)
require(len(FIXED) == 1 and PARAMS[FIXED[0]] == (0, 0, 0, 0, 0, 0),
        'there should be exactly one terminal state, the identity core')
# the odd eighth-root scalars: the family is closed under translations only modulo mu_8
ODD_SCALAR = {n: s['odd_scalar'] for n, s in STATS.items() if s['odd_scalar']}
require(set(ODD_SCALAR) == {'translate_x', 'translate_y', 'translate_z'},
        'an odd zeta_8 scalar was absorbed outside the translations')
require(sum(ODD_SCALAR.values()) == 768, 'the count of odd-scalar applications changed')
# how far the identity moves reach inside the Toffoli family
IDENT_EDGES = [(i, apply_move(mv, i)[0]) for mv in MOVES if mv['name'] in IDENTITY_MOVES
               for i in range(512) if apply_move(mv, i) is not None and apply_move(mv, i)[0] != i]
IDENT_TT = sum(1 for i, j in IDENT_EDGES if PARAMS[i][0] == 1 and PARAMS[j][0] == 1)

# --------------------------------------------------------------------------------------------------
# negative controls
# --------------------------------------------------------------------------------------------------
def control(label, moves, core_index, expect='FAIL'):
    """Run the whole check with a move or a core removed.

    `expect=None` records the outcome without demanding it.  Both a theorem-bearing failure (`CoverFailure`) and a
    plain structural assertion are caught, so a perturbation that trips either is recorded rather than aborting.
    The STRONG criterion -- thirty orbits with exactly one listed core each -- is evaluated too, because "every
    orbit meets a core" alone would not notice a move whose loss merges or splits the classes.
    """
    rec = {'control': label, 'expected': expect or 'recorded, not demanded'}
    try:
        orbs, _, _, fixed = orbits(moves, core_index)
        strong = len(orbs) == 30 and all(sum(1 for i in o if i in core_index) == 1 for o in orbs)
        rec.update(observed='PASS', detail='every orbit still meets the listed core set',
                   strong_criterion='PASS' if strong else 'FAIL', orbit_count=len(orbs),
                   terminal_states=len(fixed))
    except CoverFailure as exc:
        rec.update(observed='FAIL', detail=str(exc), strong_criterion='FAIL')
        if getattr(exc, 'coreless', None) is not None:
            rec['coreless_classes'] = exc.coreless
    except AssertionError as exc:
        rec.update(observed='FAIL', detail='structural assertion fired: %s' % (exc,), strong_criterion='FAIL')
    assert expect is None or rec['observed'] == expect, rec
    return rec


def absorb_control():
    """Enlarge the absorbed group by the non-Clifford CS_xy phase and show the structural assertions fire."""
    extra = {q.mul(d, q.mono(0, 1, 0, 0, 0, 0)) for d in ABSORB}
    non_clifford = sum(1 for d in extra if q.symp(d) is None)
    big = ABSORB | extra
    collisions = sum(1 for i in range(512) for j in range(i + 1, 512) if q.mul(INV_STATE[i], STATE[j]) in big)
    rec = {'control': 'enlarge the absorbed group by the CS_xy phase', 'expected': 'FAIL',
           'observed': 'FAIL' if non_clifford and collisions else 'PASS', 'strong_criterion': 'FAIL',
           'detail': '%d of the %d added elements are not Clifford, so the absorbed-group assertion fires, and %d '
                     'pairs of parametrized forms fall into one coset, so the transversal assertion fires'
                     % (non_clifford, len(extra), collisions)}
    assert rec['observed'] == 'FAIL', rec
    return rec


CONTROLS = []


def run_controls():
    out = []
    keep = [mv for mv in MOVES if mv['name'] not in IDENTITY_MOVES]
    out.append(control('drop-move: both target-Clifford identities (H_z and R_z conjugations)', keep, CORE_INDEX))
    # redundancy probes: each identity is dropped with its inverse, so the move relation stays symmetric
    for label, drop in (('identity_Hz', ('identity_Hz',)),
                        ('identity_Rz (with its inverse)', ('identity_Rz', 'identity_Rz_dagger')),
                        ('conj_Sz (with its inverse)', ('conj_Sz', 'conj_Sz_dagger'))):
        part = [mv for mv in MOVES if mv['name'] not in drop]
        out.append(control('drop-move: %s only (redundancy probe)' % label, part, CORE_INDEX, expect=None))
    # restrict the identities to images that leave the Toffoli family
    restricted = [dict(mv, affine_exit_only=True) if mv['name'] in IDENTITY_MOVES else mv for mv in MOVES]
    out.append(control('restrict the identity moves to Toffoli-to-affine images only (the %d Toffoli-to-Toffoli '
                       'edges dropped)' % IDENT_TT,
                       restricted, CORE_INDEX, expect='PASS'))
    assert out[-1]['strong_criterion'] == 'PASS' and out[-1]['orbit_count'] == 30, out[-1]
    lone = [row['cores'][0] for row in TABLE if len(row['cores']) == 1]
    name = lone[0]
    dropped = {i: n for i, n in CORE_INDEX.items() if n != name}
    out.append(control('drop-core: %s (every core is the unique core of its orbit, so the bijection guarantees '
                       'this control fires)' % name, MOVES, dropped))
    out.append(absorb_control())
    return out


ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument('--negative-control', action='store_true',
                help='run only the negative controls and print them; write no record')
args = ap.parse_args()
CONTROLS = run_controls()
if args.negative_control:
    print(json.dumps({'negative_controls': CONTROLS}, indent=1), flush=True)
    raise SystemExit(0)

# --------------------------------------------------------------------------------------------------
# record
# --------------------------------------------------------------------------------------------------
def describe(mv):
    d = {'move': mv['name'], 'citation': mv['citation'],
         'applies': STATS[mv['name']]['applies'], 'changes_state': STATS[mv['name']]['moves_state'],
         'applies_affine': STATS[mv['name']]['applies_t0'], 'applies_toffoli': STATS[mv['name']]['applies_t1'],
         'derived_affine_corrections': len(STATS[mv['name']]['corrections'])}
    if 'note' in mv:
        d['note'] = mv['note']
    if 'state_dep' in mv:
        d['note'] = mv['state_dep']
    if mv['name'] in ACTION:
        d['action_on_parameters'] = ACTION[mv['name']][0]
    if mv['name'] in IDENTITY_MOVES:
        d['domain'] = IDENTITY_DOMAIN
        d['action_on_parameters'] = [[list(PARAMS[i]), list(PARAMS[apply_move(mv, i)[0]])]
                                     for i in range(512) if apply_move(mv, i) is not None]
    if mv['kind'] == 'monomial' and mv['Ci'] is not None:
        d['output_clifford'] = [list(mv['Co'][0]), list(mv['Co'][1])]
        d['input_clifford'] = [list(mv['Ci'][0]), list(mv['Ci'][1])]
    if 'linear' in mv:
        d['linear_class'] = mv['linear']
        d['permutation'] = mv['perm']
        d['action_on_parameters'] = LINEAR_ACTION
    d['odd_scalar_applications'] = STATS[mv['name']]['odd_scalar']
    return d


LINEAR_ACTION = ('substitution n -> L^{-1} n in the phase polynomial and p -> L^{-T} p in the mark, where L is the '
                 'matrix of the `permutation` field -- that field is the CONJUGATOR Lambda, the move being '
                 'u -> Lambda u Lambda^{-1}, so the substitution matrix is Lambda^{-1}, not Lambda.  Both halves '
                 'are asserted on every applicable state.  The inversion does not change the orbits, because the '
                 '168 linear changes form a group closed under inverses.')


named = [describe(mv) for mv in MOVES if not mv['name'].startswith('linear_')]
lin = [describe(mv) for mv in MOVES if mv['name'].startswith('linear_')]
LIN_CLASSES = {}
for cls in ('shear', 'control', 'radical', 'general'):
    sel = [d for d in lin if d['linear_class'] == cls]
    LIN_CLASSES[cls] = len(sel)
    require(all(d['applies_affine'] == 256 for d in sel),
            'the %s linear changes do not all preserve the affine family' % cls)
    require(all(d['applies_toffoli'] == (0 if cls == 'general' else 256) for d in sel),
            'radical-line applicability of the %s linear changes is not as the proof states' % cls)
require(sum(1 for d in lin if d['applies_toffoli']) == 24,
        'exactly the 24 linear changes fixing e_z should apply to the Toffoli forms')
require(len(ORBS) == 30 and all(len(row['cores']) == 1 for row in TABLE),
        'the thirty listed cores are not in bijection with the orbits')

inputs = {}
for rel in ('src/compact/compact_core.py', 'src/completion_core.py', 'src/compact/check_core_levels.py',
            'evidence/compact/core-levels.json', 'src/compact/check_block_cover.py'):
    inputs[rel] = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
# the consumed record must have been produced against the compact_core.py that is on disk now
assert levels['inputs']['src/compact/compact_core.py'] == inputs['src/compact/compact_core.py'], \
    'core-levels.json was produced against a different compact_core.py'

r = {'status': 'PASS',
     'statement': 'Every one of the 512 rooted block forms pi_0 i^(a xy + b xz + c yz + d xyz) T_p is carried by the '
                  'moves listed below to one of the thirty listed canonical cores; the '
                  'orbits of those moves are exactly the classes recorded below.',
     'inputs_relative_to': 'package root',
     'rooted_forms': 512, 'listed_cores': 30, 'marked_non_third_level': 14, 'unmarked_non_third_level': 6,
     'third_level_cores': 10,
     'absorbed_group': {'description': 'input-side diagonal Cliffords times eighth-root scalars, exponents '
                                       's + 2*lambda.n + 4*mu(n)', 'size': len(ABSORB),
                        'every_element_clifford': True, 'family_is_transversal': True, 'is_a_group': True,
                        'is_the_full_input_side_diagonal_clifford_group': True,
                        'closure_check': 'of the 8^7 = 2,097,152 diagonal monomials normalized to phase zero at '
                                         'n = 0, exactly the 512 normalized members of this group are Clifford; '
                                         'settled on one representative of each of the %d cosets, since the group '
                                         'lies inside the Clifford group' % COSETS,
                        'generators_clifford_checked': GENERATORS_CHECKED},
     'scalar_quotient': {'odd_zeta8_applications': sum(ODD_SCALAR.values()), 'by_move': ODD_SCALAR,
                         'note': 'every application whose image is an ODD power of zeta_8 times a parametrized '
                                 'form is a translation: translating through a mark with p_j = 1 inverts the root '
                                 'orientation, T_p -> zeta_8 T_p^{-1}.  The 512-form family is therefore closed '
                                 'under translations only modulo mu_8, and this record certifies the cover modulo '
                                 'mu_8; scalars do not enter the orbit computation.  No scalar tracking is attempted.'},
     'identity_move_reach': {'edges': len(IDENT_EDGES), 'toffoli_to_toffoli': IDENT_TT,
                             'note': '%d of the %d state-changing identity applications stay inside the Toffoli family; '
                                     'the control that forbids them still leaves 30 orbits with exactly one '
                                     'listed core each.' % (IDENT_TT, len(IDENT_EDGES))},
     'effective_move_families': {'count': 5,
                                 'note': 'translations, S_z conjugation, linear changes (shears included), and the '
                                         'two target-Clifford identities.  root_Sp is vacuous on the '
                                         'parametrization -- S_p is itself absorbed -- so it generates no edge.'},
     'moves': named + [{'move': 'linear_000..linear_167', 'citation': lin[0]['citation'],
                        'count': 168, 'apply_to_affine_forms': 168, 'apply_to_toffoli_forms': 24,
                        'classes': LIN_CLASSES,
                        'action_on_parameters': LINEAR_ACTION + '  For the Toffoli forms the derived output '
                                                'correction is affine exactly for the 24 changes fixing the '
                                                'radical line.',
                        'note': 'linear coordinate changes; all 168 preserve the affine family, exactly the 24 '
                                'fixing the radical line <e_z> preserve the Toffoli family (classes: shear = the '
                                'four target shears z -> z + lambda(x,y), which include the identity; control = the '
                                'five further lifts of GL(2,2) on the controls; radical = the remaining fifteen '
                                'changes fixing e_z; general = the 144 that move the radical line)'}],
     'linear_moves_preserving_the_radical_line': [d for d in lin if d['applies_toffoli']],
     'identities_verified': ['H_z pi H_z = Q', 'R_z (pi Q) R_z^dagger = CS_xy^{-1} Q'],
     'orbit_count': len(ORBS), 'orbits': TABLE,
     'cores_reached': reached,
     'terminal_states': [list(PARAMS[i]) for i in FIXED],
     'mark_orbits': {'affine_forms': 'the 168 linear changes are transitive on the seven nonzero marks',
                     'toffoli_forms': 'the 24 radical-line preserving changes have the two orbits '
                                      '{x, y, x+y} (tangent) and {z, x+z, y+z, x+y+z} (transverse)'},
     'invariants': {'is_c3': 'constant on orbits', 'dim_gamma_A': 'constant on orbits',
                    'cubic_parity': 'constant on orbits',
                    'core_levels_agreement': 'is_c3 of all 512 forms agrees with evidence/compact/core-levels.json. '
                                             'This recomputes the same function with the same primitives, so it is '
                                             'a staleness check on compact_core.py and on that record, not an '
                                             'independent confirmation of the levels.',
                    'cover_table_rows': {n: TABLE_ROW[n] for n in sorted(TABLE_ROW)}},
     'forms': [{'parameters': list(tp), 'orbit': ORBIT_OF[i], 'dim_gamma': GAMMA[i], 'third_level': ISC3[i]}
               for i, tp in enumerate(PARAMS)],
     'negative_controls': CONTROLS,
     'negative_control_exercised': True,
     'inputs': inputs,
     'elapsed_seconds': round(time.perf_counter() - t0, 2)}
(EVID / 'block-cover.json').write_text(json.dumps(r, indent=1) + '\n', encoding='utf-8')
print(json.dumps({k: v for k, v in r.items()
                  if k not in ('forms', 'linear_moves_preserving_the_radical_line', 'moves')}, indent=1), flush=True)
