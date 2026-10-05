/*  ge4_check.c -- exact census kernel for the 1,010 four-qubit gates of (C7).

    Exact integer arithmetic throughout (Z[zeta8] as 4-tuples, zeta8^4 = -1).

    Conventions (unit): n = 4 qubits, N = 16, computational basis |x>, bit j of x is
    qubit j, qubit 0 low.  Hermitian Pauli P_{a,b} = i^{|a&b|} X^a Z^b, label v = a + 16 b,
    P_v |x> = i^{|a&b|} (-1)^{|b&x|} |x xor a>.  Products act rightmost first.
    A mu8-monomial M is M|x> = zeta8^{w[x]} |pi[x]>.

    Modes:
      pointdomain <triangular.txt>          build T = AGL-orbits of triangular involutions,
                                            verify T = I2 := {involutions pi with ord(pi T_a pi T_b)
                                            a power of two for all a,b} over all 46,206,736
                                            involutions of 16 points; write pointdomain.bin
      lock <roots.txt> <pointdomain.bin> <out.txt> [first last]
                                            for each root M1 (a mu8-monomial Hermitian involution),
                                            enumerate every skeleton (chamber L, involution pi' in T)
                                            of a class member h = C_L M' C_L^{-1} that commutes or
                                            anticommutes with M1, and intersect the set of Pauli
                                            planes met by the normalised-ray space of every such h.

    That this enumeration covers every partner, and that a surviving plane gives a returned physical
    Pauli ray, are part of the written proof in the paper; they are not checked here.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define N 16
#define NL 256
#define NLAG 2295
#define NPLANES 10795
#define PBITS ((NPLANES + 63) / 64)

typedef struct { int64_t c[4]; } z8;
static const z8 Z0 = {{0, 0, 0, 0}};
static int popc(unsigned x) { return __builtin_popcount(x); }
static z8 zunit(int e) { z8 r = Z0; e = ((e % 8) + 8) % 8; if (e < 4) r.c[e] = 1; else r.c[e - 4] = -1; return r; }
static z8 zadd(z8 a, z8 b) { z8 r; for (int k = 0; k < 4; k++) r.c[k] = a.c[k] + b.c[k]; return r; }
static z8 zmul(z8 a, z8 b) {
    int64_t t[8] = {0};
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) t[i + j] += a.c[i] * b.c[j];
    z8 r; for (int k = 0; k < 4; k++) r.c[k] = t[k] - t[k + 4]; return r;
}
static int ziszero(z8 a) { return !a.c[0] && !a.c[1] && !a.c[2] && !a.c[3]; }
static int zeq(z8 a, z8 b) { return a.c[0] == b.c[0] && a.c[1] == b.c[1] && a.c[2] == b.c[2] && a.c[3] == b.c[3]; }

/* ------------------------------------------------------------ Paulis */
static int omega(int v, int w) { int a = v & 15, b = v >> 4, c = w & 15, d = w >> 4; return (popc(a & d) + popc(b & c)) & 1; }
/* P_v P_w = i^{kap(v,w)} P_{v^w} */
static int kap(int v, int w) {
    int a = v & 15, b = v >> 4, c = w & 15, d = w >> 4;
    return ((popc(a & b) + popc(c & d) - popc((a ^ c) & (b ^ d)) + 2 * popc(b & c)) % 4 + 4) % 4;
}
/* weight of P_v at |x>: P_v|x> = zeta8^{pw(v,x)} |x ^ a| */
static int pw(int v, int x) { int a = v & 15, b = v >> 4; return (2 * popc(a & b) + 4 * (popc(b & x) & 1)) % 8; }

/* ------------------------------------------------------------ Lagrangians and planes */
static int lag[NLAG][N];      /* the 16 labels of each Lagrangian (sorted) */
static int lagbasis[NLAG][4];
static int nlag = 0;
typedef struct { uint64_t b[4]; } set256;
static void set_add(set256 *s, int v) { s->b[v >> 6] |= 1ULL << (v & 63); }
static int set_has(const set256 *s, int v) { return (s->b[v >> 6] >> (v & 63)) & 1; }
static int set_cmp(const void *x, const void *y) { return memcmp(x, y, sizeof(set256)); }

static void build_lagrangians(void) {
    /* BFS over isotropic subspaces by dimension, sets stored as 256-bit masks */
    set256 *cur = malloc(sizeof(set256) * 400000), *nxt = malloc(sizeof(set256) * 400000);
    int ncur = 1; memset(&cur[0], 0, sizeof(set256)); set_add(&cur[0], 0);
    for (int dim = 0; dim < 4; dim++) {
        int nn = 0;
        for (int i = 0; i < ncur; i++) {
            int els[N], ne = 0;
            for (int v = 0; v < NL; v++) if (set_has(&cur[i], v)) els[ne++] = v;
            for (int v = 1; v < NL; v++) {
                if (set_has(&cur[i], v)) continue;
                int ok = 1; for (int k = 0; k < ne; k++) if (omega(v, els[k])) { ok = 0; break; }
                if (!ok) continue;
                set256 s = cur[i];
                for (int k = 0; k < ne; k++) set_add(&s, els[k] ^ v);
                nxt[nn++] = s;
            }
        }
        qsort(nxt, nn, sizeof(set256), set_cmp);
        int m = 0;
        for (int i = 0; i < nn; i++) if (i == 0 || set_cmp(&nxt[i], &nxt[i - 1])) nxt[m++] = nxt[i];
        set256 *t = cur; cur = nxt; nxt = t; ncur = m;
    }
    if (ncur != NLAG) { fprintf(stderr, "lagrangian count %d != %d\n", ncur, NLAG); exit(1); }
    for (int i = 0; i < ncur; i++) {
        int ne = 0;
        for (int v = 0; v < NL; v++) if (set_has(&cur[i], v)) lag[i][ne++] = v;
        /* basis: greedy */
        int nb = 0; set256 span; memset(&span, 0, sizeof span); set_add(&span, 0);
        for (int k = 0; k < N && nb < 4; k++) {
            int v = lag[i][k];
            if (set_has(&span, v)) continue;
            lagbasis[i][nb++] = v;
            set256 s2 = span;
            for (int u = 0; u < NL; u++) if (set_has(&span, u)) set_add(&s2, u ^ v);
            span = s2;
        }
    }
    nlag = ncur; free(cur); free(nxt);
}

static int planeid[NL][NL];          /* plane through nonzero rays s1 != s2 */
static int planes[NPLANES][3];
static int nplanes = 0;
static uint64_t planes_through[NL][PBITS];  /* bitset of planes containing ray s */
static void build_planes(void) {
    for (int s1 = 1; s1 < NL; s1++) for (int s2 = 1; s2 < NL; s2++) planeid[s1][s2] = -1;
    for (int s1 = 1; s1 < NL; s1++) for (int s2 = s1 + 1; s2 < NL; s2++) {
        int s3 = s1 ^ s2; if (s3 < s2) continue;  /* canonical: s1 < s2 < s3 */
        int id = nplanes++;
        planes[id][0] = s1; planes[id][1] = s2; planes[id][2] = s3;
        planeid[s1][s2] = planeid[s2][s1] = planeid[s1][s3] = planeid[s3][s1] = planeid[s2][s3] = planeid[s3][s2] = id;
    }
    if (nplanes != NPLANES) { fprintf(stderr, "plane count %d\n", nplanes); exit(1); }
    memset(planes_through, 0, sizeof planes_through);
    for (int id = 0; id < nplanes; id++) for (int k = 0; k < 3; k++) {
        int s = planes[id][k]; planes_through[s][id >> 6] |= 1ULL << (id & 63);
    }
}

/* ------------------------------------------------------------ Cliffords as label actions */
typedef struct { int sig[NL]; int eps[NL]; } cliff;   /* C P_u C^-1 = i^{eps[u]} P_{sig[u]}, eps in {0,2} */
static void cliff_id(cliff *c) { for (int u = 0; u < NL; u++) { c->sig[u] = u; c->eps[u] = 0; } }
/* transvection Clifford R_v = (I + i P_v)/sqrt2 : P_u -> P_u if commute, else i P_v P_u = i^{1+kap(v,u)} P_{v^u} */
static void apply_transvection(cliff *c, int v) {
    for (int u = 0; u < NL; u++) {
        int w = c->sig[u];
        if (omega(v, w)) { c->eps[u] = (c->eps[u] + 1 + kap(v, w)) % 4; c->sig[u] = w ^ v; }
    }
}
static int cliff_check(const cliff *c) {  /* automorphism identity for all pairs */
    for (int u = 0; u < NL; u++) if (c->eps[u] & 1) return 0;
    for (int u = 0; u < NL; u++) for (int w = 0; w < NL; w++) {
        int lhs = (c->eps[u] + c->eps[w] + kap(u, w)) % 4;
        int rhs = (c->eps[u ^ w] + kap(c->sig[u], c->sig[w])) % 4;
        if (lhs != rhs) return 0;
    }
    return 1;
}
static cliff CL[NLAG];
static void build_chamber_cliffords(void) {
    for (int L = 0; L < nlag; L++) {
        cliff c; cliff_id(&c);
        int ell[4]; for (int j = 0; j < 4; j++) ell[j] = lagbasis[L][j];
        for (int j = 0; j < 4; j++) {
            int zj = 16 << j;
            int cj = c.sig[zj];
            if (cj == ell[j]) continue;
            if (omega(cj, ell[j])) { apply_transvection(&c, cj ^ ell[j]); }
            else {
                int w = -1;
                for (int cand = 1; cand < NL && w < 0; cand++) {
                    if (!omega(cj, cand) || !omega(ell[j], cand)) continue;
                    int ok = 1; for (int i = 0; i < j; i++) if (omega(ell[i], cand)) { ok = 0; break; }
                    if (ok) w = cand;
                }
                if (w < 0) { fprintf(stderr, "no w for L=%d j=%d\n", L, j); exit(1); }
                apply_transvection(&c, cj ^ w);
                apply_transvection(&c, w ^ ell[j]);
            }
            if (c.sig[zj] != ell[j]) { fprintf(stderr, "map fail L=%d j=%d\n", L, j); exit(1); }
            for (int i = 0; i < j; i++) if (c.sig[16 << i] != ell[i]) { fprintf(stderr, "disturbed L=%d\n", L); exit(1); }
        }
        if (!cliff_check(&c)) { fprintf(stderr, "cliff check fail L=%d\n", L); exit(1); }
        CL[L] = c;
    }
}

/* ------------------------------------------------------------ point domain */
typedef uint64_t perm16;  /* 16 nibbles, nibble x = pi(x) */
static int pget(perm16 p, int x) { return (p >> (4 * x)) & 15; }
static perm16 pset(perm16 p, int x, int v) { return (p & ~(15ULL << (4 * x))) | ((uint64_t)v << (4 * x)); }
static perm16 pcompose(perm16 p, perm16 q) { perm16 r = 0; for (int x = 0; x < N; x++) r = pset(r, x, pget(p, pget(q, x))); return r; }
static perm16 pinv(perm16 p) { perm16 r = 0; for (int x = 0; x < N; x++) r = pset(r, pget(p, x), x); return r; }
static perm16 ptrans(int a) { perm16 r = 0; for (int x = 0; x < N; x++) r = pset(r, x, x ^ a); return r; }
static int is_pow2(int k) { return k > 0 && !(k & (k - 1)); }
static int perm_order_pow2(perm16 p) {  /* order of p is a power of two iff every cycle length is */
    int seen = 0;
    for (int x = 0; x < N; x++) if (!(seen >> x & 1)) {
        int len = 0, y = x;
        do { seen |= 1 << y; y = pget(p, y); len++; } while (y != x);
        if (!is_pow2(len)) return 0;
    }
    return 1;
}
static int cmp_u64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }

static perm16 *T = NULL; static int nT = 0;
static perm16 *Talpha = NULL; static int *Tcls = NULL;   /* alpha pi alpha^{-1} = class rep, class id */
static int T_index(perm16 p) {  /* binary search, -1 if absent */
    int lo = 0, hi = nT - 1;
    while (lo <= hi) { int mid = (lo + hi) / 2; if (T[mid] == p) return mid; if (T[mid] < p) lo = mid + 1; else hi = mid - 1; }
    return -1;
}
static int in_T(perm16 p) { return T_index(p) >= 0; }

/* ------------------------------------------------------------ lock mode */
typedef struct { int level, sector, lam; int pi[N]; int f[N]; } root_t;

static z8 coefM[NL];               /* Pauli coefficient numerators (x16) of M1 */
static int use_product_filter = 0; static long long killers_seen = 0, nontraceless_seen = 0;
static z8 last_g16[NL];   /* product numerators x16 of the last product_passes call */
static void monomial_coeffs(const int *pi, const int *w8, z8 *out) {
    for (int v = 0; v < NL; v++) out[v] = Z0;
    for (int x = 0; x < N; x++) {
        int a = pi[x] ^ x;
        for (int b = 0; b < N; b++) {
            int v = a + 16 * b;
            int e = w8[x] + 2 * popc(a & b) + 4 * (popc(b & pi[x]) & 1);
            out[v] = zadd(out[v], zunit(e));
        }
    }
}
/* dense A_{y^a, y} = sum_b a_{(a,b)} i^{|a&b|} (-1)^{|b&y|}   (numerators x16) */
static void densify(const z8 *acoef, z8 A[N][N]) {
    for (int x = 0; x < N; x++) for (int y = 0; y < N; y++) A[x][y] = Z0;
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) {
        z8 c = acoef[a + 16 * b]; if (ziszero(c)) continue;
        for (int y = 0; y < N; y++) A[y ^ a][y] = zadd(A[y ^ a][y], zmul(c, zunit(2 * popc(a & b) + 4 * (popc(b & y) & 1))));
    }
}
static int rank_f2(const int *vecs, int n) {
    int basis[8], nb = 0;
    for (int i = 0; i < n; i++) { int v = vecs[i]; for (int k = 0; k < nb; k++) v = v < (v ^ basis[k]) ? v : (v ^ basis[k]); if (v) basis[nb++] = v; }
    return nb;
}
/* normalised-ray space of a support set: D = span of differences; return list of s with omega(s, D) = 0 */
static int normalized_rays(const int *supp, int ns, int *out) {
    int diffs[NL], nd = 0;
    for (int i = 1; i < ns; i++) diffs[nd++] = supp[i] ^ supp[0];
    int basis[8], nb = 0;
    for (int i = 0; i < nd; i++) { int v = diffs[i]; for (int k = 0; k < nb; k++) v = v < (v ^ basis[k]) ? v : (v ^ basis[k]); if (v) basis[nb++] = v; }
    int no = 0;
    for (int s = 1; s < NL; s++) { int ok = 1; for (int k = 0; k < nb; k++) if (omega(s, basis[k])) { ok = 0; break; } if (ok) out[no++] = s; }
    return no;
}
static void bits_or_rays(uint64_t *acc, const int *rays, int nr) {
    for (int i = 0; i < nr; i++) for (int k = 0; k < PBITS; k++) acc[k] |= planes_through[rays[i]][k];
}

/* union-find with Z/4 potentials */
static int ufp[N], ufpot[N];
static int uf_find(int x, int *pot) { int p = 0, y = x; while (ufp[y] != y) { p += ufpot[y]; y = ufp[y]; } *pot = ((p % 4) + 4) % 4; return y; }
/* impose f(x) - f(y) = r ; return 0 if inconsistent */
static int uf_union(int x, int y, int r) {
    int px, py; int rx = uf_find(x, &px), ry = uf_find(y, &py);
    /* f(x) = f(rx) + px ; f(y) = f(ry) + py ; want f(x)-f(y) = r  => f(rx) - f(ry) = r - px + py */
    int d = (((r - px + py) % 4) + 4) % 4;
    if (rx == ry) return d == 0;
    ufp[rx] = ry; ufpot[rx] = d;  /* f(rx) = f(ry) + d */
    return 1;
}

typedef struct { long long skel_checked, skel_consistent[2], skel_locked_generic[2], skel_fallback[2], fallback_solutions[2]; } stats_t;

/* Hermitian solvability and explicit enumeration of phase solutions for a skeleton.
   comp[x] = component id, delta[x] = potential wrt component root.  kappa in {0,3}.
   Returns number of solutions (0 if none); if enumerate != 0, calls the callback per solution. */
static int comp_of[N], delta_of[N], ncomp;
static int comp_root_pot[N];  /* helper */
#include "fast_helpers_v3.h"
typedef void (*sol_cb)(const int *fprime, void *ctx);
static int enum_abort;   /* defined with the enumeration context below */
static long long solve_hermitian(const int *pi, int kappa, sol_cb cb, void *ctx, long long maxenum) {
    /* variables t_c, c < ncomp. Equations: t_{comp(x)} + t_{comp(pi x)} = kappa - delta(x) - delta(pi x). */
    int alpha[N], sgn[N], assigned[N]; memset(assigned, 0, sizeof assigned);
    int part_root[N], npart = 0;
    int two_t_constraint[N]; for (int c = 0; c < N; c++) two_t_constraint[c] = -1;  /* value of 2*t_root for a part, or -1 */
    int part_of[N];
    /* BFS over components via equations */
    for (int c0 = 0; c0 < ncomp; c0++) if (!assigned[c0]) {
        int stack[N], sp = 0; stack[sp++] = c0; assigned[c0] = 1; alpha[c0] = 0; sgn[c0] = 1; part_of[c0] = npart; part_root[npart] = c0;
        while (sp) {
            int c = stack[--sp];
            for (int x = 0; x < N; x++) if (comp_of[x] == c) {
                int y = pi[x], cy = comp_of[y];
                int rhs = (((kappa - delta_of[x] - delta_of[y]) % 4) + 4) % 4;
                if (cy == c) {
                    /* 2 t_c = rhs  (self relation): t_c = alpha + sgn t0 => 2 alpha + 2 sgn t0 = rhs => 2 t0 = rhs - 2 alpha (mod 4, sgn irrelevant mod 4 for factor 2) */
                    int val = (((rhs - 2 * alpha[c]) % 4) + 4) % 4;
                    if (val & 1) return 0;
                    if (two_t_constraint[npart] < 0) two_t_constraint[npart] = val; else if (two_t_constraint[npart] != val) return 0;
                } else if (!assigned[cy]) {
                    assigned[cy] = 1; alpha[cy] = (((rhs - alpha[c]) % 4) + 4) % 4; sgn[cy] = -sgn[c]; part_of[cy] = npart; stack[sp++] = cy;
                } else {
                    /* t_c + t_cy = rhs: (alpha_c + s_c t0) + (alpha_cy + s_cy t0) = rhs */
                    int s = sgn[c] + sgn[cy];
                    int base = (((rhs - alpha[c] - alpha[cy]) % 4) + 4) % 4;
                    if (s == 0) { if (base != 0) return 0; }
                    else { /* +-2 t0 = base */ int val = base; if (val & 1) return 0;
                        if (two_t_constraint[npart] < 0) two_t_constraint[npart] = val; else if (two_t_constraint[npart] != val) return 0; }
                }
            }
        }
        npart++;
    }
    long long total = 1;
    for (int p = 0; p < npart; p++) total *= (two_t_constraint[p] < 0) ? 4 : 2;
    if (!cb) { fast_family_ready=0; if(fast_enabled && fast_counts && total>=64) fc_prepare(alpha,sgn,part_of,two_t_constraint,npart); }
    if (!cb || total > maxenum) return total;
    /* enumerate: choose t0 per part */
    int t0[N];
    long long idx;
    for (idx = 0; idx < total; idx++) {
        if (enum_abort) break;   /* the live candidate set is empty: nothing left to restrict */
        long long r = idx;
        for (int p = 0; p < npart; p++) {
            if (two_t_constraint[p] < 0) { t0[p] = r % 4; r /= 4; }
            else { t0[p] = (two_t_constraint[p] / 2) + 2 * (r % 2); r /= 2; }
        }
        int fprime[N];
        for (int x = 0; x < N; x++) { int c = comp_of[x]; int t = (alpha[c] + sgn[c] * t0[part_of[c]]) ; fprime[x] = (((t + delta_of[x]) % 4) + 4) % 4; }
        /* verify equations (defensive) */
        for (int x = 0; x < N; x++) if (((fprime[x] + fprime[pi[x]]) % 4) != kappa) { fprintf(stderr, "hermitian solve bug\n"); exit(1); }
        cb(fprime, ctx);
    }
    return total;
}

/* ---- atlas fingerprints (necessary class test for the product M1 h) */
static int nfp = 0, fplen[128]; static int64_t fpdata[128][64][4];
static int fp_npairs[128], fp_pairs[128][32][2], fp_has_pairs = 0;
static void load_fingerprints(const char *path) {
    FILE *f = fopen(path, "r"); if (!f) { perror(path); exit(1); }
    if (fscanf(f, "%d", &nfp) != 1 || nfp > 128) { fprintf(stderr, "bad fingerprint file\n"); exit(1); }
    fp_has_pairs = strstr(path, "fingerprints2") != NULL;
    for (int i = 0; i < nfp; i++) { fscanf(f, "%d", &fplen[i]); if (fplen[i] > 64) { fprintf(stderr, "fingerprint too long long\n"); exit(1); }
        for (int k = 0; k < fplen[i]; k++) for (int c = 0; c < 4; c++) { long long v; fscanf(f, "%lld", &v); fpdata[i][k][c] = v; }
        if (fp_has_pairs) { fscanf(f, "%d", &fp_npairs[i]); if (fp_npairs[i] > 32) { fprintf(stderr, "too many pairs\n"); exit(1); }
            for (int k = 0; k < fp_npairs[i]; k++) fscanf(f, "%d %d", &fp_pairs[i][k][0], &fp_pairs[i][k][1]); } }
    fclose(f);
}
static int cmp_z8(const void *a, const void *b) { const z8 *x = a, *y = b; for (int c = 0; c < 4; c++) { if (x->c[c] < y->c[c]) return -1; if (x->c[c] > y->c[c]) return 1; } return 0; }
static z8 zcanon(z8 a) { z8 b; for (int c = 0; c < 4; c++) b.c[c] = -a.c[c]; return cmp_z8(&a, &b) <= 0 ? a : b; }
static int fingerprint_match(z8 *num16 /* x16 numerators, modified */) {
    z8 fp[NL]; int n = 0;
    for (int u = 0; u < NL; u++) if (!ziszero(num16[u])) fp[n++] = zcanon(num16[u]);
    qsort(fp, n, sizeof(z8), cmp_z8);
    for (int i = 0; i < nfp; i++) {
        if (fplen[i] != n) continue;
        int ok = 1; for (int k = 0; k < n && ok; k++) for (int c = 0; c < 4; c++) if (fp[k].c[c] != fpdata[i][k][c]) { ok = 0; break; }
        if (ok) return i + 1;
    }
    return 0;
}
static long long prod_reject_pairs = 0;
static int fixed_rays_count(const int *supp, int ns);
/* (|F|,|N|) of a coefficient vector (x16 numerators) must be one of the pairs allowed for its fingerprint */
static int pairs_ok(const z8 *g16, int fpidx) {
    if (!fp_has_pairs) return 1;
    int supp[NL], ns = 0; for (int u = 0; u < NL; u++) if (!ziszero(g16[u])) supp[ns++] = u;
    int rays[NL]; int nN = normalized_rays(supp, ns, rays); int nF = fixed_rays_count(supp, ns);
    for (int k = 0; k < fp_npairs[fpidx]; k++) if (fp_pairs[fpidx][k][0] == nF && fp_pairs[fpidx][k][1] == nN) return 1;
    prod_reject_pairs++; return 0;
}
/* necessary conditions for g = M1 h (times i if anticommuting) to be a class member:
   traceless; Lemma 2.1 integrality (coefficients in Z/8 real, or sqrt2 Z/16, uniformly); atlas fingerprint */
static long long prod_reject_trace = 0, prod_reject_integral = 0, prod_reject_fp = 0;
static int msupp[NL], nmsupp = 0;   /* support of coefM, set per root */
static int product_passes(const z8 *hco16, int sign) {
    z8 g[NL];
    for(int u=0;u<NL;u++){
        z8 sum=Z0;
        for(int i=0;i<nmsupp;i++){
            int v=msupp[i];z8 t=hco16[u^v];if(ziszero(t))continue;
            sum=zadd(sum,zmul(zmul(coefM[v],t),zunit(2*kap(v,u^v))));
        }
        if(sign)sum=zmul(sum,zunit(2));
        if(u==0 && !ziszero(sum)){prod_reject_trace++;return 0;}
        for(int j=0;j<4;j++) {
            if(sum.c[j]%16){prod_reject_integral++;return 0;}
            sum.c[j]/=16;
        }
        g[u]=sum;
    }
    for(int u=0;u<NL;u++)last_g16[u]=g[u];
    int fi=fingerprint_match(g);if(!fi){prod_reject_fp++;return 0;}
    return pairs_ok(g,fi-1);
}

/* ---------------------------------------------------------------- exact tests (used only on killers) */
typedef struct { int perm[N]; int w[N]; } mono;   /* M|x> = zeta8^{w[x]} |perm[x]> */
static mono mono_mul(const mono *A, const mono *B) { mono r; for (int x = 0; x < N; x++) { r.perm[x] = A->perm[B->perm[x]]; r.w[x] = (B->w[x] + A->w[B->perm[x]]) % 8; } return r; }
static mono mono_dag(const mono *A) { mono r; for (int x = 0; x < N; x++) { r.perm[A->perm[x]] = x; r.w[A->perm[x]] = (8 - A->w[x]) % 8; } return r; }
static mono mono_pauli(int v) { mono r; int a = v & 15; for (int x = 0; x < N; x++) { r.perm[x] = x ^ a; r.w[x] = pw(v, x); } return r; }
static mono mono_corr(const mono *M, int v) { mono P = mono_pauli(v), D = mono_dag(M), t = mono_mul(M, &P); return mono_mul(&t, &D); }
static int mono_is_pauli(const mono *M) {
    int a = M->perm[0]; for (int x = 0; x < N; x++) if (M->perm[x] != (x ^ a)) return 0;
    int w0 = M->w[0], b = 0;
    for (int x = 0; x < N; x++) { int r = (M->w[x] - w0 + 8) % 8; if (r != 0 && r != 4) return 0; }
    for (int j = 0; j < 4; j++) if ((M->w[1 << j] - w0 + 8) % 8 == 4) b |= 1 << j;
    for (int x = 0; x < N; x++) if ((M->w[x] - w0 + 8) % 8 != 4 * (popc(b & x) & 1)) return 0;
    return 1;
}
static int mono_is_affine(const mono *M) {
    int c = M->perm[0], lin[4]; for (int j = 0; j < 4; j++) lin[j] = M->perm[1 << j] ^ c;
    for (int x = 0; x < N; x++) { int y = c; for (int j = 0; j < 4; j++) if ((x >> j) & 1) y ^= lin[j]; if (M->perm[x] != y) return 0; }
    return 1;
}
static int mono_is_clifford(const mono *M) {   /* all 255 corrections are Pauli rays (no generator shortcut) */
    if (!mono_is_affine(M)) return 0;
    for (int v = 1; v < NL; v++) { mono c = mono_corr(M, v); if (!mono_is_pauli(&c)) return 0; }
    return 1;
}

/* C2 is a group, so eight generator images suffice for C3.
   This shortcut is NOT used for levels above three. */
static int mono_clifford8(const mono *m) {
    for(int j=0;j<8;j++) {
        mono c=mono_corr(m,1<<j);
        if(!mono_is_pauli(&c))return 0;
    }
    return 1;
}
static int mono_low_tier(const mono *m) {
    if(mono_is_pauli(m))return 1;
    if(mono_clifford8(m))return 2;
    for(int j=0;j<8;j++) {
        mono c=mono_corr(m,1<<j);
        if(!mono_clifford8(&c))return 4;
    }
    return 3;
}
static int test_low_tiers(const char *path) {
    FILE *f=fopen(path,"r");if(!f)return 2;
    int count=0,id,level,sector,lam;
    while(fscanf(f,"%d %d %d %d",&id,&level,&sector,&lam)==4) {
        mono m;int a;
        for(int x=0;x<16;x++)if(fscanf(f,"%d",&m.perm[x])!=1)return 2;
        for(int x=0;x<16;x++){if(fscanf(f,"%d",&a)!=1)return 2;m.w[x]=(lam+2*a)&7;}
        int fast=mono_low_tier(&m), slow=4;
        if(mono_is_pauli(&m))slow=1;
        else if(mono_is_clifford(&m))slow=2;
        else {
            int ok=1;
            for(int v=1;v<256 && ok;v++) {mono c=mono_corr(&m,v);ok=mono_is_clifford(&c);}
            if(ok)slow=3;
        }
        if(fast!=slow || fast!=(level<4?level:4)){
            fprintf(stderr,"TIER_MISMATCH root=%d fast=%d slow=%d stored=%d\n",id,fast,slow,level);return 2;
        }
        count++;
    }
    fclose(f);if(count!=1020)return 2;
    printf("LOW_TIER_TEST_PASS roots=%d generator_vs_all_Paulis=equal\n",count);return 0;
}

/* exact class test for a dense Hermitian involution g given by x16 Pauli numerators:
   find a chamber L normalised by g, transport to the computational chart, require a mu8-monomial
   whose permutation lies in T.  Returns 1 and fills *gm if so. */
static void conv_pauli_right(const z8 *g, int u, z8 *out) { for (int w = 0; w < NL; w++) { z8 t = g[w ^ u]; if (ziszero(t)) { out[w] = Z0; continue; } out[w] = zmul(t, zunit(2 * kap(w ^ u, u))); } }  /* (g P_u)_w = g_{w^u} i^{kap(w^u,u)} */
static void conv_full(const z8 *a, const z8 *b, z8 *out) { for (int u = 0; u < NL; u++) { z8 s = Z0; for (int v = 0; v < NL; v++) { if (ziszero(a[v])) continue; z8 t = b[u ^ v]; if (ziszero(t)) continue; s = zadd(s, zmul(zmul(a[v], t), zunit(2 * kap(v, u ^ v)))); } out[u] = s; } }
static long long transport_noninvolution = 0;
static int class_indeterminate = 0;            /* set by exact_class_test when its self-check fails (fail closed) */
static long long transport_failures = 0;            /* killers kept because of an indeterminate class test */
static unsigned long long partner_digest = 0, partner_count = 0;   /* reference-check mode: FNV-1a digest of every traceless partner's ray list */
static unsigned long long fnv1a_rays(const int *rays, int nr) {
    unsigned long long h = 1469598103934665603ULL;
    for (int i = 0; i < nr; i++) { h ^= (unsigned long long)(rays[i] & 255); h *= 1099511628211ULL; }
    return h;
}
static int exact_class_test(const z8 *g16, mono *gm) {
    /* scan all chambers: g normalises A_L iff g is monomial in the chart C_L (4096-op relabel + densify) */
    for (int L = 0; L < nlag; L++) {
        const cliff *C = &CL[L];
        z8 acoef[NL]; for (int u = 0; u < NL; u++) acoef[u] = zmul(zunit(2 * C->eps[u]), g16[C->sig[u]]);
        z8 A[N][N]; densify(acoef, A);   /* numerators x16 x16 = x256 of g'?  densify expects x16 numerators, giving entries x16 */
        int perm[N], w[N], good = 1;
        for (int y = 0; y < N && good; y++) { int cnt = 0, xx = -1; for (int x = 0; x < N; x++) if (!ziszero(A[x][y])) { cnt++; xx = x; } if (cnt != 1) { good = 0; break; }
            /* entry must be 16 * zeta8^e */
            int e = -1; for (int t = 0; t < 8; t++) { z8 u = zunit(t); int m = 1; for (int c = 0; c < 4; c++) if (A[xx][y].c[c] != 16 * u.c[c]) { m = 0; break; } if (m) { e = t; break; } }
            if (e < 0) { good = 0; break; } perm[y] = xx; w[y] = e; }
        if (!good) continue;
        int seen = 0; for (int y = 0; y < N; y++) seen |= 1 << perm[y]; if (seen != 0xFFFF) continue;
        /* self-check: the transported matrix must be a Hermitian involution: z_y z_{perm y} = 1 and perm^2 = 1 */
        { int okinv = 1; for (int y = 0; y < N; y++) if (perm[perm[y]] != y || (w[y] + w[perm[y]]) % 8 != 0) okinv = 0;
          if (!okinv) { transport_noninvolution++; if (transport_noninvolution == 1) { fprintf(stderr, "TRANSPORT NOT INVOLUTION at L=%d: w:", L); for (int y = 0; y < N; y++) fprintf(stderr, " %d", w[y]); fprintf(stderr, " perm:"); for (int y = 0; y < N; y++) fprintf(stderr, " %d", perm[y]);
              /* dense checks on A (entries x16): A^2 == 256 I ?  A == A^dagger ? */
              int sq_ok = 1, herm_ok = 1;
              for (int x = 0; x < N; x++) for (int y = 0; y < N; y++) { z8 s = Z0; for (int k = 0; k < N; k++) s = zadd(s, zmul(A[x][k], A[k][y])); z8 want = Z0; if (x == y) want.c[0] = 256; if (!zeq(s, want)) sq_ok = 0;
                  z8 c = A[y][x]; z8 cj; cj.c[0] = c.c[0]; cj.c[1] = -c.c[3]; cj.c[2] = -c.c[2]; cj.c[3] = -c.c[1]; if (!zeq(A[x][y], cj)) herm_ok = 0; }
              fprintf(stderr, " | A^2=256I:%d hermitian:%d", sq_ok, herm_ok);
              /* also check the input g: g^2 = I via Pauli convolution of g16 with itself (x256) */
              { z8 gg[NL]; conv_full(g16, g16, gg); int gsq = 1; for (int u = 0; u < NL; u++) { z8 want = Z0; if (u == 0) want.c[0] = 256; if (!zeq(gg[u], want)) gsq = 0; } fprintf(stderr, " g^2=I:%d", gsq); }
              fprintf(stderr, "\nMISS_G"); for (int u = 0; u < NL; u++) if (!ziszero(g16[u])) fprintf(stderr, " %d:%lld,%lld,%lld,%lld", u, (long long)g16[u].c[0], (long long)g16[u].c[1], (long long)g16[u].c[2], (long long)g16[u].c[3]);
              fprintf(stderr, " | sig(z_j):"); for (int j = 0; j < 4; j++) fprintf(stderr, " %d", C->sig[16 << j]); fprintf(stderr, " lagbasis:"); for (int j = 0; j < 4; j++) fprintf(stderr, " %d", lagbasis[L][j]); fprintf(stderr, "\n"); }
              class_indeterminate = 1; return 0; } }   /* fail closed: the caller treats an indeterminate test as "not rejected" */
        perm16 p = 0; for (int y = 0; y < N; y++) p = pset(p, y, perm[y]);
        if (!in_T(p)) continue;
        for (int y = 0; y < N; y++) { gm->perm[y] = perm[y]; gm->w[y] = w[y]; }
        return 1;
    }
    return 0;
}
static long long exact_reject_partner_level = 0, exact_reject_product_class = 0, exact_reject_product_level = 0, exact_kept = 0, exact_tests = 0;
static long long level_decided_partner = 0, level_decided_product = 0;   /* rigidity rejections that the LEVEL coordinate alone decided (|N| and |F| tied with the root) */
static int only_chart = -1;   /* reference-check mode: process one chamber only and list every surviving plane */
static int root_level = 6, root_nrays = 0, root_nfixed = 0, use_rigidity = 0, candidates_inside = 0, fixplane_s1 = 0, fixplane_s2 = 0; static long long exact_reject_more_rigid = 0, table_miss_partner = 0, table_miss_product = 0;
static int fixed_rays_count(const int *supp, int ns) {   /* |supp^perp \ 0| */
    int basis[8], nb = 0;
    for (int i = 0; i < ns; i++) { int v = supp[i]; for (int k = 0; k < nb; k++) v = v < (v ^ basis[k]) ? v : (v ^ basis[k]); if (v) basis[nb++] = v; }
    int cnt = 0; for (int s = 1; s < NL; s++) { int ok = 1; for (int k = 0; k < nb; k++) if (omega(s, basis[k])) { ok = 0; break; } cnt += ok; }
    return cnt;
}
static double t_level = 0, t_class = 0;

/* The invariant of a phase family is
       S_j = S_0 ∩ ⋂_{h ∈ E_j} Met(h),   E_j = eligible solutions processed so far,
   implemented by restricting the LIVE candidate set immediately at every eligible solution
   (family_restrict); there is no accumulator, and the killer test reads the live set. A solution
   whose eligibility is unknown (table miss) restricts like an eligible one. `do_restrict` is 0 only
   in the fixplane killer-dump debug mode, which is not part of the certificate. */
typedef struct { const cliff *C; int lam; int sign; int ti; int L; const int *pi; long long count, kept, restricted; uint64_t *surv; int do_restrict;
                 int ptr[NL], npt; } enum_ctx;   /* ptr: chart rays r = (p,q) passing the point test pi(x+p) = pi(x)+p (phase-independent) */
static long long ray_mismatch = 0;   /* partners on which the direct ray test and the coefficient method disagree (must stay 0) */
static int enum_abort = 0;        /* set when the live candidate set becomes empty: the phase enumeration stops */
static int root_incomplete = 0;   /* number of phase families of the current root that were NOT enumerated (fiber too large) */
static void family_restrict(uint64_t *surv, const uint64_t *met) { for (int k = 0; k < PBITS; k++) surv[k] &= met[k]; }
static int bits_any(const uint64_t *b) { for (int k = 0; k < PBITS; k++) if (b[k]) return 1; return 0; }
static int debug_killers = 0, debug_printed = 0;
static void solution_cb(const int *fprime, void *vctx) {
    enum_ctx *ctx = (enum_ctx *)vctx;
    if(fast_enabled && fast_meets(ctx->pi,fprime)){ctx->count++;fast_quick++;return;}
    int w8[N]; for (int x = 0; x < N; x++) w8[x] = (ctx->lam + 2 * fprime[x]) % 8;
    { z8 tr = Z0; for (int x = 0; x < N; x++) if (ctx->pi[x] == x) tr = zadd(tr, zunit(w8[x])); if (!ziszero(tr)) { ctx->count++; nontraceless_seen++; return; } }  /* a frame element is traceless */
    ctx->count++;
    /* N(h) by the direct ray test: for a chart ray r = (p,q) passing the
       point test, H_f P_r H_f = +- P_r  iff  k_r(f;x) = f(x+p) - f(x) + 2 q.(pi x + x) (mod 4) is constant in x
       and even; the physical ray is sig[r]. No coefficient expansion is needed for Met(h). */
    int direct_fixed=0;
    uint64_t rmask[4] = {0, 0, 0, 0};
    for (int i = 0; i < ctx->npt; i++) {
        int r = ctx->ptr[i], p = r & 15, q = r >> 4, k0 = -1, ok = 1;
        for (int x = 0; x < N; x++) {
            int k = ((fprime[x ^ p] - fprime[x] + 2 * (popc(q & (ctx->pi[x] ^ x)) & 1)) % 4 + 4) % 4;
            if (k0 < 0) k0 = k; else if (k != k0) { ok = 0; break; }
        }
        if (ok && (k0 & 1) == 0) { if(k0==0)direct_fixed++; int s = ctx->C->sig[r]; rmask[s >> 6] |= 1ULL << (s & 63); }
    }
    int rays[NL]; int nr = 0;
    for (int s = 1; s < NL; s++) if ((rmask[s >> 6] >> (s & 63)) & 1) rays[nr++] = s;
    if (only_chart >= 0) {   /* reference-check contract: per-partner ray lists, and the coefficient method as a self-check */
        partner_digest += fnv1a_rays(rays, nr); partner_count++;
        z8 co[NL]; monomial_coeffs(ctx->pi, w8, co); z8 hco[NL]; for (int u = 0; u < NL; u++) hco[u] = Z0;
        for (int u = 0; u < NL; u++) if (!ziszero(co[u])) hco[ctx->C->sig[u]] = zmul(co[u], zunit(2 * ctx->C->eps[u]));
        int supp[NL], ns = 0; for (int u = 0; u < NL; u++) if (!ziszero(hco[u])) supp[ns++] = u;
        int rays2[NL]; int nr2 = normalized_rays(supp, ns, rays2);
        if (nr2 != nr || memcmp(rays, rays2, sizeof(int) * nr) != 0) ray_mismatch++;
    }
    uint64_t met[PBITS]; memset(met, 0, sizeof met); bits_or_rays(met, rays, nr);
    /* a solution that meets every LIVE candidate plane restricts nothing, so its eligibility is irrelevant */
    int killer = 0; for (int k = 0; k < PBITS; k++) if (ctx->surv[k] & ~met[k]) { killer = 1; break; }
    if (!killer) { ctx->kept++; return; }
    if (fast_enabled && use_product_filter && use_rigidity==3 &&
        (nr<root_nrays || (nr==root_nrays && direct_fixed<root_nfixed))) {
        exact_reject_more_rigid++;fast_rigid_skip++;return;
    }
    if (use_product_filter) {
        /* a killer: its coefficient vector is needed for the product and rigidity tests */
        z8 co[NL]; monomial_coeffs(ctx->pi, w8, co);
        z8 hco[NL]; for (int u = 0; u < NL; u++) hco[u] = Z0;
        for (int u = 0; u < NL; u++) if (!ziszero(co[u])) hco[ctx->C->sig[u]] = zmul(co[u], zunit(2 * ctx->C->eps[u]));   /* h = C M' C^-1 */
        if (!product_passes(hco, ctx->sign)) return;   /* fails a necessary class condition on the product: restricts nothing */
        int supp[NL], ns = 0;
        for (int u = 0; u < NL; u++) if (!ziszero(hco[u])) supp[ns++] = u;
        { int rays2[NL]; int nr2 = normalized_rays(supp, ns, rays2); if (nr2 != nr || memcmp(rays, rays2, sizeof(int) * nr) != 0) ray_mismatch++; }
        {   /* decide eligibility exactly; a rejected solution restricts nothing */
            exact_tests++;
            clock_t c0 = clock();
            mono hm;for(int x=0;x<N;x++){hm.perm[x]=ctx->pi[x];hm.w[x]=w8[x];}
            int plevel = (nr==root_nrays && direct_fixed==root_nfixed) ? mono_low_tier(&hm) : 4;
            t_level += (double)(clock() - c0) / CLOCKS_PER_SEC;
            if (plevel < 0) { table_miss_partner++; if (table_miss_partner == 1) fprintf(stderr, "level table miss (partner): ti=%d cls=%d sector=%d\n", ctx->ti, Tcls[ctx->ti], ctx->lam); }
            else if (use_rigidity == 3) {
                /* order O2: more rigid = fewer normalised rays, then fewer fixed rays, then LOWER level */
                int nf = fixed_rays_count(supp, ns);
                if (nr == root_nrays && nf == root_nfixed && plevel < root_level) level_decided_partner++;
                if (nr < root_nrays || (nr == root_nrays && (nf < root_nfixed || (nf == root_nfixed && plevel < root_level)))) { exact_reject_more_rigid++; return; }
            } else
            if (plevel > root_level) { exact_reject_partner_level++; return; }
            /* rigidity exclusion: a partner of the same level with a strictly smaller normalised-ray
               space is strictly more rigid than the root and would have been chosen as the frame's root */
            if (use_rigidity == 1 && plevel == root_level && nr < root_nrays) { exact_reject_more_rigid++; return; }
            if (use_rigidity == 2 && plevel == root_level) {
                /* rigidity = (level, exactly-fixed rays asc, normalised rays asc); F(h) = supp(h)^perp */
                int nf = fixed_rays_count(supp, ns);
                if (nf < root_nfixed || (nf == root_nfixed && nr < root_nrays)) { exact_reject_more_rigid++; return; }
            }
            int gsupp0[NL], gns0 = 0; for (int u = 0; u < NL; u++) if (!ziszero(last_g16[u])) gsupp0[gns0++] = u;
            int grays0[NL]; int gnr0 = normalized_rays(gsupp0, gns0, grays0); int gnf0 = fixed_rays_count(gsupp0, gns0);
            if (use_rigidity == 3 && (gnr0 < root_nrays || (gnr0 == root_nrays && gnf0 < root_nfixed))) { exact_reject_product_level++; return; }   /* product strictly more rigid: cheap reject before the chamber scan */
            c0 = clock();
            mono gm; class_indeterminate = 0; int cls = exact_class_test(last_g16, &gm);
            t_class += (double)(clock() - c0) / CLOCKS_PER_SEC;
            /* fail closed: a transport self-check failure makes the class test indeterminate, and an
               indeterminate test never rejects (the killer restricts as if eligible) */
            if (!cls && !class_indeterminate) { exact_reject_product_class++; return; }
            int glevel = -1, gti = -1, gsec = -1;
            if (cls) {
                glevel=(gnr0==root_nrays && gnf0==root_nfixed) ? mono_low_tier(&gm) : 4;
            } else transport_failures++;
            if (glevel < 0) {   /* conservative: keep the killer, record the miss */
                table_miss_product++;
                if (table_miss_product == 1) { fprintf(stderr, "level table miss (product): gti=%d gsec=%d sign=%d lam=%d L=%d w:", gti, gsec, ctx->sign, ctx->lam, ctx->L); for (int y = 0; y < N; y++) fprintf(stderr, " %d", gm.w[y]); fprintf(stderr, " perm:"); for (int y = 0; y < N; y++) fprintf(stderr, " %d", gm.perm[y]);
                    fprintf(stderr, "\nMISS_H"); for (int u = 0; u < NL; u++) if (!ziszero(hco[u])) fprintf(stderr, " %d:%lld,%lld,%lld,%lld", u, (long long)hco[u].c[0], (long long)hco[u].c[1], (long long)hco[u].c[2], (long long)hco[u].c[3]);
                    fprintf(stderr, "\nMISS_G"); for (int u = 0; u < NL; u++) if (!ziszero(last_g16[u])) fprintf(stderr, " %d:%lld,%lld,%lld,%lld", u, (long long)last_g16[u].c[0], (long long)last_g16[u].c[1], (long long)last_g16[u].c[2], (long long)last_g16[u].c[3]);
                    fprintf(stderr, "\n"); }
            } else if (use_rigidity == 3) {
                /* the product is a frame element too: its rigidity must not exceed the root's */
                int gsupp[NL], gns = 0; for (int u = 0; u < NL; u++) if (!ziszero(last_g16[u])) gsupp[gns++] = u;
                int grays[NL]; int gnr = normalized_rays(gsupp, gns, grays); int gnf = fixed_rays_count(gsupp, gns);
                if (gnr == root_nrays && gnf == root_nfixed && glevel < root_level) level_decided_product++;
                if (gnr < root_nrays || (gnr == root_nrays && (gnf < root_nfixed || (gnf == root_nfixed && glevel < root_level)))) { exact_reject_product_level++; return; }
            } else
            if (glevel > root_level) { exact_reject_product_level++; return; }
            if (debug_killers && debug_printed < 400) {
                debug_printed++;
                /* in fixplane mode only dump solutions that fail the fixed plane */
                if (fixplane_s1 > 0) { int id = planeid[fixplane_s1][fixplane_s2]; if ((met[id >> 6] >> (id & 63)) & 1) goto skip_dump; }
                fprintf(stderr, "KILLER L=%d ti=%d sign=%d kappa=%d plevel=%d glevel=%d nF=%d nN=%d rays:", ctx->L, ctx->ti, ctx->sign, ctx->lam ? 3 : 0, plevel, glevel, fixed_rays_count(supp, ns), nr);
                for (int i = 0; i < nr; i++) fprintf(stderr, " %d", rays[i]);
                fprintf(stderr, " | pi:"); for (int x = 0; x < N; x++) fprintf(stderr, " %d", ctx->pi[x]);
                fprintf(stderr, " | f:"); for (int x = 0; x < N; x++) fprintf(stderr, " %d", fprime[x]);
                fprintf(stderr, "\nHCOEF"); for (int u = 0; u < NL; u++) if (!ziszero(hco[u])) fprintf(stderr, " %d:%lld,%lld,%lld,%lld", u, (long long)hco[u].c[0], (long long)hco[u].c[1], (long long)hco[u].c[2], (long long)hco[u].c[3]);
                fprintf(stderr, "\n");
                skip_dump: ;
            }
            exact_kept++;
        }
    }
    /* an eligible killer (or one of unknown eligibility): restrict the live candidate set now */
    ctx->kept++; ctx->restricted++;
    if (nr == 0) killers_seen++;
    if (ctx->do_restrict) { family_restrict(ctx->surv, met); if (!bits_any(ctx->surv)) enum_abort = 1; }
}

static int mode_lock(const char *rootfile, const char *pdfile, const char *outfile, int first, int last) {
    FILE *f = fopen(pdfile, "rb"); if (!f) { perror(pdfile); return 1; }
    fread(&nT, sizeof(int), 1, f); T = malloc(sizeof(perm16) * nT); fread(T, sizeof(perm16), nT, f);
    Talpha = malloc(sizeof(perm16) * nT); Tcls = malloc(sizeof(int) * nT);
    if (fread(Talpha, sizeof(perm16), nT, f) != (size_t)nT || fread(Tcls, sizeof(int), nT, f) != (size_t)nT) { fprintf(stderr, "point domain file lacks conjugators (rebuild with pointdomain)\n"); return 1; }
    fclose(f);
    int *Tp = malloc(sizeof(int) * nT * N);
    for (int i = 0; i < nT; i++) for (int x = 0; x < N; x++) Tp[i * N + x] = pget(T[i], x);
    printf("point domain: %d involutions\n", nT);
    build_lagrangians(); build_planes(); build_chamber_cliffords();
    printf("lagrangians %d planes %d cliffords built and checked\n", nlag, nplanes);
    root_t *roots = malloc(sizeof(root_t) * 1100); int nr = 0;
    f = fopen(rootfile, "r"); if (!f) { perror(rootfile); return 1; }
    while (1) { int idx; root_t r; if (fscanf(f, "%d %d %d %d", &idx, &r.level, &r.sector, &r.lam) != 4) break;
        for (int x = 0; x < N; x++) fscanf(f, "%d", &r.pi[x]); for (int x = 0; x < N; x++) fscanf(f, "%d", &r.f[x]); roots[nr++] = r; }
    fclose(f);
    printf("roots: %d\n", nr);
    if(fast_enabled)fast_load_certificate();
    if (last < 0 || last >= nr) last = nr - 1;
    FILE *out = fopen(outfile, first == 0 ? "w" : "a");
    for (int ri = first; ri <= last; ri++) {
        root_t *R = &roots[ri];
        fast_early=fast_common=fast_four=fast_phase_skip=fast_quick=fast_rigid_skip=0;
        if(fast_enabled) {
            if(fast_cert_mode[ri]==3){ fprintf(out,"LARGEN_SKIPPED %d (baseline coverage retained)\n",ri);continue; }
            for(int z=0;z<3;z++)fast_plane[z]=fast_cert_plane[ri][z];
            use_product_filter=(fast_cert_mode[ri]==2);use_rigidity=3;
        }
        clock_t t0 = clock();
        root_level = 4; /* recomputed below; stored exact level is not trusted */

        exact_reject_partner_level = exact_reject_product_class = exact_reject_product_level = exact_kept = exact_tests = 0;
        int w8[N]; for (int x = 0; x < N; x++) w8[x] = (R->lam + 2 * R->f[x]) % 8;
        mono root_mono;
        for(int x=0;x<N;x++){root_mono.perm[x]=R->pi[x];root_mono.w[x]=w8[x];}
        root_level=mono_low_tier(&root_mono);
        monomial_coeffs(R->pi, w8, coefM);
        /* survivors: planes meeting N(M1) */
        uint64_t surv[PBITS]; memset(surv, 0, sizeof surv);
        { int supp[NL], ns = 0; for (int u = 0; u < NL; u++) if (!ziszero(coefM[u])) supp[ns++] = u;
          int rays[NL]; int n0 = normalized_rays(supp, ns, rays); root_nrays = n0; root_nfixed = fixed_rays_count(supp, ns);
          nmsupp = ns; for (int i = 0; i < ns; i++) msupp[i] = supp[i];
          if (fixplane_s1 > 0) {   /* a single prescribed candidate plane; killers are dumped, no early exit */
              int id = planeid[fixplane_s1][fixplane_s2]; if (id < 0) { fprintf(stderr, "bad fixplane\n"); exit(1); }
              surv[id >> 6] |= 1ULL << (id & 63);
          } else if (candidates_inside && n0 >= 3) {   /* candidate planes contained in N(M1) */
              for (int i = 0; i < n0; i++) for (int j = i + 1; j < n0; j++) { int id = planeid[rays[i]][rays[j]]; surv[id >> 6] |= 1ULL << (id & 63); }
          } else bits_or_rays(surv, rays, n0); }
        if(fast_enabled) {
            fast_id=planeid[fast_plane[0]][fast_plane[1]];
            int root_meets=0;
            for(int z=0;z<3;z++){
                int sign=-1,ok=1;
                for(int u=0;u<NL;u++)if(!ziszero(coefM[u])){
                    int v=omega(fast_plane[z],u);if(sign<0)sign=v;else if(sign!=v){ok=0;break;}
                }
                if(ok)root_meets=1;
            }
            if(!root_meets){fprintf(stderr,"proposed plane does not meet root\n");return 2;}
            memset(surv,0,sizeof surv);surv[fast_id>>6]=1ULL<<(fast_id&63);
        }
        exact_reject_more_rigid = 0; t_level = t_class = 0; root_incomplete = 0;
        long long restricting_cnt[2] = {0, 0};
        stats_t st; memset(&st, 0, sizeof st);
        for (int L = (only_chart >= 0 ? only_chart : fast_first_chart); L < (only_chart >= 0 ? only_chart + 1 : fast_last_chart+1); L++) {
            if(fast_subset && !fast_chart_selected[L])continue;
            const cliff *C = &CL[L];
            if(fast_enabled)for(int z=0;z<3;z++){
                fast_pulled[z]=-1;
                for(int u=1;u<NL;u++)if(C->sig[u]==fast_plane[z]){fast_pulled[z]=u;break;}
                if(fast_pulled[z]<0)return 2;
            }
            z8 acoef[NL]; for (int u = 0; u < NL; u++) acoef[u] = zmul(zunit(2 * C->eps[u]), coefM[C->sig[u]]);  /* a_u = i^{-e_u} m_{sig u}, e even */
            z8 A[N][N]; densify(acoef, A);
            if(fast_enabled) fast_colorize(A);
            unsigned rowmask[N]; for (int x = 0; x < N; x++) { rowmask[x] = 0; for (int y = 0; y < N; y++) if (!ziszero(A[x][y])) rowmask[x] |= 1u << y; }
            for (int ti = 0; ti < nT; ti++) {
                const int *pp = Tp + ti * N;
                if(fast_enabled && fast_generic(pp)){st.skel_checked++;fast_early++;continue;}
                /* support invariance: supp(A[pp x]) == pp(supp(A[x])) */
                int ok = 1;
                for (int x = 0; x < N && ok; x++) {
                    unsigned m = rowmask[x], pm = 0;
                    while (m) { int y = __builtin_ctz(m); m &= m - 1; pm |= 1u << pp[y]; }
                    if (pm != rowmask[pp[x]]) ok = 0;
                }
                st.skel_checked++;
                if (!ok) continue;
                /* ratios r_xy with A[pp x][pp y] = i^{r} A[x][y] */
                int rxy[N][N]; int consistent_ratio = 1;
                for (int x = 0; x < N && consistent_ratio; x++) for (int y = 0; y < N; y++) {
                    if (ziszero(A[x][y])) { rxy[x][y] = -1; continue; }
                    int found = -1;
                    if(fast_enabled) {
                        if(fast_color[x][y]==fast_color[pp[x]][pp[y]])
                            found=(fast_entry_phase[pp[x]][pp[y]]-fast_entry_phase[x][y])&3;
                    } else
                    for (int r = 0; r < 4; r++) if (zeq(A[pp[x]][pp[y]], zmul(zunit(2 * r), A[x][y]))) { found = r; break; }
                    if (found < 0) { consistent_ratio = 0; break; }
                    rxy[x][y] = found;
                }
                if (!consistent_ratio) continue;
                for (int sign = 0; sign < 2; sign++) {   /* sign 0: commute (eps=+1); sign 1: anticommute (multiply by i^2) */
                    for (int x = 0; x < N; x++) { ufp[x] = x; ufpot[x] = 0; }
                    int okc = 1;
                    for (int x = 0; x < N && okc; x++) for (int y = 0; y < N; y++) if (rxy[x][y] >= 0)
                        if (!uf_union(x, y, (rxy[x][y] + 2 * sign) % 4)) { okc = 0; break; }
                    if (!okc) continue;
                    /* components and potentials */
                    int rootid[N]; ncomp = 0;
                    for (int x = 0; x < N; x++) { int p; int rt = uf_find(x, &p); delta_of[x] = p; int c = -1; for (int k = 0; k < ncomp; k++) if (rootid[k] == rt) c = k; if (c < 0) { rootid[ncomp] = rt; c = ncomp++; } comp_of[x] = c; }
                    for (int kap8 = 0; kap8 < 2; kap8++) {   /* even sector kappa=0, lam=0; odd sector kappa=3, lam=1 */
                        int kappa = kap8 ? 3 : 0, lam = kap8 ? 1 : 0;
                        long long nsol = solve_hermitian(pp, kappa, NULL, NULL, 0);
                        if (nsol == 0) continue;
                        st.skel_consistent[sign]++;
                        if(fast_enabled && fast_family_ready && fc_cover(pp)){
                            fast_phase_skip+=nsol;continue;
                        }
                        /* generic lock: maximal support S_max = { (a,b) : F_a nonempty } mapped by sig */
                        int hasa[N]; memset(hasa, 0, sizeof hasa); for (int x = 0; x < N; x++) hasa[pp[x] ^ x] = 1;
                        int supp[NL], ns = 0; for (int a = 0; a < N; a++) if (hasa[a]) for (int b = 0; b < N; b++) supp[ns++] = C->sig[a + 16 * b];
                        int rays[NL]; int nray = normalized_rays(supp, ns, rays);
                        uint64_t met[PBITS]; memset(met, 0, sizeof met); bits_or_rays(met, rays, nray);
                        int generic_ok = 1;
                        for (int k = 0; k < PBITS; k++) if (surv[k] & ~met[k]) { generic_ok = 0; break; }
                        if (generic_ok && only_chart < 0) { st.skel_locked_generic[sign]++; continue; }   /* reference-check mode enumerates every family explicitly */
                        /* fallback: explicit phase enumeration */
                        st.skel_fallback[sign]++;
                        enum_ctx ctx; ctx.C = C; ctx.lam = lam; ctx.sign = sign; ctx.ti = ti; ctx.L = L; ctx.pi = pp; ctx.count = 0; ctx.kept = 0; ctx.restricted = 0; ctx.surv = surv; ctx.do_restrict = (fixplane_s1 == 0);
                        ctx.npt = 0;   /* chart rays passing the phase-independent point test for this permutation */
                        for (int r = 1; r < NL; r++) { int p = r & 15, okp = 1; for (int x = 0; x < N; x++) if (pp[x ^ p] != (pp[x] ^ p)) { okp = 0; break; } if (okp) ctx.ptr[ctx.npt++] = r; }
                        enum_abort = 0;
                        long long got = solve_hermitian(pp, kappa, solution_cb, &ctx, 1L << 24);
                        if (got > (1L << 24)) { root_incomplete++; fprintf(out, "ROOT %d FIBER_TOO_LARGE L=%d ti=%d sign=%d kappa=%d fiber=%lld\n", ri, L, ti, sign, kappa, got); fflush(out); continue; }
                        st.fallback_solutions[sign] += ctx.count; restricting_cnt[sign] += ctx.restricted;
                        /* the candidate set was restricted in place by every eligible killer (family_restrict) */
                    }
                }
            }
            /* early exit if no plane survives (not in fixplane mode: there we want every killer) */
            int any = 0; for (int k = 0; k < PBITS; k++) if (surv[k]) { any = 1; break; }
            if (!any && fixplane_s1 == 0) break;
        }
        int nsurv = 0, ex[8];
        for (int id = 0; id < nplanes; id++) if ((surv[id >> 6] >> (id & 63)) & 1) { if (nsurv < 8) ex[nsurv] = id; nsurv++; }
        int fast_failed = root_incomplete || transport_failures || ray_mismatch;
        /* an unenumerated family means the root is NOT decided: the survivor field says so and no count is printed */
        char survstr[16]; if (root_incomplete) strcpy(survstr, "INCOMPLETE"); else sprintf(survstr, "%d", nsurv);
        fprintf(out, "ROOT %d level %d sector %d survivors %s skel_checked %lld consistent %lld/%lld generic %lld/%lld fallback %lld/%lld fallback_solutions %lld/%lld prod_reject trace %lld integral %lld fp %lld pairs %lld killers %lld time %.6fs planes",
                ri, R->level, R->sector, survstr, st.skel_checked, st.skel_consistent[0], st.skel_consistent[1], st.skel_locked_generic[0], st.skel_locked_generic[1],
                st.skel_fallback[0], st.skel_fallback[1], st.fallback_solutions[0], st.fallback_solutions[1], prod_reject_trace, prod_reject_integral, prod_reject_fp, prod_reject_pairs, killers_seen, (double)(clock() - t0) / CLOCKS_PER_SEC);
        prod_reject_trace = prod_reject_integral = prod_reject_fp = prod_reject_pairs = killers_seen = 0;
        fprintf(out, " exact_tests %lld reject_partner_level %lld reject_more_rigid %lld reject_product_class %lld reject_product_level %lld exact_kept %lld root_nrays %d",
                exact_tests, exact_reject_partner_level, exact_reject_more_rigid, exact_reject_product_class, exact_reject_product_level, exact_kept, root_nrays);
        fprintf(out, " t_level %.0f t_class %.0f nontraceless %lld table_miss %lld/%lld", t_level, t_class, nontraceless_seen, table_miss_partner, table_miss_product); nontraceless_seen = table_miss_partner = table_miss_product = 0;
        fprintf(out, " restricting %lld/%lld incomplete_families %d transport_failures %lld kernel publication-lowtier-v1 scope %s", restricting_cnt[0], restricting_cnt[1], root_incomplete, transport_failures,
                (!fast_subset && only_chart < 0 && fast_first_chart==0 && fast_last_chart==NLAG-1 && fixplane_s1 == 0 && !debug_killers) ? "full" : "diag");
        transport_failures = 0;
        fprintf(out, " ray_mismatch %lld level_decided %lld/%lld", ray_mismatch, level_decided_partner, level_decided_product); ray_mismatch = 0; level_decided_partner = level_decided_product = 0;
        if (only_chart >= 0) { fprintf(out, " partners %llu partner_digest %016llx", partner_count, partner_digest); partner_count = 0; partner_digest = 0; }
        for (int k = 0; k < nsurv && k < 8; k++) fprintf(out, " (%d,%d,%d)", planes[ex[k]][0], planes[ex[k]][1], planes[ex[k]][2]);
        if (only_chart >= 0) { fprintf(out, " lagbasis %d %d %d %d allplanes", lagbasis[only_chart][0], lagbasis[only_chart][1], lagbasis[only_chart][2], lagbasis[only_chart][3]); for (int id = 0; id < nplanes; id++) if ((surv[id >> 6] >> (id & 63)) & 1) fprintf(out, " %d,%d,%d", planes[id][0], planes[id][1], planes[id][2]); }
        fprintf(out," guided %d charts %d:%d early %lld common %lld four %lld phase_skip %lld quick %lld",
                fast_enabled,fast_first_chart,fast_last_chart,fast_early,fast_common,fast_four,fast_phase_skip,fast_quick);
        fprintf(out," early_rigid %lld",fast_rigid_skip);
        fprintf(out, "\n"); fflush(out);
        if(fast_enabled && (nsurv!=1 || fast_failed)){
            fprintf(stderr,"guided acceptance failed\n");fclose(out);return 2;
        }
        printf("root %d done: survivors %d (%.0fs)\n", ri, nsurv, (double)(clock() - t0) / CLOCKS_PER_SEC); fflush(stdout);
    }
    fclose(out);
    return 0;
}

/* ---------------------------------------------------------------- synthetic test of the family invariant
   Cases over one or two candidate planes (bits 0 and 1) and a list of (met, eligible) solutions. The
   repaired update restricts the live set at every eligible solution; the legacy accumulator (kept here
   only as a test oracle) unions the met sets of eligible solutions and intersects once at the end. */
static void legacy_family(uint64_t *surv, const uint64_t (*mets)[PBITS], const int *elig, int n) {
    uint64_t acc[PBITS]; memset(acc, 0, sizeof acc); int kept = 0;
    for (int i = 0; i < n; i++) if (elig[i]) { for (int k = 0; k < PBITS; k++) acc[k] |= mets[i][k]; kept++; }
    if (kept > 0) for (int k = 0; k < PBITS; k++) surv[k] &= acc[k];
}
static void repaired_family(uint64_t *surv, const uint64_t (*mets)[PBITS], const int *elig, int n, int *processed) {
    *processed = 0;
    for (int i = 0; i < n; i++) { (*processed)++; if (!bits_any(surv)) break;   /* early termination: nothing left to restrict */
        if (!elig[i]) continue; family_restrict(surv, mets[i]); }
}
static int test_family(void) {
    /* met sets: A meets plane 0 only, B misses plane 0 (meets plane 1 only), Z meets nothing, AB meets both */
    uint64_t A[PBITS] = {1}, B[PBITS] = {2}, Z[PBITS] = {0}, AB[PBITS] = {3};
    struct { const char *name; uint64_t s0; int n; const uint64_t *m[3]; int e[3]; uint64_t want; int legacy_must_differ; } cases[] = {
        {"meets then misses (one plane, both eligible)",          1, 2, {A, B}, {1, 1}, 0, 1},
        {"misses then meets (reversed order)",                    1, 2, {B, A}, {1, 1}, 0, 1},
        {"two solutions meeting disjoint plane sets",             3, 2, {A, B}, {1, 1}, 0, 1},
        {"no eligible solution: candidates unchanged",            3, 2, {Z, Z}, {0, 0}, 3, 0},
        {"ineligible failure restricts nothing",                  1, 1, {Z},    {0},    1, 0},
        {"unresolved lookup treated as eligible",                 1, 1, {Z},    {1},    0, 0},
        {"early termination after the set empties",               1, 3, {B, AB, AB}, {1, 1, 1}, 0, 1},
        {"eligible solution meeting everything keeps the set",    3, 1, {AB},   {1},    3, 0},
    };
    int ncase = sizeof cases / sizeof cases[0], fails = 0;
    for (int c = 0; c < ncase; c++) {
        uint64_t mets[3][PBITS]; memset(mets, 0, sizeof mets); for (int i = 0; i < cases[c].n; i++) memcpy(mets[i], cases[c].m[i], sizeof mets[i]);
        uint64_t s1[PBITS] = {0}, s2[PBITS] = {0}; s1[0] = s2[0] = cases[c].s0; int processed;
        repaired_family(s1, mets, cases[c].e, cases[c].n, &processed);
        legacy_family(s2, mets, cases[c].e, cases[c].n);
        int ok = (s1[0] == cases[c].want) && (!cases[c].legacy_must_differ || s2[0] != cases[c].want);
        if (c == 6 && processed != 2) ok = 0;   /* the third solution must not be processed */
        printf("%-56s repaired %llu (want %llu) legacy %llu processed %d : %s\n", cases[c].name, (unsigned long long)s1[0], (unsigned long long)cases[c].want, (unsigned long long)s2[0], processed, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }
    printf(fails ? "FAMILY TEST FAILED (%d)\n" : "FAMILY TEST PASSED: repaired invariant correct on all %d cases, legacy accumulator wrong where expected\n", fails ? fails : ncase);
    return fails ? 1 : 0;
}

#ifdef FAST_TESTS
#include "fast_tests.h"
#endif
int main(int argc, char **argv) {
    if(argc>=3 && !strcmp(argv[1],"testtiers"))return test_low_tiers(argv[2]);
#ifdef FAST_TESTS
    if(argc>=2 && !strcmp(argv[1],"fasttests"))return fast_tests();
#endif
    if (argc >= 2 && !strcmp(argv[1], "testfamily")) return test_family();

    if (argc >= 5 && !strcmp(argv[1], "lock")) {
        int first = argc > 5 ? atoi(argv[5]) : 0, last = argc > 6 ? atoi(argv[6]) : -1;
        int optstart = 7;
        if (argc > 7 && strstr(argv[7], ".txt")) { load_fingerprints(argv[7]); use_product_filter = 1; printf("product filter on: %d fingerprints\n", nfp);
            optstart = 8; /* no hierarchy-level table is loaded */ }
        for (int a = optstart; a < argc; a++) {
            if(!strcmp(argv[a],"guided"))fast_enabled=1;
            if(!strncmp(argv[a],"chartlist=",10)){
                FILE *cf=fopen(argv[a]+10,"r");if(!cf)return 2;
                int x;fast_subset=1;
                while(fscanf(cf,"%d",&x)==1){if(x<0 || x>=NLAG || fast_chart_selected[x])return 2;fast_chart_selected[x]=1;}
                fclose(cf);
            }
            if(!strcmp(argv[a],"nocount"))fast_counts=0;
            if(!strncmp(argv[a],"chrange=",8)){
                if(sscanf(argv[a]+8,"%d,%d",&fast_first_chart,&fast_last_chart)!=2 ||
                    fast_first_chart<0 || fast_last_chart>=NLAG || fast_first_chart>fast_last_chart)return 2;
            }
            if (!strcmp(argv[a], "rigid")) { use_rigidity = 1; printf("rigidity exclusion on (level, N)\n"); }
            if (!strcmp(argv[a], "rigid2")) { use_rigidity = 2; printf("rigidity exclusion on (level, F, N)\n"); }
            if (!strcmp(argv[a], "rigid3")) { use_rigidity = 3; printf("rigidity order O2: (N asc, F asc, level asc)\n"); }
            if (!strcmp(argv[a], "debug")) { debug_killers = 1; }
            if (!strncmp(argv[a], "fixplane=", 9)) { sscanf(argv[a] + 9, "%d,%d", &fixplane_s1, &fixplane_s2); debug_killers = 1; debug_printed = -1000000000; printf("fixed plane through %d,%d; killers dumped\n", fixplane_s1, fixplane_s2); }
            if (!strcmp(argv[a], "inside")) { candidates_inside = 1; printf("candidate planes inside N(M1)\n"); }
            if (!strncmp(argv[a], "onlychart=", 10)) { only_chart = atoi(argv[a] + 10); printf("reference-check mode: chamber %d only, all surviving planes listed\n", only_chart); }
        }
        return mode_lock(argv[2], argv[3], argv[4], first, last);
    }
    if (argc >= 2 && !strcmp(argv[1], "cliffall")) {   /* dump the label actions of all chambers: "L a b c d" then 256 lines "u sig eps" */
        build_lagrangians(); build_planes(); build_chamber_cliffords();
        for (int L = 0; L < nlag; L++) { printf("L %d basis %d %d %d %d\n", L, lagbasis[L][0], lagbasis[L][1], lagbasis[L][2], lagbasis[L][3]);
            for (int u = 0; u < NL; u++) printf("%d %d %d\n", u, CL[L].sig[u], CL[L].eps[u]); }
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "cliff")) {   /* dump the label action of chamber L */
        build_lagrangians(); build_planes(); build_chamber_cliffords();
        int L = atoi(argv[2]); printf("L %d basis %d %d %d %d\n", L, lagbasis[L][0], lagbasis[L][1], lagbasis[L][2], lagbasis[L][3]);
        for (int u = 0; u < NL; u++) printf("%d %d %d\n", u, CL[L].sig[u], CL[L].eps[u]);
        return 0;
    }
    if (argc >= 2 && !strcmp(argv[1], "selftest")) {
        build_lagrangians(); build_planes(); build_chamber_cliffords();
        printf("selftest: lagrangians %d planes %d cliffords ok\n", nlag, nplanes);
        /* kappa formula check against explicit products on |x> for all v,w */
        for (int v = 0; v < NL; v++) for (int w = 0; w < NL; w++) {
            int a = v & 15, c = w & 15;
            for (int x = 0; x < N; x++) {
                int lhs = (pw(w, x) + pw(v, x ^ c)) % 8;               /* P_v P_w |x> weight */
                int rhs = (pw(v ^ w, x) + 2 * kap(v, w)) % 8;            /* i^kap P_{v^w}|x> */
                if (lhs != rhs || (a ^ c) != ((v ^ w) & 15)) { printf("kappa mismatch v=%d w=%d x=%d\n", v, w, x); return 1; }
            }
        }
        printf("kappa formula verified on all 65536 pairs and all x\n");
        return 0;
    }
    fprintf(stderr, "usage: ge4lock pointdomain <tri.txt> <classes.txt> <out.bin> | lock <roots.txt> <pd.bin> <out.txt> [first last] | selftest\n");
    return 1;
}
