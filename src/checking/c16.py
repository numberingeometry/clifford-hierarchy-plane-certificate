"""Fact (C1): stationary gates CS_yz^a CCZ^b and complete affine flags. Exact.

Basis order is |xyz>, x the most significant bit. Pauli labels are packed (x_bits << 3) | z_bits with
both three-bit words in this order, and the Hermitian representatives are i^(x.z) X^x Z^z. Operators act
on input columns; products compose right to left. Matrix entries are Gaussian integers and phases are
exponents in Z/4; there are no tolerances or eigenvector approximations.

Frame ordering (argued in the paper): a joint eigenbasis of independent commuting Hermitian Paulis admits
Clifford column phases, and its eigenvalue ordering gives one Clifford frame. Every other Clifford frame
for the same Pauli MASA changes this ordering by an affine permutation, and the implication tested over
all affine alpha and all complete affine flags is invariant under that change. So matching projectors
covers every Clifford frame of each witness. Internal key: C16.
"""

from itertools import combinations, permutations
import json
import time


class CheckFailure(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise CheckFailure(message)


IDENTITY = tuple(range(8))
ZERO = (0,0)
ONE = (1,0)


def rotate(value, exponent):
    r,i = value
    return ((r,i),(-i,r),(-r,-i),(i,-r))[exponent % 4]


def add(left,right):
    return left[0]+right[0],left[1]+right[1]


def multiply(left,right):
    return left[0]*right[0]-left[1]*right[1],left[0]*right[1]+left[1]*right[0]


def mono_compose(left,right):
    p,f = left
    q,g = right
    return tuple(p[q[j]] for j in range(8)),tuple((g[j]+f[q[j]])%4 for j in range(8))


def mono_inverse(operator):
    p,f = operator
    inverse = [0]*8
    for j,value in enumerate(p):
        inverse[value]=j
    return tuple(inverse),tuple(-f[inverse[j]]%4 for j in range(8))


def mono_conjugate(operator,target):
    return mono_compose(mono_compose(operator,target),mono_inverse(operator))


def pauli(label):
    x,z = label >> 3,label & 7
    return tuple(j^x for j in range(8)),tuple(((x&z).bit_count()+2*(z&j).bit_count())%4 for j in range(8))


PAULIS = tuple(pauli(label) for label in range(64))
GENERATORS = (8,16,32,1,2,4)


def symplectic(left,right):
    x,z = left>>3,left&7
    a,b = right>>3,right&7
    return ((x&b).bit_count()+(z&a).bit_count())&1


def span(vectors):
    values = {0}
    for vector in vectors:
        values |= {value^vector for value in tuple(values)}
    return frozenset(values)


def basis_of(values):
    basis = []
    covered = {0}
    for value in sorted(values):
        if value not in covered:
            basis.append(value)
            covered |= {x^value for x in tuple(covered)}
    return tuple(basis)


def lagrangians():
    out = set()
    for basis in combinations(range(1,64),3):
        if any(symplectic(a,b) for a,b in combinations(basis,2)):
            continue
        L = span(basis)
        if len(L) == 8:
            out.add(L)
    require(len(out) == 135,"Lagrangian domain is incomplete")
    require(all(len(L)==8 and len(basis_of(L))==3 for L in out),"invalid Lagrangian dimension")
    require(all(not any(symplectic(a,b) for a in L for b in L) for L in out),"nonisotropic Lagrangian")
    return sorted(out,key=lambda L:tuple(sorted(L)))


def support(operator):
    p,f = operator
    out = []
    for label,(q,g) in enumerate(PAULIS):
        total = ZERO
        for j in range(8):
            if p[j] == q[j]:
                total = add(total,rotate(ONE,f[j]-g[j]))
        if total != ZERO:
            out.append(label)
    return tuple(out)


def is_clifford(operator):
    return all(len(support(mono_conjugate(operator,PAULIS[label])))==1 for label in GENERATORS)


def is_third_level(operator):
    return all(is_clifford(mono_conjugate(operator,PAULIS[label])) for label in GENERATORS)


def witnesses(operator,lags):
    correction_supports = {label:support(mono_conjugate(operator,PAULIS[label])) for label in range(64)}
    found = []
    for L in lags:
        union = set().union(*(correction_supports[label] for label in L))
        M = span(union)
        if len(M)==8 and not any(symplectic(a,b) for a in M for b in M):
            found.append((L,M))
    return found


def matrix_identity():
    return tuple(ONE if r==c else ZERO for r in range(8) for c in range(8))


def projector_numerator(basis,eigenvalue_bits):
    # N=8*Pi is the exact product of (I+(-1)^bit P_b) for three generators.
    N = matrix_identity()
    for label,bit in zip(basis,eigenvalue_bits):
        p,f = PAULIS[label]
        sign_exponent = 2*bit
        N = tuple(add(N[8*r+c],rotate(N[8*r+p[c]],f[c]+sign_exponent))
                  for r in range(8) for c in range(8))
    require(all(N[8*r+c] == (N[8*c+r][0],-N[8*c+r][1]) for r in range(8) for c in range(8)),"stabilizer projector is not Hermitian")
    trace = ZERO
    for j in range(8):
        trace = add(trace,N[9*j])
    require(trace == (8,0),"stabilizer projector trace does not give rank one")
    for r in range(8):
        for c in range(8):
            entry = ZERO
            for k in range(8):
                entry = add(entry,multiply(N[8*r+k],N[8*k+c]))
            require(entry == (8*N[8*r+c][0],8*N[8*r+c][1]),"stabilizer projector is not idempotent")
    return N


def matrix_conjugate_monomial(operator,matrix):
    p,f = operator
    out = [ZERO]*64
    for r in range(8):
        for c in range(8):
            out[8*p[r]+p[c]] = rotate(matrix[8*r+c],f[r]-f[c])
    return tuple(out)


def frame_permutation(operator,L):
    basis = basis_of(L)
    require(len(basis)==3,"frame basis has the wrong dimension")
    projectors = tuple(projector_numerator(basis,((j>>2)&1,(j>>1)&1,j&1)) for j in range(8))
    lookup = {projector:j for j,projector in enumerate(projectors)}
    require(len(lookup)==8,"stabilizer eigenprojectors are duplicated")
    require(all(sum(projectors[j][k][axis] for j in range(8)) == 8*matrix_identity()[k][axis]
                for k in range(64) for axis in (0,1)),"stabilizer projectors do not sum to identity")
    permutation = []
    for projector in projectors:
        image = matrix_conjugate_monomial(operator,projector)
        require(image in lookup,"stationary gate does not permute its input-frame projectors")
        permutation.append(lookup[image])
    require(len(set(permutation))==8,"frame action is not a permutation")
    return tuple(permutation),basis


def compose(p,q):
    return tuple(p[q[j]] for j in range(8))


def affine_values(p):
    return all(p[j] == p[0] ^ (p[4]^p[0] if j&4 else 0) ^ (p[2]^p[0] if j&2 else 0)
               ^ (p[1]^p[0] if j&1 else 0) for j in range(8))


def affine_permutations():
    out = [p for p in permutations(range(8)) if affine_values(p)]
    require(len(out)==1344,"affine-permutation domain is incomplete")
    return out


def cosets(W):
    return frozenset(frozenset(x^w for w in W) for x in range(8))


def complete_flags():
    lines = [frozenset((0,v)) for v in range(1,8)]
    planes = [frozenset(j for j in range(8) if not ((normal&j).bit_count()&1)) for normal in range(1,8)]
    flags = tuple((cosets(line),cosets(plane)) for line in lines for plane in planes if line <= plane)
    require(len(flags)==len(set(flags))==21,"complete-affine-flag domain is incomplete")
    return flags


def flag_mask(permutation,flags):
    mask = 0
    for index,flag in enumerate(flags):
        if all(frozenset(permutation[j] for j in coset) in partition for partition in flag for coset in partition):
            mask |= 1<<index
    return mask


def validate_summary(summary):
    expected = {"input_lagrangians":135,"complete_affine_flags":21,"affine_permutations":1344,
                "witness_counts":{"I":135,"CS_yz":3,"CCZ":15,"CS_yz_CCZ":3},
                "frames_checked":156,"frame_affine_pairs_checked":209664,
                "implication_failures":0,"counterexample_witness_count":3}
    for key,value in expected.items():
        require(summary.get(key)==value,f"C16 {key}: expected {value!r}, received {summary.get(key)!r}")
    require(summary.get("all_four_gates_stationary") is True,"C16 stationarity failed")
    require(summary.get("counterexample_stationary") is True,"C16 counterexample is not stationary")
    require(summary.get("counterexample_third_level") is True,"C16 counterexample is not third-level")
    require(summary.get("counterexample_clifford") is False,"C16 counterexample unexpectedly Clifford")
    require(summary.get("counterexample_product_flag_indices"),"C16 counterexample antecedent fails")
    require(summary.get("counterexample_common_flag_indices")==[],"C16 counterexample consequent holds")


def run():
    started = time.perf_counter()
    lags = lagrangians()
    affines = affine_permutations()
    flags = complete_flags()
    masks = {}
    def mask(p):
        if p not in masks:
            masks[p] = flag_mask(p,flags)
        return masks[p]
    require(all(mask(tuple(j^translation for j in range(8)))==(1<<21)-1 for translation in range(8)),"translations do not preserve every declared affine flag")
    gates = {
        "I":(IDENTITY,(0,)*8),
        "CS_yz":(IDENTITY,tuple(((j>>1)&1)*(j&1) for j in range(8))),
        "CCZ":(IDENTITY,tuple(2*((j>>2)&1)*((j>>1)&1)*(j&1) for j in range(8))),
        "CS_yz_CCZ":(IDENTITY,tuple((((j>>1)&1)*(j&1)+2*((j>>2)&1)*((j>>1)&1)*(j&1))%4 for j in range(8))),
    }
    counts = {}
    frames = []
    pairs_checked = failures = 0
    for name,operator in gates.items():
        found = witnesses(operator,lags)
        require(all(L==M for L,M in found),f"{name} is not stationary")
        counts[name] = len(found)
        for L,M in found:
            permutation,basis = frame_permutation(operator,L)
            pmask = mask(permutation)
            hypothesis_count = 0
            for alpha in affines:
                pairs_checked += 1
                if mask(compose(alpha,permutation)):
                    hypothesis_count += 1
                    if not (mask(alpha)&pmask):
                        failures += 1
            require(failures==0,f"C16 implication failed in a {name} witness frame")
            frames.append({"gate":name,"input_lagrangian":sorted(L),"hermitian_pauli_basis":list(basis),
                           "eigenvalue_frame_permutation":list(permutation),"affine_maps_checked":len(affines),
                           "antecedent_satisfying_affine_maps":hypothesis_count})
    # Fixed paper control; no nondeterministic class representative or group cap.
    counter_permutation = tuple(4*((((j>>2)&1)^((j>>1)&1)^(j&1)^1))
                                +2*((((j>>2)&1)^((j>>1)&1)^1))
                                +(((j>>2)&1)^(j&1)^1) for j in range(8))
    counter_phases = tuple(((j>>2)&1)*(((j>>1)&1)^(j&1)) for j in range(8))
    D = counter_permutation,counter_phases
    alpha = tuple(4*(j&1)+2*(((j>>1)&1)^(j&1))+((j>>2)&1) for j in range(8))
    require(counter_permutation==(7,2,1,4,0,5,6,3),"paper counterexample coordinate convention changed")
    require(counter_phases==(0,0,0,0,0,1,1,0),"paper counterexample phase convention changed")
    require(alpha in affines,"paper counterexample alpha is not affine")
    counter_witnesses = witnesses(D,lags)
    require(all(L==M for L,M in counter_witnesses),"paper counterexample is not stationary")
    diagonal_lagrangian = frozenset(range(8))
    require(any(L==diagonal_lagrangian for L,M in counter_witnesses),"computational frame is not a counterexample witness")
    product_mask = mask(compose(alpha,counter_permutation))
    common_mask = mask(alpha)&mask(counter_permutation)
    result = {
        "fact":"stationary-cs-ccz-gates-and-complete-affine-flags",
        "route":"recompute","status":"accepted",
        "input_lagrangians":len(lags),"complete_affine_flags":len(flags),"affine_permutations":len(affines),
        "witness_counts":counts,"all_four_gates_stationary":True,"frames_checked":len(frames),
        "frame_affine_pairs_checked":pairs_checked,"implication_failures":failures,
        "frame_ordering_analytic_premise":"Clifford eigenvalue frames for one Pauli MASA differ by affine orderings; universal affine-map and flag testing is invariant under those changes",
        "frame_records":frames,
        "counterexample_witness_count":len(counter_witnesses),"counterexample_stationary":True,
        "counterexample_third_level":is_third_level(D),"counterexample_clifford":is_clifford(D),
        "counterexample_permutation":list(counter_permutation),"counterexample_phases_mu4":list(counter_phases),
        "counterexample_alpha":list(alpha),
        "counterexample_product_flag_indices":[j for j in range(21) if product_mask>>j&1],
        "counterexample_common_flag_indices":[j for j in range(21) if common_mask>>j&1],
        "distinct_permutations_screened_against_all_flags":len(masks),
        "seconds":time.perf_counter()-started,
    }
    validate_summary(result)
    return result


if __name__ == "__main__":
    print(json.dumps(run(),sort_keys=True))
