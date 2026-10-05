"""Fact (C9): normal form of the triples (rho, u, D). Exact census.

Basis order is |xyz>, x the most significant bit; permutations compose right to left. Every affine
permutation and every A tau_0 B is generated, and every admissible invariant affine plane is checked for
all 10,752 permutations. Affine normalization brings every exceptional configuration to one of the two
stated representatives for rho = tau_0. Internal key: C15.
"""

from itertools import permutations
import json
import time


class CheckFailure(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise CheckFailure(message)


IDENTITY = tuple(range(8))
TAU0 = tuple(j ^ ((j >> 2 & 1)*(j >> 1 & 1)) for j in range(8))


def compose(p, q):
    return tuple(p[q[j]] for j in range(8))


def inverse(p):
    out = [0]*8
    for j, value in enumerate(p):
        out[value] = j
    return tuple(out)


def affine_values(values):
    return all(values[j] == values[0] ^ (values[4]^values[0] if j & 4 else 0)
               ^ (values[2]^values[0] if j & 2 else 0) ^ (values[1]^values[0] if j & 1 else 0)
               for j in range(8))


def affine_domain():
    out = [p for p in permutations(range(8)) if affine_values(p)]
    require(len(out) == 1344, "affine permutation enumeration is incomplete")
    return out


def dot(c, value):
    return (c & value).bit_count() & 1


def rank(vectors):
    basis = {}
    for value in vectors:
        while value:
            pivot = value.bit_length()-1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def invariant_planes():
    planes = {frozenset(j for j in range(8) if dot(normal,j) == constant)
              for normal in range(1,8) for constant in (0,1)}
    require(len(planes) == 14, "affine-plane enumeration is incomplete")
    configs = [(u,D) for u in range(1,8) for D in sorted(planes,key=lambda s:tuple(sorted(s)))
               if frozenset(j^u for j in D) == D]
    require(len(configs) == 42, "invariant affine-plane domain is incomplete")
    return configs


def characters(rho):
    return tuple(c for c in range(1,8) if affine_values(tuple(dot(c,rho[j]) for j in range(8))))


def exceptional(rho_characters, u):
    return rank(c for c in rho_characters if dot(c,u) == 0) < 2


def action(alpha, config):
    u, D = config
    return alpha[u]^alpha[0], frozenset(alpha[j] for j in D)


def encode_config(config):
    u,D = config
    return {"u": u, "D": sorted(D)}


def validate_summary(summary):
    expected = {"affine_permutations":1344, "toffoli_double_coset_size":9408,
                "complete_permutation_domain":10752, "configurations_per_permutation":42,
                "full_configuration_domain":451584, "identity_exceptional_configurations":0,
                "tau0_exceptional_configurations":36, "stabilizer_size":192,
                "full_exceptional_configurations":338688}
    for key,value in expected.items():
        require(summary.get(key) == value, f"C15 {key}: expected {value}, received {summary.get(key)!r}")
    require(summary.get("representative_orbit_sizes") == {"D_x_eq_1_u_y":12,"D_z_eq_1_u_y":24}, "C15 representative orbit sizes changed")
    require(summary.get("two_representative_orbits_disjoint") is True, "C15 orbits overlap")
    require(summary.get("two_representative_orbits_equal_exceptional_set") is True, "C15 orbit coverage failed")
    require(summary.get("every_full_domain_exception_normalizes") is True, "C15 full-domain normalization failed")


def run():
    started = time.perf_counter()
    affines = affine_domain()
    affset = set(affines)
    inverse_affines = {a:inverse(a) for a in affines}
    # Each supplied left normalization is produced directly from the complete
    # double-coset generation, rather than obtained by an unbounded search.
    double_coset = {}
    for left in affines:
        left_tau = compose(left,TAU0)
        normalize_left = inverse_affines[left]
        for right in affines:
            rho = compose(left_tau,right)
            double_coset.setdefault(rho,normalize_left)
    require(len(double_coset) == 9408, "Toffoli double-coset size changed")
    require(not (affset & double_coset.keys()), "affine and Toffoli domains overlap")
    domain = affset | double_coset.keys()
    require(len(domain) == 10752, "complete affine/Toffoli domain is incomplete")
    configs = invariant_planes()
    for u,D in configs:
        tau = tuple(j^u if j in D else j for j in range(8))
        require(tau in affset, "an invariant-plane shear is not affine")
    identity_exceptions = {config for config in configs if exceptional(characters(IDENTITY),config[0])}
    reduced_exceptions = {config for config in configs if exceptional(characters(TAU0),config[0])}
    stabilizer = [a for a in affines if compose(TAU0,compose(a,TAU0)) in affset]
    rep_x = 2,frozenset(j for j in range(8) if j & 4)
    rep_z = 2,frozenset(j for j in range(8) if j & 1)
    orbit_x = {action(a,rep_x) for a in stabilizer}
    orbit_z = {action(a,rep_z) for a in stabilizer}
    require(rep_x in reduced_exceptions and rep_z in reduced_exceptions, "advertised representatives are not exceptional")
    require(orbit_x <= reduced_exceptions and orbit_z <= reduced_exceptions, "stabilizer does not preserve the exceptional set")
    full_checked = 0
    full_exceptions = 0
    for rho in sorted(domain):
        rho_characters = characters(rho)
        if rho in affset:
            require(len(rho_characters) == 7, "an affine permutation lost an affine character")
            normalizer = None
        else:
            normalizer = double_coset[rho]
            require(compose(TAU0,compose(normalizer,rho)) in affset, "stored affine normalization is invalid")
        for config in configs:
            full_checked += 1
            if exceptional(rho_characters,config[0]):
                full_exceptions += 1
                require(normalizer is not None, "affine rho has an exceptional configuration")
                normalized = action(normalizer,config)
                require(normalized in reduced_exceptions, "full-domain exceptional configuration does not normalize")
    result = {
        "fact":"normal-form-of-exceptional-balanced-configurations",
        "route":"recompute",
        "status":"accepted",
        "affine_permutations":len(affines),
        "toffoli_double_coset_size":len(double_coset),
        "complete_permutation_domain":len(domain),
        "configurations_per_permutation":len(configs),
        "full_configuration_domain":full_checked,
        "identity_exceptional_configurations":len(identity_exceptions),
        "tau0_exceptional_configurations":len(reduced_exceptions),
        "stabilizer_size":len(stabilizer),
        "full_exceptional_configurations":full_exceptions,
        "representative_orbit_sizes":{"D_x_eq_1_u_y":len(orbit_x),"D_z_eq_1_u_y":len(orbit_z)},
        "two_representative_orbits_disjoint":not (orbit_x & orbit_z),
        "two_representative_orbits_equal_exceptional_set":orbit_x | orbit_z == reduced_exceptions,
        "every_full_domain_exception_normalizes":True,
        "reduced_exceptional_configurations":[encode_config(config) for config in sorted(reduced_exceptions,key=lambda x:(x[0],tuple(sorted(x[1]))))],
        "seconds":time.perf_counter()-started,
    }
    validate_summary(result)
    return result


if __name__ == "__main__":
    print(json.dumps(run(),sort_keys=True))
