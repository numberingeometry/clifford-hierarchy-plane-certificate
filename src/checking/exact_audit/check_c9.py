"""Fact (C8): Pauli terms and the cubic coefficient, checked exactly.

(a) d odd: for every non-diagonal Pauli P, D^dagger P D has at least eight Pauli terms (the minimum is
    exactly eight).
(b) d even: for every p, {v : D T_p P_v T_p^dagger D^dagger in C_2(3)} is a subspace of dimension at
    least five.
D = pi_0 i^(a xy + b xz + c yz + d xyz), pi_0 in {I, pi}, a, b, c in F_2, d in Z/4. Key: C9.
"""
import itertools, time
from .mono3 import *


if not __debug__:
    raise RuntimeError("Run without -O: exact audit assertions must be enabled.")

def recompute():
    t0 = time.perf_counter()
    out = {"fact": "C9", "manuscript": "327e9e7"}

    # (a)
    counts = []
    gates_a = 0
    for t, a, b, c, d in itertools.product((0, 1), (0, 1), (0, 1), (0, 1), (1, 3)):
        D = mono(t, a, b, c, d, 0)
        gates_a += 1
        for v in range(64):
            if v & 7 == 0:
                continue  # diagonal
            counts.append(len(pauli_support(mul(mul(inv(D), PAULIS[v]), D))))
    out["a"] = {"gates": gates_a, "paulis_per_gate": 56, "cases": len(counts),
                "minimum_terms": min(counts), "histogram": {str(k): counts.count(k) for k in sorted(set(counts))},
                "holds": min(counts) >= 8, "minimum_exactly_eight": min(counts) == 8}

    # (b)
    forms = 0
    dims = []
    nonsub = []
    for t, a, b, c, d, p in itertools.product((0, 1), (0, 1), (0, 1), (0, 1), (0, 2), range(8)):
        D = mono(t, a, b, c, d, 0)
        T = mono(0, 0, 0, 0, 0, p)
        DT = mul(D, T)
        forms += 1
        good = [v for v in range(64) if is_clifford(conj(DT, PAULIS[v]))]
        if not is_subspace(good):
            nonsub.append([t, a, b, c, d, p])
        dims.append(span_dim(good))
    out["b"] = {"forms": forms, "labels_per_form": 64, "non_subspace_forms": nonsub,
                "minimum_dimension": min(dims), "dimension_histogram": {str(k): dims.count(k) for k in sorted(set(dims))},
                "holds": not nonsub and min(dims) >= 5}
    out["holds"] = out["a"]["holds"] and out["b"]["holds"]
    out["seconds"] = round(time.perf_counter() - t0, 2)
    return out
