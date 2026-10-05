"""Exact recomputation of (C19): the order of A and the characteristic polynomial of A^2 B.

A and B are the symplectic actions of Appendix G. The binary roots p_0, ..., p_5 have as polar form the
graph form of the Dynkin diagram of E_6 (edges 01, 12, 23, 34, 25), and A = t_{p_1} ... t_{p_5},
B = t_{p_0} t_{p_2} with t_p(v) = v + omega(v, p) p. The polar form has rank 6, so the roots are a basis of
V_3, and the order and the characteristic polynomial computed in the root basis are those of A and A^2 B.
Products compose right to left. The binary root p_i is bit i. Key: C18 = (C19).
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
from . import binary as f

KEYS = ("C18",)
E6_EDGES = ((0, 1), (1, 2), (2, 3), (3, 4), (2, 5))


def e6_case():
    n = 6
    polar = f.graph_polar(n, E6_EDGES)
    f.need(f.rank(polar) == n, "E6 polar form is singular")
    roots = tuple(1 << j for j in range(n))
    T = [f.transvection(v, polar) for v in roots]
    for t in T:
        f.need(f.multiply(t, t) == f.identity(n), "E6 root transvection is not involutory")
    A, B = f.product(*T[1:6]), f.product(T[0], T[2])
    order_A = f.order(A, 8)
    f.need(order_A == 8, "E6 A order changed")
    polynomial = f.characteristic_polynomial(f.product(A, A, B))
    f.need(polynomial == (1 << 6) | (1 << 3) | 1, "E6 characteristic polynomial changed")
    irreducible = f.polynomial_irreducible(polynomial)
    f.need(irreducible, "E6 characteristic polynomial is reducible")
    return {"name": "E6", "dimension": n, "edges": [list(e) for e in E6_EDGES], "polar_rank": n,
            "order_A": order_A, "characteristic_polynomial_bits": polynomial,
            "characteristic_polynomial": "z^6+z^3+1", "polynomial_irreducible": irreducible}


def recompute():
    return {"case": e6_case(), "arithmetic": "exact F2 packed integer matrices and polynomials"}


def validate(result):
    case = result['case']
    f.need(case['name'] == 'E6' and case['dimension'] == case['polar_rank'] == 6 and
           case['edges'] == [list(e) for e in E6_EDGES], "C19 root domain missing or changed")
    f.need(case['order_A'] == 8, "C19 order of A failed")
    f.need(case['characteristic_polynomial_bits'] == 73 and case['characteristic_polynomial'] == 'z^6+z^3+1'
           and case['polynomial_irreducible'] is True, "C19 characteristic polynomial failed")


def run(fact, output_dir=None):
    f.need(fact in KEYS, "Unsupported root-word fact: " + str(fact))
    started = time.perf_counter()
    result = recompute()
    validate(result)
    owner = Path(__file__).parent
    report = {'status': 'accepted', 'route': 'recompute', 'fact': fact, 'navigation': fact,
              'exact': True, 'result': result, 'elapsed_seconds': time.perf_counter() - started,
              'source_sha256': {name: hashlib.sha256((owner / name).read_bytes()).hexdigest()
                                for name in ('runner.py', 'binary.py')}}
    if output_dir is not None:
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        with (output / 'result.json').open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fact', choices=KEYS)
    parser.add_argument('--output-dir')
    args = parser.parse_args()
    print(json.dumps(run(args.fact, args.output_dir), indent=2))


if __name__ == '__main__':
    main()
