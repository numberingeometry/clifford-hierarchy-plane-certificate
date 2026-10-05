"""Fact (C16): a 2-group criterion in S_6. Exact recomputation.

Permutation products act rightmost first. A generated subgroup is a 2-group exactly when its generators
lie in a common Sylow 2-subgroup; the 45 Sylow 2-subgroups are enumerated. Internal key: C8.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import permutations
import json
from pathlib import Path
import time


def need(condition, message):
    if not condition:
        raise ValueError(message)


def run():
    started = time.perf_counter()
    perms = list(permutations(range(6)))
    ident = tuple(range(6))

    def mul(a, b):
        return tuple(a[b[j]] for j in range(6))

    def inv(a):
        return tuple(a.index(j) for j in range(6))

    edge = frozenset((4, 5))
    matching = {frozenset((0, 1)), frozenset((2, 3))}
    base = frozenset(g for g in perms if frozenset(g[j] for j in edge) == edge
                     and {frozenset(g[j] for j in e) for e in matching} == matching)
    need(len(perms) == 720 and len(set(perms)) == 720, "complete S6 domain")
    need(len(base) == 16 and ident in base, "order-sixteen base Sylow")
    need(all(mul(a, b) in base for a in base for b in base), "base subgroup closure")
    # The largest power of two dividing 6! is sixteen. Sylow conjugacy covers
    # every Sylow two-subgroup with the complete 720-element conjugating domain.
    sylows = set()
    for g in perms:
        gi = inv(g)
        sylows.add(frozenset(mul(mul(g, h), gi) for h in base))
    sylows = sorted(sylows, key=lambda h: sorted(h))
    need(len(sylows) == 45 and all(len(h) == 16 for h in sylows), "complete Sylow domain")
    masks = {g: sum(1 << j for j, subgroup in enumerate(sylows) if g in subgroup)
             for g in perms}
    two = [g for g in perms if masks[g]]
    need(len(two) == 256, "two-element domain size")
    for g in perms:
        g2 = mul(g, g)
        need(bool(masks[g]) == (mul(g2, g2) == ident), "Sylow cover and order-four test disagree")
    membership = [masks[g] for g in two]
    inverses = [inv(g) for g in two]
    quotient = [[masks[mul(c, ai)] for ai in inverses] for c in two]
    checked = passed = bad = 0
    first = None
    for ai, ma in enumerate(membership):
        for bi, mb in enumerate(membership):
            qba = quotient[bi][ai]
            for ci, mc in enumerate(membership):
                checked += 1
                if not (ma & quotient[ci][bi] and mb & quotient[ci][ai] and mc & qba):
                    continue
                passed += 1
                if not (ma & mb & mc):
                    bad += 1
                    if first is None:
                        first = [two[ai], two[bi], two[ci]]
    need(checked == 256 ** 3 == 16777216, "complete triple census")
    need(passed == 153616, "directional hypothesis count changed")
    need(bad == 0 and first is None, "two-group conclusion failed: " + str(first))
    return {"status": "PASS", "route": "recompute", "exact": True,
            "fact": "s6-two-group-criterion", "permutations": len(perms),
            "base_sylow_order": len(base), "sylows": len(sylows), "two_elements": len(two),
            "triples_checked": checked, "directional_hypothesis_pass": passed,
            "counterexamples": bad, "first_counterexample": first,
            "seconds": time.perf_counter() - started,
            "conditions": "<a,c b^-1>, <b,c a^-1>, <c,b a^-1> are two-groups",
            "conclusion": "<a,b,c> is a two-group",
            "arithmetic": "exact permutations and integer Sylow-membership masks",
            "independent_kernel": False,
            "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a new execution record")
    args = parser.parse_args()
    record = run()
    encoded = json.dumps(record, indent=2) + "\n"
    if args.output is not None:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
