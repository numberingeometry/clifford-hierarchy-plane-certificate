from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/checking'))
"""Completeness controls and fail-closed tests for the C17 recomputation."""
from copy import deepcopy
from itertools import combinations
import unittest
from unittest.mock import patch

try:
    from checking import c17
except ImportError:
    import c17


def brute_maximal(adjacency):
    n = len(adjacency)
    cliques = []
    for subset in range(1 << n):
        vertices = [i for i in range(n) if subset >> i & 1]
        if not all(adjacency[a] >> b & 1 for a, b in combinations(vertices, 2)):
            continue
        if any(all(adjacency[other] >> vertex & 1 for vertex in vertices)
               for other in range(n) if not subset >> other & 1):
            continue
        cliques.append(subset)
    return sorted(cliques)


class C17Controls(unittest.TestCase):
    def test_complete_bron_kerbosch_on_every_graph_through_five_vertices(self):
        # An independent exhaustive subset algorithm tests clique completeness,
        # including isolated vertices, disconnected graphs and the empty graph.
        for n in range(6):
            edges = list(combinations(range(n), 2))
            for graph in range(1 << len(edges)):
                adjacency = [0] * n
                for j, (a, b) in enumerate(edges):
                    if graph >> j & 1:
                        adjacency[a] |= 1 << b
                        adjacency[b] |= 1 << a
                actual, _ = c17.maximal_cliques(adjacency)
                self.assertEqual(actual, brute_maximal(adjacency), (n, graph))

    def test_malformed_graph_rejected(self):
        for adjacency in ([1], [2], [2, 0]):
            with self.subTest(adjacency=adjacency), self.assertRaises(ValueError):
                c17.maximal_cliques(adjacency)

    def test_omitted_computed_clique_rejected(self):
        original = c17.maximal_cliques

        def omitted(adjacency):
            cliques, nodes = original(adjacency)
            return cliques[:-1], nodes

        with patch.object(c17, "maximal_cliques", omitted):
            with self.assertRaisesRegex(ValueError, "maximal-clique histogram"):
                c17.run()

    def test_changed_counts_or_unresolved_candidate_rejected(self):
        record = c17.run()
        changes = {
            "permutations": 719,
            "sylows": 44,
            "two_elements": 255,
            "maximal_cliques": 224,
            "deficient_subsets_checked": 12239,
            "distinct_bad_triples": 1439,
            "unresolved_deficient_subsets": 1,
        }
        for key, wrong in changes.items():
            with self.subTest(key=key):
                changed = deepcopy(record)
                changed[key] = wrong
                with self.assertRaises(ValueError):
                    c17.validate_summary(changed)

    def test_false_or_out_of_domain_bad_triple_rejected(self):
        # All these quotient elements share the same Sylow; the proposed
        # triple is therefore not an obstruction to right-coset containment.
        quotient_masks = [[1] * 3 for _ in range(3)]
        with self.assertRaisesRegex(ValueError, "lies in a Sylow right coset"):
            c17.validate_bad_triple((0, 1, 2), quotient_masks, (0, 1, 2))
        with self.assertRaisesRegex(ValueError, "not in its candidate set"):
            c17.validate_bad_triple((0, 1, 2), quotient_masks, (0, 1))
        with self.assertRaisesRegex(ValueError, "not in its candidate set"):
            c17.validate_bad_triple((0, 1, 1), quotient_masks, (0, 1))


if __name__ == "__main__":
    unittest.main()
