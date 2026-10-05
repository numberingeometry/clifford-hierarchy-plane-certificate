# (C17) Sylow cosets in $S_6$ and three-element subsets

Fact (C17) of the paper, label `cert:s6-sylow-cosets-and-three-element-subsets`, stated in a box at its first use (used in `lem:fib-sharp-small`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $F\subseteq\Sp(4,2)\cong S_6$. If every subset of $F$ with at most three elements lies in a right coset $Pg$ of a Sylow $2$-subgroup $P$, then so does $F$.
```

## Check

```sh
python check.py math --fact s6-sylow-cosets-and-three-element-subsets
```

- Finite domain: The graph on the 255 nonidentity elements of 2-power order; its 45 Sylow and 180 other maximal cliques; the 17,820 subsets with at least three elements of the 180 cliques.
- Checked: Each of the 12,240 subsets in no Sylow subgroup (10,320 distinct) contains a triple a, b, c with ba^{-1} and ca^{-1} in no common Sylow subgroup; 1,440 such triples occur.
- Runtime: under 1 s.
- Implementation: `src/checking/c17.py`.
- Release run: accepted in under 1 s; record in [results/s6-sylow-cosets-and-three-element-subsets.json](../../results/s6-sylow-cosets-and-three-element-subsets.json).
