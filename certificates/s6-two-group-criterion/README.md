# (C16) A $2$-group criterion in $S_6$

Fact (C16) of the paper, label `cert:s6-two-group-criterion`, stated in a box at its first use (used in `lem:fib-sharp-small`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $a,b,c\in\Sp(4,2)\cong S_6$. If $\langle a,cb^{-1}\rangle$, $\langle b,ca^{-1}\rangle$ and $\langle c,ba^{-1}\rangle$ are $2$-groups, then so is $\langle a,b,c\rangle$.
```

## Check

```sh
python check.py math --fact s6-two-group-criterion
```

- Finite domain: The 256^3 ordered triples of elements of 2-power order in Sp(4,2), isomorphic to S_6; its 45 Sylow 2-subgroups.
- Checked: 153,616 triples satisfy the hypotheses, and all of them satisfy the conclusion.
- Runtime: 2 s.
- Implementation: `src/checking/c8.py`.
- Release run: accepted in 2 s; record in [results/s6-two-group-criterion.json](../../results/s6-two-group-criterion.json).
