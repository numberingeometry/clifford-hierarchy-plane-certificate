# (C19) Root words for $E_6$

Fact (C19) of the paper, label `cert:e6-root-words`, stated in a box at its first use (used in `thm:exact-cells`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $A$ and $B$ be as in Appendix~\ref{subsec:orthogonal-fibers}. Then $A$ has order $8$ and $\chi_{A^2B}(z)=z^6+z^3+1$.
```

## Check

```sh
python check.py math --fact e6-root-words
```

- Finite domain: The binary root matrices A and B of the E_6 example on the six-dimensional root space with the graph form of the Dynkin diagram.
- Checked: The graph form has rank 6; A has order 8; A^2 B has characteristic polynomial z^6+z^3+1, which is irreducible over F_2.
- Runtime: under 1 s.
- Implementation: `src/checking/root_words/`.
- Release run: accepted in under 1 s; record in [results/e6-root-words.json](../../results/e6-root-words.json).
