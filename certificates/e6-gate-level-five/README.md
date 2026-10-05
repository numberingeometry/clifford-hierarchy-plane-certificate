# (C18) The gate $U_{E_6}$ lies in $\mathcal{C}_5(5)$

Fact (C18) of the paper, label `cert:e6-gate-level-five`, stated in a box at its first use (used in `thm:exact-cells`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
The gate $U_{E_6}$ of~\eqref{eq:fib-e6-family} lies in $\C_5(5)$.
```

## Check

```sh
python check.py math --fact e6-gate-level-five
```

- Finite domain: The gate U_{E6} as a 32 x 32 matrix over Z[i,1/2]; all 1,024 Paulis at levels 4 and 5, ten generators at levels 2 and 3.
- Checked: U_{E6} lies in C_5(5) (and not in C_4(5)), decided exactly by the hierarchy recursion; every Clifford block it uses is audited against direct matrix arithmetic.
- Runtime: 31 s.
- Implementation: `src/checking/e6_level.py`, `src/checking/e2_hierarchy.py`.
- Release run: accepted in 31 s; record in [results/e6-gate-level-five.json](../../results/e6-gate-level-five.json).
