# (C11) Mixed blocks with two input Lagrangians

Fact (C11) of the paper, label `cert:mixed-blocks-with-two-input-lagrangians`, stated in a box at its first use (used in `app:return-multiple`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $D=\pi_0\,i^{axy+bxz+cyz+dxyz}$ with $\pi_0\in\{I,\pi\}$, $a,b,c\in\F_2$ and $d\in\Z/4$, a block form of Lemma~\ref{app:return-block}\textup{(4)} with $p'=0$. If $|X(D)|\ge2$, then $D\in\C_3(3)$.
```

## Check

```sh
python check.py math --fact mixed-blocks-with-two-input-lagrangians
```

- Finite domain: The 64 forms pi_0 i^(axy+bxz+cyz+dxyz) and the 135 Lagrangians of V_3.
- Checked: The input Lagrangians of each form: 20 forms have at least two (3, 15 or 135), and all 20 are third-level.
- Runtime: 12 s.
- Implementation: `src/checking/conventional.py, src/checking/conventional_audit.py`.
- Release run: accepted in 12 s; record in [results/mixed-blocks-with-two-input-lagrangians.json](../../results/mixed-blocks-with-two-input-lagrangians.json).
