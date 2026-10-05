# (C13) Invariant sets of input Lagrangians

Fact (C13) of the paper, label `cert:invariant-sets-of-input-lagrangians`, stated in a box at its first use (used in `app:return-multiple`). Facts it depends on: (C11), (C12).

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $A$ be a non-third-level normal form of \textup{(\ref{cert:block-form-clifford-classes})}, and let $D$ be a Clifford conjugate of one of the four gates $\mathrm{CS}_{yz}^a\CCZ^b$. If $X(D)$ is invariant under $J_A$, then $L_Z\in X(D)$.
```

## Check

```sh
python check.py math --fact invariant-sets-of-input-lagrangians
```

- Finite domain: The 1,396 groups E_D of the Clifford conjugates D of CS_yz^a CCZ^b, with 451 distinct sets X(D); the 20 non-third-level cores.
- Checked: The 1,396 groups are exactly the Sp(6,2)-conjugates of the groups of the four gates; for each core, J_A (of order 8 or 64) and the 620 pairs with X(D) invariant under J_A, each with L_Z in X(D).
- Runtime: 18 s.
- Implementation: `src/checking/conventional.py, src/checking/conventional_audit.py`.
- Release run: accepted in 18 s; record in [results/invariant-sets-of-input-lagrangians.json](../../results/invariant-sets-of-input-lagrangians.json).
