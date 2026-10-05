# (C14) Reduction of $B_0$ to $256$ gates

Fact (C14) of the paper, label `cert:reduction-to-common-input-pair-domain`, stated in a box at its first use (used in `app:return-multiple`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $A$ be a non-third-level normal form with $T$-factor $T_{p_A}$, and $B_0=\Pi_\rho\,i^{g}T_{p_A}$ with $g:\F_2^3\to\Z/4$ and $\rho$ affine or affinely equivalent to $\pi$ of Lemma~\ref{app:return-block}\textup{(4)}. Let $\rho_0$ be the identity, and for $0\ne k\in\F_2^3$ let $\rho_k=\lambda_k^{-1}\pi\lambda_k$ with $\lambda_k$ linear and $\lambda_k(k)=(0,0,1)$. Modulo left multiplication by monomial Clifford gates, $B_0$ is one of the $256$ gates $\Pi_{\rho_k}\,i^{axy+bxz+cyz+dxyz}\,T_{p_A}$ with $k\in\F_2^3$, $a,b,c\in\F_2$ and $d\in\Z/4$, or $\rho(x+k)+\rho(x)$ does not depend on $x$ exactly for $k\in\{0,(1,1,1)\}$, and then some $B_0P_vA^\dagger\notin\C_3(3)$.
```

## Check

```sh
python check.py math --fact reduction-to-common-input-pair-domain
```

- Finite domain: The 10,752 permutations affine or affinely equivalent to pi; all g: F_2^3 -> Z/4 modulo diagonal Clifford phases; the 20 non-third-level cores.
- Checked: Every permutation lies in a left affine coset of one of eight listed ones; the 32 listed phases represent all classes for the identity and the six other radical directions and 16 for (1,1,1), where a non-third-level B_0 P_v A^dagger is found; the 360 common-input pairs are recomputed.
- Runtime: 11 s.
- Implementation: `src/checking/exact_audit/`.
- Release run: accepted in 11 s; record in [results/reduction-to-common-input-pair-domain.json](../../results/reduction-to-common-input-pair-domain.json).
