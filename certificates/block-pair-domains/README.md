# (C15) Pair domains for blocks

Fact (C15) of the paper, label `cert:block-pair-domains`, stated in a box at its first use (used in `app:return-multiple`). Facts it depends on: (C12), (C14).

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
(a) Let $A$ be a non-third-level normal form of \textup{(\ref{cert:block-form-clifford-classes})} and $B_0$ one of the $256$ gates of \textup{(\ref{cert:reduction-to-common-input-pair-domain})} for this $A$ ($5{,}120$ pairs). For the $360$ pairs with every $D_v=B_0P_vA^\dagger\in\C_3(3)$, each $D_v=G_vY_v$ with $G_v\in\C_2(3)$ and $Y_v$ a Clifford conjugate of some $\mathrm{CS}_{yz}^a\CCZ^b$, and $\mathcal B(A,FG_0^{-1}B_0)$ holds for every $F\in\C_2(3)$ with $\gamma_F\gamma_{G_0^{-1}G_v}e$ unipotent for all $v\in V_3$ and $e\in E_{Y_v}$.

(b) Let $A$ be a third-level normal form and $D$ a Clifford conjugate of some $\mathrm{CS}_{yz}^a\CCZ^b$ with $X(D)$ invariant under $E_A$ ($13{,}960$ pairs, reduced to $2{,}076$). Then $\mathcal B(A,FDA)$ holds for every $F\in\C_2(3)$ with $\gamma_Fef$ unipotent for all $e\in E_D$ and $f\in E_A$.
```

## Check

```sh
python check.py math --fact block-pair-domains
```

- Finite domain: (a) 5,120 pairs, of which 360 have every D_v third-level; (b) 13,960 pairs, reduced to 2,076; all 1,451,520 elements of Sp(6,2).
- Checked: The exact C kernels are compiled, the 360 and 2,076 records regenerated, and all 572,416 and 18,331,136 qualifying combinations satisfy B(A,B).
- Runtime: 2 min.
- Implementation: `src/checking/conventional_large.py`.
- Release run: accepted in 2 min; record in [results/block-pair-domains.json](../../results/block-pair-domains.json).
