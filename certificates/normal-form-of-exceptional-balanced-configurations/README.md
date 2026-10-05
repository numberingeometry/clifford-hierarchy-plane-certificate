# (C9) Normal form of the triples $(\rho,u,D)$

Fact (C9) of the paper, label `cert:normal-form-of-exceptional-balanced-configurations`, stated in a box at its first use (used in `app:return-balanced`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $\tau_0(x,y,z)=(x,y,z+xy)$, and let $\rho$ be a permutation of $\F_2^3$ that is affine or affinely equivalent to $\tau_0$. Let $u\ne0$ and let $D$ be an affine plane with $D+u=D$ such that $\tau\rho$, where $\tau(x)=x+u\,[x\in D]$, is again affine or affinely equivalent to $\tau_0$. If fewer than two independent $c$ with $c\cdot u=0$ make $c\cdot\rho$ affine, then after replacing $(\rho,u,D)$ by $(\alpha\rho\beta,\bar\alpha u,\alpha(D))$ for affine $\alpha,\beta$, with $\bar\alpha$ the linear part of $\alpha$, in coordinates $(s,y,z)$ we have $\rho(s,y,z)=(s,y,z+sy)$, $u=e_y$, and $D=\{s=1\}$ or $D=\{z=1\}$.
```

## Check

```sh
python check.py math --fact normal-form-of-exceptional-balanced-configurations
```

- Finite domain: The 1,344 affine maps; the permutations rho in {1, tau_0}; all u != 0 and all affine planes D with D+u = D.
- Checked: Exactly 36 exceptional pairs (u,D) for rho = tau_0, forming two orbits, of sizes 12 and 24, under the 192 affine alpha with tau_0^{-1} alpha tau_0 affine.
- Runtime: 2 s.
- Implementation: `src/checking/c15.py`.
- Release run: accepted in 2 s; record in [results/normal-form-of-exceptional-balanced-configurations.json](../../results/normal-form-of-exceptional-balanced-configurations.json).
