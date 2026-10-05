# (C10) The case $\mathrm{rank}\,E=2$

Fact (C10) of the paper, label `cert:balanced-rank-two-branch`, stated in a box at its first use (used in `app:return-balanced`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
In Case~2b of the proof of Lemma~\ref{app:return-balanced}: for $\ell=z$ and the correction $h=WX_sW^\dagger$, the monomial gate $Z_1h$, in the basis $\ket\pm\otimes\ket x$, has a permutation of cycle type $2^26^2$, so this case does not occur. For $\ell=s$, the Pauli-MASA condition on the off-diagonal block makes every $\operatorname{diag}\theta_t$ Clifford, and among the $64$ choices of $a,b,c,d,e,f\in\F_2$ the three unipotence conditions hold exactly when $a=b=c=e=0$ and $d=f$.
```

## Check

```sh
python check.py math --fact balanced-rank-two-branch
```

- Finite domain: The permutation of Z_1 h on 16 points; the 512 two-target phase tables; the 64 parameter choices a, b, c, d, e, f.
- Checked: Z_1 h has cycle type 2^2 6^2; every stationary phase table is Clifford; the unipotence conditions hold exactly when a = b = c = e = 0 and d = f, and then a common witness exists.
- Runtime: 2 s.
- Implementation: `src/checking/c7.py`.
- Release run: accepted in 2 s; record in [results/balanced-rank-two-branch.json](../../results/balanced-rank-two-branch.json).
