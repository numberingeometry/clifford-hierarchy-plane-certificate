# (C12) Clifford classes of the block forms

Fact (C12) of the paper, label `cert:block-form-clifford-classes`, stated in a box at its first use (used in `app:return-multiple`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Under $D\mapsto C_1DC_2$ with $C_1,C_2\in\C_2(3)$, the $512$ block forms $D$ of Lemma~\ref{app:return-block}\textup{(4)} fall into $30$ classes. Ten classes are third-level, and every correction $DP_vD^\dagger$ of the other twenty is third-level. We call the listed representative of a class its normal form. The third-level classes containing a form with $p'=0$ are those of the four gates $\mathrm{CS}_{yz}^a\CCZ^b$.
```

## Check

```sh
python check.py math --fact block-form-clifford-classes
```

- Finite domain: The 512 block forms pi_0 i^(axy+bxz+cyz+dxyz) T_p' and the Clifford moves used in the paper.
- Checked: 30 orbits, each containing exactly one listed core; 10 third-level cores; all 1,280 corrections and right derivatives of the other 20 are third-level; 52,320 normalization steps checked exactly; the third-level classes meeting p' = 0 are exactly the four classes of CS_yz^a CCZ^b.
- Runtime: 47 s.
- Implementation: `src/checking/conventional.py, src/checking/conventional_audit.py`.
- Release run: accepted in 47 s; record in [results/block-form-clifford-classes.json](../../results/block-form-clifford-classes.json).
