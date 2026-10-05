# (C6) Representatives of $\mathscr R_4$ up to Clifford conjugacy

Fact (C6) of the paper, label `cert:r4-representatives-up-to-clifford-conjugacy`, stated in a box at its first use (used in `app:finite-endpoint`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
A list of $1{,}020$ gates of the form~\eqref{core:eq:row} on four qubits contains a Clifford conjugate of every $h\in\mathscr R_4$. Exactly ten listed $M$ have $|N(M)|-1\ge63$, and they lie in $\C_2(4)$.
```

## Check

```sh
python check.py math --fact r4-representatives-up-to-clifford-conjugacy
```

- Finite domain: 32,768 triangular seeds, 2,064 involutions, 163,216 affine-closed point maps, 43 phase forests and 2,818,048 phase rows; the 1,020 listed gates.
- Checked: Complete coverage of the forms (core:eq:row) up to Clifford conjugacy by the 1,020 listed gates; N(M) for each of them; the ten with |N(M)|-1 >= 63 lie in C_2(4).
- Runtime: 5 s.
- Implementation: `four-qubit/src/check_coverage.py`.
- Release run: accepted in 5 s; record in [results/r4-representatives-up-to-clifford-conjugacy.json](../../results/r4-representatives-up-to-clifford-conjugacy.json).
