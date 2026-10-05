# (C7) A plane meeting the fixed Pauli labels

Fact (C7) of the paper, label `cert:plane-meeting-fixed-pauli-labels`, stated in a box at its first use (used in `app:finite-endpoint`). Facts it depends on: (C6).

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
For each of the other $1{,}010$ listed $M$ there is a plane $\Pi_M\le V_4$ with
\[
N(h)\cap\Pi_M\ne\{0\}\qquad\text{for all }h\in\{M\}\cup\mathcal K(M),
\]
where $N(h)=\{s:hP_sh^\dagger=\pm P_s\}$.
```

## Check

```sh
python check.py math --fact plane-meeting-fixed-pauli-labels --allow-expensive --workers 8
```

- Finite domain: The 1,010 listed gates with |N(M)|-1 < 63; the 2,295 Lagrangians of V_4; the 163,216 eligible point maps; both scalar sectors and every solution of the phase conditions.
- Checked: The exact C kernel is compiled, passes its regression tests, and is run afresh on all 1,010 gates; for each gate the computed plane, level and sector must be those listed in data/lock-n4.txt.
- Runtime: 8.2 h with one worker; `--workers N` distributes the gates over N processes.
- Implementation: `src/checking/conventional_large.py`.
- Release run: accepted in 84 min with 6 workers; record in [results/plane-meeting-fixed-pauli-labels.json](../../results/plane-meeting-fixed-pauli-labels.json); census summary in [results/plane-meeting-census/](../../results/plane-meeting-census/).
