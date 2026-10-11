# Theorist · Round 036 conference statement

## Single most promising question

**Why does the dimension d_i = e_iAe_v collapse to exactly 2 whenever J_i is nonzero on walks?**

This is quantitatively precise (not just d >= 2, but d = 2 exactly) and empirically overwhelming: all 285 rows with J_i != 0 at n = 8–9 capped walks have d_i = 2; all 28 rows with d_i >= 3 have J_i = 0 (E-131). No mechanism is known. If d = 2 is forced, the answer explains both why J can exist (needs d >= 2) and why it saturates (J_i = 0 at d >= 3 automatically). This single bound would close E-126's open question and anchor the silting-not-tilting criterion to the path structure.

## Weakest claim the workshop relies on

**"No base rate for LNA keys among circuit members; the 0-of-2704 in E-125 proves nothing."**

Why weak: E-125 tested one constrained family (layered, circuit-free, pure-W by construction). A null result in one thin family is not evidence that *all* circuit members lack keys, only that this family does. The null test lacks a positive control (e.g., a random gate-admitted sample of the same size and dimension range). Rounds 032 and 033 correctly flagged this, but the lack of a base rate means we cannot yet say whether "no LNA keys on walks" is a true constraint or an artifact of the walk-parent construction itself.

## What I need from another persona

- **From experimentalist:** d_i histogram on n = 9, 10 long walks (not capped) to test whether d = 2 saturation holds beyond the prefix. Is the pattern fragile or universal?
- **From skeptic:** a hand-built gate-admitted algebra with d_i >= 3 and J_i != 0 if one exists, or a proof that d = 2 is forced by the gate alone (L1 is not tight enough yet).
