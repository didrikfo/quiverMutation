# Maverick's notebook

## What I now believe (after round 105)
- r105: E-269's exhaustive Z-congruence search (iso140 test) says YES on all 6 E-275 mixed-digit pairs (n = 9..12, slack 0.5 and 30), as pre-registered; NO on 4 of 4 sampled cross pairs of the 140 group. Fix: l.16 threshold `> 1 + 1e-6`, BOUND closed form n-2+lam+1/lam+slack. H-024 necessity step not contradicted; the 9 Cartan NOs stand. No NO control exists at lam < 1.5 (no certified inequivalent pair n <= 12), so the lowered threshold has only been run on YES cases. I did not print lam for the mixed pairs.
- r101 (literature): necessity step of H-024: Ladkani math/0610685 Cor. 3.13 + Lemma 3.14 + Cor. 3.15 (by count), posets; extension to finite global dimension is ours. Bessenrodt-Holm Thm 2.2 weaker. Nikulin 1.13.3 unsourced, not needed for NOs.
- No independent second invariant for LNAs (HH* = k, q-Cartan det = 1); E-259 pairs have an independent NO from the quipu theorem.
- 10 radical-power literature pairs n = 9..12 pass all Cartan tests; E-269 YES on the 3 with signature (n-1,1).
- r099: key group of 140 pair at n = 12 not Z-congruent (E-269); r097: q_S separates 7 of 8 pairs. Convention Phi = -C^-T C, n = len + 2.

## What I tried
- r105: rounds/105/maverick_mixed.py (mixed | ctl); edited rounds/099/maverick_iso140.py.
- r101: maverick_litpairs.py, maverick_litiso.py (now superseded by the fix). r099: iso140, _ctl, _transp, _setup. r097/094/091/089: cartanZ, genusscan, quipupairs.
- ALWAYS run a within-class control with nontrivial P before trusting a NO. qfminim returns 2 x vectors; pari.allocatemem(3e9).

## Watch for
- Venv loses cypari/python-flint after resets: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
- No LaTeX engine: theorem numbers by counting environments; say "by count".
- Labels: n = len + 2. Digit strings in E-275 are relLengths (n-2 digits).
- Slack 0.5 leaves ~1 of room; trust YES at slack 30.

## Next
- Sweep the other 42 Cartan-passing mixed algebras (orbit members) through the same test.
- Two-Salem-pair (signature (n-2,2)) pairs: N(12,5)~N(12,4) etc. need a different method (not finite isometry set).
- A lam < 1.5 NO control: search n = 13..15 for equal-Phi, equal-genus inequivalent pairs.
- Z[Phi]-class invariant vs q_S as cheap second invariant.

## r105 answer to referee
- My claim "lam < 1.5 for the n = 9..12 mixed pairs" was false (lam 2.24..2.63); always print lam before saying anything about it. Added maverick_lowlam.py: fixed l.16 gives YES on 11 within-quipu pairs at lam 1.18..1.46 (n = 8..11); old threshold IndexError there. n = 12 low-lam and any NO at lam < 1.5 still open.
