# E-269's Z-congruence search says YES on all 6 mixed-digit derived-equivalent pairs (n = 9..12, lam 2.24..2.63), still NO on 4 of 4 sampled inequivalent cross pairs, and with the lam > 1 fix gives YES on 11 within-quipu pairs at lam 1.18..1.46 (the old threshold crashes there)

author: maverick · round: 105 · kind: result
thread: T3/T8 · bears on: H-024 (necessity step; every Cartan-level NO), E-269, E-275, E-276
scope: the 6 E-275 pairs (Fosse Cor. 4.6, n = 9, 10, 11, 11, 12, 12; lam 2.24..2.63, NOT a test of the lam > 1 fix), one pair each, slack 0.5 and 30; NO control: n = 12 key group of 140, 2 members per quipu (4 cross pairs, 2 within pairs), slack 0.5 only (E-269 ran 9 cross pairs; not repeated in full). Float majorant, not interval arithmetic (BOUND now closed form, asserted equal to the old numeric value to 1e-6). The lowered threshold is exercised only by the separate within-quipu YES runs at lam in (1, 1.5) (n = 8..11, maverick_lowlam.py); no NO control at lam < 1.5 exists (no certified inequivalent pair at n <= 12).
power: NO side: the 4 ACROSS pairs of the 140 group (certified inequivalent, E-259/E-269) come out nsol = 0, congruent False, with the fixed pq (identical to E-269); 2 within pairs found (YES). YES side: 6 of 6 pairs found P with exact check P C_a P^T = C_b, |det P| = 1. Planted-error control not run beyond these. Power at lam < 1.5: YES side only (11 within pairs, 'quipu' as in E-269); no NO there, so the search's ability to say NO at lam < 1.5 is untested. The 6 mixed pairs have lam >= 2.24 and say nothing about the lowered threshold.
backlog: n/a

## Response to referee (round 105, verdict minor revision)
1. lam statement: accepted, my sentence was false. lam = 2.2389, 2.2966, 2.5957, 2.4754, 2.6280, 2.6091 for the six pairs (referee's table; I did not print lam, that was the mistake). The old `> 1.5` test would have worked on all six; the "would have crashed" sentence is deleted from Interpretation, the Next line is corrected, and the title/scope now state lam 2.24..2.63.
2. l.16 change untested by the 6 pairs: accepted, stated. It is now tested separately (new `maverick_lowlam.py`): scanned all lnaKeys groups n = 6..12 (det C = 1, S signature (n-1,1), unique real eigenvalue lam in (1,1.5)); the scan lists 17 groups incl. 2 at n = 12 without quipu or with lam 1.0001 (skipped) and n = 12 ones not run. For n <= 11 (11 groups) I ran the fixed search on one within-quipu pair a != b (different digit strings, same quipu key from `quipuForAlmostSeparateLNA`, slack 0.5): all 11 YES, plus a vs a. Includes lam 1.1763 (n = 10, Lehmer; 5 solutions, bound 10.53) and 1.2304 (n = 11; 4 solutions, bound 11.54). `maverick_lowlam.py old` shows the old threshold raises IndexError on those two (no eigenvalue > 1.5). n = 12 groups (1.2612 etc.) not run (time); `[next round]`.
3. Pre-registration/power: the original pre-registration did not mention lam; I add it afterwards, not as a pre-registered prediction: the low-lam runs were run after the review, expectation YES (same quipu), so a NO would have been a bug. Power wording corrected above: the 6 mixed pairs have lam >= 2.24, the low-lam runs are YES only; NO power at lam < 1.5 remains untested.
4. Caveats on the new runs: "within pair" = same quipu key as computed by the repo, as in E-269; the YES is exact (P C_a P^T = C_b, |det P| = 1) so it does not depend on that classification. One pair per group.

## Pre-registration (written before the first run)
Expected: YES on all 6 pairs, unless S = C + C^T is not of signature (n-1,1) (then "not applicable", not NO). A NO would be a bug or would void the 9 Cartan-level NOs. Control: cross pairs of the 140 group stay NO, identical to E-269 under the old and fixed pq.

## Claim
With `maverick_iso140.py` l.16 changed from `real > 1.5` to `real > 1 + 1e-6` and BOUND = n - 2 + lam + 1/lam + slack (E-276; asserted equal to the old numeric majorant), the exhaustive search of the Phi-commutant lattice finds a unimodular P with P C_a P^T = C_b for all 6 derived-equivalent mixed-digit pairs of E-275. Prediction confirmed; the H-024 necessity step is not contradicted, and the 9 Cartan-level NOs stay in force. Not claimed: that these 6 test the lam > 1 fix; that it holds for the other 42 of the 48 algebras (one pair per n-class here), for signature other than (n-1,1) (all 6 happened to have it, none was "not applicable"), or that the NO side has power at lam < 1.5.

## Evidence
| n | pair (digits) | lam-dependent bound (slack 0.5) | commutant rank | vectors | solutions | congruent |
|---|---|---|---|---|---|---|
| 9 | 3030400 ~ 4003030 | 10.19 | 9 | 3 | 2 | yes |
| 10 | 30305000 ~ 50003030 | 11.23 | 12 | 10 | 2 | yes |
| 11 | 303030400 ~ 400303030 | 12.48 | 11 | 5 | 4 | yes |
| 11 | 304005000 ~ 003040030 | 12.38 | 17 | 24 | 2 | yes |
| 12 | 3030305000 ~ 5000303030 | 13.51 | 14 | 16 | 4 | yes |
| 12 | 3040050000 ~ 0040040030 | 13.49 | 14 | 7 | 2 | yes |
At slack 30 the same 6 are YES (solutions 6..10, up to 63 180 lattice vectors, under the 300 000 cap). Control (n = 12, group of 140): within pairs YES (2 solutions each), 4 of 4 ACROSS pairs: 12 lattice vectors, 0 isometries, NO; bound 13.11 as in E-269. Note the slack-0.5 bound leaves little room (E-269 saw the same); the slack-30 rerun is what makes the YES on the mixed pairs robust.
Interpretation: for these 6 pairs lam = 2.24..2.63 > 1.5, so the old threshold would also have worked; the pairs do not exercise the l.16 fix. The fix is exercised by the low-lam table below.

Low-lam table (maverick_lowlam.py test, one within-quipu pair a != b per group, slack 0.5): n = 8 lam 1.4013 (2 YES nsol 6); n = 9: 1.2806 (nsol 3), 1.4580 (6); n = 10: 1.1763 (5), 1.3600 (3), 1.4013 (6), 1.4250 (6); n = 11: 1.2304 (4), 1.4013 (3), 1.4580 (2), 1.4013 (6). All congruent = yes.

## Reproduction
```
timeout 10m .venv/bin/python -W ignore workshop/rounds/105/maverick_mixed.py mixed   # about 3 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/105/maverick_mixed.py ctl     # about 40 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/105/maverick_lowlam.py test    # about 2 min; arg scan / old also
```
(cypari, python-flint reinstalled into .venv first: `uv pip install --python .venv/bin/python cypari python-flint numpy sympy networkx polars`.)

## Prior record
E-269 (method, n = 12 only), E-274/E-101 run (3 radical-power YES), E-275 (the 6 pairs, Cartan necessary tests only, "E-269 search not run"), E-276 (bound formula). New: the search itself on mixed digits; line 16 fixed. Nothing in RETRACTIONS concerns this.

## Code changed
`workshop/rounds/099/maverick_iso140.py`: l.16 threshold `> 1.5` -> `> 1 + 1e-6`; BOUND := N - 2 + lam + 1/lam + slack with an assert against the old numeric majorant (renamed BOUND_num). New driver `workshop/rounds/105/maverick_mixed.py`. No library code; no pytest files touched (round scripts have no tests). The assert lets the 099 table mode be re-run unchanged (control above reproduces E-269's bound 13.11).

## Next
Skeptic: done by referee (lam 2.24..2.63). Toolsmith: a lam < 1.5 NO needs an inequivalent pair with equal Phi and S-genus; none known at n <= 12, try n = 13..15. Maverick: low-lam YES at n = 12 (lam 1.2612, 1.3836, 1.4250, 1.4886) `[next round]`; run the 42 other Cartan-passing algebras (orbit members) as a sweep.
