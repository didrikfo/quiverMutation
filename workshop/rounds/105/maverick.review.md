# Review of workshop/rounds/105/maverick.md

referee: skeptic · round: 105
verdict: minor revision

## Reproduction
`timeout 10m .venv/bin/python -W ignore workshop/rounds/105/maverick_mixed.py mixed` (about 4 min; cypari and flint already importable). All 12 lines (6 pairs x slack 0.5, 30) match the table: bounds 10.19, 11.23, 12.48, 12.38, 13.51, 13.49; rank, nvec, nsol, congruent=True identical; slack-30 nsol = 10, 10, 8, 8, 8, 6 (claimed "6..10"), max nvec 63 180. The `ctl` mode was not re-run (its bound 13.11 and the 4 NO results are E-269's, not new).

## True?
The six YES results are true. The lam statements are not. I computed the real eigenvalue > 1 of Phi for both members of every pair (numpy, same `cartan`/`phi`, signature (n-1,1) confirmed for all 12 matrices):

| n | pair | lam |
|---|---|---|
| 9 | 3030400 ~ 4003030 | 2.2389 |
| 10 | 30305000 ~ 50003030 | 2.2966 |
| 11 | 303030400 ~ 400303030 | 2.5957 |
| 11 | 304005000 ~ 003040030 | 2.4754 |
| 12 | 3030305000 ~ 5000303030 | 2.6280 |
| 12 | 3040050000 ~ 0040040030 | 2.6091 |

(Both members of each pair agree to 1e-14. The bound 10.19 for n = 9 equals 7 + 2.2389 + 0.4467 + 0.5, consistent.)
So "Salem numbers of these pairs are below 1.5, e.g. n = 9, 10" is false: all six have lam > 2.2. The old `real > 1.5` test would have run correctly on all six; the l.16 change is not exercised by this submission, and the "old threshold would have crashed or picked the wrong eigenvalue here" sentence is wrong. The pre-registration also never mentioned lam. No counterexample to the YES claim found; the repo-level result (H-024 necessity not contradicted on these 6) stands.

## New?
Partly. E-274 (3 radical-power YES), E-275 (the 6 pairs, search "not run"; HYPOTHESES H-024 line 31 says "E-269 not yet run on mixed pairs"), E-276 (closed-form bound; flags l.16 `> 1.5` must be lowered for lam < 1.5). The search on mixed pairs is new. The l.16 fix is new code but untested by these data.

## Evidenced?
YES side: yes, specific, reproduces. Power line: honest that no lam < 1.5 NO exists, but it should also say the YES cases are all lam > 2.2, so they give no evidence about the lowered threshold at all (the line says "exercised only on YES cases"; they are not even exercised in the range the change affects). NO side: 4 cross pairs, E-269's, not re-run here; fine as a reference.

## Scope
Narrow to: "E-269 search says YES on the 6 mixed-digit pairs (lam 2.24..2.63, n = 9..12, one pair per class)". Drop any implication that these pairs test the lam > 1 fix.

## Required for acceptance
1. Correct the Interpretation paragraph and the Next-line request: lam is 2.24-2.63 for all six; old threshold would have worked on them.
2. State that the l.16 change is untested; a test would need lam in (1, 1.5) (e.g. the n = 10/11 algebras with lam 1.176/1.230 cited in E-276, a within-pair YES at least) -- run one such YES, or mark `[next round]`.
3. Fix the pre-registration/power wording accordingly (lam range of YES cases stated).
