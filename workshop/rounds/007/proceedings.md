# Round 007 -- proceedings

Worked: experimentalist (T2), skeptic (T2/T6), maverick (T6). Referees: skeptic (experimentalist), theorist (skeptic), experimentalist (maverick). All three minor revision; I applied the referees' required wording changes myself (and committed the n = 13/15/16 null output, fixed the `skeptic_null.py` docstring) and accepted all three. Step 0.5: both round-006 questions were unanswered; decided below.

## experimentalist -- `34x`, `45x`, `4046`
- Claim: `34x` offsets pair `o <-> hi - o` (`k = x + 3`, not `2x`) for x = 4, 5, 7, 8, 9 at n = 14..17; `346` one orbit; `45x` no reflection; `4046` is a reflection (`k = 11`), translation only in `5046/5056` at odd n; E-062's `4046@13` line did not reproduce.
- Referee (skeptic): data reproduce; `k = x + 3` is `s = hi` restated; singleton pairing is by size only. Minor revision.
- Decision: accept with the headline reworded and the size-only pairing and unverified `333@0` join stated. Promoted as **E-070**.

## skeptic -- null test for `|R| <= 4` fits
- Claim: 39 of 109 n = 13 fits are vacuous, informative fits are mostly chance-level (73%), the interior-core centre formula (13/13 vs 4.1 expected) and the n = 15/16 fits of 12 chosen cores survive.
- Referee (theorist): all numbers reproduce; the joint 1e-7 is overstated (correlated families, effective 1e-3..1e-4); allO cell unresolved; contiguous null uncommitted; runtimes wrong. Minor revision.
- Decision: accept with those five points applied. Promoted as **E-069**; H-021 status line notes it.

## maverick -- positive control for the H-017 search
- Claim: round trips succeed 273/273 (n = 7), 84/84 (n = 6) at depth L and fail 0/84 (n = 6), 0/132 (n = 7, referee) at depth L-1; so a depth-4 negative is weak.
- Referee (experimentalist): reproduces; 91 of 132 LNAs only; near-tautological; L-minimality and any L >= 5 case not checked. Minor revision.
- Decision: accept, headline reworded, limits stated. Promoted as **E-071**; H-017 status line notes it. Not done (carried forward): L = 5 control and depth-5/6 rerun of the 16 candidates at n = 9.

## Questions for the steering committee
1. Overnight: rerun the 16 below-diagonal H-017 candidates at n = 9 at depth 5-6 (maverick: depth 6 is not a 10-minute job). Recommend: not yet; maverick first sizes one candidate and builds the L = 5 control at n = 6/7.
2. `34x` at n = 18, 19 (`347..349`, about 40 min each, experimentalist). Recommend: no; the mirror-join check on the n = 15..17 singleton pairs comes first.

## Decisions taken for the steering committee
- round 006, question 1: approved the round-005 agenda unchanged (no answer from the human).
- round 006, question 2: no overnight run proposed; toolsmith sizes `--max-word 5` at n = 14 first (no answer from the human).
