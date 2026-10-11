# Experimentalist notebook (rewritten each round)

## What I now believe (after round 057)
- Power control for the J = 0 join test (T10, agenda 2) at n = 9: F-010 pair 3060000 / 3304000 (equal polynomial, F-047 profile differs, certified inequivalent) has 0 J = 0 joins at depth 2..6 per side (balls 4080 / 2304); certified-equivalent pairs do join (3060000~3030000 total 5; 3060000~6000030 total 9, 6+6; 2223030 and the second quipu pair missed at 4+4). So the test can say "no" and "joined" is not produced for an inequivalent pair, to depth 6.
- BUT in all those balls the J = 0 ball equals the gate+key ball and all 1551 edges up to depth 4 pass the Hom test: no J != 0 step occurs. The control does not test the premise; it tests that key-guarded tilting steps do not cross a Z-conjugacy invariant (nearly tautological).
- Only n = 9 F-010 is certified (plus 3 n = 10 profile-separated groups, rows not collected; the Phi^18 group is uncertified, E-172).
- Earlier (054): group A (05040330 -> 33460000) has 5 total-7 J = 0 paths, non-labelled meetings (E-169); F-037 5 paths; n = 10 merges.py depths 3-5 finds no link (E-164).

## What I tried
- 057: `experimentalist_powerjoin.py` (J0 / ALL walks, WL + VF2 join, reuses amerge), `_powerhom.py` (skeptic_tilt.stepTest on every edge). Costs: depth 4/5/6 per row 19/53/155 s per mode.
- 054: amerge fwd/back/join; 051: merges.py, f037replay; earlier deepreplay, blocks, BFS, W walks.
- Lessons: background + until-loop polling; never pkill -f; a miss at fixed depth is not a verdict (2 of 4 equivalent controls missed at 4+4); relabelling-aware join needs hash + exact iso.

## What I would do next
1. OVERNIGHT proposal (written in the submission): F-010 pair at depth 7-8 per side, report first J != 0 edge in either ball and its Hom verdict (growth ~2.7x/level).
2. Collect n = 10 rows of the 3 separated groups (maverick) and run the same script with depth 4.
3. Same join for F-037 (34504030 -> 50505000) with amerge parametrised; depth 9 n = 7 c1, c2 (049).
- Watch: depth < 7 never meets a J != 0 step at n = 9; do not quote this as support for the premise.
