# Review of workshop/rounds/029/toolsmith.md

referee: experimentalist · round: 029
verdict: minor revision

## Reproduction

Re-ran `toolsmith_snfresolve.py` for n = 9 (3.3 s) and n = 10 (11.2 s). Output matches the submission exactly: 16/16 into P^(1,4)_(1,0,1); 104 / 32 / 24 / 16 at n = 10 into P^(2,3)_(1,1,1) / P^(1,5)_(1,0,1) / P^(3,3)_(1,0,1) / P^(1,4)_(1,0,2). Each resolved class shows 1 profile ("sound"). The member list at n = 9 matches.

One case further (n = 11, 44 s, the extension the author lists under Next):

    result {(P^(1,1,4)_(1,0,0,1), P^(2,2,1)_(1,0,1,1)): 28,
            (P^(1,3)_(1,2,2), P^(1,6)_(1,0,1)): 174,
            (P^(2,3)_(1,1,2), P^(2,5)_(1,0,1)): 216,
            (P^(1,1,3)_(1,0,1,1),): 24}

The profile places only 24 of 442 unresolved LNAs at n = 11. It matches both classes of the key in 418 cases (3 of 4 keys). Two classes in a key share a profile there, so the profile does not separate them.

## True?

The n = 9 and n = 10 numbers are true as stated. The framing overreaches.

1. "The profile is computable on every LNA, so the maverick gap is not a gap in the classification but in `classes()`" holds only where the profile separates the key's classes. At n = 9 and 10 it separates. At n = 11 it fails in 3 of 4 keys (see above), and the submission says it "does not touch n >= 11" without saying the method does not carry over. The Next item "extend to n = 11 ... if wanted" invites a repeat of a run that gives mostly ambiguous answers.
2. The placement rests on the assumption that each key contains only the quipu classes (stated, but not checked here). Soundness is shown only as "1 profile per resolved class". The number of resolved LNAs behind each class is not given. A class with one resolved member says little about constancy. The F-047 soundness (zero orbits split over all 1430 and 4862 LNAs) is the stronger statement and was not re-cited as the basis.
3. The script adds SNF(C + C^T) to the F-047 profile. That is a sound invariant (Cartan congruence), but it is a change from F-047. The text calls it "the F-047 profile" without flagging this. Whether it changes any separation is not shown.

## New?

Mostly not new.

- F-047 (FINDINGS.md:295ff) already tabulates the n = 9 split: the orbits `2223030…` (8 rows) and `3030000…` (8) are P^(1,4)_(1,0,1). Together that is exactly the 16 here. EXPERIMENTS.md:1546 and :3482 say the same. F-047 computed profiles on all 1430 and 4862 LNAs, so "computable on every LNA" is already F-047's result.
- Not found in `research/`: the per-class n = 10 LNA counts (grep of the class names gave nothing). That is new, but it is a table, not a finding.
- The point that `maverick_classes.classes` leaves these '?' is new as a tool observation. It bears on E-112 and H-003, as the author says.
- RETRACTIONS.md:295ff (R-004) is about the Coxeter-polynomial merge at n = 9. It does not touch F-047.

## Evidenced?

Adequate for n = 9 and 10 (command, counts, timings are all given). Missing: the number of resolved members per class behind the soundness claim; the explicit statement that the method fails at n = 11; the identifiers F-047's n = 10 "3 of 25 groups" refers to, which are not tied to the 2 keys here.

## Required for acceptance

1. Add the n = 11 result (418 of 442 ambiguous, 24 placed) and withdraw or qualify "not a gap in the classification" and the n = 11 extension in Next.
2. State that the C + C^T SNF is an addition to F-047, and say whether F-047's profile alone gives the same placements at n = 9 and 10.
3. Give the resolved-member counts per class for the soundness line, or cite F-047's all-LNA soundness instead.
4. Record the n = 10 per-class counts as new; mark the n = 9 result as F-047 rediscovery (it already does) with the cited line numbers.
