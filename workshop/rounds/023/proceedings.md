# Round 023 -- proceedings

Ordinary. Worked: skeptic, scholar, maverick. Referees: theorist (on skeptic), experimentalist (on scholar), skeptic (on maverick). Questions of round 022 were unanswered; I took the recommended options (see below).

## skeptic -- "reject iff long square" off the walks (T5)
Claim: off the walks a genuine long relation always makes the step reject (362 apparent counterexamples were redundant presentations); 26 gate-admitted rejections have no square as `hasLongSquare` tests it (19 two-out-arrow, 7 one-out-arrow); two examples are unreachable by the key-preserving walk.
Referee (theorist): minor revision. Counts reproduce (n = 6). The 19 two-out cases are False by construction of `hasLongSquare`; reachability not re-run; "none reachable" in the title overreaches; "genuine" needs a definition.
Decision: **accept**, scoped as the referee asked. Promoted as **E-105**, with the straw-man point and the untested reachability in Limits.

## scholar -- derivation of "reject iff long square" (T5)
Claim: square => reject needs only a minimal relation `c·alpha` with `c` a combination of >= 2 paths; reject => square additionally needs out-degree 1; on class-0 walks at n = 5..7, J != 0 iff out-degree 1 and long square (1 139 steps at n = 6, 162 at n = 7).
Referee (experimentalist): minor revision, but decisive: n = 8 class 0 (200 s) has 42 gate-admitted J != 0 steps with out-degree 2 and no long square (and 2 with out-degree 1), so "walks reach no D-type" and the step-by-step iff hold for n <= 7 only. Four fixes required (n = 8 in the table, identify the 42 and run the Cartan test, reword Next, state class/caps in the claim).
Decision: **revise**. Not promoted. The n = 8 observation is recorded inside E-105's Limits as the referee's run.

## maverick -- two-big-relation LNAs (T6)
Claim: D1 holds on all 236 / 981 LNAs with >= 2 big relations (n = 8, 9); E-103's peeling formula fails on 1/6/24 (n = 8/9/10), all depth 2; a fitted mirror-chain rule has 0 mismatches at n = 8..10 and 8 of 8 out-of-sample predictions.
Referee (skeptic): **accept**. All numbers reproduce; ran n = 12 (5 of 5 exact, including a depth-4 case `2223030222` and shapes other than `303`); nothing prior beyond E-103.
Decision: **accept**. Promoted as **E-106**, including the referee's n = 12 run. A fit, not a proof; the data at n = 8, 9 cannot separate the chain from "blocked implies depth >= 2".

## Questions for the steering committee
1. **Agenda:** keep the round-020 agenda (recommend yes). Round 024 is a conference.
2. **Overnight:** none (recommend). The scholar's n = 8 class-0 walk gives the out-degree-2 rejects in 200 s; nothing needs more than 10 minutes.

## Decisions taken for the steering committee
- Round 022, question 1 (agenda): kept the round-020 agenda -- decided by the chair of round 023; no answer from the human.
- Round 022, question 2 (overnight): none -- decided by the chair of round 023; no answer from the human.
