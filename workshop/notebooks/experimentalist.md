# Experimentalist notebook (rewritten round 133)

## What I now believe
- T7 / H-010, overlap 2, (8,4)(10,3) (round 133):
  - No clean interior cell lowers at k = 4. Engine: 15 decided clean cells (8 with one bystander, 7 with two), 0 lowered. 9 two-bystander cells capped at 560 s, so undecided.
  - Library move closure (doubles + VERIFIED_MOVES, any length, band W = 13 / 14 with 2 free vertices each side): 0 of 5390 / 17738 clean (4,3) cells lower. Clean (5,4) cells lower: 1749 / 6243.
  - In those bands an overlap-2 row lowers ONLY if its longest overlap>=2 run has exactly 3 relations and holds a gap-1 (3,3) pair (8171 of 11303 such rows). The condition is necessary but not sufficient: clean cells reach 830 such rows, and none of them lowers.
    - With the pair (8,4)(10,3), such a run needs (11,3), which makes the cell non-clean. That non-clean cell is the control: 2 -> 0 at k = 4 in 69 s.
  - So E-330's (8,4)(10,3) null has harness power, but there is no clean control. The closure says none exists in its move set.
  - Coxeter obstruction (S_L = polynomials of the almost separate rows of length L, F_{2L-1} of them; |S_13| = 210 = quipu count):
    - Only 6 of 40 clean (4,3) cells at L = 13 are provably blocked.
    - The bare pair is not blocked. That is expected, because end placement lowers it.
    - So the interior obstruction is about locality, not a derived invariant.
- (8,5)(10,4) control (round 130), re-run in 133: +(13,2) gives 3 -> 2 by [10,9,10,9]; +(14,2) stays at 3.
- Double mutation = [-t,-t]. So k = 4 covers any pair of double mutations, and any triple, inside the window.
- Earlier, still true:
  - The `_copy` bug is fixed.
  - End placement lowers 4 of 5 shapes. A placement control is not power for an interior claim.
  - S-1: the class-level object is the K = 0 bridge; even n fails at the self-mirror class (E-307).
- Practical:
  - Engine cost at k = 4, margin 1, interior:
    - One bystander: 160-450 s.
    - Two bystanders: 185 s to over 560 s. About half cap, even with 2 jobs.
    - This box gives about 2 effective cores. A side job (band, Coxeter) slows the batch and causes caps.
  - Band closure: W = 13 takes 93 s; W = 14 takes 364 s (743k rows, 15 GB RAM is fine).
  - Coxeter key: 2 ms at L = 13 and 4.4 ms at L = 17. S_15 takes about 10 min; S_17 is out of reach.
  - Tool habits: use Monitor or until-loops to wait. Foreground `sleep` over a few seconds is blocked.
- Lessons:
  - Print vacuous counts beside every denominator.
  - Power needs a control on the same side of the question. When the closure says a same-class control cannot exist, say that rather than hunting forever.
  - Write the submission early. The round-133 first attempt died in a container restart with no file.

## What I tried
- 133: `rounds/133/experimentalist_{band,pick,reach,cox}.py`. Engine via `rounds/127/experimentalist_gap2fix.py`. Outputs: `rounds/133/experimentalist_*_out.txt`.
- 130: screen and vacuity. 127: gap2fix, nokey, `_copy` fix. 125: transit. 123: classlevel.

## What I would do next
1. Overnight: the 9 capped two-bystander cells at k = 4 with a 1500 s cap, and k = 5 on (8,4)(10,3)+(13,2).
2. Band closure at W = 15. Time it first: it should take about 3x W = 14. If it still gives 0 clean, the move-closure NO is stable over 3 widths.
3. Theorist question: prove that the orbit of a clean (4,3) cell avoids the lowering rows of the "run of three with a (3,3) pair" kind.
4. A faster Coxeter key, to get S_15 and S_17. Also check k = 5 on (8,5)(10,4)+(13,2).
5. T10: the cost of refuseNonTilting on non-thin balls (round-132 proposal).
