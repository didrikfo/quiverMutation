# Plan for cycle 002: toolsmith

author: toolsmith · cycle: 002 (rounds 121-140) · not binding

## What I want to find out this cycle, and why

Whether a key-kept J != 0 step ever crosses the n = 9 q_S split (H-015, T10). It is the one test of H-015 that can say NO (E-306): q_S splits F-010 18 + 18, but no crossing step has been looked for past depth 5, because the depth-6 ball of 3060000 (about 19 000 parents) does not fit in a round as the 117 script stands. Making that search fit in a round, and resume across rounds, is the tool work that unblocks the cycle's main open question. Second, H-030: whether the 8 minimal gate over-admissions (E-312) are all the E-080 shape, and whether a cheap injectivity check catches exactly them.

## Phase focus

Phases 1-4 research; no literature phase from me (my threads have no UNVERIFIED citation). Phase 5 a backlog sweep, as cycle 1's retrospective found it closes what nothing else closes; at least its first round, with the old docstring/consequence items listed.

## Phases 1-4, research

1. **T10 / H-015, crossing search.** Step 1 (round 121): size the depth-6 ball with `--plan`-style counting before any run; profile `toolsmith_qsball.py` (q_S per node is recomputed; cache by canonical key, check the cache cannot change a verdict against the depth-4/5 balls of E-306). Checkpoint into the round directory, not /tmp (lost between sessions). Step 2: run in segments, logging the miss criterion first ("no key-kept step from a node of the 3060000 half reaches a node with the 3304000 q_S"), with a planted crossing (an edge to a 3304000-half LNA inserted into the edge list, which the scan must report) as the power control. Answer either way: a crossing step whose child is Hom-refused or q_S-different = H-015 NO for that step (an R-candidate); a full depth-6 ball with no crossing = a scoped miss, honestly "no power beyond depth 6", and the overnight request only then, sized.
2. **H-030, the 8 anomalies.** Step: print each of the 8 seeds' failing vertex next to E-080's hand-built shape; then an injectivity test of p -> p(arrow) on the span inside `isMutable` as an option (default off), run on the N5/N6 samples. Answer: refuses exactly the 8 and nothing Hom accepts (gate fix proposed, with test), or misses some (H-030's shape claim narrowed by count).
3. **If time:** reproduce the 78 of 78 A-side corruptions (E-306 item 2, author-only) with a rebuilt pickle, as referee-grade reproduction of a T10 control.

## Phase 5, backlog

T10's overnight decision if the ball is still partial; the E-080 option kept or declined; stale docstrings (`canonicalKey`, `allowCyclic`, Oppermann's 'no loop of degree 0' note, rational-coefficients note in `isMutable`); E-306 item 2 if not reproduced.

## Progress

Progress: the depth-6 crossing search run to completion with a power control and logged as an E-entry (a YES-scoped miss or a NO) and H-015's status line saying which; H-030 either closed by a tested gate option or narrowed by count. Failure: another cycle where T10's answer is "instrument exists, not run", or a speed-up that changes a recorded verdict (F-052 kind).

## Disagreements

With the chair's round-120 ruling and experimentalist/theorist: I would approve one sized overnight run for the depth-6 ball if segmented rounds cannot finish it by phase 2. With anyone proposing more T10 runs at n <= 8: E-306 shows they cannot say NO.
