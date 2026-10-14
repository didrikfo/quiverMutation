# Round 064 -- conference proceedings

All six personas wrote position statements (`rounds/064/<id>.md`); nothing was run. Five of six (theorist and maverick excepted) name the same single most promising question: **a J != 0 step on a join path, with a join test that can say NO across it** (experimentalist, skeptic, scholar, toolsmith; the theorist names G2 itself, the maverick the n = 11 class count). All six name the same weakest claim: "25 of 25 failing n = 7 children are in the derived class", which rests on the J = 0 premise (End(T) cannot see J, E-181; Thm 1.1 does not give J = 0, E-189; G2 only as a span statement, E-183; no join test has ever returned NO across a J != 0 step). The skeptic asks that every use of "25 of 25" carry "under the J = 0 premise"; STATE does so.

## Ledger
E-entries since the last conference (E-183 .. E-191):
- **H-015** (guard sufficient): stays OPEN. E-183, E-184, E-186, E-189, E-191 add power controls and narrow the premise to G2 beyond span plus the left/right convention; none of them tests J. Key-guard law remains class 0 only.
- **H-010** (overlap reducible only at an end): stays SUPPORTED. E-187 proves the single-mutation rule for right mutation relative to the repo's `procedure`; left mutation, two bystanders, k >= 4 and the reading against arXiv:2112.08129 open. No status change (the proof is of a rule about `procedure`, not of H-010 as worded).
- **H-020**: stays SUPPORTED. E-188 closes the width 6..8 length gaps (0 failures in 8 248); widths 9..11 and lengths above 13 untested; H3 (completeness) still the weakest part.
- **H-021, H-017**: stay OPEN; untouched by E-183..E-191 apart from E-185 (cord-bearing control passes), which does not change H-017's status.
- **T3/T8**: E-190 corrects "3 candidates" to "2 classes" at n = 10; no H- or F- text depends on it. Candidate C is open at n = 11 only.
- No new finding (F-), no retraction (R-). No thread closed: T7 and T4 stay as breadth slots; T9 stays parked.

## Proposed agenda (ranked)
1. **T10/T5: give the J = 0 join test a J != 0 step.** Toolsmith: size with `--plan` a depth-limited ball search at n = 8 or 9 that records every J != 0 edge on a path between two LNAs known equivalent, or states that none exists; any path found is replayed independently by the skeptic. Parallel (theorist, scholar): G2 beyond span on the parallel-arrow and cyclic cases (generator count, sign, transpose), Sect. 8 and the cancellation step of 1504.02617. Skeptic: left/right convention by hand on the APR example (Oppermann l. 403).
2. **T3/T8, Candidate C at n = 11 (maverick, toolsmith).** Count derived classes inside the four periodic groups (535, 64, 24, 15) with an invariant finer than signature; tie to n = 10 by deletion. Level: idea; one script, about 40 s. Scholar first checks whether 0911.5137 / 1310.1557 already give the count.
3. **S-1 / vertex deletion at n = 9, 10 (experimentalist, toolsmith).** Pure counting on known classes (F-011, E-175/E-178 pairs): how often do same-class pairs map to the same class under a fixed rule for the deleted vertex. A negative result is welcome. Breadth slot (dormant since 055).
4. **T7 (theorist, scholar).** Left mutation by duality; the induced rule on LNA rows (theorist's new line: after a right slide with (r_{v-1}, r_v) = (2, 0), which v' give an LNA again, first test on H-010's frozen pair (8:3)(9:3)); reading against arXiv:2112.08129. Experimentalist's cheap T4 closure (widths 9..11 at w+3..w+4, `rounds/062/experimentalist_gaps.py --plan`) is a filler slot.
5. **Housekeeping:** `compare2` / docstring rewords and `canonicalKey` cap (toolsmith); the Smith form of C + C^T for the F-010 pair vs E-154 (one line, maverick; H-017 stays parked).

Disagreement: none on priorities; the theorist ranks G2 above the J != 0 search, the other four the reverse. The agenda runs both in parallel (item 1). The maverick alone puts Candidate C first; it takes item 2.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- No unanswered questions from round 063. No overnight run proposed (nothing above exceeds the 10-minute command limit).
- arxiv.org: the scholar asked whether it is reachable; STEERING (2026-10-09 note by the human) says yes, and STATE's "blocked (403)" standing fact was stale; it is corrected in STATE. Scholar may fetch arXiv:2112.08129 under the standing permission.
- Proposed agenda set as the open threads; ordinary rounds work from it until the human changes it in `STEERING.md`.
