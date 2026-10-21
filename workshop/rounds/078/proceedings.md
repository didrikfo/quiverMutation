# Round 078 -- proceedings
kind: ordinary. Merge of origin/main: already up to date. Round 077 left no questions; nothing to settle.

## toolsmith (T10) -- refereed by skeptic: minor revision, response done
Claim: the 8 E-216 chains have recorded F/R move lists and verify edge by edge (Hom + Cartan, sign-aware iso at the meets); child and parent end at the same seed in 3 of 4 reps; `mutateAtVertex` refuses a cyclic D(i).
Decision: **accept**. The referee reproduced 8 of 8 and the suite; required items (broken recorded scripts, left/multi-vertex and loop tests, seed indices, title) done; I checked the script diffs are `allowCyclic=True` only. Promoted as E-222. Verifier power is thin (E-212); five argument-needing caller scripts were not rerun (on the board).
Decision on the code (the author asked the chair): **keep** the refusal. The search never reaches it (E-204, E-223); direct callers that hit it were wrong in the E-211 sense.

## theorist (T5) -- refereed by scholar: minor revision, response done
Claim: the guard does not force D(i) acyclic; inside the search step 7 runs only at acyclic nodes.
Decision: **note** (rediscovery of E-211 plus a one-line reading of search.py l. 420/437). Recorded as E-223. T5's item "does the guard force D(i) acyclic" is answered NO; the search is covered by E-199.

## experimentalist (S-1) -- refereed by theorist: minor revision, response done
Claim: no labelled-key meet for 90000000 / 00000030 at 5 + 4 under J = 0 F/R moves (positive control 5 of 5); n = 15 scan shard plan.
Decision: **accept, narrowed** (a bound, labelled key, move set not Hom-verified at n = 10; the control tests the join, not completeness). Recorded as E-224. Levels 4-5 and the killed `run 5 5` were not rerun.

## Consequences
None contradict a finding or hypothesis. E-216's "not Hom-verified" is now Hom-verified for its 8 chains (E-222); "25 of 25 / 70 joined" still rests on the J = 0 premise; the rep 2 seed link (1 vs 8) is open on T10.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- No unanswered question from round 077; nothing recorded in STEERING.
- Overnight: none approved (n = 15 scan, 7 + 7 join are proposals only; no `--budget-hours` script).
- `CyclicDescendantsError` kept in `procedure.mutateAtVertex`.
