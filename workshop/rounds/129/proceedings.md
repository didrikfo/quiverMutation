# Round 129 -- proceedings
Cycle 2, phase 3 (round 1 of 4), ordinary. Step 0.5: round 128 left no questions; nothing to settle.

## Submissions
- **maverick (revision of 127)**: a pure-F key-guarded walk crosses at n = 11 (dual 303305000 -> 403033030, one J != 0 step, SNF(S) (1,2,4) -> (1,2,2)); the shipped `mutationSearchDepthFirst` reaches it, `refuseNonTilting=True` removes it; the "3 edges with no R inverse" are an ideal-vs-generator artefact. Referee skeptic: **accept** (independent Cartan/SNF replay of the 10 steps; withdrew the round-127 objection; extra run `REFUSE=1` on the dual root). Decision: **accept**; promoted as **E-332**. Wording kept to "a pure-F walk", not "the library search finds" 403033030.
- **theorist**: all 50 of E-323's J != 0 steps are non-tilting irreducible silting mutations by an independent Hom computation; the equivalence is already on record (E-130, E-171). Referee experimentalist: **accept** (0 failures on 7 500 further steps). Decision: **accept**; promoted as **E-333** (new: n = 11, 12 scope and the agenda correction).
- **skeptic**: Hom-level non-tilting of E-323's F3 step (Hom(T,T[-1]) = 2); pair separated by Phi-rank invariants that share the SNF(S) necessity chain; the NO needs only E-308/E-272. Referee toolsmith: **minor revision**; response did the three required items (checked: title and claim reworded, moduli stated, opinion marked). Decision: **accept**; promoted as **E-334**.

## Consequences
- **E-323 (4) "every crossing needs an R move"** and the round-127 reading of E-331 do not generalise: E-332 is a pure-F crossing at n = 11. H-032 status line updated. Docstring of `mutationSearchDepthFirst` ("guard keeps a walk in one derived class") is false for the walk as shipped (backlog seed, now with a counterexample).
- **Agenda corollary of round 128** "option-on walks cannot leave a class (relative to AI + E-327)": rests on H-023, not E-327; non-tilting does not imply out of class (E-316). Thread T10 updated.
- E-323's blocker is half done: the Hom-level check is supplied (E-333, E-334); no non-SNF argument for the pair exists independent of the necessity chain (a second computation only). With E-332's independent SNF replay, the F-candidate (E-323/E-332 crossings) is written at conference 132 at its scope or retired.

## Questions for the steering committee
1. Should `refuseNonTilting` become the default of the library walk, given E-332 (default walk leaves the class at n = 11) against E-316 (19 J != 0 steps stay in class, so refusing may cost connectivity)? **Recommended: no, keep it off**; document the crossing in the docstring and have toolsmith run F-041 at n = 10, 11 with the option on first. If unanswered, the chair of round 130 takes this.

## Decisions taken for the steering committee
(none from step 0). Chair: no overnight run proposed (toolsmith's depth-7 n = 11 ball is 10-25 min per root and fits the per-command limit). Breadth deviation recorded in call.md: next ordinary round takes a dormant thread (T4/T6).
