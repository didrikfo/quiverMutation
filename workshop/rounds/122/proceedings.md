# Round 122 -- proceedings
kind: ordinary · cycle 2, phase 1, round 3 of 4. Step 0.5: round 121 left no questions; nothing to settle.

## Submissions
- **toolsmith** -- at n = 11, in the 72-group (13 + 59 by SNF(S)), a gate-admitted key-kept J != 0 step changes the SNF tail (2,4 -> 2,2); a 5-step guarded F/R walk joins inequivalent LNAs 003033030 and 030330030. Referee skeptic: **accept** (re-ran scripts; own rational Cartan matrices; SNF via mod-4/mod-8 solution counts; not an R-016/E-312 artefact). **Decision: accept.** Promoted E-318; H-015 -> REFUTED in the step reading (referee's wording: the library's F-only walk is not shown to reach it).
- **theorist** -- H-031 as three lemmas under (LF)+(H1)+(H3) on 13 formal models; M5 passes directly; M9 new. Referee experimentalist: minor revision (title count, M3 weight, missing printed checks, pivot coefficient, M9 overclaim); response did all items, M9 narrowed. **Decision: accept, narrowed** to the formal models, Lemma 3 relying on E-240 and a diamond lemma cited from memory. Promoted E-319; H-031 status line updated (still OPEN).
- **scholar** -- ENS math/0101145 proofs give a non-completed pair-cancellation lemma; BHL a second source for H-024; M4 correction. Referee theorist: minor revision (mis-citation of Cor. 3.11, letter counts vs "equals", scope of "only step left", plants, Q-invertibility, E-265 conflation); all seven done. **Decision: accept** at the narrowed scope. Promoted E-320. Literature file for ENS/BHL not yet written (owed).

## Consequences
- H-015's sufficiency reading is refuted at n = 11 by an R-move J != 0 step. The R-candidate "a key group is one derived class" is now broken at two groups (E-314, E-318). The `mutationSearchDepthFirst` docstring ("the guard keeps a walk in one derived class") is false for walks with R moves; it was already flagged unsupported. Whether any recorded F-only merge uses a J != 0 step is untested: thread T10 stays open with that replay as its next check.
- E-320 corrects round-121 `scholar.md` item 3 (M4 keeps 3 Theta letters), not E-315's recorded text.
- Ledger hint for conference 124: strongest R-candidate = "key group = one derived class" (E-314, E-318); F-candidate for H-031 would need (LF) on Construction-2.2 parents.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
(none from step 0). Chair: no overnight run proposed; the depth-6 ball stays unapproved (E-318 found the crossing with a depth-4 search).
