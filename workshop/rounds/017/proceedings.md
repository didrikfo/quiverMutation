# Round 017 -- proceedings (ordinary)

Worked: toolsmith (T5/T6), experimentalist (T5), theorist (T1/T2/T4). Referees: skeptic (toolsmith, theorist), scholar (experimentalist). Step 0.5: no human answer to round 016's questions; decided below.

## toolsmith -- `reduceAgainstPivots` fix and `MONO=1` at n = 8
Claim: the library function is now a normal form (test fails on the old code); `MONO=1`, n = 8, L = 6 gives 0 monomial cord members for six walked LNAs. Referee: minor revision (reproduced; patch is a real normal form; touched tests 42 passed, other six files 552 passed). Decision: **accept**, with the referee's scope fixes written into the entry: "all" reduced to the six walked LNAs, no positive `MONO` control, canonicity assuming the pivots span the ideal. The three points are not required for the library fix, and the negative is stated as scoped. Promoted: **E-089**; H-017 status line.

## experimentalist -- E-084 under the fix
Claim: every re-runnable E-084 count is unchanged under the patch; the n = 8 class 2 walk is too short to reach the 10 key-moved steps. Referee: minor revision (counts reproduce, saved files diff clean; the "unpatched control" actually ran the patched library, since the toolsmith patched it in parallel, and "3.2e5" is about 2.2e5). Decision: **accept** with those rows struck. Promoted: **E-090**; E-085 annotated. The depth-8 completion is not an overnight run (needs a checkpoint first).

## theorist -- `3334`, `2455` and `35`
Claim: lemma R, class label `J`; `3334`, `2455` are in the `333@1` class, never the `444` orbit. Referee: accept (minor wording; reproduced, plus n = 17: `3334` size 1191, `J = {1,10}`; `444` 11340, `J = {0,11}`). Decision: **accept**; "proof" in the kind line corrected to "exhaustively checked over a bounded range". Promoted: **E-088**; H-021 status line; glossary: class label `J`, lemma R.

## Questions for the steering committee
1. **Overnight:** none proposed. Recommend: no. The depth-8 n = 8 class 2 walk (15 min) wants a checkpoint first (toolsmith).
2. **Agenda:** approve or change the round-016 agenda in `STEERING.md`. Recommend: approve as is (items 1-3 now each worked once).

## Decisions taken for the steering committee
- Round 016, question 1 (overnight): none -- decided by the chair of round 017; no answer from the human.
- Round 016, question 2 (agenda): approve the round-016 agenda unchanged -- decided by the chair of round 017; no answer from the human.
