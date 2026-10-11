# Round 013 -- proceedings

Ordinary. Worked: theorist, toolsmith, skeptic. Referees: skeptic (on theorist), theorist (on toolsmith), experimentalist (on skeptic). All three reviews were "minor revision"; I applied the fixes in the promoted entries (below) and accepted all three. Agenda items 1, 2 and 4 of the round-012 agenda were worked.

## theorist -- staircase between the parity orbits (T1/T3)
Claim: `3a@o -> 3a@(o+2)` takes a - 1 moves (a = 5..9); `4@0 -> 3@1` is one table rule; `5046`/`5056` have two orbits of different sizes at n = 17; no invariant for the parity found. Referee: reproduced every path, extended to a = 10..12; did not finish the n = 17 run. Decision: **accept** (E-082). Fixes carried: a = 10..12 added, n = 12 stated for neighbour counts, the n = 17 sizes marked unreproduced with no saved output, and the "valid move" caveat kept. Not met: a saved n = 17 output file (re-run is a ~3.5 min job for round 014).

## toolsmith -- node counts for H-017 (T6)
Claim: n = 9 depth-6 negatives are 5e4-6e4-node walks (4 of 16 measured); an n = 7 depth-6 control finds 16 of 16, 0 of 4 at depth 5. Referee: reproduced the depth-4 counts, found the depth-7 factor mislabelled (5.2 per step, not 27x) and the cross-n "3-4 times" comparison to carry no weight on coverage. Decision: **accept** (E-083) with those corrections; the title claims only the node counts. The first control run still has no saved output (stated in Limits).

## skeptic -- E-077 counted by orbit (T2/T4)
Claim: the letter-4 contrast is one orbit per n (the `444` orbit) and that orbit holds both `34`-words and 4-no-`34` words, so the data cannot separate "letter 4" from "collapse to `34`". Referee: re-ran all four scans, identical; "refutes" to "sharpens"; E-077's "20 of 25" at n = 14 does not match the scan's 11 of 25. Decision: **accept** (E-081) as a sharpening. The 20-vs-11 discrepancy is stated in the entry and left unreconciled; E-077 itself is not edited.

Promoted: E-081, E-082, E-083; status lines of H-021 and H-017; glossary entry "Staircase".

## Questions for the steering committee
1. **Overnight: n = 17 key-coarser lists** (`toolsmith_orbitclass.py`, 2+ h). The theorist's run covers only `5046`/`5056`. Recommend: not yet; a skeptic or experimentalist first re-runs and saves the n = 17 `5046` output, and the theorist tries an invariant on the mutated-vertex multiset.
2. **H-017 depth 7 / n = 9 control with a non-hereditary member** (overnight, about 25-30 min per candidate, about 2-3 h for 4). Recommend: no; the toolsmith first sizes an n = 8 control with a non-hereditary source.

## Decisions taken for the steering committee
- Round 012, question 1: n = 17 lists not run overnight (see STEERING).
- Round 012, question 2: round-012 agenda approved unchanged (see STEERING).
