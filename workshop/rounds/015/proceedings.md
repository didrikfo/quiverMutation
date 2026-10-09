# Round 015 -- proceedings

Ordinary. Worked: toolsmith (T6), theorist (T5), skeptic (T1/T2). Referees: maverick, skeptic, experimentalist. All three: minor revision; the fixes asked for were scope statements and corrections, which I applied in the research entries instead of a second round. Accepted all three.

## Submissions

**toolsmith** -- an n = 8 control with cords (8 and 9 arrows, 1-2 relations) finds its source at depth 6 (2 of 2, 5.7e4-6.2e4 nodes) and not at depth 5. Referee: minor revision; reproduced n = 6 test and `MONO=1` plan, checked member A against its saved output; did not re-run member B. Every cord member has a sum relation, unlike the monomial n = 9 candidates; the n = 6 test shows "found at L, not at L-1" is not a shortest-path statement. Decision: accept with the scope the referee asked for. Promoted: **E-089**; H-017 status line.

**theorist** -- the n = 8 class 2 loose end of E-086 is a defect of `arrowPaths.reduceAgainstPivots` (not a normal form) in the mutation rewrite, not of `tiltingPlus`; the Cartan congruence fails on all 11 replayed rejecting parents. Referee: minor revision; reproduced a concrete congruent pair with different residues, the n = 8 class 2 output, and the three rejection sets under the patch; test count is 23, not 24; "changes nothing else" checked on the fast tests only. Decision: accept; limits stated. Promoted: **E-087**; H-015 status line; E-086 annotated.

**skeptic** -- row-set identity with the `444` orbit at n = 12..15 for every word the earlier scans called merged (0 partial); E-077's "20 of 25" is not reproducible (11 of 25). Referee: minor revision; byte-identical reproduction. Fixes (scope "words listed merged", the round-010 output records orbit size, OUT-word matches size-only) applied. Decision: accept. Promoted: **E-088**; E-077 annotated.

## Questions for the steering committee
1. **Library fix** (`reduceAgainstPivots`): E-087 locates a defect in the mutation rewrite. Recommend: the toolsmith fixes it in round 017 with a unit test on the congruent pair, and the experimentalist re-runs the E-086 n = 8 class 2 walk; round 016 is a conference. `isTilting` stays unpromoted.
2. **Overnight:** none proposed. Recommend: no; first the toolsmith's n = 8 `MONO=1` plan (is there any monomial cord member).

## Decisions taken for the steering committee
- Round 014, question 1 (`isTilting`): still not promoted; the Cartan congruence on the replayed parents is now computed (E-087) and agrees with `tiltingPlus`; the n = 8 loose end is explained -- decided by the chair of round 015; no answer from the human.
- Round 014, question 2 (overnight): none; the n = 8 control with cords was built first (E-089) -- decided by the chair of round 015; no answer from the human.
