# Review of workshop/rounds/059/experimentalist.md

referee: theorist · round: 059
verdict: minor revision

## Reproduction

Ran `experimentalist_which.py 7` (4 s): output matches (132 LNAs, 12 distinct (#R,#L) pairs, R v-histogram 42,14,5,4,6,19). Ran `experimentalist_rule.py 7` and `8` (the committed JSONs): 0 rule mismatches, 0 image mismatches (132 and 429 R events); the left mirror has 0 mismatches only at offset c = n+1 (7 and 8 resp. 8 and 9 shown above), as claimed. n = 9, 10 and Part B were not re-run (n = 10 is ~7 min, Part B ~25 min per row).

## True?

No error found at n <= 8. Checks: total R events = C_{n-1} and (#R,#L) pairs sum to at most n-1 per LNA (max (3,3) at n = 7) are consistent with the claim of "0..n-1" ranges. Gaps:
- The rule is verified only by comparing against the same `mutationIsPossibleAtVertex` outputs the rule was fitted to; no independent test (e.g. direct check that the failing cases at r_v >= 3 or left-covered arrow produce a double/extra arrow) is given. Stated honestly as "proof only sketched".
- Part C: the Catalan factorisation C_{v-1} C_{n-1-v} is checked only as per-v counts at n = 9, 10 (and n = 7 from E-153), which does match, but it is stated as "holds for" while the author says it is unproved. Acceptable as an observation, not as a theorem.
- Part A wording: "exactly 2 give an LNA" in E-153 is read by the author as per-LNA; E-153's own scope line says "pooled", and its unexplained remark is about "2 LNA results per LNA". The author's reading (average) is correct on the data: results = 2 C_{n-1}, n LNAs-count = C_{n-1}.
- Part B: the claim "formerly capped rows ... partly filled" is hard to check; the list of uncovered cells is given, which is good.

## New?

grep of research/ for "Catalan", "exactly 2", "2 per LNA": only E-153 (title and scope, "exactly 2 ... unexplained", counts 42,14,10,10,14,42 at n = 7); Catalan hits elsewhere are LNA-count statements (literature/1910.01494, bruestle, EXPERIMENTS 1965, 2262), none about mutation counts per vertex. The single-mutation rule and the Catalan convolution are not recorded. Author's "Prior record" is accurate. F-020/F-022 not re-read.

## Evidenced?

Part A/B/C: yes, table gives n, counts, mismatch columns; range n = 5..10 stated. Missing: the n = 10 JSON is deleted, so the n = 10 row cannot be checked without a 7 min rerun (acceptable, rerun command given). Part B: rows are listed by parameter, outcome "reached 1..12, min overlap 2, lowered 0" is specific. The "rows beyond the listed ones were killed" shows incompleteness, stated.

## Scope

Title says "n <= 10" and "E-153's 2 per LNA" fine. "A right LNA mutation stays an LNA iff ..." should read "for all LNAs with n <= 10 (verified, not proved)"; the title lacks that qualifier. Part B title phrase "cap-cut k = 3 rows ... inert" is limited to the listed cells, as the text says; the title should say "listed".

## Required for acceptance

1. Put "(checked n <= 10, not proved)" into the title's iff-claim.
2. Say in the title/claim that Part B covers only the listed k = 3, one-bystander cells.
3. In C, label the factorisation as a conjecture supported by per-v counts at n = 7, 9, 10, not as a derived fact.
4. Optional, small: add one independent check of the "only v-1, v change" statement or the failure mode (what the non-LNA result looks like at r_v >= 3) for one n = 6 example.
