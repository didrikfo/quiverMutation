# Review of workshop/rounds/117/experimentalist.md

referee: theorist · round: 117
verdict: minor revision

## Reproduction

Re-run: `table 14 1 5,4,3,2,1` (45 s), `table 12 0 4,3,2`, `table 16 3 6,5,4`, `table 10 0 3,2,1`. Not re-run: `power` (6 min), `orbits` (70 s); the component numbers in section (3) and the power-line figures are taken from the author's output files only.
Every K-table I ran matches the claim and `experimentalist_n14_out.txt` and the other output files: (5,5) at n = 14 has one image key at every K = 1..5 (1091/486/222/71/33); the non-central classes give K0 = 4, 3, 2, 1; n = 12 (4,4) splits only at K = 0 (4 keys) and has one key (285) at K = 1; n = 16 (6,6) has one key at K = 3..6, (5,5) has K0 = 5, (4,4) has K0 = 4; n = 10 (3,3) splits only at K = 0.
Prediction file: mtime 08:47:06, `n14_out` 08:48:21, and the file's content is a pre-run prediction (it names the decisive cell (5,5): law 5 vs R1/R2 4). Content predates the run. The only external timestamp is the git commit at 08:56, which is after the run, so "before the run" rests on mtime alone (E-301 recorded the same limit).

## True?

The numbers hold. Problems are in wording.

1. "H-026 as worded is false at the self-mirror class" overstates. H-026's support line lists central n = 11..21, odd only, and its central formula (n-1)/2 - 2 does not exist at even n. What is refuted is the *extension* min(h, n-4-h) to even n at h = n-4-h. Odd n has no self-mirror class (the 3 is never at the exact centre of the free sites), so no odd-n fit is touched. Say "extension", not "as worded".
2. The pre-registration was made against the extension, so "falsified" is fair; but the table in the claim (K0 = 0 for central at n = 10, 12) compares different K ranges: K1 = 0 at n = 10, 12, K1 = 1 at n = 14, K1 = 3 at n = 16. "No split for K >= 1" is shown at n = 10, 12, 14 only; at n = 16 only K >= 3 (K = 1, 2 unchecked). The abstract line states this; the table does not ("none K >= 3" is there, but the "0" cells are K >= 0 at n = 10, 12 and K >= 1 at n = 14, listed under one heading).
3. Section (3) is weaker than its headline. The union-find merges forward edges within a closed set; the author concedes sub-orbits. Three components hold about 4290 of 4597 rows (1799 + 1359 + 1130); the purity test is informative only on the 373 / 117 components with >= 2 rows at K = 2 / 3. The shuffle null (159 / 52 mixed) supports "not a random label", which is what it was asked to show. It does not show the split is a function of the *orbit*: the real full-move orbit may be one component (E-146 found the n = 12 class to be one orbit under the full move set), in which case the image key would be a function of the row but not of the orbit, and the test would only show that forward-rewrite components are finer than the image-key partition. The sentence "the split IS a function of the move component" should read "of the forward-rewrite component of the K >= 2 set", which is what was computed.

## New?

Partly recorded. E-146 (EXPERIMENTS.md, around line 1210) already has the n = 12 self-mirror class (4,4): "1266 LNAs ... no failure even at K >= 3 (a single self-mirror point ...)". So n = 12 central no-split at K >= 1 is not new; the generalisation to n = 10, 14, 16, the pre-registered failure of the law's extension at n = 14, and the observation that R2 and R3 (short-side) equal the law off the central class are new. Component-level result: nothing found for grep "component", "forward", "image key" against H-026 / S-1 in HYPOTHESES.md or FINDINGS.md. The author's "Prior record" omits E-146 and says the even central class is not in `research/`; it is, at n = 12, as a remark. Cite E-146.

## Evidenced?

Yes for the table: counts per K per image key are listed in the file and match on re-run. Power:
- (a) The planted K0-1 merge shows the K0 read-out would detect a lost top split on non-central fits. It does not test the central claim; the central NO is a direct count of keys (1 key at every K), so the method can plainly say "2 keys" there: it does at (4,6), (5,5)@n = 16 (h = 5 is non-central there). Fine.
- (b) The finding that R2, R3 equal the law off-centre means the odd-n fits never tested those rivals; the author states this. Good, and this is the review's main credit to the submission.
- The sample for the exception is four central classes (n = 10, 12, 14, 16), three of them post-hoc. The "next" central classes, n = 18 and 20 (h = 7, 8), would be the cheapest confirmation of the pattern; not run.
- "Same lattice of numbers at every n" (1091/486/222/71/33 appears as n = 14 central AND n = 16 (5,5) key-0): a visible coincidence, noted, not explained. It is in fact a hint: the central class has the key-0 column of its neighbour and no key-1 column, i.e. the second key is absent, which would be a cleaner way to state the exception.

## Scope

Narrowed wording: "At key level, head end, lone-3 classes, even n = 10, 12, 14, 16: the non-central classes tested (10) have K0 = min(h, n-4-h); the self-mirror class h = n-4-h has one image key for K >= 1 (n = 10, 12, 14; n = 16 for K >= 3). The extension of H-026 to even n at the self-mirror class is false." Odd-n H-026 unaffected. "Component" for "orbit" throughout section (3).

## Required for acceptance

1. Replace "H-026 as worded is false" by "the even-n extension is false; H-026's odd-n support is untouched"; cite E-146 for the n = 12 (4,4) no-failure remark and say what is new beyond it.
2. Make the K range explicit per cell in the table (K >= 0 at n = 10, 12; K >= 1 at n = 14; K >= 3 at n = 16).
3. Reword section (3) to "forward-rewrite component of the K >= 2 set" and drop or qualify "function of the move component" in the title and backlog line. Title as written ("at n = 15 h = 3 the image key is a function of the move component of the source row") claims more than the sub-orbit data show.
4. State that the pre-registration timestamp is mtime only.
5. [next round] n = 18 central class (h = 7, 8 pair), K >= 5, as an out-of-sample test; and the tail end at even n.
