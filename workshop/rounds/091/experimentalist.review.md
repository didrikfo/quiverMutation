# Review of workshop/rounds/091/experimentalist.md

referee: toolsmith · round: 091
verdict: minor revision

## Reproduction

I could not re-run anything: the Bash tool returned "auto mode classifier gave no verdict" on every attempt (four tries, including a 6-deep rescan of `/tmp/e91/st8.pkl` into `/tmp/e91r`). So this review rests on reading code and the author's saved outputs, and the numbers are not independently reproduced. I read `workshop/rounds/091/experimentalist_jscan_out.txt`, the logs `/tmp/e91/{b2,build,C,R,A}.log`, `experimentalist_jscan.py` and `toolsmith_d6scan.py` (`jcands`, `base`).

Consistent with the claim:
- The per-depth node counts in the scan output (1, 8, 31, 90, 239, 632, 1690, 4570, 12949) are the differences of the cumulative ball sizes in the claim, and they sum to 20 210.
- J:moved 654 + J:kept 1 = 655, as stated.
- The ball logs match the Evidence table: C ball 1, 14, 100, 543, 2692; 3304000 ball to depth 5 is 4046; 2400230 ball to depth 6 is 2691.
- The scan covers every ball node. It is not restricted to a depth-8 front: the `d > MAXD` filter only cuts above 8.

Not reproduced: the meet results ("0 WL tests", "1 test, total 9, Hom-verified"). There is no meet output saved for any of the three meets, only the side logs.

## True?

Nothing wrong found. Points the author did not check or state:

1. **The `base` key is never shown to equal the 2400230 key.** `jcands` takes "kept" to mean the child key equals `base`, the key of `KEYROOT`, which defaults to 3060000. The report applies this to a 2400230 ball without saying that key(2400230) = key(3060000) = key(3304000). The result is consistent: C meets the 2400230 ball by isomorphism, so C has that key, and the 1 "kept" and 654 "moved" counts would be meaningless otherwise. But this is a silent default that a reader must reconstruct. State it, or give the one-line key comparison.
2. **"None at depth <= 5" is not in any saved output.** The scan prints only totals by node depth and by status, with no per-depth flagged counts. The claim and the comparison to 354 flagged of 2958 at 3304000 depth 5 rest on an unrecorded breakdown.
3. **The "7" in "1 + 5 of 7" in the Structural paragraph is unexplained.** It is E-239's 1 plus E-244's 6. E-244's scope says depth-6 parents only, so these are disjoint and the arithmetic holds, but the 7 appears nowhere in the text. Also, "1 + 5 of 7" leaves 8427 as the only recorded unresolved candidate, which is consistent with the text.
4. **The structural argument is partly vacuous.** A child that joins its own component by J = 0 cannot bridge. That is true, but it holds for the 2400230 child C too, as the measured part shows. The conclusion "no bridge" therefore does not depend on the scan at all for any J != 0 step whose child lies in its own J = 0 component, which is every one found so far (E-234, E-239, E-244, and now C). The only informative cases are 8427 and unscanned exits, as the author says.
5. **Weak control.** The only evidence that the test could detect a bridge is E-226's earlier 5 + 5 and 6 + 8. There is no positive control for the C(4) vs 3304000(5) meet, such as C(4) against a ball it is known to meet. The C(4) vs 2400230(6) meet does serve as one (it found 9 = 3 + 6), which should be said explicitly.

## New?

Nothing found for 2400230 as a J != 0 scan in `research/`. A grep for `2400230` gives FINDINGS.md:417 (SNF orbit table), E-238, E-237, E-226 and the F-047 mentions. E-234, E-239 and E-244 are the 3060000 and 3304000 sides, as cited. The claim is new as a scan of the 2400230 ball. Its content is a null partial result and adds no finding.

## Evidenced?

Mostly. The sizes, timings, depths and the single kept path are specific. Missing:
- the saved meet outputs (the three `skeptic_join.py meet` lines);
- the per-depth flagged-count breakdown behind "none at depth <= 5";
- a statement that `base` = the 2400230 class key (see True? 1);
- the C Hom-verification for the 9 = 3 + 6 meet, which is asserted but not logged (the `j0.pkl` child is available).

## Scope

The title says "does not bridge 2400230 to 3304000". The explored sets show only that C(4) and 3304000(5) do not meet, and neither do 2400230(6) and 3304000(5). Narrowed wording: "no bridged total <= 13 exists among the explored balls (C radius 4, 3304000 radius 5, 2400230 radius 6)". The text already says this in the body, and the title "(null, partial coverage)" is acceptable. The claim that "the move-set gap is not closed by it" should read "this one step does not close it within the sizes built". The "Next" section's hope that (c) will give a J = 0 baseline at 14 is a plan, not a claim.

## Required for acceptance

1. Add the sentence that key(2400230) = base key (3060000 default `KEYROOT`), with the command or log line that shows it. Doable in one sitting.
2. Save the three meet outputs and the per-depth flagged counts (`jcands` count by node depth) to a file next to `experimentalist_jscan_out.txt`, or drop "none at depth <= 5". Doable in one sitting (about 4 min plus 2 s each).
3. Explain the 7 in "1 + 5 of 7" (E-239's 1 plus E-244's 6, disjoint).
4. State that the C(4) vs 2400230(6) meet is the positive control for the meet tool.
5. Narrow the title and claim to the explored radii as above.
6. I did not reproduce the numbers. The author, or another referee with a working shell, should rerun the 49 s C-ball and the `meet` commands to confirm before the claim is cited elsewhere. `[next round]` if the shell is unavailable.
