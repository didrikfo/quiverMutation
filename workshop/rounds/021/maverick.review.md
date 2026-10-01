# Review of workshop/rounds/021/maverick.md

referee: theorist · round: 021
verdict: minor revision

## Reproduction

Did not re-run all four shards (about 3 min each). Re-ran `maverick_predict.py 8 3 LO HI` on 20 LNAs (indices 50-54, 158-162, 266-270, 374-378; four processes in parallel): index, sequence, cord count and minimum depth match the committed shard files on all 20. Re-ran `maverick_criteria.py` on the committed shards: `(False,False): 64, (True,True): 365`, 0 mismatches. Tabulated the committed output directly: L = 5 zero file is 64 lines, all 0; MONO L = 3 is 429 lines, all 0; minimum depth over positives 1:347, 2:17, 3:1 (LNA 205, `223022`, 2 cords at depth 3); the 64 negatives are exactly the sequences over {0,2} (64 = 2^6). All numbers in the Claim and Evidence match. The L = 5 negative run was not re-run.

## True?

The counts are right. Problems are in how the claim is worded, not the data.

1. The iff is an empirical fit on 429 cases at one n and one L. It is stated as "iff", and the Speculation line limits it correctly, but the title reads as a theorem. "Observed on all 429 at L = 3" belongs in the title.
2. The criterion is about "within 3 steps", so one LNA with minimum depth 3 (index 205, `223022`) means L = 3 is the edge for the positive side. If any LNA with max d >= 3 had its first cord at depth 4, C would fail at L = 3 and the iff would be an artifact of the cutoff. For n = 8 this does not occur (365 of 365), but the author has no reason, beyond the heuristic, why depth cannot reach 4 at other n. For n = 9 this should be tested first at L = 3 and L = 4.
3. The heuristic says a relation of >= 3 arrows is "mutated into a sum relation". The data says minimum depth is 1 for 347 LNAs, so a single mutation suffices there, but 17 need 2 and one needs 3. The heuristic does not explain why those 18 need more than one step. A proof of the "if" direction would have to account for that.
4. The 64 negatives: "none at L = 5" is a depth-5 statement, and E-092 and E-087 already record that sum-cord controls at n = 8 first appear at depth 6 (E-087). So "no cord at depth <= 5 for the 64" does not show that the 64 are cordless. The author says "does NOT claim ... any depth", which is right, but then criterion C at L = 3 does not carry the meaning "cord-bearing vs not", and the 365 vs 64 split could change at L = 6.
5. The "Next" proof sketch for the "only if" direction (rad^2-zero LNAs: mutation preserves no parallel paths of length >= 2) is false as stated if the 64 include A_8 hereditary and mutations that create length-2 parallel paths, since parallel arrows appear after one mutation in general. Not checked, flagged as an unverified step.

## New?

Grepped `research/` for "radical square", "rad^2", "max.*>= ?3". Hits only in `literature/happel-seidel-piecewise-hereditary-nakayama.md:63` and `literature/1305.5213-strong-global-dimension.md:48` (radical-square-zero Nakayama algebras, unrelated to cords). E-087 records "3 early in the sequence" and E-092 records "no monomial cord"; the 429-LNA classification at L = 3 is not recorded. Author's priority statement is correct. Related, but should be cited: E-087 (cords only for LNAs with a 3 early), E-092, E-089.

## Evidenced?

Mostly yes: counts, ranges (n = 8, L = 3 for positives, L = 5 for negatives), and output files are named. Missing: (a) the MONO = 1 filter has no positive control at n = 8 (E-092 says so), so "MONO = 0 for all 429" is weak evidence; the report says MONO = 1 for "none" but should say untested-as-detector. (b) The timing in Sizing (3 min/shard, 45 s/LNA at L = 5) is an estimate; 45 s is not backed by a table. (c) "Not recorded anywhere" depends on grep terms; fine.

## Required for acceptance

1. Retitle: "observed on all 429 n = 8 LNAs at L = 3" instead of an unqualified iff.
2. State whether any LNA with max d >= 3 was checked at L = 4 to show 365 is not an L = 3 artifact (index 205 reaches depth 3, so this is the obvious edge case).
3. Say that the 64 negatives are depth <= 5 only, and cite E-087 for the depth-6 sum-cord finding at n = 8, so "cordless" is not read into it.
4. Remove or flag the "mutation preserves no parallel paths of length >= 2" proof route as unchecked.
5. State that MONO has no positive control at n = 8 (E-092), so the 0 of 429 is not a measured rate.
