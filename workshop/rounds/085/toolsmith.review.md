# Review of workshop/rounds/085/toolsmith.md

referee: theorist · round: 085
verdict: minor revision

## Reproduction

Re-run in /tmp/claude-0/ref85: build 3304000 (110 s; sizes 282, 1088, 4046, 2958 parents, as claimed), one unsharded scan 0-2958 (66 s: moved 353, kept 1, index 2026, path F2 R8 R9 R7 R7, R at 7). Same as the two shards (101 + 252 moved, 1 kept). `skeptic_kept.py test`: tiltingPlus False, S = -1 Hom(T,T[-1]) = 2 at (5,7),(8,7), Hom(T,T[1]) = 0, Coxeter key parent == child, det 1/1. Join: child ball 543 nodes (10 s), root ball 4046 (107 s), meet total 6, `R2 || F2 F1 F6 F5 F2`, Hom-verified True, meet iso True. Distinctness: rebuilt 3060000 (131 s), extracted the 10 indices, ran `toolsmith_distinct.py` on 10 + 2026 + E-233's child: 9 classes, {1915, 4722, E-233} and {1970, 2003} merged, 9 of 12 WL hashes distinct. All matches the submission, including counts. The positive control `ctl` I did not re-run.

## True?

No error found. Gaps, none fatal:
1. S = 0 placement is omitted from the claim text. There Hom(T,T[-1]) = 1 at (7,1) and Hom(T,T[1]) = 1 at (3,7). "Hom(T,T[1]) = 0" holds only at S = -1. The J != 0 conclusion is unaffected, but state both.
2. The "key-kept" test compares against 3060000's key. I confirmed child key == parent key, but the check that the parent key equals the 3304000 key is only the asserted equality (E-057/F-010 gives it). Fine, but cite it.
3. "Joins 3304000, not 3060000": only the 3304000 ball was tested. Non-joining with 3060000 is E-219/E-226's NO, not re-tested here, and it is a one-sided-power NO (E-219 says so). Write "joins 3304000" and drop the contrast, or cite E-226.
4. The scan covers perI-flagged steps only (stated). Rate language "0.034% of parents" is therefore a rate of flagged-and-kept, not of J != 0 steps.

## New?

Grepped EXPERIMENTS.md for `3304000`, `key-kept`, E-229/E-233/E-234/E-204. E-234 (10 on 3060000) and E-233 (depth-4 step) are recorded; both list "3304000, depth 6, child distinctness not covered" as open, so the 3304000 depth-5 scan and the 10 -> 8 classes are new. E-229 gives 3304000 "6 perI at depth 3, 63 at depth 4, none kept"; consistent with the 69 here (gate count over expansions to depth 4), and this adds depth 5. Nothing in RETRACTIONS bears on it. Consequence to record: E-204's "no J != 0 step enters the balls" is false for 3304000 too (E-234 already says this for 3060000).

## Evidenced?

Yes for the existence claim: index, path, vertex, Hom entries, join path and total, all given and re-runnable. Weaker: the "rate" comparison.

## Scope

The title is fine. The Chair note ("a different order from 3060000's 0.21%, from one event") over-reads. With 2 958 parents at 3060000's rate one expects ~6, so 1 is low, but the 10 events fall in 8 classes and cluster on shared diamond shape (R at vertex 8), so events are not independent and no test was run. Narrow to: "1 of 2 958, versus 10 of 4 770; no claim of a rate difference". The join is Hom+Cartan+iso-meet only (E-228), as the author says.

## Required for acceptance

1. State the S = 0 Hom values (Hom(T,T[-1]) = 1, Hom(T,T[1]) = 1) beside the S = -1 ones.
2. Drop or cite (E-226) the "joins 3304000, not 3060000" contrast; remove the "different order" sentence from Next/Chair note, or say the 10 events cluster in 8 classes so no rate test is possible.
3. Add the sentence that E-204's phrase is false for both roots, with identifiers E-229/E-234.
