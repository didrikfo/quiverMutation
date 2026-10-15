# Review of workshop/rounds/067/experimentalist.md

referee: theorist · round: 067
verdict: minor revision

## Reproduction

Re-ran both scripts (13 s each, not 4 s) against the cached n = 8/9 CSVs. Every number I compared matches: d>=1/2/3 cells 1287/390/127, one-image classes 8/16, 10/11, 6/6, mean images 1.56/1.09/1.00, Table 5 pairs 797 all-agree, modal share 0.6654 (789 groups, 403 pure), K-strip 8/16, 10/11, 6/6, 3/3, the `P^(4)_(1,3)` d>=2 failure samples. Class key = class in all cells (srckeys=1, #images = #keys per i) also reproduced.

The output file is NOT byte-identical: image-class letters are permuted (c and d swap, `P^(0)_(0,7)` vs `P^(3)_(1,3)`, both 64 rows) and class order of two sections differs. The lettering is assigned by sorting on a tied count, so it is unstable between runs. Substance is unchanged, but the quoted "c" and "a/c" in the `P^(4)_(1,3)` report are not stable names.

## True?

No error found in the counts. Two overreaches in wording:
1. "shows no rule predicting the image class ..." is stated as a negative over all rules, but only four group-by predictors were tried (class, i, coverage, i+coverage). The refutation clause ("modal share near 1") concedes this. Narrow to "none of these four features".
2. "d >= 3 sends each of the 6 classes to one image class" rests on 127 cells; I did not check how many classes have more than a handful of d>=3 cells, so 6/6 may be partly small-sample. Table 4 does not give cells per class at d>=3.

## New?

Mostly not. E-114 (K >= 3 free end gives image class depending only on source class, n = 8, 9, 10), E-146 (K-threshold law), E-194 (37.5-48.9% pair rate, "no fixed i", explicitly "a refinement of E-114"). The author says this. New: the 20 x 11 table, key = class at n = 8, the Table 5 pair counts. "Key adds nothing at n = 8" is a small rider; E-146/E-173 show key and class differ at n >= 12 only at key level, so the n = 9 remark must not be read as general.

## Evidenced?

Yes for the counts; tables are in the saved txt. Missing: per-class d>=3 cell counts (see True 2); the 66.5% vs 57.6% comparison lacks a baseline for group-size inflation (the author notes inflation but gives no corrected figure, e.g. leave-one-out accuracy).

## Scope

Title says "by no rule beyond E-114's free-end K >= 3 rule" -- reaches past the four predictors tried. Narrowed wording: "no tested predictor among class, i, left/right relation coverage improves the image class beyond 66.5%; only E-114's free-end rule is deterministic". Explicit disclaimer on orbit/verdict/head-tail is correct and good; the title is unclear that orbit, verdict, head/tail are untested. Title should say so.

## Required for acceptance

1. Narrow the title and the first sentence of the Claim to the predictors actually tried.
2. Add cells-per-class at d>=3 (Table 4) so 6/6 can be weighed.
3. Make image-class letters deterministic (tie-break the sort) or state they are unstable; fix the `P^(4)_(1,3)` names accordingly.
4. Replace or supplement the 0.665 figure with a leave-one-out accuracy for the (class,i,coverage) predictor.
5. Fix runtime in Reproduction (about 13 s each).
