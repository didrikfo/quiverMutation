# Theorist notebook (rewritten round 019)

## What I believe now
- Round 019 (T1/T2): the one in-S placement of a split four-letter word is fixed by the right gap g = n - end of last relation, g depends on the word
  only (n = 12..16): g = 0 for 2224, 2334, 4556/4667/4778/4889; g = 1 for 224x (x>=5), 344x (x>=4); 3344 is the left-gap-1 exception.
- R (applied to the end, each step checked against `doubleMutation.rewritesOf`) never reaches 333@0: it stops at a boundary shape
  (4, 4y, 334, 3 b b+2); the rest is the `34 <-> 403` shuttle (4@0 -> 3@1 slid along). 0 of 45 reduce by R alone.
- R is valid only when the shortened interval (a-1) does not swallow a left neighbour (3344 case). E-088's test was on isolated runs.
- Short shapes are in S only at the boundary: `4` at L0/g0, `4y` at L0/g1, `334 335 357` at g0, `44` everywhere (n = 12..16).
- Round 017: orbits labelled by J = offsets of 333@o' held (pairs {j, n-6-j}); `444`, `34` j = 0; `35 55 455 2455 3334` j = 1; `3x` j = x-4.
  k(33x) = 2x from the class table, E-065 upper bound not proved. Lemma R tested only by exhaustion.
- Nulls from 013 stand: no invariant separates classes; "class" = closed orbit under this move set.
- 015: reduceAgainstPivots defect (n = 8 class 2); fix in rounds/015/theorist_fix.py, toolsmith owns the library.

## What I tried
- 019: `theorist_{gaps,rchain,shapes,split,jlabel}.py`. Fastest: membership-only gap table, then R-closure terminals. My first R closure was wrong
  (accepted a non-neighbour); filtering against rewritesOf caught it.
- 017: label every placement by (orbit size, J); hand-apply the double mutation to intervals to read the lemma.

## Next
- Derive "4y in S only at L0 / g1" (the shuttle), then g of any word from its R-terminal; predict g for 5-letter words, test at n = 13, 15.
- Test words outside S (3334, 2455, 3335) for their terminal and g.
- Why is the shadow orbit (no 333) a parity class? Same question as P/Q.
- Blind spots: tables not proofs; letters <= 8 for shapes; I never ran n = 17; no inequivalence proof anywhere.
