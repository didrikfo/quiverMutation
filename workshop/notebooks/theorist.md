# Theorist notebook (rewritten round 003)

## What I believe now
- H-021 without the mirror clause: "a core pairs at n" = orbits closed under o -> s-o (some nontrivial pair), s = n-k(c),
  shortfall d = hi - s, n-independent (12 cores 13->14 fine). Statement is about pairing cores only.
- d = tail - head (H-020) exactly when the slide has an outside block strictly inside (13/13 at n=13, 3/3 at 14): then s = first o + last o,
  the reflection folds the outside block onto itself. Near-tautology given orbit verdicts are constant; the content is that the fold happens.
- Exceptions (4045 3556 4506 4556; 4046 5046 5056 all-outside) are all cases with the block touching an end, where several centres
  are slide-consistent and the orbits pick one off by 1.
- All-inside slides (62 of 109 fit cores) are silent on d; d>0 there (33x: d = x-3) is not explained. That is where T2 is still open.

## What I tried
- Slides per orbit via batch._verdictFor(free=REDUCED): 395 verdicts at n=13, ~5 min on 4 procs (scripts theorist_slides.py etc.).
- Checked s == hi+h-t; then refined to first+last outside offset; 12 cores at n=14.
- Did not read research/ entries beyond H-021, H-020, F-053 headers.

## Next
- Test interior blocks with large m at n>=15 (skeptic request); n=15 shift check.
- Explain d for 33x / 4033 (all-inside): look at which orbit merges; maybe Coxeter key data (skeptic's scan n=24) gives k(c) directly.
- T4: proof that the block folds (rule table acts uniformly in the interior).
- Blind spot to watch: "first+last o" is almost a tautology; do not sell it as mechanism.
