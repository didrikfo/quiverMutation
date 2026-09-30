# Theorist notebook (rewritten round 004)

## What I believe now
- H-021' (pairing cores, s = n - k(c), d n-independent) stands. New: for the family 33x, k = 2x, hi = n-x-3, d = x-3, exact for
  n = 13..17 (x = 3..6, 20 fits); the d top offsets s+1..hi are singleton orbits, the rest exact reflection pairs.
- Interior/end-touch split now committed (rounds/004/theorist_shortfall.py): int 13/13, end0 7/10, endhi 10/11, allO 10/13, allI 45/62.
  It explains none of the 7 failures. Slide-consistency (s in cons) holds for 109/109: it is necessary, never selective.
- 4 end-touch failures: orbits merge one adjacent pair only, at the smaller consistent centre (4045 3556 4556) or the larger (4506).
  3 allO failures 4046 5046 5056: orbits {0,2},{1,3} look like a parity/translation-by-2, not a reflection. Fit has little power at |R| <= 4.
- Fitted centre is unique for most cores (theorist_k_n13.txt); ambiguous for 36 405 5004 5006 5046 5056 5066.
- Families: 3344/3355/3366 have d=-2 constant, 3444.. d=2 constant; only 33x has d moving with the digit.

## What I tried
- Column + slide-consistent centre set; k table for all fit cores at n=13; census for 33..336 at n=14..17 (20 s each, cheap).
- Did not try: x >= 7, other 3-digit families at n > 13, a null model for weak fits, the rule-table reason for k = 2x.

## Next
- Explain k = 2x from the rule table / T4 (why do the top x-3 offsets have no partner?). Check 34x, 44x, 45x for linear k in last digit at n=14..17.
- Ask skeptic for a null test of the fit at |R| <= 4; check 5046/5056 orbits at n=14 for a period-2 translation.
- Blind spot: "s in cons" is tautological given constant verdicts; do not sell it. Fits and descriptions (k=2x) are not mechanism.
