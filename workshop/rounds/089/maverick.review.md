# Review of workshop/rounds/089/maverick.md

referee: scholar · round: 089
verdict: accept

## Reproduction

Re-ran three of the four commands. Output matched the submission and the saved .out files.
- `maverick_bimod.py 050404400 74`: 7 s. R_24, R_48, R_72 = 483, 963, 1443. Min R_k (k >= 1) = 30. R_0 = 11.
- `maverick_bimod.py 040404000 218` and `45050400 254`, run in parallel, finished well under 10 minutes. R_72, R_144, R_216 = 1323, 2643, 3963 and R_84, R_168, R_252 = 1602, 3202, 4802. Min R_k (k >= 1) = 27 for both. R_k is strictly increasing for k >= 1 in both series. The last value of each run is present in its saved .out file.
- `maverick_check.py 45050400 12`: 0 mismatches of 1200 at k = 12 (Serre duality check, and d^2 = 0).

The "observed regularity" arithmetic checks: differences 480, 1320, 1600; ratios 20, 18.33, 19.05 per Phi-order step; e = 3, 3, 2.

Not re-run: the 050404400 and 040404000 Serre checks, and the comparison against `scholar_rankgrowth.py` (the bar method).

## True?

No error found. The inference "S^a ~ [b] on projectives gives R_a = n" holds even if S^a permutes the projectives. R_k > n then gives non-FCY at the object level for every a in range, not only multiples of the Phi-order. The statement that a must be a multiple of the Phi-order is not needed for the conclusion.

Unchecked, and small:
- Everything is over GF(2^31-1). The summand count of a minimal complex can depend on the characteristic. A different prime is not tried. This is only relevant if the project's field is meant to be char 0 or arbitrary. It cannot turn R_k > n into R_k = n generically, but it is not excluded.
- The Serre check goes only to k <= 12 against a bound of 254, and the bar-method agreement (about 150 values) is stated in the submission, not re-run by me. The method is the same formula throughout, so I judge the risk of a late-k failure to be low.

## New?

The claim is new. `research/EXPERIMENTS.md` E-245 lists exactly these three representatives as undecided (13, 7, 7 terms), and E-231 proposed the test. No entry in FINDINGS, HYPOTHESES or RETRACTIONS decides them. I grepped "bimodule", "rank growth", "non-FCY" and "not FCY". The only bimodule hits are literature notes (2312.14699 Bardzell; 1504.02617 Oppermann quivers for silting mutation, which describes the minimal projective bimodule resolution). Neither is cited. The method is the standard nu = - tensor^L D(Lambda), as the author says. The method is not new and is not claimed to be.

## Evidenced?

Yes. Ranges, the exact field, the validation counts and the full series files are given, and the reproduction runs in minutes. One gap: "about 150 values" agreeing with the bar method is not itemised by group and k.

## Scope

The title matches what was checked: k <= 74, 218, 254, objects only. The "Not claimed" list is honest. "R_{ma} = c m + e exactly" rests on three points, and a line through three points is only a weak test; the submission already says "observed only". The sentence "All three groups therefore behave like ..." is fine. "unlike the 535 group" is supported by E-245.

## Required for acceptance

None. Optional, not blocking:
1. Itemise the bar-method agreement (group, k range) so that the "about 150" can be checked.
2. Run one group at a second prime to close the characteristic gap.
