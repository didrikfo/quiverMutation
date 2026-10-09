# Review of rounds/050/skeptic.md

referee: experimentalist · round: 050
verdict: minor revision

## Reproduction

Private dir /tmp/exp50. Re-ran: selftest (10 s); `skeptic_agree.py 6` (n=6 cross-tab: 252 gate-yes/tilting, 168 gate-no/J!=0, no off-diagonal cell: matches); both collectors (c1 366 s, c2 273 s: 16 and 9 failing steps, matches); `skeptic_replay.py` paths and parents, classes 1 and 2: edges 56+79 = 135 (child), 124+65 = 189 (parent), all `tilt`, `cart`, J=0; pathsOK 7+8 = 15, 16+9 = 25; BAD []. `skeptic_failsteps.py` on c1 and c2 (tail shown): failers have Hom(T,T[-1]) = 1 (tuple field), cartan-congruent True. Rebuilt paths agree with rounds/049/toolsmith_paths_logs.txt (spot check: key 31a15e, moves R2 R4 R7 R1 | F1 F5 F4 R7 F3 identical). Same numbers as claimed. Did not re-run the 200 control steps separately or n=5.

## True?

No counterexample found. Two small points:
- Selftest prints `Hom(A,A)==Cartan False True`: H equals Cartan^T, not Cartan. The submission says "equals the Cartan matrix" in the self-test sentence; the body correctly says Cartan^T. Wording fix only.
- The test is not run on anything outside the printed paths/walk (e.g. n = 8, class 3+); the claim does not say it was.

## New?

Nothing found for "Hom(T,T[" / "tilting complex" / Okuyama beyond E-126 (Hom(T,T[-1]) = sum J_i, theory) and E-223 / literature notes (rickard, 2509.12983 CHZ). No recorded computation of Hom(T,T[m]) on E-155 edges. New as a computation; the author states that it adds no information over `tiltingPlus` on 1081 tests except code/formulation independence. I agree: the novelty is independence, not a result.

## Evidenced?

Mostly yes: counts per set, cross-tab, reproduction commands with timings, limits (a)-(f) stated. Gaps: (1) the 25 failers rejected is the only power check at n=7; the cross-tab with "= gate" at n=5,6 is LNAs only. (2) No test that H = Cartan^T could catch a wrong End(T) (author admits Cartan-level only). (3) Generation not tested, so "genuine tilting step" in the title is slightly stronger than "Hom-vanishing + Cartan match".

## Scope

Title says "genuine tilting steps" and "premise holds". Narrowed wording: "all 324 printed edges at n = 7, classes 1 and 2, pass Hom(T,T[±1]) = 0 and the Cartan match of End(T) with the child, by an independent implementation; generation assumed (Okuyama-Rickard)". Title otherwise matches the sample (40 paths, 15+25).

## Required for acceptance

1. Retitle/reword "genuine tilting steps" to "pass an independent Hom-vanishing + Cartan test (generation assumed)".
2. Fix the selftest sentence (Cartan vs Cartan^T).
3. [next round] Quiver-level check of End(T) (already listed by the author), and the 5 misses / child 12.
