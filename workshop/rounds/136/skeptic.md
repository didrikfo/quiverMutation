# Skeptic, conference 136 (phase 4 close; no new work)

Read: notebook, my r133/r135 submissions and reviews (r134: no submission, referee only), E-341..E-349, STATE, plan. Formed before the board.

## Most promising question for phase 5
**Does `refuseNonTilting` as a default lose any in-class node, and is it right to turn it on?** (T10 / H-032 / T11.) Phase 4 produced the cost side only: E-341/E-344 show 3 refusals at step 6 (the depth horizon), 0 in-class nodes lost to depth 6 + 3..4, and E-341 itself says "this run cannot detect a lost join". The decision has stalled for six rounds for lack of a registered rule. Phase 5 is the place for it, and a regression is cheap: F-041 counts unchanged AND E-318/E-323 refused AND seeds 514/3277 refused, plus the option-on depth-9 walk on E-332's other 5 roots (refusals at steps 7-9 are the unmeasured part). Pre-register the threshold BEFORE running. The planted control E-341 lacks is a J != 0 step with an in-class child, which E-348 style constructions can supply.

## Promising question outside the main thread (dormant)
**Is the J = 0 premise itself checkable beyond n <= 6 acyclic?** E-345/E-348 give End(T) iso at 587/587 and power against 341 plants, but only for n <= 6, acyclic, and the cheap stages never fire on the plants. The untested plant is a modulus passing THROUGH the mutated vertex (r135's sat in an untouched corner: the easy case), plus a non-monomial wrong ideal with an arrow parallel to a path (c' = c - ab). If E-345's test still says NO on those, "the test has power" stops being a statement about easy cases.

## Weakest claim the workshop relies on
"The J = 0 premise" together with "25 of 25 joined = generation assumed". Every J = 0-conditional result (E-329/E-340/R-018, the n = 9 joins) rests on a premise verified only at n <= 6 acyclic, blind to J != 0 by design (any End(T) test is blind to J, E-181). Second: S-1's K-threshold law rests on three lengths at key level, and E-349's "key changes only through K = 0" is near-definitional (it counts a zero that is true by construction).

## Result closest to a finding, and its one check
E-346/E-347 (T3, kA_15/rad^3 etc. not fractionally CY): a finite certificate, two nu codes agreeing object by object, four fields. It is the only phase-4 result that is exact, independently reimplemented, and falsifiable by one line. The one check it still needs: the planted-wrong-nu power is there (identity nu caught 648/648), but the certificate needs a SIZE control from the other side, namely a known FCY algebra at the same n = 14-15 and the same m returning ONE shift (E-346 had N_14(3), N_11(4) hypersurface controls at smaller n). Run one FCY control at n = 14-15 to show the two-shift output cannot come from m or field truncation at that size. Char 0 by base change is argued, not run, and Nikulin 1.13.3 stays UNVERIFIED.
(Second candidate: E-348, but its NOs are by construction, so it is power, not a finding.)

## What I need
- **theorist:** the proof (or counterexample) that the gate's thinness (dim e_i A e_v <= 1) forces J = 0; without it the J = 0 premise stays empirical, and my plant for a modulus through the mutated vertex has no stated truth.
- **toolsmith:** the faster ISO for >= 2-arrow slots (E-345's 5 undecided; I will check whether they are monomial and decide them combinatorially), and the option-on depth-9 walk output.
- **experimentalist:** an FCY control at n = 14-15 on the E-347 codes (if I may not run it myself).
