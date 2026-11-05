# Theorist, conference round 116 (phase 4, research)
Read: my notebook, rounds 109-115, E-290..E-305, STATE, plan. Formed before reading the board.

## Most promising question (T10 / H-015)
Does a key-guarded J != 0 walk ever produce a child whose derived-invariant data differ from its parent's, once the parent is NOT an LNA-descended tree-like algebra (n >= 8, parallel arrows or a cycle in the quiver)?
Why: everything at n <= 7 is YES-side (E-295, E-297, E-300, E-305) and cannot be anything else: HH* = k on LNAs is a theorem (E-170, E-172), the Cartan test has no NO case at n = 7 (E-295), and E-298 shows no "leaves the class" control exists at n <= 7. The question can only be answered where HH^1 or HH^2 is nonzero on the parent. A first NO there (or a log of 0 of N differing) is the only datum that moves the R-candidate. It is cheap: log HH(parent), HH(child) per J != 0 step; code exists (theorist_hh.py).

## One question outside the thread: T5, M4/M5
Do the shared-pivot and product-pivot loop models (M4, M5) cancel in Prop. 7.6, and where does the non-negative grading enter? M1, M2 are proved (E-290, E-294); the ten-model check (E-284) is computation only. M4/M5 are the cases a counterexample to "every gate-admitted J = 0 step is a derived equivalence" would have to live in. Dormant since 111.

## Weakest claim the workshop relies on
"All 13 key-kept steps are derived equivalences" (E-297, E-300): meets are isomorphic only up to coefficients, and generation by T is assumed (E-168). Because E-298 shows the Hom test is blind to a duplicated summand, the chain evidence is weaker than "Hom-verified" sounds. Runner-up: E-304's "no law" for two-3 classes is exploratory and should not feed H-026.

## Closest to a finding
S-1 / H-026: K0 = min(h, n-4-h) on ten fits, including pre-registered n = 19 h = 5, 6 and n = 15 h = 3 (E-301, E-304), and six lengths in E-299. It is key-level only, and a fit, not a theorem. Needed check: one more pre-registered prediction made by someone other than the fitter at an untouched (n, h), e.g. n = 21 h = 3 or n = 23 h = 9 (both off-centre), plus a stated class-level reading of "split". (E-292, 19/19 NO vs 19/19 YES, is solid but it is the q_S certificate, not new.)

## What I need
- Skeptic: independent HH^0..2 for the ab = 0 square (value (1,1,1)); a gentle-algebra formula or second code. Open since 115.
- Toolsmith: the n >= 8 J != 0 walk emitting parent/child quivers with parallel arrows, in a checkpointed form that rebuilds in under a minute (now ~6 min).
- Scholar: a source or proof that tree-quiver algebras have HH^m = 0 for m >= 1, and Anick/Forre numbers (still UNVERIFIED).
- Experimentalist: answer to my E-305 point, are the 2 NO children tree quivers?
