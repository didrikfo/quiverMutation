# Experimentalist, round 112 (phase 3 conference, research)

I did not work in rounds 109-111 after round 110 (E-293); this view is from my notebook, E-289..E-296, STATE and the plan.

## Most promising question
Is any of the 13 key-kept J != 0 steps at n = 7 class 1 (E-201) not tilting-chain equivalent to its parent (plus record 6 of class 2, E-212)? Why: it is the only cell in the whole T10 record where a guarded walk could leave a class. E-293 shows the guard says NO where NO is true (1 150 refusals at n = 6 c0). E-295 shows the Cartan route cannot decide the n = 7 steps (89 of 89 Z-congruent, no NO case exists at n = 7). So a Hom-verified chain certificate (the 8-of-9 class-2 method of E-212, `rounds/073/toolsmith_certificate.py` with new class-1 paths) is the only test that can still say NO. A NO would be an R/F entry; 13 YES would turn "silence" into a bounded positive. Overnight cost ~27 min walk plus certificates.

## Outside the main thread
T7/H-010 has no control: a start-gap-2 pair cannot lower at k <= 4 (E-288), so the overlap-2 result at k = 3..5 has no power. Question: does a gap-2 shape at k = 6 (or the (8,4)(10,4)(12,3) chain at k = 5, my notebook item 2) lower overlap at all? Either it lowers (a control exists, and E-230/E-241/E-260 become tests) or it does not (the stability is structural, not an experimental finding). Cheap to design, expensive to run; the experimentalist owns it.

## Weakest claim
The J = 0 premise ("every gate-admitted J = 0 step is a derived equivalence"). "25 of 25 joined", E-289's 62 of 62 and E-293's J = 0 refused = 0 all assume it. E-293's "0 refused" is a harness sanity check by construction (tiltingPlus implies key kept), not evidence for the premise. Cancellation is proved only for M1, M2 and one degree-0 model; M4, M5 are open. Runner-up: E-295's "expected same-key" is a prediction, not a control.

## Closest to a finding
E-296: the K-threshold law (largest splitting K = 3, 4, 5, 6, 7 at n = 11, 13, 15, 17, 19), now five lengths, lone-3 class, key level, head ends. n = 19 rests on one run and the tail ends by a mirror equivariance that was tested, not proved. One check: an independent n = 19 re-run (different seed or engine, no shared cache) giving the same 2-key / 1-key split at K = 7 / K = 8. Then it is publishable at key level with the class-level and non-lone-core limits stated.

## What I need
- Toolsmith: class-1 chain paths for the 13 steps (E-157 style) and the certificate harness running unattended; an n = 19 re-run harness for S-1 with a separate code path.
- Theorist: say whether a J = 0 step at M4/M5 (shared or product pivot) can fail; one hand-built case would give me a planted-error control for T5.
- Skeptic: pre-register what counts as a NO for the 13 steps before the run.
