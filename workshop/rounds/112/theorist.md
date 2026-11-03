# Conference 112, theorist (phase 3, research)

Not run: reading only (notebook, 109-111 submissions and reviews, E-289..E-296, STATE, plan, STEERING).

## Most promising question for the next phase
Does the T5 loop-model cancellation ever bear on a library step? i.e. state ONE theorem: every dg parent that a gate-admitted J = 0 step can produce is "pivot-generator or Lemma-A shape" (M1 proved E-290, M2 proved E-294; M4, M5 are the same g-filtration, one sitting), AND decide whether the premise of H-023 is about models the library can reach at all. Why: E-211/E-279 say library steps never meet loops (cyclic D(i) is refused), so M1/M2 are proved but are not the hypothesis H-023 needs. The real question is the converse: for the finite-dimensional, non-loop parents the gate admits, is Lemma A's input (E-235 step 4, E-256 a) satisfied by a proof rather than by 25 joins? That is the plan's question 2, and it is still not closed.

## One question outside the phase's main thread
T10 / R-candidate: a Hom-verified chain for the 13 key-kept J != 0 steps of n = 7 class 1 (E-201, E-295). Cartan congruence cannot separate them (same key implies same genus), so the NO side has no case at n = 7. Cheapest route in my range: compute Hom(T,T[m]) for the one-step cone on 2 or 3 of the 13, as I did for M1 (planted-zero-map control included). If any is not tilting, the guard leaves a class (R-candidate dies); if all are, the guard has no counterexample at n = 7 with power.

## Weakest claim the workshop relies on
That "every gate-admitted J = 0 step is a derived equivalence" is supported by cancellation proofs (E-235, E-290, E-294). Those prove a quasi-isomorphism kM -> N (Prop. 7.6); the passage from that to a tilting complex with End(T) = A rests on Oppermann's Thm 7.5 and the E-202 sign set, both taken as printed (never re-derived), and on generation (Bruning) UNVERIFIED. Also E-269's float bound / the unsourced necessity step behind every Cartan-level NO (STATE).

## Closest to a finding
E-294 (M2 cancels, proved by the g-filtration reduction to Lemma A, machine-checked to weight 12 with a planted wrong-g control that is refused, star-ideal control nonzero in 37 of 66 cells). One check it needs: Lemma A (E-240) was tested for arrow degrees 0..2 only; the degree-0 loop has u in degree -1 (checked only to weight 6, 300 instances). Rerun E-240's check for degree -1 arrows to weight 12, or prove the Koszul-sign independence; this is the referee's item 1. The degree-0 claim stays "one parent, w(g) = 3, weight <= 12" until then.
Runner-up: E-296 (K-threshold (n-1)/2 - 2 at n = 15, 17, 19), but n = 19 rests on one run and there is no planted control.

## What I need
- skeptic: a pre-registered planted-error for the M4/M5 filtration and the degree -1 Lemma A rerun; agreement on whether M1/M2 count toward H-023 at all.
- toolsmith: a Construction-2.2 algebra with a loop at i (or a Hom/End(T) harness that accepts the cyclic D(i)); the 13 n = 7 key-kept steps exported as (parent, vertex) for a Hom check.
- scholar: where non-negative grading is used after Prop. 7.6 (to say what a degree-0 loop breaks); Anick primary.
