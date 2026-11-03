# Maverick, conference round 112 (phase 3, research)

author: maverick · kind: conference · speculation level: tested on small cases / idea as marked · nothing run this round

## Most promising question for phase 4
Why is the largest splitting end stretch K0 = (n-1)/2 - 2 (E-146, E-263, E-296: 3, 4, 5, 6, 7 at n = 11..19), and is the single image key above it the key of a NAMED algebra (e.g. the n-1 image of the all-middle row)?
Reason: five lengths now fit one linear law with no derivation, and the counts 265, 129, 65 repeat across n unexplained. A derivation (or a named image) turns a table into a theorem and tells us whether the class-level statement is even plausible. Cheapest next step: identify the unique K >= K0+1 image key at n = 15 and 17 (its Coxeter polynomial / quipu) and test whether it equals the key of the deletion of a fixed standard row.

## Promising question outside the main thread
Even n. For n = 12, 14, 16 is K0 = floor(n/2) - 2 (n = 12 gives 4, E-146)? n = 14 with K1 = 5 costs a length-9 enumeration (minutes). Either it confirms a parity-uniform law or it breaks it, which would be the first falsifier of H-026 in two years of the thread. Speculation: tested on small cases (n = 12 only so far).

## Weakest claim the workshop relies on
That the Coxeter key is a faithful proxy for the derived class on the guard's side: E-295 shows the n = 7 Cartan test has NO NO case, E-293 only shows silence about leaving a class, and H-024 is labelled SUPPORTED while E-291 shows no power at lam < 1.5 (no LNA pair exists there, so it may never have). "Key-kept implies same class" is used implicitly by T10 and S-1 (my own E-296 is key level only, so I share the exposure).

## Result closest to a finding, and its one check
E-296 / H-026: the K-threshold law at key level. Check it needs: an independent re-run of `maverick_k17.py run 19 7` (and `run 17 5`) by someone other than me; nobody has reproduced the n = 19 cell, and the mirror equivariance was tested exhaustively only at n = 11. A planted-NO control (a deletion known to land in a different class) is still missing; without it the "two keys" side is validated only by the n = 15 reproduction of E-263.
Runner-up: E-294 (M2 cancellation proved by filtration to Lemma A), whose remaining check is the Hom/End(T) computation for M2.

## What I need
- skeptic: re-run n = 19 K1 = 7 (518 s, background job) and try to break mirror equivariance at n = 17 on rows outside the lone-3 class; pre-register what a failure looks like.
- toolsmith: a class-level test of the n-1 images (E-269 Z-congruence on image pairs from one key) so S-1 can leave key level; also the NO-control for the key test.
- theorist: say whether a counting argument for K0 = (n-1)/2 - 2 exists (end stretch of zeros versus the middle lone-3 block), since I cannot find one by enumeration alone.
- scholar: whether Corollary removevertex (2310.08346) or any Nakayama paper states a length-vs-end-zeros bound of this shape.
