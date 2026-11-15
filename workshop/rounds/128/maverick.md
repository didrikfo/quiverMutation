# Maverick, round 128 (phase 2 conference)

Read: notebook; rounds 125-127 (mine and reviews); E-323, E-324, E-328, E-330, E-331 and STATE. Not run: anything.

## Most promising question
Is the repo's F step exactly the irreducible silting mutation of End(T) at the vertex, and so is a J = 0 pure-F walk a descending chain in the silting poset of the start LNA?
Why: it turns T10 from a search whose every NO is capped by depth (E-331, 129 813 steps, zero hits, and a deeper ball does not fit in one command) into an interval question that can be settled by argument. E-328 already certifies the any-slot irreducible-mutation property on every step of the n = 5-7 A_n/B_n paths, so the claim has a first test on real steps. What is missing is the general step and the direction. If F is a left or a right mutation only, "every crossing needs an R move" becomes a structural statement, not a count. Also it explains why 3 of 397 J = 0 F edges have no R inverse (E-327 shows the gate and the tilting test disagree, so the library F step need not be the mutation).

## One promising question outside the thread
H-025 / backlog item 8: search for a periodic-not-FCY Coxeter transformation Phi on LNAs, the converse of "Phi periodic => FCY". T3/T8 has sat at F-055 for several rounds with nothing run on the converse, and a counterexample would break the Cartan-level necessity that every NO rests on. It is cheap: Phi is computable at n <= 12 from the existing key tooling, and we only need one LNA with finite-order Phi and infinite global dimension or non-FCY.

## Weakest claim
Every Cartan-level NO (E-292, E-324, the X_a separations) rests on one unsourced necessity step: Nikulin 1.13.3 is still secondary-only, and Ladkani 3.5 + Rickard + Aihara-Iyama gives necessity only for the form P = r^-T as stated. For the single n = 11 pair, all Cartan witnesses use the same H-024 step, so they are one witness and not several. Second weakest: "no pure-F route" (E-331 plus R-only balls). The R-only part is a lower bound, because 3 of 397 J = 0 F edges have no J = 0 R inverse, so backward-by-R is not a complete reverse search. The review is right that the claim must be narrowed.

## Closest to a finding
E-331: zero key-kept J != 0 steps in 129 813 gate-admitted F-only steps at n = 12, depth 5, 20 roots (one per SNF class of the 10 split groups). The one check it needs: a power control. The same detector must fire on a known F-only J != 0 step in the same pipeline, in practice the n = 7 parents at depth 7-8 of E-151/E-204, or a planted step. Until it fires, the zero is "not found", and the roots are one per SNF class, not per orbit.

## What I need
- theorist: the direction and exactness of the F step against Aihara-Iyama (mutation at one summand, left or right), and whether J = 0 forces a descending chain. One page is enough.
- toolsmith: speed for depth 6-7 at n = 12 (q_S cache) and the exhaustive edge-inversion count the review asked for (J = 0 F edges with and without a J = 0 R inverse, n = 11 and 12), with the three non-invertible edges recorded.
- skeptic: ruling on whether forward F-only closure from the parents' SNF class with iso-tests is the acceptance bar for my revision.

Most promising question: is the repo's F step exactly the irreducible silting mutation, so that J = 0 pure-F reachability is an interval question in the silting poset?
