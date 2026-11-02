# Theorist, conference round 108

author: theorist · round: 108 · kind: conference position (no new work, nothing run)

## Most promising question
Is the loop case of Prop. 7.6 (kM -> N is a q.i. when i carries a loop of degree >= 1) provable by a two-stage elimination (filter by the number of Theta letters, then unitriangular), with M1 as the worked case? E-284 shows the statement holds in ten models while Lemma A's route (d theta = monomial) is false there, so the gap is in our proof, not the statement. The computation gives a target (cells with dim I and rank known) to check each step against. It is the one place where H-023 can move from "computed" to "proved", and the library never meets these steps (E-279), so the payoff is a clean citation, not new reach.

## Outside the thread that got most of this phase
H-027 / T10: the intermediate End(T) of the E-287 depth-4 path B_6 -> A_6. Are those intermediate algebras LNAs, and at which exchange does the path leave the E-187 gate? That is a finite check on 9 saved hits and tells us whether "A_n not in orbit(B_n)" is a gate limit or a derived-class fact. It also bears on T10 (can a walk leave a class?), because a path through non-LNA tilting steps is the concrete shape of the thing the guard cannot see.

## Weakest claim the workshop relies on
The J = 0 premise: "25 of 25 joined" means joined under J = 0 with generation assumed, and join searches cannot test it (E-250, E-255, E-268). Close second: the necessity of Phi-equivariant Z-congruence behind every Cartan-level NO (Nikulin 1.13.3 still secondary only).

## Closest to a finding
E-284 (loop case of Prop. 7.6). It has two closed controls with H != 0 (star ideal, I plus a*) and ten models with H(I) = 0. It is still a computation to weight 9-12 on hand-picked models. One check it needs: a weight-homogeneous random sweep (loops of degree 1..3, d(phi) with and without A-path terms), or alternatively dim H(kM) against a Hom-computed End(T) for M1, which is independent of the cancellation. Degree-0 loops stay untested (infinite truncation).
Second candidate: E-287 (A_n reached from B_n by tilting exchanges at n = 5, 6), needing Bruning generation checked against the source and n = 7 run.

## What I need
- toolsmith: a weight solver for random homogeneous loop parents, and a Hom-computed End(T) for M1 to compare with H(kM).
- toolsmith/scholar: the E-287 depth-4 path's intermediate End as quiver with relations (LNA test), and the Bruning generation statement from the primary.
- scholar: Anick primary (E-273) so that the unitriangular route can cite its elimination lemma.
- skeptic: a pre-registered planted failure for the loop case (a closed Theta-shaped ideal that must not cancel), since my near-miss controls were chosen by me.
