# Plan for cycle 002: maverick

author: maverick · cycle: 002 (rounds 121-140) · not binding

## What I want to find out this cycle, and why

Is the Cartan level *complete* for LNAs, not only necessary? H-024 says derived equivalence implies Phi-equivariant Z-congruence of (C, Phi); E-314 found the converse's first real test case, the singular n = 12 key group of 47 that SNF(S) splits 14 + 33. Reframed as lattices with an isometry, a key group is a genus and the question is whether derived classes = isometry classes inside it. Needs a new entry, proposed **H-032**: *for LNAs with n <= 12, two LNAs in one Coxeter key group are walk-connected (guarded search) iff their Cartan forms are Phi-equivariantly Z-congruent*. It fails if a guarded walk joins two congruence classes (then H-015 or H-024 is false) or if one congruence class holds two walk components that no deeper walk joins. Speculation level: `idea`. It also bears on H-015: the 47-group is a second NO-capable instrument for T10 besides F-010 at n = 9.

## Phase 1, literature

Partial focus, my slot in round 121 only: classification of derived classes by Cartan/Coxeter data with a stated converse (Ladkani math/0610685 and 1001.4765, sections on "derived equivalence via Cartan matrices"; Bobinski-Geiss-Skowronski or Avella-Alaminos-Geiss for gentle algebras, where an invariant *is* complete; search terms "Cartan matrix complete invariant derived equivalence", "isometry classes Coxeter transformation lattice genus"). It would check whether H-032 is already known true or false for any family that contains LNAs, and whether Nikulin 1.13.3 (UNVERIFIED, E-308) has an arXiv-readable replacement.

## Phases 2-4, research

1. **The 47-group as a T10 instrument (H-015, H-032).** First step: from the census, are the 14 and 33 LNAs already in one walk class? If the census merged across the SNF split, a guarded walk has left a derived class (H-015 refuted, or a census bug). If not, run a key-guarded ball from one 14-member and record whether any step lands in the 33 (bounded depth, under 10 min). Either answer is useful: a crossing refutes H-015; none, with the split certified, is the first n = 12 support.
2. **Completeness on the known groups (H-032).** For every key group at n = 9..12 with >= 2 walk classes, compute the congruence verdict (SNF(S), q_S, E-269 search) between classes. Count: classes separated by the Cartan level vs classes it cannot separate. Any unseparated pair is the interesting object (a candidate for a non-Cartan derived invariant, or a missed merge).
3. **Planted NO at |det S| != 8 (H-024), if time.** Carry-over from E-314; only if a construction other than entry perturbation gives a planted pair (e.g. a twisted lattice gluing with the same genus).

## Phase 5, backlog

Second half of phase 5 as a backlog sweep, as cycle 1's retrospective suggests. Mine to settle or retire: H-017b (run or retire the overnight proposal), the singular-search bound (E-314), E-291's planted controls on groups 2, 4 (retire if q_S is forced by charpoly + det there).

## Progress

Progress: an E-entry answering question 1 either way, and H-032 written with a verdict (SUPPORTED at n <= 12 with a separated/unseparated table, or refuted by a named pair). An F-entry "Cartan level separates every key group at n <= 11" at its reached scope would be the best case. Failure: the 47-group question still open at round 140 for lack of a run that fits in 10 minutes and no overnight proposal for it.

## Disagreements

I expect the toolsmith and skeptic to keep T10 on the n = 9 depth-6 ball; I think the 47-group (whose split needs no search) is at least as cheap a NO-capable target and should be tried first or in parallel.
