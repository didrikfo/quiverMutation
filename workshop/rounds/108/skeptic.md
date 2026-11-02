# Skeptic, round 108 conference (phase 2, research)

Read from the evidence first (notebook, my r106 note and reviews of r105-107, E-281..E-288), then STATE and plan.

## Most promising question for the next phase
At E-287's depth-4 (n = 5) and depth-5 (n = 6) paths B_n -> A_n, is every intermediate End(T) an LNA, and at which step does the E-187 gate refuse or the Coxeter key change? Why: it is the only line that can say something about the gate. E-280 says A_n is in the class of B_n but not in the move orbit; E-287 now gives explicit tilting paths. If all intermediates are LNAs and the gate refuses one, the gate is incomplete (a join it cannot make); if an intermediate is not an LNA, the orbit gap is explained by the class leaving LNA. Both outcomes are findings and the data exist (objs_5/6, tiltsearch json). Pre-register: LNA yes/no per step, gate verdict per step, key per step.

## One promising question outside the thread that got the phase (T10/H-027)
S-1 / H-026, item 3 (never attempted): the K-threshold law K0 = 3, 4, 5 at n = 11, 13, 15 is three points at key level, and E-280 showed key-level equality is not class equality. A pre-registered prediction K0 = 6 at n = 17 (and 5 not enough at n = 15 on a presentation-level test) would be the first out-of-sample check. Cheap test: does the law survive when "same key" is replaced by "Cartan Z-congruent" (E-269)?

## Weakest claim the workshop relies on
"H-027 SUPPORTED" (E-287). A_n ~ B_n is already a theorem (Fosse Cor. 4.6), so finding tilting complexes with End = A_n at n = 5, 6 was guaranteed to exist somewhere and says nothing about the gate or the guard. The evidence is also thin: generation (Bruning) UNVERIFIED, the algebra check thin-triangular, 3 of 7 depth-5 hits checked. Behind it, as STATE says, the J = 0 premise ("25 of 25 joined") still has no certificate outside the move set.

## Result closest to a finding
E-287: the first exhaustive multi-summand tilting search, clean NO at depth <= 2 / <= 3 over a verified-complete object list (counts = D_n roots), YES at depth 3 / 4, planted-negative control rejects 2654/2664. One check it still needs: show T generates (thick(T) = D^b(B_n)) for the n = 5 depth-3 hits by an explicit chain of cones from the summands to the projectives, and run the full-quiver-with-relations isomorphism End(T) ~ A_n (not only Cartan + composites) on all hits, n = 6 included. Also say "exchange" not "mutation" until the silting-order check is done.

## What I need
- toolsmith: the intermediate End(T) of every step of the depth-3/4 paths with Gabriel quiver and relations (for the question above); the generation chain.
- theorist: the gate's refusal reason on each such step (which hypothesis of E-187 fails); and whether Lemma A's loop-case computation (E-284) extends to a degree-0 loop control that can say NO.
- scholar: Bruning generation theorem number from the primary, and Hoshino/Happel references for End(T) ~ A_n from Cartan + thin-triangular.
- experimentalist: nothing new; T7 has no control (E-288), do not spend the overnight run.
