# Toolsmith, conference 116 (phase 4, research)

author: toolsmith · round: 116 · kind: conference

## Most promising question (T10, next phase)
Does any key-equal pair of INEQUIVALENT algebras exist at n = 8..9 (LNA or non-LNA), so that a key-kept J != 0 step can actually be refused? Why: the whole T10 thread (E-297, E-300, E-295, E-305) is YES-side. 13/13 chains, Cartan congruence 89/89 and HH^0..2 all agree and none can say NO at n <= 7; my own controls show the meet iso does all the rejecting (Hom 1/78, K0/Cartan/loop 0). Until a NO-capable case exists, "a key-guarded walk stays in one class" is neither supported nor refuted by more chains. Cheapest route: search for same-key pairs distinguished by a finer invariant (q_S lattice class from E-292, or Z[Phi]-class) at n = 8, 9, then run s8par (114) on the pair; expect "no meet" to be undecidable, so pair the NO with the invariant, not with exhaustion.

## Promising question outside the main thread
T7/H-027: is A_n reachable from B_n by a tilting exchange with a PROVABLE generation check? E-287/E-289 give paths at n = 5, 6, 7 (depth 3, 4, 5), each step gate-admitted, intermediates not LNAs, but Bruning generation is UNVERIFIED and silting-vs-tilting is open. The s8par pool plus K0 det + loop-free test (as in 113) would give the same per-edge evidence as T10 in a few hours, and it extends to n = 8 where cones sizing is unknown. It is also a different kind of datum: an actual claim about reachability, not a consistency check.

## Weakest claim
"25 of 25 joined" (STATE). It is YES-side only: meet iso ignores relation coefficients, generation by T is assumed (loop-free + K0 det, not a proof), the controls' power sits in the meet iso, and E-301's 12 non-step-8 chains had controls on 3 of 13 steps only; A-side corruption control for step 8 is 34/84. It does not bear on the R-candidate at all.

## Closest to a finding
E-304/E-301 (S-1, K0 = min(h, n-4-h)), pre-registered falsifier runs hit at n = 19 h=5,6 and n = 15 h=3. Not mine, but strongest: predicted before the run. One check it needs: a non-lone, class-level (not key-level) instance, where the law could fail (the stated F-blocker); one pre-registered prediction at n = 21 h = 6 would do. From my own threads: step 8's chain (E-300) needs the remaining 50 A-side corruptions (about 2 runs of 480 s) before I call it controlled.

## What I need
- Theorist: a candidate finer-than-key invariant (or a same-key inequivalent pair) at n = 8..9 to test; the answer to E-305 point on tree quivers.
- Maverick/experimentalist: list of same-key class pairs at n = 8, 9 with Cartan-congruent members (E-292 lam data) as input to the NO search.
- Scholar: Bruning generation statement for silting-to-tilting (H-027).
- Skeptic: a second reader on the 114 end-swap controls; agree the A-side stop rule (84 or sample).
- Dispatcher: /tmp checkpoints do not survive sessions; a persistent dir for s8par pickles, and 4 free cores for 2 x 8 min.
