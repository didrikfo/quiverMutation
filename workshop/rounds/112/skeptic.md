# Skeptic, conference 112 (phase 3, research)

Read: notebook, r110 submission and review, E-289..E-296, STATE, plan, STEERING S-1. I have no r109/r111 submissions (I worked r110, refereed theorist in r111). I did not read the board first. Nothing was run.

## Most promising question
Is any of the 13 key-kept J != 0, D = 0 n = 7 class-1 steps (E-201) NOT a tilting step, by a Hom test (round-050 chain-map `homdim`: Hom(T,T[m]) = 0 for m != 0, K0 unimodular) that does not use the Cartan matrix or the key?
Why: it is the one blocking check for the R-candidate "a key-guarded walk stays in one class". E-293 shows the guard refuses 1 150 of 1 150 J != 0 steps at n = 6, so the guard's NO side works. E-295 shows the Cartan route has no NO case at n = 7 (89 of 89 congruent, expected from same key), so it cannot decide. Only a test independent of the key can. Either answer is a result: a failing step is a counterexample to the guard, 13 passes is the first real evidence it holds.

## Promising question outside the main thread
T6 / H-017: a statistic with power for "defect per key class". E-283 found defect NOT constant at n = 9 depth 5; it is constant only on the defect <= 1 filter, which I suspect is a selection effect. It has had no work since r105. Question: does defect separate classes at all, or just track depth/length? Cheap control: shuffle the class labels within equal (n, depth) and see whether the "constant on the filter" result survives.

## Weakest claim
The J = 0 premise (every gate-admitted J = 0 step is a derived equivalence), carrying "25 of 25 joined" and the guard's whole story. E-294 proves M2, but M2 is a formal dg algebra, not a Construction-2.2 model of a finite-dimensional algebra (kM has H^1 != 0). M4/M5 and Hom/End(T) for M2 are open, and the degree-0 case is one parent at weight 3. "Proved in models" is not "proved for the steps the search takes". Runner-up: H-024's necessity step is sourced only inside the proof of Ladkani Lemma 3.14 and Nikulin 1.13.3 is secondary; E-292 now softens it (q_S gives the NO independently on 19 of 19).

## Closest to a finding
E-296: the K-threshold law, largest splitting K = (n-1)/2 - 2, at n = 11, 13, 15, 17, 19 (5, 6, 7 at 15, 17, 19), key level, lone-3 class. Check it needs: re-run of the cells nobody reproduced (`run 19 7` and `run 17 5`, 518 s), plus a pre-registered prediction for n = 21 (K0 = 8; K = 8 two keys, K = 9 one) fixed BEFORE the run. Without that it is a fit to data already seen. Note the tail-end mirror shortcut was tested only on n = 11, random n = 13 and n = 17, not n = 19.

## What I need
- toolsmith: the Hom-verified tilting test on the 13 E-201 steps (and the planted-wrong control: 4 of 6 wrong cone placements failing, as in E-290).
- maverick: someone other than the author to re-run n = 19 K = 7, and the n = 21 run with prediction logged first.
- scholar: still the primary for "lists complete to degree 44" and Nikulin 1.13.3 (open since r102).
