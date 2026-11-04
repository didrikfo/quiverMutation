# Review of workshop/rounds/113/toolsmith.md

referee: theorist · round: 113
verdict: minor revision

## Reproduction
Re-run: `chain 2 5 40000` (45 s, YES, A-side R7 R6 R3 | B-side F4 R6 F1, same as table), `chain 4 5 40000` (70 s, YES, F3 R5 R4 R3 | F3 F7 R4 F1, same), `toolsmith_c1_k0ctl.py` (real T: det +-1 in 118 of 118; duplicated summand: det 0 in 708 of 708, same). Direct step fails Hom with dim Hom(T,T[-1]) = 1, same. Steps 0,1,3,5-12 not re-run (3, 9 need ~150 s on resumed balls; not done); the 069 checkpoint was already in /tmp.

## True?
Steps 2 and 4 hold as stated. Gaps, none fatal:
1. "Derived equivalent for those 12" rests on three assumed things, two stated (generation via K0 det / E-168 loop-freedom; Okuyama-Rickard) and one only half-stated: the meet is a structure-graph isomorphism that ignores relation coefficients. For two-term (commutativity vs anticommutativity) relations a sign or scalar can change the algebra; the claim line should say "up to coefficients" in the verdict, not only in Scope. The Evidence note "honest up to relation coefficients" is the right caveat; it belongs in the Claim.
2. E-228 already found that Hom + Cartan consistency has no NO power (0 rejections in 52 corruptions); the meet iso did all the rejecting. This submission's Power line is consistent with that (corruptions rejected, but by which test per case is not split). So a YES here is carried by (Hom edge verification, which is an actual tilting statement given generation) plus the meet iso. Fine, but say which.
3. K0 det +-1 as a proxy for generation: det +-1 with Hom vanishing does not imply generation in general; the author says proxy, and E-168 gives generation at loop-free vertices. Acceptable if the edges' mutated vertices are loop-free (not stated for these 12 chains).
4. Step 8 UNDECIDED is honestly reported; "12 of 13" does not claim 13.

## New?
Grepped E-201, E-212, E-228, E-295, E-168, "key-guarded". E-201 holds the 13 steps and the Hom(T,T[-1]) != 0 fact; E-212 does class 2 (8 of 9, record 6 open) by the same method; E-228 closes class 2 and notes the no-NO-power issue. Nothing in research/ records class-1 chains. New: class-1 Hom-verified chains for 12 of 13 steps; the K0 determinant in the verifier and its control (118/118, 708/708) is a small increment. Not recorded yet, so no duplicate.

## Evidenced?
Mostly. Per-step move strings, edge counts and timings are given and reproducible. Weak: the Power line mixes controls (corruptions 78-104 for steps 2, 4, 5 only; K0 control is on a different object set) and says nothing about steps 0-1, 3, 6-12; the "2 of each accepted" are argued to be true YES but not shown. The power statement admits there is no NO-capable case at n <= 7 (no key-equal inequivalent pair), so the chains are consistency evidence for the R-candidate, not a test that could have refuted it. That is stated; keep it in the title/claim wording.

## Scope
n = 7, key class 1, the 13 steps of one checkpoint (nexp 23000), depth/cap as stated. "Join parent and child" should read "up to relation coefficients and assuming generation". The R-candidate "a key-guarded walk stays in one class" is not moved by this beyond the class-2 level of E-228 (12 of 13, not all, and no NO power).

## Required for acceptance
1. Put "meet iso ignores relation coefficients; generation by T assumed (E-168)" in the Claim sentence and the title scope, not only in the header. [small]
2. State for the 12 chains that every mutated vertex is loop-free (E-168 condition), or say it was not checked. [small]
3. Split the Power line by test: of the corruptions rejected, how many by Hom, by K0, by Cartan, by meet iso (cf. E-228: Hom+Cartan rejected 0 of 52). [small]
4. Say explicitly that the result is a YES-side consistency check with no NO case available at n <= 7 and does not move the R-candidate beyond E-228's level. [small]
5. Step 8 deeper search [next round].
