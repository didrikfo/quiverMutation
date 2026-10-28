# Conditional on the omega-reduction of E-235, d_1 on E has the shape of Lemma A of E-240 (600 parents, no loop at i); and Hom(T,T[-1]) at S = -1 is J != 0 restated (E-128), with one new datum (E-254)

author: theorist · round: 093 · kind: proof (a) + rediscovery with a small table (b)
scope: (a) CONDITIONAL on the omega-reduction gr_omega d = d_1 and "d lowers ||.|| by 1" of E-235 (printed degrees, granted, not proved): the star==1 part d_1 of Def. 3.2 (arXiv:1504.02617v1), checked on 600 random parents with no loop at i (`sf.random_parent(rng, False)`, seeds 0..599, arrow degrees 0, 1, 2; 6 347 letters). Not checked: gr_G of the full D as an independent statement; loops at i; d^2 = 0. (b) n = 9, three replayed paths (E-233, E-239, E-254), 17 J = 0 steps (path prefixes, not a control for the J = 0 premise) and 3 J != 0 steps; E-234 (10 steps) and E-244 (6 steps) NOT replayed, taken from the record.
scope: (a) printed d_1 of Def. 3.2 (arXiv:1504.02617v1) on 600 random parents with no loop at i (round-083 generator, 6 347 letters); parents the search reaches are in this family. (b) n = 9, three replayed paths (E-233, E-239, E-254), 17 J = 0 steps and 3 J != 0 steps; E-234 (10 steps) and E-244 (6 steps) NOT replayed (need the pickles of rounds 083 / 087), taken from the record.

## Response to referee

Item 1 (retitle, make (a) conditional, say H4 tests star==1). Done. Title, Scope and Claim now say: after the omega-reduction of E-235 (gr_omega d = d_1, from "d lowers ||.|| by 1"; granted from printed degrees, not proved), d_1 on E has the hypotheses of Lemma A of E-240. The Evidence paragraph is corrected: `theorist_step4.py` H4 tests the nstar == 1 part of `m.partial`, which is d_1, not gr_G of the full D. New run `theorist_graded.py 600` (output `theorist_graded_out.txt`) looks at all terms of the full D by (star count, W-count change): star count >= 2: 0 terms; star count 1: 572 theta terms (W-count change 0, exactly the monomials of Lemma A) and 2 347 psi terms (change -1); star count 0: 2 354 terms from the parent's own d (1 435 + 333 + 497 with change 0, 11 + 78 with change +1). The star-0 terms are NOT removed by the W-count filtration; they are removed only by the omega-reduction. So "gr_G D = D_theta" is true only for d_1 after that reduction, and the earlier headline "gr_G D on E is exactly D_theta" is withdrawn as worded. The script's own "delta-0 part of the full D" counter (1 295 mismatches) is exactly these star-0 terms and is not a failure of the claim; it shows why the omega-reduction cannot be skipped.

Item 2 (degree range). `sf.random_parent(rng, False)` = `opp.random_instance` minus loops at i; it is graded, not degree 0. Over the 600 parents: arrow degrees 0: 3 303, 1: 1 141, 2: 181; letter degrees 0..3 (3 211 / 2 548 / 527 / 61); max arrow degree per parent 0: 12, 1: 429, 2: 159 parents. The 572 thetas have deg(phi) 0 (440) or 1 (132). So the claim is about non-negatively graded parents with arrow degrees up to 2, no loop at i, and a graded run was already what the sample was; no narrowing to degree 0 is needed. Loops of degree >= 1 elsewhere (not at i) occur in the sample.

Item 3 (signs). The sign of d_1(theta_{phi,alpha}) is (-1)^deg(phi) (the Koszul sign of moving phi past the star letter, which has degree 0), and "rescale theta" means theta -> (-1)^deg(phi) theta, an invertible change of the free generators, which does not change the dg algebra up to isomorphism. This is now checked, not asserted: `theorist_graded.py` compares d_1 theta to (-1)^deg(phi) phi alpha^* with sign included, 572 of 572 exact, 0 mismatches, including 132 with deg(phi) = 1. E-240's own check was in degrees 0..2, but E-240 is the free monomial lemma with D = 0 on X; its statement is sign-insensitive after rescaling, and I make no claim that E-240's cells covered signs of this shape beyond that rescaling.

Item 4 (title and counts in (b)). Done: the "3 of 3 / 0 of 17" is removed from title and Claim; (b) is stated as "Hom(T,T[-1]) at S = -1 is J != 0 restated (E-128); E-254's value is added". The table now carries record identifiers: E-233 (3 / 1 recorded), E-239 (2 / 0 and 1 / 1 recorded), E-234 (Hom[-1] total 2 to 4 recorded, not replayed). Only new datum: E-254, F3 step, 2 / 0 at S = -1 and 1 / 1 at S = 0. The 17 J = 0 rows add nothing beyond E-128 and are path prefixes, not a control of the J = 0 premise (said in Evidence).

Referee points not in the required list: terms of star count 0 or >= 2 are now inspected (above). Reproduction note: both original commands re-ran identically for the referee.

## Claim

(a) CONDITIONAL YES. Granting the omega-reduction of E-235 (gr_omega d = d_1, from "d lowers ||.|| by 1", read off printed degrees, not proved), the W-count graded part of d_1 on E = T(U,W,c,Y,Theta) is exactly D_theta (this is not gr_G of the full D: the parent's star-0 terms are removed only by the omega-reduction): D_theta(theta_{phi,alpha}) = +-phi alpha^* and D_theta = 0 on every other arrow letter. Take X1 = U_i (non-A arrows out of i, not ending at i), Y = {alpha^*}, X = all non-theta letters. These are the hypotheses of Lemma A of E-240 (free path algebra over k^n, X1 and Y disjoint subsets of X, D = 0 on X, one theta for every composable pair (x,y) in X1 x Y, D theta = +-xy). So H(E, D_theta) = T(U,W,c,Y)/(phi alpha^*) = Nbar, and step 4 goes through with Lemma A discharged. The check of d_1 does not fail anywhere; where it would fail is exactly a loop at i (see below). What this does NOT remove: the omega-reduction itself and "d lowers ||.|| by 1" (step 2) and d^2 = 0 are still not proved from Def. 3.2; the printed-v1 caveat stays.

(b) The recorded J != 0 key-kept steps have Hom(T,T[-1]) != 0 at the S = -1 placement and Hom(T,T[1]) = 0 there; the replayed J = 0 steps have both 0 at S = -1. That is a theorem, not a finding: Hom(T,T[-1]) = sum_i J_i for acyclic parents with T silting (E-128 L2, AI 2.31 read in E-171), and E-161 already recorded "the test is Ladkani's criterion restated". So it is a detector of J != 0 by construction and says nothing about the J = 0 premise. At S = 0 (the wrong placement) Hom(T,T[-1]) and Hom(T,T[1]) are nonzero at J = 0 steps (up to 5 and 10), so the "detector" must be stated at S = -1. Not claimed: anything about E-234 / E-244 beyond the record, or about cyclic parents.

## Evidence

(a) Where each hypothesis is used. Letters (Def. 3.2): U = non-A arrows not ending at i, W = psi ending at i, c = alpha psi, Y = alpha^*, Theta = phi alpha^{-1} for phi in U_i, alpha in A.
1. D = 0 on X after gr_G. d_1(psi) = sum_alpha alpha^* c_{alpha psi}: W-count drops by one (c has W-count 0), so it vanishes in gr_G. d_1 = 0 on U, c, Y (phi in U does not end at i: "no loop at i" is used here and only here, plus in the next item). So gr_G D = 0 on all non-theta arrows. D_theta keeps the W-count (phi alpha^* has none, theta has none).
2. D_theta(theta_{phi,alpha}) = +-phi alpha^*, a single monomial. With a loop phi at i there is a second term sum alpha'^*(alpha' phi alpha^{-1}) (E-235 control: wrong shape on 6 968 of 10 206 loop letters); that is the failure mode, and it is the E-248 / E-254 case "(c) loops", not the 25 failing n = 7 parents (acyclic Gabriel quiver, no loop at i, round 090).
3. Disjointness and completeness: X1 = U_i letters, Y = star letters, disjoint as letter sets. Thetas exist for all pairs U_i x A (no pair missing, none extra), so every composable (x,y) in X1 x Y carries one theta (alpha^* ends at i = source of phi, so all pairs compose). Sign (-1)^deg(phi) per theta (checked 572 / 572 in `theorist_graded.py`): rescale theta -> (-1)^deg(phi) theta, an isomorphism of the free generators.
4. The weight splitting of E-240 handles the G-filtration: G_q is exhaustive, G_{-1} = 0, and within a fixed weight W-count <= weight, so there are finitely many q and no convergence issue; gr_G pi is the map of Lemma A; five lemma gives G_q E -> G_q Nbar quasi-iso by induction, and a colimit finishes.
Machine check (`theorist_step4.py`; NOTE H4 tests the nstar == 1 part of `m.partial`, i.e. d_1, not gr_G of the full D; see `theorist_graded.py` for the full-D term table; 600 parents, 6 347 letters, 572 theta, 1 187 psi with 2 347 terms, 453 parents with |A| > 1): H1 (theta set = U_i x A) 0 failures, H2 (d_1 theta a single +-phi alpha^*) 0, H3 (X1 and Y disjoint) 0, H4 (every non-theta letter has d_1 = 0 or is psi with W-count down by exactly 1, theta keeps W-count) 0. This re-checks the shape of E-235 part A, now phrased as the hypotheses of Lemma A. Weakest step: item 4 uses E-240's proof as written (it was checked line by line by its referee for the free algebra with arbitrary extra spectator letters; E-240's machine check had at most 6 arrows, so the spectator-heavy shape is covered by the proof, not the cells).

(b) Table: dimension totals (sum over all (T_i, T_j)) of Hom(T,T[-1]) / Hom(T,T[1]), `stepTest` (round 050, code independent of tiltingPlus), replay of the recorded paths by the repo's mutation. Full output `theorist_homtable_out.txt`.

| case (root) | step | tiltingPlus | S = -1: [-1] / [1] | S = 0: [-1] / [1] |
|---|---|---|---|---|
| E-233 (3060000; recorded 3 / 1) | 4 J = 0 steps F1 R9 R8 R8 | True | 0 / 0 each | 5/1, 3/1, 3/2, 1/10 |
| E-233 | R8 (J != 0, perI {4,5,6}) | False | **3 / 0** | 0 / 1 |
| E-239 (3304000; recorded 2 / 0, 1 / 1) | 5 J = 0 steps F2 R8 R9 R7 R7 | True | 0 / 0 each | 3/1, 3/1, 4/1, 3/2, 1/8 |
| E-239 | R7 (J != 0, perI {5,8}) | False | **2 / 0** | 1 / 1 |
| E-254 (2400230) | 8 J = 0 steps F1 F1 F2 F3 R9 R7 R9 F3 | True | 0 / 0 each | 3/1 ... 1/8 (all 8 nonzero) |
| E-254 | F3 (J != 0, perI {1,5}) | False | **2 / 0** | 1 / 1 |

At S = -1 (identity E-128, not a finding): Hom[-1] total equals sum_i J_i (3 = 1+1+1, 2 = 1+1, 2 = 1+1), as E-128 says. E-233 and E-239 reproduce the recorded Hom values (3 / 1; 2 / 0 and 1 / 1): positive control. E-254's Hom was not in the record; it is the only new datum (2 / 0 at S = -1, 1 / 1 at S = 0). From the record, not replayed: E-234 (10 steps, Hom[-1] total 2 to 4, all fail), E-244 (6 steps, all Hom[-1] != 0).
Note the 17 J = 0 steps are the steps on the path to the J != 0 parent, a convenience sample of the recorded walks, not a sample of all J = 0 steps; E-161 / E-050 has the broad control (200 J = 0 steps accepted, 1 081 step tests agree with tiltingPlus; 25 of 25 failing steps rejected).

What a NO with power would have to be. Hom(T,T[-1]) != 0 at S = -1 says only that this one-step complex T is silting and not tilting, so End(T) is not obtained by a tilting step from the parent; it is silent on whether the repo's child (or any later algebra) is derived equivalent to the parent by some other T', and every recorded J != 0 child is in fact joined to the class by J = 0 moves. A NO needs an object that is not a property of T: (i) an invariant of the child as an algebra (derived invariant not seen by Cartan, e.g. Euler-form defect / Hochschild data) that differs from every member of the class, or (ii) a proof that no tilting complex has End equal to the child (a statement about all of D^b, e.g. via a numerical invariant of K_0 with the Euler form: the Cartan congruence class is already equal at J != 0 steps, E-161 "Cartan matrices stay congruent"), or (iii) a search whose miss means something (E-238 power table: the verifier has thin NO power, E-212: 1 rejection in 195 corruptions). Hom(T,T[-1]) is (by E-128) the J-test and can only enter as the label "J != 0".

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/093/theorist_step4.py 600      # a few seconds; theorist_step4_out.txt
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/093/theorist_homtable.py   # under 1 min; theorist_homtable_out.txt
timeout 10m .venv/bin/python workshop/rounds/093/theorist_graded.py 600   # graded / full-D term table; theorist_graded_out.txt
```
Both import round-083 / 049 / 050 code (`theorist_starfilt.py`, `toolsmith_tiltpath.py`, `skeptic_tilt.py`).

## Prior record

E-240 states the gap this closes ("until that check is done ... conditionally"). E-235 step 4 and E-248 / E-254 on loops. (b): E-128 (Hom(T,T[-1]) = sum J_i, acyclic), E-161 (Hom test accepts all 324 printed edges, rejects all 25, "Ladkani's criterion restated"), E-168 (Hom(T,T[-1]) != 0 on the 16 failing steps), E-233 / E-239 / E-234 / E-244 (the Hom values), E-228 (Hom+Cartan has no NO power). New: the E-254 value, the S = 0 versus S = -1 placement contrast across J = 0 steps, and the statement that the detector question is settled by E-128 as an identity, not an empirical regularity.

## Code changed

None in `quivermutation/`. New: `theorist_step4.py`, `theorist_homtable.py`, `theorist_graded.py` and their outputs. No tests touched.

## Next

- Theorist: line-by-line derivation of "d lowers ||.|| by 1" from Def. 3.2, then d^2 = 0 (the remaining inputs of E-235).
- Skeptic / scholar: option (i) above, is there any derived invariant of the child (not of T) that is computable at n = 9 and would reject a J != 0 child? E-253's non-Cartan invariant search (agenda 2) is the same question.
- Whoever replays E-234 / E-244: add the S = -1 and S = 0 columns for all 16 steps and J = 0 controls from the same parents (flagged-but-J = 0 steps), not only the path prefix.
