# The dg cancellation of Prop. 7.6 (arXiv:1504.02617v1) when there is no loop at i, of any degree: a proof modulo the cited Anick/Halperin-Lemaire lemma (Lemma A) and the printed-degree fact that d lowers ||.|| by 1, with machine checks of the star-count shape and omega-monotonicity on about 16k letters each, and toys with one and two (phi, alpha) pairs

author: theorist · round: 083 · kind: proof
thread: T5 · bears on: H-015 (J = 0 premise, link "repo step 7 = End(T)"), E-214, E-217, E-227
scope: v1 LaTeX only (journal version UNVERIFIED). A proof, NOT fully rigorous: it rests on (i) Lemma A, cited (Anick strongly free sets), (ii) "d lowers ||.|| by exactly 1 on each printed term", read off the printed degrees and machine-checked, not proved line by line from Def. 3.2. Machine checks: star-count shape on 4 249 letters (Part A), omega-monotonicity on 15 966 letters of 1 500 parents (Part C), exact-rank toys with 1 and 2 (phi, alpha) pairs to weight 12 / 15 (Parts B, D). NOT covered: a loop of degree >= 1 at i (Step 4 letters); the journal version; a real kM beyond letters.

## Response to referee (round 083)

Verdict was minor revision. Point by point:

1. "Rigorous" overstated. ACCEPTED. Title and scope narrowed to the referee's wording (done above). Omega-monotonicity is not re-derived from Def. 3.2 this sitting: it needs the printed terms of d on every letter type (psi, c, phi, theta, products), which I have only as code; I state it as an input (Step 2) and as machine-checked (item 3).
2. Lemma A named. DONE. Cited statement: if the quadratic monomials xy (x in X_1, y in Y, X_1 and Y disjoint) have no overlap, the set {xy} is strongly free (Anick 1982/1986; Halperin-Lemaire 1987), the Anick resolution of T/(xy) has no chains of length >= 3, and T<theta_{xy}> with D theta = xy is a free resolution-model of T/(xy), so T<theta> -> T/(xy) is a quasi-isomorphism. The Morse-matching acyclicity is a re-proof of that, not an extra input. Referee's one-line argument (position of leftmost bad spot decreases) I did NOT verify as a full gradient-path argument: a path may re-create a bad spot to the left only through an xy pair, which needs a y next to an x, impossible when no overlap, and that is the idea of the cited lemma; I leave it cited. Scholar to vet the exact theorem number.
3. Omega check added. DONE. `theorist_omega.py` Part C (my own merge of `experimentalist_omega.py`, other seeds 5000-6499): 1 500 parents, 15 966 letters (1 137 parents with |A| > 1, 1 474 theta letters, 7 400 star terms): ||.|| >= 0 on all letters and all terms of d(L): 0 violations; no term of d(L) with omega above omega(L): 0; no star-free term of d(L) with omega equal to omega(L): 0. Agrees with the referee's 6 454-letter run. These are the two facts Steps 2-3 need.
4. Cycles wording. DONE. Coverage paragraph now says: hypothesis (H) allows cycles through i and the generator (Part A, B, D) produces them; parents the search actually reaches have no cycles at all, so for those only the loop-free, cycle-free case is used.
5. Part B with two phi, two alpha (asked for next round). DONE early: Part D of `theorist_omega.py`: 13 letters, vertices i, j1, j2, k, m1, m2, letters as1, as2, psi (d1 = c1 as1 + c2 as2), phi1, phi2, four thetas t_pq (d1 = as_q phi_p), cycles i -> m_p -> k -> i. Kernel K of E -> Nbar acyclic in all 377 cells (weight <= 15, 64 608 words), 0 bad; control (thetas not killed) gives 360 nonzero cells. This exercises several phi, several alpha, and the W-count interplay (psi has two terms), at the level of the associated graded complex d1 only; it is still a toy, not kM.
6. Random-instance ranges. Part A/C use `random_instance` of round 069: 4 vertices (i = 0), |A| 1 to 3 arrows of degree 0 out of i (targets 1, 2), 1 to 2 degree-0 arrows into i, arrows 1->3, 3->2, 1 to 3 degree-1 loops at i (removed), 1 to 3 further degree-1 arrows between random vertices (loops and 2-cycles allowed), up to 2 degree-2 arrows from products of degree-1 pairs, d of a degree-1 arrow = 0 to 3 random degree-0 paths of length <= 4 with coefficients in {1, -1, 2}. Note only degree-1 loops at i are generated; a degree-0 non-A loop at i is excluded by the generator and by (H).
7. Not done: nothing required is left over. Not done and flagged: line-by-line derivation of ||.|| decrease; journal version; loops of degree >= 1 at i.

Files: `theorist_omega.py`, `theorist_omega_out.txt` (new); `theorist_starfilt.py` unchanged.

## Claim

Hypotheses (H): the graded dg quiver (kQ, d) is non-negatively graded, i is a vertex with NO loop at i (of degree 0 as in Thm 1.1, and also none of degree >= 1). Nothing else: not finite dimension, not acyclicity of D(i) or of Q, not J, parallel arrows and cycles anywhere (even through i) allowed. Then kM -> N := kM/(theta, d theta) is a quasi-isomorphism, where theta runs over the letters phi alpha^{-1} (Step 3; Step 4 letters do not exist under (H)), i.e. Prop. 7.6 holds with its first relation only.
Coverage: (H) allows cycles through i and the generator produces them, but every parent the search reaches (`noCycles`, search.py l. 377/420, `simple_cycles` sees self-loops) has no cycle at all and satisfies (H), because its Construction 2.2 model Q^(inf) has every arrow s -> t with s < t strictly in the reachability order of the acyclic parent (relations and H_k classes live in e_t(.)e_s with a path s ~> t of length >= 1). So no loop, no cycle, D(i) acyclic. E-211 also satisfies (H) (no loop at 1): its failure is step 7, not cancellation, as E-214/E-220 said.
Not claimed: loops of degree >= 1 at i; the rest of Thm 1.1 (d^2 = 0 is E-202, numerical; W1, W2 of E-217; Thm 7.5); that T is tilting (J).

## Evidence (proof)

Notation. Letters of kM under (H): U = non-A arrows not ending at i; W = arrows psi ending at i (none starts at i); c = compositions alpha psi; Y = stars alpha^*; Theta = phi alpha^{-1} with phi in U_i (U starting at i, non-A). Degrees ||.|| >= 0 on all letters (||alpha psi|| = ||psi|| - 1 >= 0 since psi in W has ||psi|| = |psi| + 1; phi out of i has |phi| >= 1 as it is not in A). s(w) = number of star letters; omega(w) = ||w|| + s(w).

1. Shape of d (Def. 3.2, l. 383-395). On a letter, d has terms with s = 0 (d_0) and s = 1 (d_1), never more: only the term sum_alpha alpha^* (alpha psi) in d(psi), psi in W, and the term +-phi alpha^* in d(phi alpha^{-1}) contain a star; d(phi) (phi in U), d(alpha psi), d(alpha^*) = 0 and the rest of d(theta) are star-free (the Delta-insertions give products theta . c, no stars). Hence d = d_0 + d_1 as derivations, d_1^2 = 0. MACHINE: Part A below.
   d_1(psi) = sum_alpha alpha^* c_{alpha psi}; d_1(theta_{phi,alpha}) = +-phi alpha^* (phi in U_i, so d_1 phi = 0 because phi does not end at i: this is where "no loop at i" is used); d_1 = 0 on U, c, Y.
2. Filtration F_p = span{ w : omega(w) <= p }, increasing, exhaustive, F_{-1} = 0. d lowers ||.|| by 1 and raises s by 0 or 1, so d F_p subset F_p and gr d = d_1. So gr kM = E := (T(U,W,c,Y,Theta), d_1).
3. N has a basis of normal words. N = T(U,W,c,Y)/(rho_{x,alpha}), rho = x alpha^* - e, x in U_i, e = (d theta mod theta part, star-free, omega(e) < omega(x alpha^*)). The leading words x alpha^* have no overlap (a leading word is (letter of U)(letter of Y); U and Y are disjoint letter sets, so no suffix of one is a prefix of another, in either reading order), and rewriting x alpha^* -> e lowers omega, hence terminates: by the diamond lemma the normal words (no x alpha^* factor) are a basis of N, F_p N has the normal words of omega <= p, and gr N = T(U,W,c,Y)/(x alpha^*) =: Nbar with the differential d_1 induced (d_1 psi = sum alpha^* c). gr(pi): E -> Nbar kills Theta and the words with a factor x alpha^*; its kernel is K = ideal generated by theta and x alpha^*.
4. gr(pi) is a quasi-isomorphism. Filter E and Nbar by G_q = {W-count <= q}. On generators d_1 = D_theta + D_W, D_theta(theta) = +-x alpha^* keeps the W-count (x in U), D_W(psi) = sum alpha^* c lowers it by 1. So gr_G d = D_theta on E and 0 on Nbar. Lemma A (below) gives H(E, D_theta) = Nbar (W, c, rest of U are spectator letters). G is exhaustive and bounded below, so induction on q and passage to the colimit give: E -> Nbar is a quasi-isomorphism.
5. By induction on p (F_{-1} = 0) and exactness of colimits, gr_p(pi) quasi-iso for all p gives pi: kM -> N quasi-iso. QED.
Lemma A. T = T_R(X) free path algebra, X_1, Y subset X disjoint, all composable pairs (x,y) in X_1 x Y relations, T<theta_{xy}> with D theta_{xy} = xy and D = 0 on X: then T<theta> -> T/(xy) is a quasi-isomorphism. Proof: algebraic Morse matching on words in X + Theta: match a word whose leftmost bad spot is theta_{xy} with the word where theta is replaced by xy (a bad spot is a theta or an adjacent pair x,y). Well defined because x not in Y and y not in X_1 (no overlap, neither (a,x) nor (y,b) can be bad), coefficient +-1; critical cells = words without theta and without a bad pair = normal words, D = 0 between them. Acyclicity of the matching: cited (statement: {xy} with no overlaps is strongly free, Anick resolution has no chains of length >= 3; theorem number to be vetted by scholar), this is the Anick chain complex with no chains of length >= 3 (Anick 1982 "strongly free", Halperin-Lemaire 1987; Sköldberg 2006, Jöllenbeck-Welker 2009); I did not rederive it (LEAST SURE STEP; toy check below).
Why the printed words "never overlap" is not the whole proof (E-227): non-overlap gives only Lemma A for the star-part d_1; the content of steps 2-3 is that d_0 (which contains theta in d(phi) through Delta) is of lower filtration, i.e. the star count.

## Machine checks (Part A, B of the script)
A. For random parents (round-069 instance generator, arrows of degree 0, 1, 2, loops and 2-cycles anywhere, loops at i removed): 4 249 M-letters, 412 of type theta: 0 letters with a term of star count >= 2, 0 letters whose s = 1 part differs from step 1, 0 with a star-count-2 part in d^2 (that is d_1^2 = 0). Control (loops at i allowed): 14 403 letters, star count <= 1 still, but d_1 on 6 968 of 10 206 loop letters (psi, theta, a-psi-b^{-1}) is NOT the step-1 shape: the d_1 there contains alpha^* (alpha phi) alpha^{-1}-type terms, so the monomial Lemma A does not apply; this is the open case.
C. Omega-monotonicity (`theorist_omega.py`, 1 500 parents, 15 966 letters): 0 negative ||.||, 0 terms of d(L) with omega > omega(L), 0 star-free terms with omega = omega(L).
D. Two phi, two alpha toy: H(K) = 0 in 377 cells to weight 15 (64 608 words); control nonzero in 360.
B. Toy E (cycle i -> m -> k -> i present, letters as above, 8 letters): exact ranks mod 1 000 003 of the kernel K = (theta, phi alpha^*): H(K) = 0 in all 755 (weight, source, target, degree) cells with weight <= 12 (17 998 words). Control (theta not killed): H(K) != 0, so the test can say no.

## Reproduction

```
.venv/bin/python workshop/rounds/083/theorist_starfilt.py 400 12     # 1 s; output theorist_starfilt_out.txt
```
It imports `workshop/rounds/069/scholar_oppermann_sign.py` (executable Def. 3.2 with the printed signs).

## Prior record

E-214, E-196, E-220: hypotheses of Thm 1.1 only; E-227: no overlap lemma printed, "not an Anick/Groebner hypothesis" without an order. New: the order (omega = degree + star count) that makes it one, the reduction to Lemma A through a second filtration by W-count, the observation that the cancellation needs only "no loop at i of any degree", so acyclicity of D(i) is not used, and the coverage of every `noCycles` parent. Consistent with E-217 (Gamma_red = N), E-211.

## Code changed

None in `quivermutation/`. New: `workshop/rounds/083/theorist_starfilt.py` (+ `_out.txt`). No tests touched.

## Next

- Loops of degree >= 1 at i (Step 4, Theta2 letters): d_1(phi alpha^{-1}) = +-phi alpha^* +- sum alpha'^* (alpha' phi alpha^{-1}) is not monomial; try a third filtration by loop count, or the matrix structure psi, theta1, c, theta2. Needed only for direct callers, not the search.
- scholar: check the Morse-matching acyclicity against Sköldberg/Anick statement and the journal version of Prop. 7.6.
- toolsmith: Part B on a real kM of a small cyclic parent (needs the full sign-checked d on products).
- Remaining links of "repo child = End(T)": d^2 = 0 (Sect. 9, E-202 numerical only), W1/W2 of E-217, Thm 7.5, Constr. 2.2 for the repo parents.
