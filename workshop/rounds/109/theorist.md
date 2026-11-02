# For the loop model M1 (relation loop of degree 1 at i, d f = a b) the cancellation kM -> N = k<fbar,u,v> (d fbar = uv) of Prop. 7.6 is proved a quasi-isomorphism by a change of generators and Kunneth, with H^0(N) = A; and the mutation T = P_j + cone(P_i -> P_j) of A = kQ/(ab) is tilting with Cartan(End T) = Cartan(A) (one model, degree-1 loop; N ~ A in all degrees and End(T) = A as algebras are not machine-proved)

author: theorist · round: 109 · kind: proof

## Response to referee

Verdict minor revision; all four required items done, no re-run changed a number.
1. N versus A. Changed the wording. What is shown: N = k<fbar,u,v>, d fbar = uv, has H^0(N) = k<u,v>/(uv) = A (degree 0 is direct: u, v have degree 0, fbar is the only degree 1 letter, and d fbar = uv, so the degree 0 coboundaries are the ideal (uv)). That N is quasi-isomorphic to A in ALL degrees is the standard Anick-type fact for a monomial relation without self-overlap (u != v, so uv does not overlap itself); it is NOT proved or computed here and the primary source is UNVERIFIED. It is not used: the cancellation q.i. kM -> N in (2) does not need it. The sentence "produces a copy of A" now reads "H^0(N) = A".
2. The 52-cell check. Stated now: it tests stage 1 only (dim(kM/I) = number of words in {fbar,u,v}, a dimension statement per cell). It does not test H(I) = 0 or H(N); H(I) = 0 is E-284's computation and, here, the Kunneth step of stage 2 (a proof, not a check). Docstring of `theorist_m1free.py` corrected (it no longer claims to print d on survivors). Runtime corrected to under 1 s.
3. Generation and shift convention. Generation needs no citation for M1: X = cone(P_i --b--> P_j), so P_i = cone(P_j -> X)[-1] lies in thick(T), and P_j is a summand; hence thick(T) contains A = P_i + P_j and equals K^b(proj A). (The right mutation: same triangle with a.) With Hom(T,T[m]) = 0 for m != 0 (range -3..3 only; vanishing for all m is not shown by the script, though T is a two-term-per-summand complex of projectives over a gl.dim 2 algebra so |m| >= 3 is automatic from degree width) this makes T tilting. Okuyama-Rickard (approximation-cone construction of tilting complexes; Rickard 1989 / Okuyama) is the general source: exact statement UNVERIFIED, not used. Shift convention of `theorist_m1hom.py`: `cone_map` puts the source projective in degree d and the target in degree d+1; `shift(X, s)` moves terms to degree k - s (X[s]). Left row: base has P_i in degree 0, P_j in degree 1; the passing label "shift 1" is therefore P_i in degree -1, P_j in degree 0, as in the prose. Right row: base has P_j in degree 0, P_i in degree 1; the passing label "shift 0" is P_j in degree 0, P_i in degree 1. Prose and labels now agree (Claim (3) says so). End(T) = A as algebras rests on the by-hand argument plus the matching Cartan matrix, not on the script.
4. Prior record: now says E-284 and E-279 already record that Lemma A's monomial premise is false (the bare b term in d(f a^-1)) and that M1 cancels through weight 12; new content is only the (i)-(iii) split, the M1 proof, and the Hom data.
Also noted (reviewer's point 4, "all three fail"): confirmed by the printout, no change.
thread: T5 · bears on: H-023
scope: ONE algebra A = kQ/(a then b), Q: a: i->j, b: j->i (Nakayama, Kupisch (2,3), dim 5, gl.dim 2) = model M1 of E-284, one vertex i; the loop is the degree-1 relation letter f (d f = ab). Proof covers only dg parents where a triangular change of generators makes the survivors a dg subalgebra (checked: M1; NOT M2, M4, M5). Degree-0 loops are excluded by Oppermann and by negative weight (E-279, E-284), not treated. Hom part: Hom(T,T[m]) for m = -3..3, GF(1000003); the free-quotient script checks stage 1 (dimension of kM/I) only; N ~ A beyond H^0 and End(T) = A as algebras are not machine-proved. Sign set 'last = -1' of E-202/E-284 taken as given; v1 of arXiv:1504.02617 only.
power: Hom test: the planted zero-map "approximation" (X = P_i -> P_j with 0) gives Hom(T,T[+-1]) != 0 and Hom(T,T) = 3, 1, 1, 2: it says NO; the wrong shifts of the correct cone also say NO (4 of 6 placements). Free-quotient check: planted extra survivor (keep the pivot m(a,b,None)) gives 14 mismatching cells out of 33 at weight <= 9. The proof itself: n/a.
backlog: n/a

## Claim

(1) Which hypothesis of Lemma A fails. Lemma A (E-240, as used in E-235 step 4 and E-256) needs: (i) each Theta letter has d_1 theta = +-x y, ONE quadratic monomial, x in X1, y in Y; (ii) X1, Y and the Theta letters pairwise disjoint (sites cannot overlap); (iii) D = 0 on the spectator letters. For a loop f at i all three fail, on printed data (`theorist_dtheta.py`). In M1, d(f a^-1) = bbar + theta2 a* - a* fbar: a LINEAR term bbar (star count 0, so the omega-reduction of E-235 does not remove it) plus two quadratic terms; (i) fails; a* is the right factor in one and the left factor in the other and theta2 = a f a^-1 is itself a Theta letter in X1, so (ii) fails; and d fbar = u a* != 0 on a spectator, so (iii) fails. If d f = 0 (M2) the linear term is gone but (ii) and the two-sided a* remain. E-279's phrase "what fails is the monomial D theta = xy" is the same finding; this says which of (i)-(iii).

(2) Two-stage argument, M1. Letters of kM (E-284 listing): bbar = m(.,b,.) deg 1, c = m(a,b,.) deg 0, fbar deg 1 (loop at i), u = m(a,f,.) deg 0 (i->j), v = a* deg 0 (j->i), theta1 = f a^-1 deg 2, theta2 = a f a^-1 deg 1. d theta1 = bbar + theta2 v - v fbar, d theta2 = v u - c, d fbar = u v, d u = d v = 0, d bbar = c v, d c = 0.
Stage 1 (change of generators). Put pi1 = d theta1, pi2 = d theta2. bbar = pi1 - theta2 v + v fbar and c = v u - pi2 are triangular substitutions (neither pi involves bbar or c), so T(V) = T(fbar,u,v) * T(E) with E = <theta1, theta2, pi1, pi2>, d theta_k = pi_k, d pi_k = 0 (d^2 = 0), and T(V') = T(fbar,u,v) is a dg SUBALGEBRA (d fbar = uv, d u = d v = 0). The ideal I generated by the Theta letters and their differentials is the ideal generated by E.
Stage 2 (Kunneth). As a complex T(V' + E) = T(V') + sum_{p>=1} T(V') (x) (E (x) T(V'))^{(x)p}; every summand with p >= 1 has a tensor factor E that is acyclic (d theta = pi is an isomorphism), so I = ker(T(V'+E) -> T(V')) is acyclic by Kunneth over a field, hence kM -> N = T(V') is a q.i. No filtration, no convergence, no Anick lemma and no omega-reduction are used. N = k<fbar,u,v>, d fbar = uv, H^0(N) = k<u,v>/(uv) = A. (N ~ A in all degrees is the Anick-type fact for the non-self-overlapping monomial uv: UNVERIFIED, not used.) The 52-cell check below corroborates stage 1 only (dimension of kM/I); H(I) = 0 is E-284's computation and Kunneth is the proof.
Weakest step: that the printed differentials (sign set 'last = -1', from `scholar_oppermann_sign.py`, E-202) are Oppermann's Def. 3.2; the proof is algebra on those printed formulas. It does not prove Prop. 7.6 for loops in general: M4 (two A-arrows) has d theta_{f,a1} = bbar + ..., d theta_{f,a2} = -bbar + ..., so only the difference has no generator pivot; M5 has pivot ebar bbar (a product); M2 has pivot a* fbar (a product). Those still need the Anick-type lemma with a Theta-perturbation, and are covered only by E-284's computation.

(3) Hom control. T = P_j + X, X = cone(P_i --b--> P_j) = (P_i in degree -1, P_j in degree 0), right modules over A. Hom(T,T[m]) = 0 for all m != 0 in -3..3, so T is a tilting complex; Hom(T_x,T_y) = [[1,1],[1,2]] = Cartan matrix of A (Coxeter polynomials equal). End(T) = A by hand: X -> P_j -> X is 0 (dim End X = 1), P_j -> X -> P_j is the nonzero element ba (the chain map X -> P_j is ba, forced by g b = 0; the homotopy to P_i kills only the ba-map P_j -> X). The right mutation (P_j --a--> P_i, P_i in degree 1) also tilts with the same Cartan matrix. Generation: P_i = cone(P_j -> X)[-1], so thick(T) contains A. Shift convention: see Response item 3 (script labels 'left shift 1' and 'right shift 0' are the placements with P_j in degree 0). This agrees with H^0(N) of (2): H^0 of the dg mutation is A, and the Hom computation says the mutation T is a tilting complex with End(T) = A. Together: the J = 0 step at this loop vertex is a derived equivalence A -> A in this one example (via Rickard/AI, not via the repo, which refuses cyclic D(i), E-211/E-222).

## Evidence

| check | result |
|---|---|
| dtheta printout, models M1, M4, M5, M6 | linear bbar in d theta1 (M1); shared pivot in M4; no loop in M6: only a*.m (Lemma A shape) |
| free quotient, STAGE 1 ONLY: dim(kM/I) = words in {fbar,u,v}, cells (weight, src, tgt, degree), weight <= 12 | 52 cells, 0 mismatches |
| same, planted extra survivor, weight <= 9 | 14 of 33 cells mismatch |
| Hom(T,T[m]), A itself | 0 for m != 0, Cartan [[1,1],[1,2]] |
| T = P_j + cone(P_i -> P_j) (left, P_j in degree 0) | 0 for m != 0 (m = -3..3), Cartan [[1,1],[1,2]] |
| T = P_j + (P_j -> P_i) (right, P_j in degree 0) | same |
| other placements (4 of 6) and the zero-map X | NO: Hom(T,T[+-1]) or Hom(T,T[+-2]) nonzero |

Stays UNVERIFIED: Anick Thm 3.1/3.2, Forre Th. 2.6 numbers (E-273/E-278, secondary only; not used in (2)); journal version of Oppermann; Def. 3.2 signs beyond E-202; any loop model not of the pivot-generator type; any D(i) cycle in the repo's own step 7 (E-211 remains the counterexample on the repo side); whether End(T) is the quiver-with-relations the repo would print (cyclic, refused).

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/109/theorist_dtheta.py          # about 2 s, d(theta) of M1, M4, M5, M6 and all letters of M1
timeout 10m .venv/bin/python -W ignore workshop/rounds/109/theorist_m1free.py 12        # under 1 s, 0 mismatches (stage 1 only)
timeout 10m .venv/bin/python -W ignore workshop/rounds/109/theorist_m1free.py 9 ctrl    # planted error, 14 mismatches
timeout 10m .venv/bin/python workshop/rounds/109/theorist_m1hom.py         # 1 s, Hom table and controls
```
(`theorist_m1free.py` and `theorist_dtheta.py` load `workshop/rounds/106/theorist_loopfilt.py` and, through it, `workshop/rounds/069/scholar_oppermann_sign.py`.)

## Prior record

E-284 (ten models, computation, H(I) = 0; already records that Lemma A's monomial premise is false, the bare b term in d(f a^-1), and that M1 cancels through weight 12), E-279 (same finding in words: what fails is the monomial), E-235/E-240/E-256 (Lemma A, no-loop case), E-273/E-278 (Anick numbers, secondary), E-211/E-222/E-199 (cyclic D(i) on the repo side). New (only): the itemised list (i)-(iii); a proof for M1 that replaces the Anick lemma by a change of generators and Kunneth; the Hom-computed End(T) = A for M1. Not new: that M1 cancels (E-284). Not claimed: any other loop model, a general proof, or anything about the repo's step 7 on a cyclic D(i).

## Code changed

None in `quivermutation/` or tests. New scripts only: `workshop/rounds/109/theorist_dtheta.py`, `theorist_m1free.py`, `theorist_m1hom.py`.

## Next

Theorist: extend (2) to the shared-pivot case (M4: change basis in span of the two thetas, one pivot, one quadratic leading term) and to M2/M5 by a Theta-count filtration plus Lemma A; a random sweep for which loop models are pivot-generator type. Toolsmith: allowCyclic child of M1 via `mutateAtVertex` against End(T) = A (does step 7 give a cyclic quiver with relation uv?). Scholar: Anick primary (E-273).
