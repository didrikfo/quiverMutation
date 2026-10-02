# Scholar's notebook (after round 027)

## Believe
- Ladkani 2.3(c) = AI 2.32(b) = `tiltingPlus`: one map g_i : p |-> (p b)_b. Authors' statement, not derived here.
- E-097: Cartan fails iff some J_i = ker g_i != 0. J is intrinsic to the parent A; step 7 never exhibits an element of J, it only supports the Cartan link.
  r027: on 4 hand cases the real rewrite (checkCartan=True) fails exactly when J != 0, including the non-W kernels.
- r027 lemma (proved for monomial + two-term relations): gate-admitted, J_i != 0 iff graph Gamma_i (vertices = nonzero classes of e_iAe_{t b}; one edge per class of e_iAe_v joining [pb1],[pb2], zero = pendant to ground) has a circuit (cycle or ground-ground path); J = flows.
  W (E-107) = length-2 ground path. Length-2 2-cycle ("nn", p1b1=p2b1 and p1b2=p2b2 nonzero) is the "cancels" branch: gate-admitted, J != 0, Cartan fails, W false (= E-103's 19 two-out hand-built). Length >= 3 circuits exist (case H). Out-degree 1: J != 0 => p1 b = p2 b with p1 != p2 (derived, not necessarily a generator; case G).
- On walks (n=7 c0, n=8 c0, c1; 155 J != 0 rows) J is always a length-2 ground path (nz/zn at out-deg 2, n at out-deg 1): 0 nn, 0 longer. Why: OPEN.
- Sum-type relations (p + q) are real; `alg.rels` loses signs, use `relationsFrom`. A test comparing normal forms must scale-normalise (my first version missed 55 rows because of this).
- Rejecting parents (distinct): n=6 class 0 1123, n=7 156. Out-deg 2 rejects only from n=8.
- CHZ Cor 3.6 "monomial?" UNVERIFIED; arxiv.org blocked (403). Do not retry.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same.

## Did
- R001..R025 (rounds/*/scholar_*.py). R027: scholar_pairtest.py (equal-signature pair test, --hand cases D/G/H, walks n=7,8).
- Do not run two walks plus other persona jobs and `pkill -f` in the same shell (killed my shell once); background with nohup, poll with sleep < 120 s.

## Next
- Prove/refute: LNA-derived algebras have no nn 2-cycle or length >= 3 circuit (relation shape: one relation per start vertex?). Try to build H or D inside an LNA walk by hand.
- n = 9 class 0 pairtest (overnight, needs checkpointing); tally >= 3-term relations among J != 0 rows.
- Classes 1-3 at n = 6..8 for out-degree 1 unrun. Diagonal and i,j != k Cartan entries underived.
- Lesson: a test that "agrees with the data" can be wrong in the other direction; check the sign/coefficient model first (arrowRels vs rels).
