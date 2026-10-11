# On gate-admitted algebras dim J_i <= d_i - 1 (d_i = dim e_iAe_v), so d_i = 2 forces dim J_i = 1; "d_i <= 2 whenever J_i != 0" on walks is NOT implied and stays empirical; Hom(N,N[-1]) = 0 on acyclic algebras, so silting-not-tilting is exactly J != 0

author: theorist · round: 034 · kind: proof
thread: T5 · bears on: E-124, E-126, E-112, E-116

## Claim

Setting: A = kQ/I finite dimensional (right modules, paths left to right), v a vertex, out-arrows b of v, d_i = dim e_iAe_v, J_i = {x in e_iAe_v : x b = 0 in A for all out-arrows b} (= Hom(S_v, e_iA), E-124).

(L1) If v passes `procedure.isMutable` then J_i != e_iAe_v for every i != v with d_i >= 1; hence dim J_i <= d_i - 1, J_i != 0 forces d_i >= 2, and d_i = 2 forces dim J_i = 1. No monomial, scalar or length hypothesis.
(L2) If Q has no oriented path from an out-neighbour t of v back to v (e.g. Q acyclic), then Hom_K(N,N[-1]) = 0 and Hom_K(D,N[-1]) = 0, where T = D + N is the 2-term complex of the mutation (N = [D' -g-> P_v] in degrees 0,1, D' = sum_b P_{t(b)}, D = sum_{j != v} P_j). Hence Hom(T,T[-1]) = Hom(N,D[-1]) = sum_i J_i, and (given that T is silting, AI 2.31, taken from the literature, not re-derived) "silting not tilting iff some J_i != 0" holds for acyclic A.
Not claimed: that d_i <= 2 at J_i != 0. That is the entire empirical content of E-126's "always inside a 2-dim e_iAe_v".

## Evidence

Proof of L1. isMutable (procedure.py:172-183) rejects v iff some basis path p into v, nonzero in A, has p b in I for every out-arrow b. So admission says: no nonzero path image lies in J_i. Paths span e_iAe_v, so if J_i were all of e_iAe_v (d_i >= 1) it would contain a path image, which would be nonzero for some path. Contradiction; so J_i is a proper subspace. d_i = 1: the single path image spans, J_i = 0. d_i = 2: proper subspace of dimension <= 1, and J_i != 0 gives 1. Weakest step: none; this is the observation of E-116 ("the gate tests single paths") read as a dimension bound. It says why dim J_i = 1 is the shape of the exceptions rather than why d_i = 2: the kernel element is a difference x = p1 - p2 (no single path in J), so two paths are needed, and with exactly two there is room for exactly one.

Proof of L2. Chain maps f: N -> N[-1] are f^1: P_v -> D' alone (f^0 lands in N^{-1} = 0, f^2 from N^2 = 0), subject to f^1 g = 0 and g f^1 = 0; there are no homotopies (they would be N^2 -> N^0). So Hom_K(N,N[-1]) is a subspace of Hom(P_v, D') = sum_b e_{t(b)} A e_v, i.e. of paths from t(b) back to v modulo I; with the arrow b: v -> t(b) these form an oriented cycle. None exists in acyclic Q, so the space is 0. Hom(D,N[-1]): D sits in degree 0, N[-1]^0 = N^{-1} = 0, so 0. Hom(N,D[-1]): f^1: P_v -> D^0 = D with f^1 g = 0, i.e. x in e_iAe_v with x b = 0 for all b, which is J_i. Weakest step: identifying g with the b-multiplications and that this N is AI's mutation of D + P_v for the repo's mutation (E-124 is the same reading; AI text not compared). The silting property (Hom(T,T[>0]) = 0, generation) is cited, not proved here.

Example (E-080 long square, hand-built, n = 5): a->b, a->c, b->d, c->d, d->e, relation abde = acde (v = d, i = a). d_a = 2 (abd, acd), x = abd - acd, x(de) = 0, so J_a = k x; d_b = d_c = 1, J = 0. Hom(N,N[-1]) = Hom(P_d,P_e)-part = e_eAe_d = 0 (no path e -> d). So Hom(T,T[-1]) = k, concentrated in Hom(N,P_a[-1]); T is silting and not tilting, and Hom(T,T[-1]) != 0 detects it, with the same one dimension as the Cartan defect at (row v, i = a) (E-097, E-124).

Checks (script output): hand examples: E-080 gate True, (d_i, dim J_i) = (2,1),(1,0),(1,0); layered m = 3 (i -> a1,a2,a3 -> v -> t1,t2; p1 b = p2 b = p3 b for both b) gate True, (3,2): so the bound d - 1 is attained at d = 3, and dim J = 1 is NOT forced in general. Walk prefixes (gate-admitted v, all i with a path to v): n = 7 c0, 1 500 expansions: (d,J) counts (0,0) 2 107, (1,0) 5 054, (2,1) 8; n = 8 c0, 600 expansions: (0,0) 828, (1,0) 2 870, (2,0) 57, (2,1) 2; violations of J <= d - 1: 0; rows with a path from an out-neighbour back to v: 0 (walks are acyclic). No (d >= 3) pair occurs in these prefixes at all; E-126's 429 J != 0 rows (d = 2 each) are consistent.

Missing hypothesis for "dim J_i = 1 on walks": a bound d_i <= 2 at gate-admitted (v, i) on walk algebras, or at least when J_i != 0. E-118's cone estimate gives d up to 8, so it does not supply it. By L1 the dimension of J_i then is d_i - 1 at most; walks would need an argument that e_iAe_v has two paths and no more, e.g. from the 2-term structure of the earlier mutation that created the second path (not attempted).

## Reproduction

```
.venv/bin/python workshop/rounds/034/theorist_dimji.py hand                    # 2 s
timeout 10m .venv/bin/python workshop/rounds/034/theorist_dimji.py walk 7 0 1500   # 31 s
timeout 10m .venv/bin/python workshop/rounds/034/theorist_dimji.py walk 8 0 600    # 31 s
```

## Prior record

E-126 records dim J_i = 1 as an observation; E-116 states the single-path nature of the gate; E-124 gives J_i = Hom(N, P_i[-1]); the scholar's round-033 note (point 6) states Hom(T,T[m<0]) goes through J with Hom(D,N[<0]) = 0 but not Hom(N,N[-1]) (E-124 referee: "not addressed"). L1 is a corollary of the gate definition and probably elementary; L2 closes the referee's gap for acyclic algebras. grep of research/ for "d_i - 1" / "Hom(N,N" found nothing else. Nothing in RETRACTIONS bears on it.

## Code changed

None in the library. New: workshop/rounds/034/theorist_dimji.py (no tests apply).

## Next

Skeptic: break L2 on a cyclic quiver (does a cycle through v ever occur and give Hom(N,N[-1]) != 0 with J = 0?). Experimentalist: d_i histogram at gate-admitted (v, i) on long walks, and any (d >= 3, J != 0) member of a derived class (would refute "d <= 2 at J != 0 on walks"). Scholar: compare L2 with AI 2.31/2.32 text.
