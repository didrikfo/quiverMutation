# At n = 7, the End(T) comparison rejects 1392 of 1392 equal-dimension variants in the line-point x bilinear-form family at step 12 of class 1, and 18 of 18 monomial/binomial one-or-two-relation replacements (distinct ideals) at the 9 class-2 failing steps, and accepts all 594 isomorphic ones

author: skeptic · round: 062 · kind: result
thread: T10 (and T5) · bears on: E-174, E-177, E-180, H-015
scope: n = 7; class 1 step 12 only (c1 numbering of E-174, the child with two parallel pairs (5,7), (7,4)); class 2 all 9 failing steps (truth over the full torus, via GF(13)); label-preserving iso only; wrong algebras = the child's quiver with relations changed (family below, not exhaustive); capped pool for class 2 (one or two relations replaced); the 370 path edges and n = 10 edges not rerun.

## Response to referee

1. Title: done; both families named (class 1: line-point x bilinear-form family; class 2: one/two-relation replacement by monomial or +-1 binomial from the path pool).
2. Parallel arrows: done. `no_parallel_arrows` (new) asserts all (source,target) pairs of arrows are distinct and no loops; it passes at all 9 steps, so with `no_parallel_longer` Aut = torus holds as stated.
3. Distinct ideals and hashes: singles, pairs and sign flips now share one deduplication set per step (ideals compared as subspaces). Class 2 has 36 distinct equal-dimension ideals (18 non-iso, 18 iso), not 60; the old 60 had repeats. Skipped by dimension: 9882 (not deduplicated). sha256 of c1.pkl and c2.pkl are in Reproduction. Step 12 is `recs[12]` of the failing-only records of c1.pkl, so pickle order matters.
4. NOT labels: done beyond mu_4. At each class-2 step the child ideal is monomials plus exactly one binomial c1 p + c2 q; its torus image is the same with q scaled by lam = t^q/t^p. `torus_iso_full` tries every lam in GF(13)* (a superset of what the torus reaches), so NOT is a proof of NOT for the whole torus. Full-torus truth and the old mu_4 truth agree on all 36 ideals; 0 disagreements. Residual caveat: arithmetic is mod 13, not over C.
Class 1 re-run unchanged (1392/1392, 576/576, 0 disagreements).

## Claim

E-177 left two cases out of the equal-dimension wrong-algebra control: step 12 of class 1 and class 2. Both are now done. The assertion of E-180 (`dimcheck`) cannot reject any of these variants, since dim K Q/I' equals dim End(T) pair by pair for each of them (filtered with the own linear algebra `quotient_dims`). The remaining test, `symcheck2` (E-174), then gives the verdict. Against independent truth (below) it rejects every non-isomorphic variant (step 12: 1392 of 1392; class 2: 18 of 18 distinct ideals) and accepts every isomorphic one (576 of 576; 18 of 18). No disagreement. So the premise "End(T) iso to the child" is not vacuous at these ten steps against wrong algebras of equal dimension and equal dimension matrix.
It does NOT claim: power against wrong algebras outside these families (e.g. different quiver, relations of higher length mixtures with other coefficients), or anything about J (E-181: End(T) = H_0 cannot see J). Refuted by a single equal-dimension non-isomorphic variant with verdict "iso", or an isomorphic variant with verdict "NO".

## Evidence

Step 12 (child: arrows 1->5, 2->6, 3->5, 6->5, P1 = {b4,b5}: 5->7, P2 = {c0,c1}: 7->4; 7 relations). Family: the three line relations x*(l0 b4 + l1 b5) = 0 for x in {1->5, 3->5, 6->5}, each re-chosen among 4 points of P^1 (64 choices), times the kernel of a bilinear form phi on P1 (x) P2 replacing the three relations on 5->7->4 (40 forms phi with entries in {-1,0,1} up to sign; rank 2 and rank 1). 2560 variants; 592 have a different dimension (skipped by the dimension assertion anyway); 1968 have equal dimension and equal dimension matrix.
Truth (independent of both tests, proved by orbit counting; assertion in the script that no arrow has a longer parallel path, so Aut over the quiver is GL2(P1) x GL2(P2) x scalars): phi of rank 2: iso to the child iff the coincidence pattern of the three points equals the true one (all distinct), since GL2(P2) is transitive on nondegenerate phi and GL2(P1) is simply transitive on ordered triples of distinct points up to scalars; phi of rank 1: never iso (the rank of phi is an invariant). The true child has all points distinct and phi of rank 2.

| class | variants (equal dim) | truth iso: test accepts | truth not iso: test rejects | disagreements |
|---|---|---|---|---|
| 1, step 12 | 1968 (1392 non-iso, 576 iso) | 576 / 576 | 1392 / 1392 | 0 |
| 2, 9 steps | 36 distinct ideals (18 non-iso, 18 iso) | 18 / 18 | 18 / 18 | 0 |

Step 12 breakdown of equal-dimension variants by (rank phi, point pattern) is in the `.out` file; non-iso ones are 960 rank-2 variants with a coincident pair or triple of points and 432 rank-1 variants (all patterns). All 1392 NO verdicts are the string "NO label-preserving iso (1 in the ideal ...)", not a timeout or an arrow-count mismatch.
Class 2 (no parallel arrows, 7 arrows, one commutativity-type binomial and 3-4 monomial relations per step): the pool is every path of length >= 2 as a monomial and every difference of two such paths with the same ends; variants replace one relation (i), or two (i, j), or flip the sign of the binomial (control). Truth: Aut is the torus (checked: no arrow has a longer parallel path), so iso iff some t in mu_4^arrows (computed in GF(13), first tried t in {+-1}) carries the ideal of the child onto the variant ideal, compared as subspaces in every path space. Most replacements change the dimension (9882 replacements skipped: dimension rejects them without help from symcheck2); the equal-dimension survivors are few (2 non-iso and 2 iso per step, deduplicated by ideal across singles, pairs and flips).
Truth for class 2 uses the full torus (lam search over GF(13)*, response 4); it agrees with the earlier mu_4 search on every variant.
Caveat: equal-dimension non-iso variants at class 2 are scarce because zero relations fix the Cartan matrix; the control there is weaker than at step 12.

## Reproduction

Pickles (about 8 min and 6 min each, run each twice: the first stops at DEADLINE, the second resumes; the environment first needed `uv pip install -e '.[test]'`):
```
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 2 20000 /tmp/tsm/c2.pkl 100
CLS=1 timeout 10m .venv/bin/python workshop/rounds/062/skeptic_wrongalg2.py   # 37 s, step 12
CLS=2 timeout 10m .venv/bin/python workshop/rounds/062/skeptic_wrongalg2.py   # 80 s, 9 steps (PAIRS=0 skips the two-relation replacements; step 12 = recs[12] of the failing records of c1.pkl)
```
Outputs: `workshop/rounds/062/skeptic_wrongalg2_c1.out`, `..._c2.out`. Pickles: c1 231406 bytes sha256 79d498c04b5227e65fb1abe578dfbca66b00e59d565b30a135f48061853c1192; c2 212071 bytes sha256 469d9134d0de9bdc2bdc844b7d564ea53901899ab10c85ff47c1da67289c80b8 (16 and 9 failing steps, as in E-180).

## Prior record

E-177 (equal-dim control at 7 one-parallel-pair c1 steps; listed step 12 and class 2 as not done), E-180 (dropped-relation control 125/125, incl. step 12 and class 2; stated equal-dimension at those steps untried), E-181 (End(T) = H_0, so End(T) does not see J), E-184. This closes the E-177/E-180 open item for the child's relation set; it does not change what the comparison can say about J.

## Code changed

None to existing files. New: `workshop/rounds/062/skeptic_wrongalg2.py` (reuses `rounds/057/toolsmith_endt2.py`'s `symcheck2`, `quotient_dims`, `dimcheck`). No tests run (nothing under `quivermutation/` touched).

## Next

- toolsmith/skeptic: the same equal-dim control on the 370 n = 7 path edges and the 35 n = 10 edges (E-180 left them); wrong algebras there with a different quiver (the verdict "arrow counts differ" is then trivial).
- theorist: an orbit-count proof for class 1 steps 0-11 would remove the need for the invariant-based truth; the E-181 question (Thm 1.1 and J = 0) remains the live premise, not this one.
