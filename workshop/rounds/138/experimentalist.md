# F-057 (a) at n = 11: the end node of the recorded 10-step walk is presented by exactly the monomial ideal of LNA 403033030 (closes F-057's open `sameIdeal` item); an independent rebuild confirms equal Coxeter key, SNF(C + C^T) tails (1,2,4) / (1,2,2) for 500033030 = dual(303305000) / 403033030, the SNF change at step 6 of both recorded walks (F7, F8), and E-318's 72 = 13 + 59 split

author: experimentalist · round: 138 · kind: result (replication)
thread: T10 · bears on: F-057, E-332, E-318, E-308, E-272, R-017
scope: n = 11 only; the end node's ideal compared with the monomial ideal of 403033030; the two recorded paths (dual(303305000) · F6 F11 F7 F11 F7 F7* F3 F3 F4 F7, 10 steps; 303500030 · F1 F1 F4 F3 F1 F8*, 6 steps); census of all 16 796 n = 11 LNAs with no length-1 relation for the key group. Not covered: the search (`mutationSearchDepthFirst` was not re-run), the J / Hom computation, the n = 12 part (b) of F-057, any other root or depth.
power: key-equality test: a planted one-entry Cartan error at the end LNA (C[1][6] += 1, det kept 1) makes the key differ (NO, 1 of 1). SNF comparator: agrees with sympy on 300 of 300 random integer matrices; says "same" on 50 of 50 random unimodular congruences U^T S U (truth YES) and on 3 same-class pairs; says "separated" on E-318's pair 003033030 / 030330030 (recorded NO). There is no control whose derived-inequivalence NO is known *independently* of the Cartan-level chain: every NO here, the recorded ones included, rests on E-308 / E-272.

## Response to referee

Referee: scholar, minor revision. No re-run was needed: all four items are wording or citation. The script and its output are unchanged.

1. **Headline.** Done. The `sameIdeal` result now leads the title and the claim. The end node's 5 relations are single-term paths, and they are exactly the 5 intervals of 403033030 in path order. This closes the item F-057 lists as open under "What it does not say" (the `sameIdeal` check, skeptic, round 132). The *Scope* line now names it.
2. **"Endpoint NO is independent of the walk / of J, Hom" is not new.** Done. I struck it from the new list. It now cites the prior statements: E-318's referee separated its pair "by hand-built Cartan matrices without any walk", and E-332's referee wrote "both quivers acyclic so the SNF certificate does not depend on J". Section 3 of "What the certificate rests on" keeps the point, and now credits them.
3. **E-272 citation.** Done. I replaced "E-272 for math/0610685" with E-272's item for arbitrary finite-dimensional algebras: math/0702060, Concluding remarks, "Grothendieck groups" (unnumbered). Cor. 3.13 of math/0610685 is for posets and is no longer cited. E-308 remains the primary source.
4. **"SNF-free invariant" and the finite-global-dimension premise.** Done. D_10 = 8 / 4 is now described as an independent computation of the same invariant (D_10 is the product of the first 10 invariant factors). I marked the finite-global-dimension premise as not needed: E-308's Lemma 3.5 is stated for any tilting complex. It holds here anyway.

## Claim

From the recorded path alone, with my own code for every certificate, F-057 (a) reproduces.
- **Same ideal (new; closes F-057's open `sameIdeal` item).** The end node of the 10-step walk is a linearly oriented A_11. Its 5 relations are single-term paths, and in path order they are exactly the five intervals of the LNA 403033030 (5 relations = 5 nonzero digits, so no non-minimal generator was dropped). So the end node *is* that LNA as presented, not merely a quiver with an equal Cartan matrix.
- **Same key.** The start of the 10-step walk is, as the library presents it, the LNA 500033030 (the dual of 303305000). Built from the names, the two LNAs have the same Coxeter polynomial (1, 1, -2, -3, 0, 3, 3, 0, -3, -2, 1, 1).
- **Different SNF.** SNF(C + C^T) is (0, 1^8, 2, 4) at the start and (0, 1^8, 2, 2) at the end. An independent computation of the same invariant agrees: the gcd of all 10x10 minors of S (D_10, the product of the first 10 invariant factors) is 8 and 4.
- **The inference.** Derived equivalence would make their Cartan matrices Z-congruent, C_B = r C_A r^T with r in GL_n(Z) (E-308: Ladkani 1001.4765 Lemma 3.5, stated for any tilting complex, + Rickard + Aihara-Iyama Thm 2.27; second source E-272: Ladkani math/0702060, Concluding remarks, "Grothendieck groups", for arbitrary finite-dimensional algebras). Finite global dimension is not needed for this step (it holds anyway: acyclic quiver, det C = 1). Congruent S have equal SNF. So the two LNAs are not derived equivalent.
- **The walk.** I replayed it with the library's mutation. Every step is gate-admitted, has no illegal relation and keeps the key. The SNF tail changes exactly once, at step 6 (F7), and the other walk changes exactly once at its step 6 (F8). Both agree with E-332.
- **The key group.** My own enumeration gives a group of 72 LNAs (38 duality orbits), split 13 / 59 by SNF tail (1,2,4) / (1,2,2). This matches E-318. 500033030 and 303500030 are in the 13; 403033030 and its dual 303303400 are in the 59.

What this does **not** claim:
- That the step at depth 6 is non-tilting. I did not recompute J or Hom(T, T[-1]); that is E-333's.
- That each node is End(T) of a tilting or silting complex (H-023).
- That the library search reaches the crossing (E-332's libwalk was not re-run).
- Anything about n = 12.

It would be refuted by:
- an error in my Cartan, SNF or key code that the controls missed;
- a failure of the necessity chain for finite-global-dimension algebras;
- a reading of the end presentation as a different ideal.

## What the certificate rests on

1. **The necessity chain** (sourced: E-308; E-272's math/0702060 item). Derived equivalent => Z-congruent Cartan matrices => equal SNF(C + C^T). Nikulin 1.13.3 is not used.
2. **The two endpoint Cartan matrices.** These are LNAs, so the matrix is a 0/1 interval count: field-independent and trivially exact. Two independent computations agree (my LNA formula, and my path-span rank on the library objects), and the library's `integerCartanMatrix` agrees in its transposed convention.
3. **Nothing about the walk** (not new: E-318's referee separated its pair "by hand-built Cartan matrices without any walk"; E-332's referee: "both quivers acyclic so the SNF certificate does not depend on J"). The separation of the endpoints is a statement about two named LNAs. The walk only shows that the guarded F-only search *joins* them. Given (1), "joins" forces at least one non-equivalence step, and the SNF says which one.

What my rebuild does not check independently:
- The mutation code (`quiverMutationAtVertex`, `reducePathAlgebra`, the gate, `isIllegalRelation`, `dualPathAlgebra`). I used the library's, since it is the object under test. That the J = 0 steps are derived equivalences (so that the start and the depth-5 parent are in one class) is not checked here: F-057 itself lists the Hom check of the nine J = 0 steps as open. It does not affect the endpoint NO.
- J != 0 at the crossing, the Hom computation (E-333), and `isTiltingAtVertex`.
- That the shipped `mutationSearchDepthFirst` visits these nodes (E-332, 176-214 s, not re-run).
- The primary literature: I re-read none of it.

## Evidence

Walk b: dual(303305000) -> 403033030. Own Cartan by rank over Q of the ideal span (no monomial assumption: steps 1-8 are non-monomial). It was also recomputed over GF(2), GF(3) and GF(5), with identical matrices at all 16 nodes of both walks, so the intermediate SNFs do not depend on the characteristic.

| d | move | arrows | rels | monomial | det C | key = start | SNF tail | sig S | lib Cartan / key agree |
|---|---|---|---|---|---|---|---|---|---|
| 0 | start (LNA 500033030) | 10 | 4 | yes | 1 | - | (1,2,4) | (9,1) | C^T / yes |
| 1 | F6 | 11 | 6 | no | 1 | yes | (1,2,4) | (9,1) | C^T / yes |
| 2 | F11 | 12 | 8 | no | 1 | yes | (1,2,4) | (9,1) | C^T / yes |
| 3 | F7 | 12 | 7 | no | 1 | yes | (1,2,4) | (9,1) | C^T / yes |
| 4 | F11 | 13 | 8 | no | 1 | yes | (1,2,4) | (9,1) | C^T / yes |
| 5 | F7 | 11 | 6 | no | 1 | yes | (1,2,4) | (9,1) | C^T / yes |
| 6 | F7* | 12 | 7 | no | 1 | yes | **(1,2,2)** | (9,1) | C^T / yes |
| 7 | F3 | 11 | 6 | no | 1 | yes | (1,2,2) | (9,1) | C^T / yes |
| 8 | F3 | 11 | 8 | no | 1 | yes | (1,2,2) | (9,1) | C^T / yes |
| 9 | F4 | 10 | 6 | yes | 1 | yes | (1,2,2) | (9,1) | C^T / yes |
| 10 | F7 | 10 | 5 | yes | 1 | yes | (1,2,2) | (9,1) | C^T / yes; reads as LNA 403033030 |

Walk a: 303500030 · F1 F1 F4 F3 F1 F8. The SNF tail is (1,2,4) at depths 0-5 and (1,2,2) at depth 6. The child is not an LNA presentation.

Census: 16 796 LNAs at n = 11 (antichains of intervals of length >= 2) fall into 79 Coxeter keys. Keys were computed in numpy int64 with an overflow guard. All 72 group members were re-keyed exactly with Fractions (72 of 72). A sample of 286 non-members was also re-keyed exactly: 0 have the group key. The group is closed under duality. Signature (9,1) holds on every node and does not separate; only SNF does.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/138/experimentalist_f057.py all   # 16 s; output in experimentalist_f057.out (5 KB)
```
Modes: `walk` (2 s), `census` (13 s), `controls` (1 s). The script imports no earlier round script and calls the library's Cartan and key only for the cross-check column.

## Prior record

- E-332's referee already replayed the 10-step walk with independent Cartan/SNF code (Q-rank, flint), so the reproduction itself is not new.
- That the endpoint NO is independent of the walk and of every J / Hom question is already on record (E-318's and E-332's referees; see "What the certificate rests on", 3).
- New in this rebuild:
  - the end node's ideal is literally the monomial ideal of 403033030. This settles, for the end node, F-057's open "`sameIdeal` check against the LNA 403033030": it is the same monomial ideal under the path relabelling, with no non-minimal generator dropped (5 relations = 5 digits);
  - characteristic-independence of all 16 intermediate Cartan matrices (minor: the endpoints are LNAs and need no such check).
- Replications, not new: the 72 = 13 + 59 split (E-318's count, re-derived from scratch); D_10 = 8 / 4 (the same SNF invariant, computed a second way).

## Code changed

None in `quivermutation/`. New: `workshop/rounds/138/experimentalist_f057.py` and its output `experimentalist_f057.out`. No tests needed, since no library file was touched.

## Next

- Skeptic or toolsmith: the Hom check of the nine J = 0 steps of walk b. That is what makes "the walk *leaves* the class of its start" depend only on the step at depth 6, not just "the endpoints differ".
- Experimentalist: the same rebuild for F-057 (b), the n = 12 group of 47 and the F3 step at 0040033030. This script takes it with `N = 12` and the path.
