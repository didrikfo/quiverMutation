# At n = 9 the quipu theorem plus Hom-verified lnaMoves certifies 2400230 ~ 3304000 without Cartan data, so E-226's non-join is a move-set gap; power table 1002 / 32 / 0 over 380 components (LNAs only)

author: skeptic · round: 085 · kind: result
thread: T10 · bears on: F-010, F-047, E-219, E-226
scope: n = 9 (soundness count repeated at n = 8), LNAs only (not tilted non-LNA algebras); certificate = relation-position combinatorics (thm:QuipuToAn, 2305.06642) carried over one-step `lnaMoves` edges; 380 components of the move graph (the forward-closure `orbitOf`, unioned), not F-047's 48 orbits; 7 149 key-equal orbit pairs stay unnamed.

## Response to referee

1. Title and Claim narrowed. The inequivalence of 3060000 and 3304000 is not new: F-010, E-057 (Ladkani Cor 3.15, EXPERIMENTS lines 452-454), F-047 (SNF split); the quipu-name route is F-003's inversion. Only (a) the certified equivalence 2400230 ~ 3304000 (theorem name directly on one side, theorem name through moves on the other) and (b) the 1002 / 32 / 0 power table are claimed as new.
2. "False negative of the E-212 search" replaced by "move-set gap". I do not claim E-212 was meant to join every equivalent key-equal LNA; what is shown is that a certified equivalent pair is not joined by the J = 0 key-kept set at 6 + 8.
3. Per-edge status of the 3304000 component (7 rows: 3304000, 3004000, 3024000, 3004002, 3034400, 3040400, 3044030; 8 directed edges). Every edge is a `lnaMoves.movesByRule` rewrite (rule-level: `verifyMove` legal + result + Coxeter, theory-side) and now also Hom-verified step by step (Hom(T,T[-1]) = Hom(T,T[1]) = 0, End(T) = Cartan(child)^T, S = -1, as `replay` of round 079) with the final algebra standardising to the predicted class: 8 / 8 pass (`skeptic_hom.py comp`, `skeptic_hom_comp_out.txt`). Wider: all 2 935 move edges of all 1 430 n = 9 rows, 0 failures (`skeptic_hom.py all`, 88 s, `skeptic_hom_all_out.txt`). So the move graph is Hom-verified at n = 9, not only the one edge of E-226. Only the naming of the three theorem-named members rests on the theorem (one code path with the tree check).
4. Tree-isomorphism output saved: `skeptic_tree.py`, `skeptic_tree_out.txt` (canonicalTreeForm differs, adjacency spectra and degree sequences equal, networkx is_isomorphic False). Names and trees are one fact (same code path), as stated in Evidence; "different names" and "non-isomorphic trees" are not two confirmations.
5. Extra soundness width: `skeptic_cert.py 8` (1 m): 429 rows, 69 components, 0 two-name components, 10 unnamed, key-equal pairs 165 same / 0 different / 73 unnamed (`skeptic_cert8_out.txt`). n = 10 not run (cert cost grows about 5x per n; over the 10-minute limit).

## Claim

There is a certificate that uses no Cartan, Coxeter or Hom data and no search. Row `3060000` names to the quipu P^(1,0,1)_(1,4) (theorem directly). Row `3304000` has overlapping relations, so the theorem is silent on it, but it is joined by `lnaMoves` to `3004000` and `3024000`, which name to P^(1,1,2)_(1,2). The two trees are not isomorphic (canonicalTreeForm differs; same degree sequence and same adjacency spectrum, which is F-010). By Happel (derived-equivalent path algebras of trees have the same tree) the pair is inequivalent. Power: in the same run, 1 002 key-equal pairs of move components carry the same name (equivalent) and 32 carry different names (inequivalent), and no component holds two different theorem names (0 conflicts: the move graph and the theorem never disagree on 1 430 rows). The key-equal equivalent control for the 3304000 side is `2400230` (names to P^(1,1,2)_(1,2) directly): so 2400230 ~ 3304000 is certified, not "an inference from naming" as E-226 says. Consequence: the 6 + 8 no-meet of E-219/E-226 is a gap in the E-212 J = 0 key-kept move set on a certified equivalent key-equal pair (not a claim that the search was meant to join all such pairs), so that search has NO demonstrated power for a pair that has 3304000 as an end beyond its own move component and its walk ends. Not claimed: that the certificate says anything about the J = 0 premise (it lives on LNAs; it is blind to what a J != 0 step does), or that it works for the 7 149 unnamed pairs.

Refutation: a pair with equal names that is provably inequivalent, a move edge that changes the theorem name, or a tree-isomorphism between Q1 and Q2.

## Evidence

- Names by `quipuForAlmostSeparateLNA(9, row)`: 3060000 ((1,0,1),(1,4)); 6000030, 2223030, 3030000 the same; 3304000 None; 3004000, 3024000, 3004002, 2400230 ((1,1,2),(1,2)).
- Move component of 3304000 has 7 rows, theorem-named members 3024000, 3004000, 3004002, all one name; 2400230 is a one-row component (no move touches it), named directly.
- Whole table, n = 9, 1 430 rows: components 380, unnamed 208, conflicts 0; key-equal component pairs: same name 1 002, different name 32, unnamed 7 149 (pairs with a side that no member of the component can name; these carry no certificate, so "7 149" is the honest limit of this method).
- Non-Cartan: input is only the start/length of relations. It is not independent of the theorem's inversion in `quipuForms` (a single code path for all names): soundness evidence is the 0 conflicts, not a second implementation. The move edges are `lnaMoves` rewrites; all 8 edges of the component, and all 2 935 edges at n = 9, are now Hom-verified (Response 3).
- Other candidate non-Cartan invariants are out: HH is k on all LNAs n <= 10 (E-172), Serre-periodicity is blind on the F-010 pair (E-172), AAG closed (R-008). So the theorem name is the only non-Cartan certificate the repo has; I could build no certificate for non-LNA (tilted) algebras and none for the J != 0 children of E-233/E-234.
- On the NO: the E-226 sentence "7 + 7 not needed" and the 3060000-side power (joins 6000030, same-profile, certified equivalent here too) survive; the 3304000-side power is limited to its 7-row component and walk ends. A NO is only as strong as the weakest side; with 2400230 un-joined at 6 + 8 the 3304000 ball is not a closure of its derived class, so the inequivalence of the pair is established by the theorem and the J = 0 NO adds nothing to it. A real NO test would first need the search to join 3304000 to 2400230 (then the 3060000 vs 3304000 meet would have power both sides).

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/085/skeptic_cert.py 9     # 5 min 18 s; output in skeptic_cert_out.txt
```
Hom: `.venv/bin/python -W ignore workshop/rounds/085/skeptic_hom.py comp` (2 s) and `... all` (88 s). Tree: `.venv/bin/python workshop/rounds/085/skeptic_tree.py`. Tree check: `quipuForms.graphFromQuipuParameters((1,0,1),(1,4))` vs `((1,1,2),(1,2))` through `canonicalTreeForm` (False), `adjacencySpectrumPolynomial` (equal).

## Prior record

F-010 already states the pair is inequivalent (two cospectral quipus; certified E-057), F-003 the inversion, F-047 the SNF split and the two orbits, E-226 the 2400230 mismatch as unexplained. New: (1) the explanation is a move-set gap, not a naming error: 2400230 is theorem-named and 3304000 is theorem-named through a verified move; (2) the soundness count (0 conflicts over 1 430 rows) and the power table; (3) the 48 orbits of F-047 and my 380 components differ in granularity (I used `orbitOf` unions, no mirror, no search merges); the same-name / different-name counts are for components, not F-047 orbits. Nothing here is in RETRACTIONS.

## Code changed

None in the library. New: `workshop/rounds/085/skeptic_cert.py`, `skeptic_cert_out.txt`. No tests run (no library file touched).

## Next

Toolsmith: why does the E-212 J = 0 key-kept move set not reach 2400230 from 3304000 (2400230 is isolated under `lnaMoves`; a join needs a non-LNA-move chain: ask for the shortest Hom-verified chain via the theorem for 3004000 -> 2400230). Experimentalist: a union-find with the mirror to compare components with F-047's 48 orbits. Skeptic: a non-Cartan certificate for tilted non-LNA algebras does not exist in the repo; the J = 0 premise cannot be tested by a name.
