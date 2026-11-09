# The proofs of ENS (math/0101145v3) Cor. 3.11 and Thm. 3.10 (Move III, Sec. 6.3) give a non-completed pair-cancellation lemma for semifree dg path algebras over k^n under (T)+(P); on the seven formal loop models it reproduces E-290's remainder and E-309's letter counts and shows M4 keeps 3 Theta letters; BHL 0906.3422v4 Prop. 2.1 is a second source for H-024's congruence step

author: scholar · round: 122 · kind: literature
thread: T5 (H-031, H-023); T3/T8 (H-024) · bears on: H-031, H-023, H-024, E-290, E-294, E-309, E-315, E-272, E-308, E-265, E-267
scope: sources read in LaTeX: ENS math/0101145v3 (Sec. 3.3, 6.3), Kalck-Yang 1803.08192v1 Sec. 2, Van den Bergh 1008.0599v3 App. A, Bastian-Holm-Ladkani 0906.3422v4 Sec. 2.1 and 1012.4661v2 Sec. 1-2, Dolgachev alg-geom/9502005 Sec. 5; skimmed by grep: Keller 0908.3499v5, Kalck-Yang 1205.1008v4, Manolescu-Rozenblyum 2605.08122v1 Sec. 2. Chekanov math/9709233 is PostScript only (no gs here): not read. The ENS hypotheses are tested only on the formal models M1-M6 (incl. M3, M3b) of `workshop/rounds/106/theorist_loopfilt.py`, not on a Construction-2.2 parent. ENS works over Z[t,t^-1] with one vertex; the transfer to K = k^n is ours (argued below, not machine-checked). Nikulin 1.13.3: one more secondary source, no replacement.
power: (T) triangularity says NO on a planted cycle in M1 (d a* := f/a: refused); (P) says NO on a planted pair d t = b + xz with b below x (refused). YES controls: the E-309 substituted M4 and M5 algebras (known direct sums) are triangular and cancel exactly their 3 trivial pairs, leaving 8 and 9 letters. Theorem numbers by `workshop/rounds/110/scholar_thmnumbers.py` (ENS, BHL use [section] counters); planted "9.9" prints nothing on each file.

## Claim

(1) **H-031, non-completed source.** In the theory of Legendrian contact homology, Etnyre-Ng-Sabloff (math/0101145v3, after Chekanov) work with semifree dg algebras (free on finitely many generators; no completion). Def. 3.8 defines elementary and tame isomorphisms: one generator goes to ±itself + u, with u in the other generators. Def. 3.9 defines stabilisation S(A) = A * E_i, with d e1 = e2. Cor. 3.11 says that projecting S(A) to A is a chain homotopy equivalence; the homotopy H is written out explicitly. The proof of invariance under Move III (Sec. 6.3, Lemmas 6.3, 6.4 and the maps g_i) is purely algebraic: if the generators are ordered so that d of each lies in the subalgebra of earlier ones (*triangularity* (T)), and d a = ±b + v with v in generators below b (*pair condition* (P)), then finitely many elementary automorphisms g_i(a_i) = a_i + H(d'a_i - Phi_{i-1} d a_i) give a tame chain isomorphism A ≅ S(A'). Then A -> A' = A/(a, da) (substitute a -> 0, b -> -v) is a quasi-isomorphism. This is the non-completed form of H-031's steps 1 and 2. The substitution is *constructed*, not assumed, and the positive-weight-per-cell worry of E-315 is replaced by the finite combinatorial condition (T). (2) **On the models.** Raw kM satisfies (T) in all of M1, M2, M3, M3b, M4, M5 and M6. Iterated ENS cancellation, with d^2 = 0 checked at each step, gives the following. M1 cancels 2 pairs and leaves {a*, af, f} with no Theta (= E-290). M4 cancels 3 pairs and leaves 8 letters, 3 of them Theta. M5 cancels 3 pairs and leaves 9 letters, 3 of them Theta (f/m, af/m, mf/m). The 8 and 9 equal the reduced parts of E-309's substituted algebras. M2, M3 and M6 have no pair. M3b cancels 2 pairs and leaves the Thetas g/a and ag/a, with a cubic term a*·f·f. So the only step H-031 still needs is Lemma A on the reduced part. **Correction to E-315 / round 121:** M4 is *not* "all Theta in C" like M1. Only 3 of its 6 linear Theta letters can be cancelled (3 targets), and 3 stay in the reduced part, as E-309's Lemma-A step for M4 already implied. (3) **H-024.** Bastian-Holm-Ladkani 0906.3422v4 Prop. 2.1 states: for finite-dimensional derived-equivalent A and B there is P in GL_n(Z) with P C_A P^T = C_B. It needs no global-dimension hypothesis: the proof uses per A, Rickard and the Euler form. Their Prop. 2.2 gives the chain (a) congruence => (b) Z-similarity of S = C C^-T => (c) Q-similarity => (d) equal characteristic polynomial. Prior art for a Cartan-level NO: 1012.4661v2 Example 2.31. There, two cluster-tilted algebras of type D_15 have the same associated polynomial, but their asymmetry matrices are not similar over F_3, "so the two algebras are not derived equivalent". That is the same certificate shape as E-265/E-292. **Not claimed:** that ENS's statement as printed covers kM (base ring Z[t,t^-1], one vertex; the transfer is ours); anything about Lemma A, Hom/End(T) or a Construction-2.2 parent; a Nikulin replacement. Refuted by: a loop model meeting (T) and (P) whose ENS quotient has H different from E-284's H(N); or an error in the k^n transfer below.

## Evidence

Mechanical test (`scholar_ens.py`: (T) checks that the digraph g -> letters of d g is acyclic; (P) as above; the cancellation is the substitution theta -> 0, b -> -v/c, repeated while a pair exists):

| model | letters | (T) | pairs cancelled (theta -> b) | left | Theta left | d^2 = 0 each step | record |
|---|---|---|---|---|---|---|---|
| M1 | 7 | yes | af/a->ab, f/a->b | 3 {a*, af, f} | none | yes | E-290: same remainder |
| M2 | 7 | yes | none | 7 | f/a, af/a | - | E-294 (Lemma A) |
| M3 | 11 | yes | none | 11 | 4 | - | open (M3) |
| M3b | 11 | yes | af/a->ab, f/a->b | 7 | g/a, ag/a (cubic word a*·f·f) | yes | open (M3b) |
| M4 | 14 | yes | a1f/a1, a2f/a1, f/a1 | 8 | a1f/a2, a2f/a2, f/a2 | yes | E-309 reduced part: 8 letters |
| M5 | 15 | yes | af/a, f/a (v: 2 words), mf/a | 9 | f/m, af/m, mf/m | yes | E-309 reduced part: 9 letters, 3 Theta |
| M6 | 8 | yes | none | 8 | m/a | - | E-235 |
| E-309 M4 / M5 (YES control) | 14 / 15 | yes | 3 / 3 | 8 / 9 | - | - | direct sums |
| plant: cycle in M1 | | NO | | | | | refused |
| plant: d t = b + xz, b below x | | | NO pair | | | | refused |

**The transfer to K = k^n, the step I am least sure of.** ENS's algebra is free over the commutative ring Z[t,t^-1]. Ours is T_K(V), a path algebra of a graded quiver over K = k^n. Every ingredient is word-level. The homotopy of Cor. 3.11 is H(w) = ±w1 e1 w2 at the first e2 of w, where e1 and e2 are parallel arrows. The elementary maps send an arrow to ± itself plus parallel paths. The induction uses only (T), d^2 = 0, tau H = 0 and tau g_i = tau. None of these uses commutativity of the base or a single vertex. Kalck-Yang 1803.08192v1 Lemma 2.1 (`lem:deleting-a-contractible-vertex`) is a non-completed dg quiver statement with a homotopy of the same word-level kind (when d x = e_i, killing vertex i is a q.i.); it supports the transfer but is not itself a pair cancellation. ENS's Lemma 6.4 is geometric in their setting; for us it holds by definition, since d' := sigma∘d is the differential of A/(a, da). The sign of H ("it is not hard to check", unproved in ENS) is not checked by me.

Other sources, no statement found: Van den Bergh 1008.0599v3 App. A. Its minimal models (Cor. A.5.5) are for pseudo-compact algebras, so it is completed. Prop. 7.2.1 uses the non-completed cofibrant form, T_l V with d(F_{i+1}V) in T_l F_i V, which is condition (T). Kalck-Yang 1205.1008v4 (complete path algebras, A_infinity minimal models) and Keller 0908.3499v5 had no non-completed cancellation lemma that grep could find. Manolescu-Rozenblyum 2605.08122v1 Sec. 2 cites "stable tame isomorphic ⇒ quasi-isomorphic" as [ENS, Cor. 3.11]; that is how I found it.

H-024 numbers: 0906.3422v4 Prop. 2.1 and Prop. 2.2 (no labels; counter); Def. 2.3 (associated polynomial). 1012.4661v2 Prop. 1.3 (same statement, citing the proof of Bocian-Skowroński Prop. 1.5 (UNVERIFIED, not on arXiv) and "[BHL09, Prop. 2.6]"). Prop. 1.4 `prop-asymmetry`, Ex. 2.31 `ex:D15` (NO via F_3), Ex. 2.21 `ex:D15_op` (Z-congruence found by Magma), Ex. 2.30 `ex:D19` (equivalent forms, derived equivalence unknown). **Discrepancy:** 1012.4661 cites the E paper's proposition as 2.6, but in 0906.3422v4 it counts as 2.1. Probably the journal version (Algebr. Represent. Theory 16 (2013) 527-551) is numbered differently; not checked. Cite "v4 Prop. 2.1". Journal data are from the arXiv API: J. Algebra 410 (2014) 277-332 for 1012.4661.

Nikulin: Dolgachev alg-geom/9502005 Prop. (5.1) cites Nikulin Cor. 1.12.3 and Thms 1.12.4 and 1.14.4 for embeddings. It has no class = genus statement. Cor. 1.13.3 stays UNVERIFIED; E-267 (28 e-prints) plus this one make 29.

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/122/scholar_ens.py      # 0.04 s: table rows, plants, controls
.venv/bin/python -I workshop/rounds/110/scholar_thmnumbers.py <math/0101145 oriented-chvII-rev.tex, \mylabel -> \label> defn:elem defn:stabilization cor:homology lem:height lem:restr-chain lem:proj-chain 9.9   # 3.8 3.9 3.11 6.2 6.3 6.4
.venv/bin/python -I workshop/rounds/110/scholar_thmnumbers.py <1012.4661 ctaD_final.tex> prop-asymmetry ex:D15 ex:D15_op ex:D19 9.9   # 1.4 2.31 2.21 2.30
.venv/bin/python -I workshop/rounds/110/scholar_thmnumbers.py <0906.3422 ctaE.tex>   # rows 2.1, 2.2 prop, 2.3 dfn
```
Fetched under the STEERING arXiv permission into the scratchpad (not committed): e-prints 1008.0599, 1205.1008, 1803.08192, 0908.3499 (round-121 leads), 2605.08122 and math/0101145 (found by an API search for "stable tame isomorphism"), math/9709233 (PostScript, unreadable here), 0906.3422 and 1012.4661 (H-024 prior art), alg-geom/9502005 (Nikulin). For the chair to commit if the record cites them: math/0101145, 0906.3422, 1012.4661.

## Prior record

A grep of research/, GLOSSARY and the notebooks for Etnyre, Chekanov, 0101145, tame isomorph, Legendrian, Bastian, 0906.3422 and 1012.4661 finds no hit. E-315 (no completed source covers H-031) stands. This adds a non-completed source with different hypotheses, and corrects E-315's M4 classification. E-308/E-272 sourced H-024's necessity through Ladkani + Aihara-Iyama. BHL is a one-line, single-source version with no global-dimension caveat (math/0702060's note says "needs finite global dimension, unchecked"), plus published prior art for the NO direction. The equivariance with the same P is BHL's (a) => (b) computation: S_B = P S_A P^-1.

## Proposed literature summaries (for the chair)

`research/literature/math-0101145-ens-stable-tame-isomorphism.md`: **Etnyre-Ng-Sabloff, Invariants of Legendrian knots and coherent orientations** (v3; J. Symplectic Geom. 1 (2002) no. 2, 321-367, as printed in 2605.08122's bibliography). Takeaways: Def. 3.8 (elementary/tame), Def. 3.9 (stabilisation), Cor. 3.11 (S(A) -> A homotopy equivalence, explicit H), and Sec. 6.3, which bootstraps to a tame isomorphism A ≅ S(A') under height order and d a = ±b + v with v below b. Caveats: Z[t,t^-1], one vertex; the k^n transfer is ours; Chekanov's original (Invent. Math. 150 (2002) 441-483) is unread.
`research/literature/0906.3422-bhl-cluster-tilted-type-e.md` (with 1012.4661): derived equivalence ⇒ Z-congruent Cartan forms for any finite-dimensional algebras (Prop. 2.1); the invariant chain (Prop. 2.2); a published F_3 NO (1012.4661 Ex. 2.31); the number discrepancy above.

## Code changed

New: `workshop/rounds/122/scholar_ens.py`. It imports the round-106 engine and `workshop/rounds/118/theorist_m45.py` read-only. No library change; no tests run.

## Next

- **Theorist:** H-031 can now be written as two steps, ENS cancellation under (T)+(P) (cite math/0101145 Cor. 3.11 + Sec. 6.3, with the k^n transfer as a lemma) followed by Lemma A on what remains. Check that (T) holds for every Construction-2.2 parent with a pivot-reducible loop. It held in all 7 models; one candidate reason is that d raises E-284 weight on non-star letters, but that is not checked. Check the sign of the ENS homotopy H. M3b leaves a cubic word a*·f·f: is that inside Lemma A?
- **Scholar (123):** read Chekanov's original via a PDF if the human can supply one; H-025 novelty retry; owed literature files (0906.0761, 0704.0649, plus the two above).
- **Chair:** amend E-315's M4 line (3 Theta letters remain in the reduced part; see the table).
