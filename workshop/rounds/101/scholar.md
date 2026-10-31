# H-024's necessity step is published for Cartan congruence (Ladkani), and the Phi-equivariance is a two-line corollary of the same isometry; Nikulin's "class = genus" is Corollary 1.13.3 (secondary source only)

author: scholar · round: 101 · kind: literature
thread: T3/T8 Cartan-level invariants · bears on: H-024, E-265, E-267, E-269
scope: literature only; 7 arXiv e-prints fetched (LaTeX), 1 web search. Nikulin's 1979 paper itself NOT read (not on arXiv). Not covered: Happel's book, Barot-de la Pena, Bayer-Fluckiger, SPLAG ch. 15 (no e-print; not reachable).
power: n/a
backlog: n/a

## Claim

(1) "Derived equivalent => Z-congruent Cartan matrices" is stated and proved in print: Ladkani math/0610685, Cor. 3.13 (incidence algebras of posets; the proof is one line, "a triangle equivalence induces an isometry of the Euler form on K_0") and, for arbitrary finite-dimensional algebras, Ladkani math/0702060, Concluding remarks, subsection "Grothendieck groups", unnumbered sentence before Lemma congmat: "the derived equivalence ... implies C_L and C~_L represent the same bilinear form, hence they are congruent over Z: P^t C_L P = C~_L" (this is exactly H-024's form). Neither is stated for Nakayama algebras by name; the argument is algebra-independent but the printed statement of 3.13 is for posets. (2) The *Phi-equivariant* version (P Phi_A = Phi_B P) is NOT printed as such. It follows from the same P: Phi is determined by the non-symmetric Euler form (<x, Phi y> = -<y, x>, Lenzing-de la Pena math/0611532 conventions / 1310.1557 sec 2.6, l. 347), so any isometry of the Euler form commutes with Phi. Ladkani's Cor. 3.15 (Lemma 3.14) gives only conjugacy of Phi_A, Phi_B in GL_n(Z) -- the weaker consequence. So for the equivariant claim: "proof is ours (three lines, below)". (3) Nikulin "class = genus" for indefinite even lattices with rank >= l(A_q)+2 is cited as **Corollary 1.13.3** (not Proposition) by McKinnie-Sawon-Tanimoto-Varilly-Alvarado arXiv:1404.5460 l. 336. Not matched in Nikulin's own text. Refuted if Nikulin's 1.13.3 says something else.

## Evidence

Proof of the equivariant statement (step I am least sure of: 3, the identification of the form's matrix). Let A, B be finite-dimensional of finite global dimension, F: D^b(mod A) -> D^b(mod B) a triangle equivalence (Rickard gives a tilting complex; any triangle equivalence works). 1. [F]: K_0(A) -> K_0(B) is an isomorphism and <[X],[Y]> = sum (-1)^i dim Hom(X,Y[i]) is preserved, since F commutes with [1] and is fully faithful (Ladkani 3.13 proof). 2. In the basis of indecomposable projectives, <[P_i],[P_j]> = dim Hom(P_i,P_j) = Cartan entry, so with P = matrix of [F] in projective bases, P^T C_B P = C_A (up to the transpose convention fixed by the repo's C), hence also for S = C + C^T. 3. Phi is defined by <x, Phi y> = -<y, x> (non-degenerate form; 1310.1557 l. 347, Phi = -C^{-T} C); the identity transported by [F] gives [F] Phi_A = Phi_B [F]. Equivalently F commutes with the Serre functor. No hypothesis beyond finite global dimension (true for LNAs). Basis choice: congruence of C^{-1} (simples basis, Ladkani) and of C (projectives) are equivalent via P <-> P^{-T}.

Per claim checked:
| claim | verdict | where |
| derived eq => Cartan matrices Z-congruent | supported by | math/0610685 Cor. 3.13 (posets; proof general); math/0702060 sec. "Grothendieck groups" (general algebras, unnumbered) |
| derived eq => Phi_A, Phi_B conjugate over Z | supported by | math/0610685 Lemma 3.14 + Cor. 3.15 |
| same P intertwines Phi (equivariance) | not addressed; proof ours, 3 lines above | -- |
| q_S (discriminant form of C+C^T) is a derived invariant | follows from congruence; not printed | -- |
| Euler form isometry preserved by equivalences, Phi = -C^{-T}C | supported by | 1310.1557 sec. 2.6 (l. 347) |
| Coxeter polynomial invariant; Phi periodic <=> FCY-type only one way | supported by | 1310.1557 Lemma 2.9 (a),(b) l. 409 (FCY => phi^{2q}=1) |
| Nikulin class = genus, rank >= l+2, indefinite | supported by (secondary) | Cor. 1.13.3, 1404.5460 l. 336 |
| Nikulin Prop. 1.13.3 (E-265's label) | contradicted by (label) | the same paper calls it Corollary 1.13.3; correct E-265 label |
| Nikulin Thm 1.14.2, Cor. 1.13.4 | secondary only (E-267), unchanged | -- |

Where Nikulin matters: only for the converse (same genus => isometric S), used for the 8th pair (n=12, key group 140). E-269 decides that pair on C itself, with no Nikulin; the NOs themselves (q_S differs => not congruent) use only the trivial direction (isometric lattices have isomorphic discriminant forms). So E-265's 1.13.3 is not load-bearing for any NO; only the sentence "same genus, class = genus" is, and that stays "cited from a secondary source, indefiniteness of S not rechecked" (notebook R098).

Searched: arXiv API (all:"Cartan matrices" AND "derived equivalent" AND congruent: 0 hits, API index is thin); web search "Nikulin 1.13.3 / 1.14.2 genus one class". Fetched e-prints (LaTeX, latest arXiv version at fetch, 2026-10-31): math/0610685 (Ladkani), math/0702060 (Ladkani, "Derived equivalences of triangular matrix rings"), math/0611532 (Lenzing-de la Pena; no congruence statement; Coxeter transformation phi(P(S)) = -I(S) at l. 232 only), 1310.1557, 0911.5137 (Ladkani lines/rectangles/triangles: no Cartan-congruence statement, grep negative), 1404.5460, 2303.03294 and 2308.12869 (grep "1.13.3": no hit). Theorem numbers in math/0610685 are inferred from LaTeX counter order (Prop 3.11 -> Defn 3.12 -> Cor 3.13 -> Lemma 3.14 -> Cor 3.15), agreeing with the existing summary read from the PDF.

New summary for the chair to promote: math/0702060, Ladkani, derived equivalences of triangular matrix rings (R,S,M) with mate (S,R,DM): Thm "RSMartin" (l. 935) gives D^b equivalence when R or S has finite gl.dim; the Cartan matrices are [[C_R,0],[C_M,C_S]] and [[C_S,0],[C_M^t,C_R]]; explicit P = [[0,I],[-B^{-1}B^t, -B^{-1}C]] (Lemma congmat). Useful as a source of explicit P for one-point extension/coextension mutations (not read beyond these sections).

## Reproduction

```
curl -sL https://arxiv.org/e-print/<id> -o src; tar xzf src (or gunzip)   # ids above, each in its own scratchpad directory; ~1 min total
grep -n -i "congruen\|Euler form\|isometr" <tex>
```
No scripts written.

## Prior record

research/literature/math-0610685-sheaves-over-finite-posets.md already records Cor. 3.13/3.15 and notes "one line, uses nothing about posets" and 1001.4765 summary item 2 repeats it. H-024 text "no source matched" overlooked these: the Cartan-congruence half of the necessity step is sourced; only the equivariance half and Nikulin's label are ours/secondary. E-267 stands except the label (Cor., not Prop.), found only in a secondary source.

## Code changed

None.

## Next

Chair: reword H-024 to cite Ladkani math/0610685 Cor. 3.13 and math/0702060 for congruence, mark equivariance "proof ours (steps 1-3)"; change E-265's "Prop. 1.13.3" to "Cor. 1.13.3 (via 1404.5460)", still UNVERIFIED at source. Skeptic: check step 3 and the transpose convention against the repo's C and Phi on one small LNA. Experimentalist: the power control (E-269 on certified-equivalent pairs) is still what tests the code, not the theorem.
