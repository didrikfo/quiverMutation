# Oppermann's v1 Prop. 7.6 has no Lemma A and no loop restriction; the Anick monomial lemma is Thm 3.1 and the leading-term criterion Thm 3.2 (lex) / Forré Th. 2.6 (any multiplicative order) of J. Algebra 78, all secondary via Gartner; the loop case needs a triangular, not diagonal, cancellation

author: theorist · round: 101 · kind: literature
thread: T5 · bears on: H-023, E-235, E-237, E-240, E-256, E-214, E-248
scope: arXiv:1504.02617 v1 LaTeX only (the arXiv listing shows a single version, published 2015-04-10, byte-identical to `research/literature/sources/1504.02617v1-...tex`; the Adv. Math. 307 version was not seen). Anick: not reachable; only two arXiv papers citing him read (1212.2118, math/9903030). No computation.
power: n/a

## Response to referee

Verdict was minor revision; all five required items done, no re-run needed (no computation claimed). The Gartner lines were re-fetched and re-read (1212.2118, paper.tex l. 852-857; bbl: DA = Anick, J. of Algebra 78:120-140, 1982; PF2 = Forré).

1. Anick's criterion number: corrected. Gartner l. 857: "shown in \cite{DA}, Th.3.2 in the case of lexicographic orders; Forré [PF2], Th.2.6 carries over to arbitrary multiplicative orders". So: monomial lemma Th. 3.1, leading-term criterion Th. 3.2 (lex) plus Forré Th. 2.6 (arbitrary multiplicative). All secondary via Gartner; the primary was not seen. Step (c) of the loop case needs the Forré extension, since omega is not lexicographic as stated. Gartner states it for homogeneous polynomials in a free algebra over a field.
2. Wording: "tail(phi) = i" now reads "target of phi = i" (source `\ta{φ} = i`; start is `\st`); for a loop both agree. "beta phi is again an arrow starting at i" now reads "the Step-2 composition arrow beta phi of M".
3. E-211 row repaired (split into one verdict plus a separate not-addressed row).
4. Claim (3) last sentence reworded to the v1-only scope.
5. Added a line: the repo code does not exercise loops of degree >= 1 at i (the parents are acyclic; `quivermutation/` only warns on loops, arrowPaths.py l. 107-127, procedure.py l. 265-284), so any code check is vacuous here. Also stated that the 6 968 / 10 206 figure is E-235's own control (EXPERIMENTS.md l. 299), not E-256's.

Title, claim (2) and the Anick paragraph changed accordingly; scope line unchanged except that it now names the Gartner lines read.

## Claim

(1) Oppermann's paper contains no "Lemma A": it cites neither Anick nor Halperin-Lemaire (grep of the source: 0 hits for anick, halperin, koszul, quasi-free). The cancellation step is Prop. 7.6, whose printed proof is two sentences ("up to terms psi beta^{-1} with psi of strictly smaller degree. Thus, cancelling ...", "similar" for the second family). Lemma A is therefore the workshop's own bridge, not a step of the source. (2) The number of the Anick statement behind Lemma A: **Anick, J. Algebra 78 (1982) 120-140, Theorem 3.1** (monomial sequence strongly free iff combinatorially free), as cited in Gartner, arXiv:1212.2118, l. 819-822, which uses it as "[DA], Th. 3.1"; the leading-term criterion is **Th. 3.2** (lexicographic orders) with **Forré, Th. 2.6** extending it to arbitrary multiplicative orders (Gartner l. 857). All secondary citations; the primaries were not seen, so the numbers stay "per Gartner, unchecked". (3) Oppermann's hypotheses on i are only "no loop of degree 0" (Thm 1.1; Setup l. 283; Sect. 3 l. 315), so loops of degree >= 1 at i are *inside* Thm 1.1, Def. 3.2 (Step 4, l. 354-358) and Thm 7.5; E-235's "no loop at i of any degree" is stronger than the source needs and was a limit of our own proof of Lemma A's applicability, not of the paper. This is a claim about v1 only; the journal version was not seen, so nothing here says what it states.

## Evidence (each claim checked; numbers are the source's auto-numbering, `numberwithin{theorem}{section}`, recounted from the LaTeX)

| Claim checked | Verdict | Where in the source |
|---|---|---|
| Hypotheses of the mutation theorem: non-negative grading, free on degree-wise finite generators, no degree-0 loop at i; nothing on cycles, finite dimension, J, parallel arrows (E-214, E-199) | supported by | Thm 1.1 (l. 263-265); Setup l. 283; Sect. 3 l. 315 and l. 321 ("none of these [A] is a loop") |
| A loop of degree >= 1 at i is allowed and handled | supported by | Step 1 l. 328 ("all loops in i keep their degree"); Step 4 l. 354-358 (phi: i -> i in Q_1, any degree; ||alpha phi beta^{-1}|| = ||phi||); Def. 3.2 last line; Sect. 7 l. 913 (phibar = M phi E for a loop) |
| Def. 3.2 is a differential, d^2 = 0 (E-235 "d^2 = 0 beyond E-202") | supported by (printed, I read the arrow cases phi, alpha phi, alpha^*, phi alpha^{-1}; "the proof for alpha phi beta^{-1} is the same calculation") | Thm 3.3 = Thm 8.1 (l. 1048), Prop. 8.3 (l. 1109-1140) |
| Endomorphism ring of the mutated silting object is Gamma, generated as printed; Gamma_red subset Gamma is a dg subalgebra and a quasi-isomorphism (E-217 "Gamma_red = N", W1/W2) | supported by | Constr. 7.3 (hat p), Lemma 7.4 (d hat p = hat red(dp)), Thm 7.5 (l. 991, proof l. 997-1004); presentation of Gamma_red l. 1009-1025 |
| kM (Sect. 3) is quasi-isomorphic to the quotient by phi alpha^* = -(-1)^{|phi|}(d phi)/alpha, alpha phi beta^* = ... (E-235 conclusion) | supported by as a statement, **not addressed** as a proof | Prop. 7.6 (l. 1027-1042); no hypothesis beyond Thm 1.1 |
| "d lowers ||.|| by 1" (E-235 weakest step) | supported by (one clause, no argument) | Prop. 7.6 proof: terms psi beta^{-1} have psi "of strictly smaller degree". My reading: d is homogeneous of degree -1 and every letter of a non-negatively graded word has degree <= the word's, so ||psi beta^{-1}|| < ||phi alpha^{-1}||; not stated in the source |
| Lemma A / any no-overlap lemma (E-237, E-240, E-227) | **not addressed** | Prop. 7.6 gives no cancellation lemma; no citation of Anick, Halperin-Lemaire, Koszul, Groebner |
| E-235 step 4 lands in the hypotheses of Lemma A (E-256 a) | **not addressed** by Oppermann (he never states Lemma A). Against the Anick statement: see below | -- |
| Cyclic D(i), E-211: does Thm 1.1/7.5 cover it? | supported by (no hypothesis on cycles, so End(T) = H_0(kM) holds for E-211 and its repo failure is step 7, as E-214/E-220 say) | Thm 1.1 |
| What the repo's step 7 does on E-211 | **not addressed** (source has no step 7) | -- |
| J = 0 premise / derived equivalence of H_0(kM) with the parent | **not addressed** (silting vs tilting is Aihara-Iyama, E-220) | -- |

**Anick, matched to the statement E-235/E-240 use.** In Gartner (1212.2118): Def. "combinatorially free" l. 804-815 (i) no sequence member is a proper submonomial of another, (ii) for any i, j (i = j allowed) and factorisations rho_i = x_i y_i, rho_j = x_j y_j, x_i != y_j (the beginning of one is never the ending of another); Thm l. 819 "strongly free iff combinatorially free", proof "[DA], Th. 3.1" with [DA] = Anick, J. Algebra 78 (1982), 120-140 (the bibliography entry). Then **"Anick's criterion"** l. 852-857 (Th. 3.2 for lexicographic orders, Forré Th. 2.6 for arbitrary multiplicative orders, both secondary): if the high terms (w.r.t. a multiplicative order) of homogeneous polynomials rho_i are combinatorially free, the sequence rho_i is strongly free. Two consequences for us:
- Condition (ii) is the "common beginning and ending" of E-237; it is the real hypothesis, and it is stated for free algebras over a field. Over R = k^n with composable paths, "beginning = ending" is automatic only for the monomial reading with matching vertices; nothing in either source addresses a semisimple base (E-237 mismatch 1 stays open).
- The criterion is *leading-term*: the relations need not be monomial. So a non-monomial cancellation such as the loop case is exactly what the criterion is designed for, provided an order exists in which the leading words are the no-overlap ones. Not cited by E-235/E-240, which proved only the monomial case directly.

**What the "loops of degree >= 1 at i" case needs** (my reading of Def. 3.2, l. 383-395; the weakest point). For a loop phi at i, target(phi) = i, so partial(phi) contains sum_beta beta^*(beta phi) (and beta phi is the Step-2 composition arrow of M). Then partial(phi alpha^{-1}) = partial(phi) alpha^{-1} + ... contains, besides the relation word phi alpha^*, the terms beta^* . (beta phi alpha^{-1}) with ||.|| = ||phi alpha^{-1}|| - 1 of the Step-4 generator beta phi alpha^{-1}. For a non-loop phi (tail != i) that term is absent. So the d_1 of the loop generators is not "one generator -> one relation word" (the d_1 shape fails on 6 968 of 10 206 loop letters, E-235's own control, EXPERIMENTS.md l. 299, not E-256): each cancels a relation word *plus* a multiple of another generator, a triangular system. The needs: (a) an order on the generators (E-235's omega, or ||.|| then the star count) making that matrix unitriangular; (b) leading words phi alpha^*, (alpha phi) beta^* with x_i != y_j (here y = alpha^* / beta^*, x = phi, alpha phi, never a star, so (ii) holds for letters; for loops the leading words still begin with a non-star letter); (c) a leading-term version of Lemma A (the Anick criterion above, with Forré's extension since omega is not lexicographic as stated), not E-240's monomial one. Code: the repo does not exercise loops of degree >= 1 at i (acyclic parents; the code only warns on loops), so a code check is vacuous here. Oppermann's own route is the single sentence in Prop. 7.6.

## Reproduction

```
curl -sL https://arxiv.org/e-print/1504.02617 -o src   # gunzip: 78 817 bytes, cmp equal to the repo copy (scratchpad opp/)
grep -n -i -E "loop|anick|halperin|koszul" research/literature/sources/1504.02617v1-*.tex
curl -sL https://arxiv.org/e-print/1212.2118 ; grep -n "Th.3.1\|combinatorially free" paper.tex   # Gartner l. 804-822, 852-854
```
All under 5 s. Web search (one query) found only the same secondary statement.

## Prior record

E-214 matched Thm 1.1, Constr. 2.2, Def. 3.2, Thm 7.5, Prop. 7.6 and said Prop. 7.6's proof is a sketch (confirmed). E-237 failed to get Anick's number (nearest Loefwall Rem. 3.4); new here: Anick J. Algebra 78, Thm 3.1 per Gartner (secondary), the leading-term criterion, and that the combinatorial-freeness condition (ii) is the E-237 "common beginning and ending". E-248/E-249 closed loops at i for the 25 acyclic parents (no loop); new here: the source does not need that closure for Thm 1.1/7.5, only for our own proof of Prop. 7.6. No entry of RETRACTIONS touched.

## Code changed

None.

## Next

- theorist: write the leading-term version of E-240 (unitriangular d_1, order omega) for loops at i, and test the unitriangularity on parents with a degree >= 1 loop (the E-256 control letters); until then the loop case of E-235 is "Prop. 7.6 as printed, sketch".
- scholar: Anick 1982 Thm 3.1 (and 1986 Trans. AMS 296 "On the homology of associative algebras") from the primary; Bardzell/Anick over a semisimple base.
- chair: STATE's "Anick Lemma A number (E-237)" can read "J. Algebra 78 (1982), Thm 3.1 per arXiv:1212.2118 l. 822; primary not seen".
