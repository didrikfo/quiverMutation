# Periodic Coxeter transformation without (twisted) fractional Calabi-Yau: not stated in the arXiv sources searched; the merely cyclotomic case is classical (affine, gentle)

author: scholar · round: 090 · kind: negative
thread: T3/T8 · bears on: E-231, E-245, E-246 (Candidate C)
scope: arXiv API (export.arxiv.org), 21 title/abstract queries run 2026-10-27, max_results 10 (counts below). Read in LaTeX source on the question: 2512.19873v1 and 2012.11927 (API v4). Grepped in source only: 2606.15595v1, 2201.04446 (version not checked), 2407.04817v3 (Thm C), 0908.3510. Abstract only: 2502.14852v1, 1310.1910v1. Nakayama: no relevant hit. Not covered: Google Scholar, MathSciNet, non-arXiv papers, full reading of any paper.

## Response to referee

Done: all five required items.
1. Periodic vs cyclotomic separated (title, claim). Pfeifer l.142-143: "cyclotomic" = some non-zero power of Phi unipotent; periodic = finite order. Cyclotomic-but-not-FCY is known and classical: affine path algebras (Pfeifer intro l.143-148; Kronecker is gentle). Periodic-but-not-(twisted)-FCY, the E-245/E-246 situation, is what was searched for and not found. Also 2407.04817v3 Thm C (`mainthm:C`, tex l.~446): psi_A(z) = (z-1)^{|Q_1|-|Q_0|} prod (z^n-(-1)^{n+m})^{phi_A(n,m)}, so every gentle algebra has a product-of-cyclotomics Coxeter polynomial (the paper says so, l.449); this does not say Phi is periodic, and Pfeifer's `prop:graded_gentle` shows most are not Serre cyclotomic.
2. Query list below now has hit counts, date 2026-10-27, max_results 10. The earlier "about 22" was wrong: 21 queries run and listed.
3. 2407.04817v3 opened (abstract and Thm C in source) and 2502.14852v1, 1310.1910v1 abstracts read: none treats periodic Phi versus FCY (2502.14852 is about infinite-dimensional gentle orders; 1310.1910 about Mahler measure). "Gentle" stays in the title only through the cyclotomic result above.
4. Depth stated in Scope: two papers in source on the question, others grepped or abstract only. The entropy remark is marked idea.
5. F-048 link added in Claim (UNVERIFIED by me against F-048: taken from the referee and from `research/literature/math-0611201...md` Thm 3.4 lines 31-34, 98).
Point 3 of the review (separable A/rad A for the closure remark): added. Re-quoted 19 of 82 from 2606.15595 source (tex l.1944, "only 19 of these 82 distributive lattices have periodic Coxeter matrix, which is a necessary condition"); it is that paper's text, not 2201.04446's, and 2201.04446 l.184 says periodic Phi "include[s]" FCY algebras. Not done: the T(A) complexity test (review point 5, optional).

## Claim

Two cases. (a) Cyclotomic (Phi^N unipotent) but not FCY: known, e.g. every affine path algebra (Pfeifer intro; Kronecker is gentle), and all gentle algebras have cyclotomic Coxeter polynomials (2407.04817v3 Thm C). (b) Periodic (finite order) Phi but not FCY, the E-245/E-246 case: not found. In the papers searched, no statement says that an algebra of finite global dimension (Nakayama, gentle, or derived equivalent to one) can have periodic or cyclotomic Coxeter transformation without being (twisted) fractionally Calabi-Yau, and none says the converse "Phi periodic => FCY" holds. The literature states only the necessary direction (FCY => Phi periodic). So E-245/E-246 are not rediscoveries of a stated result; they are also not contradicted. Three published results do frame them, all checked in source:

1. Pfeifer, "Serre cyclotomic algebras", arXiv:2512.19873v1 (Dec 2025). Def `def:serre_cyclotomic` (tex l. 544): A is (l,m,n)-Serre cyclotomic if Iwanaga-Gorenstein and the inverse dualizing bimodule N~ admits f_1..f_l in Hom(A, N~^n[m]) with Cone(f_l) (x) ... (x) Cone(f_1) = 0. It states that (1,m,n)-Serre cyclotomic = fractionally CY of dimension (m,n) (l. 548-550). `cor:serre_cyclotomic_cyclotomic_type` (l. 555): Phi_A satisfies (Phi^{2n} - 1)^l = 0. `cor:serre_cyclotomic_entropy` (l. 568): polynomial entropy of nu is <= l-1 and h_t(nu) = (m/n) t. So a periodic Phi is consistent with l = 2, which is "one step away" from FCY exactly as affine A~ is (intro, l. 145-148). The paper calls this notion the categorification of de la Pena's cyclotomic type. `prop:graded_gentle` (l. 722): a homologically smooth proper graded gentle algebra is Serre cyclotomic iff derived equivalent to hereditary A or A~. That is a gentle-class classification, not a Nakayama one, and says nothing about periodic Phi with polynomial growth. Example at l. 747: a gentle algebra of infinite global dimension, (2,2,2)-Serre cyclotomic.
2. Chan-Darpo-Iyama-Marczinzik, arXiv:2012.11927 (latest). Thm `main theorem` (tex l. 268): for A with A/rad A separable, T(A) periodic iff A has finite global dimension and is FCY. Thm `main theorem 2 again` (l. 837): T(A) twisted periodic iff every T(A)-module has complexity <= 1 iff A has finite gl.dim and is twisted FCY. Cor `derived closed` (l. 852): twisted FCY is closed under derived equivalence (finite gl.dim, perfect field). Consequence for us: S^a P_l not isomorphic to P_sigma(l)[b] for any a is a derived-invariant property (it rules out twisted FCY), and "Phi periodic, not twisted FCY" is equivalent to "T(A) has a module of complexity >= 2 although Phi is periodic". Our R_a = n test is necessary for the twisted case too (a twist permutes the projectives), so E-246's NO is a NO for twisted FCY on projectives, not only for FCY.
3. Marczinzik-Thomas-Yildirim, arXiv:2201.04446 (intro l. 184-190) and Chan et al. l. 1944 (in 2606.15595): "periodic Coxeter matrix" is stated as a necessary condition for FCY (finite gl.dim), e.g. only 19 of 82 distributive lattices on 11 points have periodic Phi. Again only the necessary direction.

Also read: Pena 1310.1557 (already `research/literature/`); Ladkani math/0611201 Thm 3.4 (already recorded; periodic + piecewise hereditary => non-negative form, so for piecewise hereditary algebras the question is vacuous). Hochschild-calculus paper 2606.15595v1 Cor `cor:fcy` (l. 1230): FCY => sigma^q = (-1)^{p-q} id on HH_0, so Cox^q = (-1)^{p-q}; recovers Pena's Thm A(a); the Serre cyclotomic case is claimed only for quasi-unipotence (l. 1237).

Why (b) is the right region: by Ladkani Thm 3.4 (record: `research/literature/math-0611201-...md`, and F-048 per the referee) periodic Phi with an indefinite form excludes piecewise hereditary algebras, so for the n = 10, 11 representatives the literature on hereditary/tubular/canonical algebras says nothing, and none of the papers read covers the remainder. The closure remark needs A/rad A separable (CDIM) and a perfect field; E-246 works over GF(2^31-1), which is perfect.

What this does not claim: that no such example exists in print. arXiv full-text search is not available; my queries match titles and abstracts only.

## Evidence

Queries (export.arxiv.org, search_query, run 2026-10-27, max_results 10; totalResults first). Counts: all:FCY AND Coxeter polynomial 2; FCY AND Nakayama 5; periodic Coxeter AND CY 1; FCY AND gentle 2; cyclotomic AND Coxeter polynomial AND finite gl.dim 0; CY AND Coxeter transformation AND finite order 0; abs Serre cyclotomic 1; ti CY AND abs Coxeter AND periodic 0; ti FCY 10; abs FCY AND Coxeter 3; Coxeter polynomial AND Nakayama 2; periodic AND Coxeter transformation AND derived 1; Serre functor AND Coxeter polynomial 0; periodic Coxeter AND fractionally 1; Coxeter transformation AND periodic AND CY 1; rowmotion AND Coxeter 0; Coxeter polynomial AND CY AND counterexample 0; periodic Coxeter matrix 0; Nakayama algebras AND derived equivalence AND CY 0; gentle AND Coxeter polynomial 1; cyclotomic type AND algebras 6 (21 queries). Hits as before:
fractionally CY AND Coxeter polynomial [2606.15595, 1310.1557]; fractionally CY AND Nakayama [none relevant]; periodic Coxeter AND Calabi-Yau [1310.1557]; fractionally CY AND gentle [2502.14852 (gentle orders, abstract read), 2512.19873]; cyclotomic AND Coxeter polynomial AND finite global dimension [none]; Calabi-Yau AND Coxeter transformation AND finite order [none]; "Serre cyclotomic" [2512.19873]; ti:"fractionally Calabi-Yau" [10 papers, all about FCY examples: 1008.1245, 2604.19582, 2012.11927, 2603.23354, 0908.3510, 2303.11656, 1807.08503, 2406.09148, 1509.07657, 1610.04918; titles only except 2012.11927 and 2603.23354 abstract]; Coxeter polynomial AND Nakayama [2203.15735 (in record), 2606.15595]; periodic AND Coxeter transformation AND derived [1310.1557]; Serre functor AND Coxeter polynomial [none]; periodic Coxeter matrix [none]; rowmotion AND Coxeter [none (2201.04446 found by citation)]; Coxeter polynomial AND Calabi-Yau AND counterexample [none]; Nakayama AND derived equivalence AND Calabi-Yau [none]; gentle AND Coxeter polynomial [2407.04817, opened]; cyclotomic type AND algebras [1310.1557, 2512.19873, 1310.1910 (Mahler measure, abstract read), 2110.15329]. Total 21 queries (earlier "about 22" corrected). Hit 2407.04817 is now opened: see Response 1.
Fetched: LaTeX e-print of 2512.19873, 2606.15595, 2012.11927, 2201.04446, 0908.3510, plus API abstracts of 2512.19873, 2606.15595, 2012.11927, 2603.23354. Why: these four are the only hits that mention both Coxeter data and FCY; 0908.3510 for Herschend-Iyama's question on derived invariance (answered by CDIM). Files are in my scratchpad (/tmp/claude-0/-home-user-quiverMutation/61ae55b4-f248-525f-bc74-a4acc5de8006/scratchpad/sch90/), not in the repo; nothing is cited from memory.

Relation to the data (idea level, not support): E-245/E-246 report periodic Phi with R_k growing about linearly in k. Total rank R_k of the minimal complexes of nu^k P_l is not the same as the DHKK dimension quantity in Pfeifer's entropy, but linear growth with h_0 = 0 is what polynomial entropy 1 looks like, i.e. consistent with (2,m,n)-Serre cyclotomic, and would be incompatible with FCY (l = 1, entropy 0). Whether these Nakayama algebras are Serre cyclotomic is open: it requires f_1, f_2 with Cone(f_2) (x) Cone(f_1) = 0 on bimodules, which we have not tested. Pfeifer's classification for gentle (prop:graded_gentle) excludes proper graded gentle algebras beyond A, A~; Nakayama algebras of the E-246 kind are not gentle in general, so it does not apply.

## Reproduction

```
curl -s -G "https://export.arxiv.org/api/query" --data-urlencode "search_query=ti:\"fractionally Calabi-Yau\"" --data-urlencode max_results=10
curl -sL https://arxiv.org/e-print/2512.19873 -o s.src   # likewise 2606.15595, 2012.11927, 2201.04446, 0908.3510
```
Each seconds. No script in the repo.

## Prior record

E-231: of 0911.5137 and 1310.1557 only, necessary conditions. E-245/E-246: the examples. This round adds: no stated counterpart in print; the Serre cyclotomic frame (2512.19873) and derived invariance of twisted FCY (2012.11927 Cor `derived closed`). Neither is in `research/literature/`.

For the chair, two files to create (text above is the summary; shape as the existing ones):
- `research/literature/2512.19873-serre-cyclotomic-algebras.md` (v1; items 1 above, with the l. numbers)
- `research/literature/2012.11927-periodic-trivial-extensions-fcy.md` (item 2 above; keep the statement that the twisted/untwisted hypotheses differ, and that A/rad A separable is assumed)
Source LaTeX is not committed; say "tex line numbers refer to the arXiv e-print as of 2026-10-27".

## Code changed

None.

## Next

- Referee/mathematician: in E-246's wording say "not twisted FCY on projectives" (the test is permutation-invariant), and cite CDIM Cor `derived closed` for derived invariance.
- Mathematician: test Serre cyclotomic of order 2 on the n = 10 representatives (does N~ ^n[m] have a bimodule map with cone square zero?); even the Grothendieck-level prediction (Phi^{2n} - 1)^2 = 0 is trivially met, so a stronger invariant is needed, e.g. h_t(nu) = (m/n) t: compute the slope from the shift of nu^k P_l and compare with -m/n.
- Scholar: fetch 2502.14852, 2407.04817, 1310.1910 abstracts-in-full if the gentle direction is wanted; ask the human for Google Scholar / MathSciNet if "not in print" must be stronger.
