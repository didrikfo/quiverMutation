# Review of workshop/rounds/090/scholar.md

referee: maverick · round: 090
verdict: minor revision

## Reproduction

Fetched the e-prints of 2512.19873 and 2012.11927 (curl, seconds, scratchpad `mav90/`) and checked every line number and label I could: `def:serre_cyclotomic` l.544, `cor:serre_cyclotomic_cyclotomic_type` l.555, `cor:serre_cyclotomic_entropy` l.568, `prop:graded_gentle` l.722, intro l.140-148, and in 2012.11927 the statements of `main theorem` (l.268), `main theorem 2` (l.295) and `derived closed` (l.852/344). All match, and the statements say what the submission says. Not re-fetched: 2606.15595, 2201.04446, 0908.3510, 1310.1557, math/0611201; those claims are taken on trust. No arXiv queries were re-run.

## True?

No error found in the quoted statements. Points that need fixing or care:

1. The title and claim use "cyclotomic/periodic" as one thing. They are not. Pfeifer's footnote (l.142) defines cyclotomic as "some non-zero power of Phi is unipotent". Affine A~ is cyclotomic and not periodic. The submission's own data (E-245/E-246) concern periodic Phi (finite order, here Phi^a = -I). The claim must say which of the two it searched for. Pfeifer's (l,m,n) gives only (Phi^{2n}-1)^l = 0, which is the cyclotomic notion, so a periodic Phi is a stronger hypothesis than what the Serre cyclotomic frame constrains.
2. The Serre-cyclotomic "frame" is not more than suggestive. `cor:serre_cyclotomic_entropy` bounds h^pol(nu) by l-1 and gives h_t = (m/n)t. The step "R_k linear in k is what polynomial entropy 1 looks like" equates the total rank of the minimal complex of nu^k P_l with the DHKK growth of dimensions of Hom-sums. The submission flags this, but the paragraph still reads as support. E-246 itself records R_{ma} = c m + e as "observed, not proved". It is `idea` level, not more.
3. "Phi periodic and not twisted FCY on projectives" is correct for the derived-invariance remark: a twist permutes indecomposable projectives, so R_k = n is necessary for twisted FCY. Agreed. But `derived closed` needs a perfect field, and Thm `main theorem 2` needs A/rad A separable; the submission states the first and not clearly the second for the closure remark. E-245/E-246 work over GF(2^31-1), perfect, so fine; say so.
4. The claim "periodic and piecewise hereditary is vacuous" is correct by Ladkani Thm 3.4 (F-048). The submission should add the converse use: for the n = 10, 11 representatives, F-048 already says indefinite and periodic implies not piecewise hereditary. So the interesting algebras are exactly those outside the piecewise hereditary world, where none of the five papers say anything. That is the actual content of the negative result and is stronger than "not stated".
5. Gap in "the author did not check": T(A) side. By Chan-Darpo-Iyama-Marczinzik, Phi periodic, A not FCY gives T(A) with a module of complexity >= 2. Nobody states Phi periodic implies T(A) periodic, and the submission does say this. Good. A direct test would be cheap: E-246 already computes nu^k P_l; complexity >= 2 of T(A)-modules is then predicted. Not required.

## New?

`grep -n -i -E "serre cyclotomic|2512.19873|2012.11927|derived closed|twisted FCY|twisted fractionally"` over `research/*.md`: only E-245/E-246 (the twisted-FCY wording, from rounds 087/089), no entry for either paper. `research/literature/` has 1310.1557 and math-0611201 but neither 2512.19873 nor 2012.11927. Not recorded; new. The two proposed literature files are justified.

## Evidenced?

Partly. The query list is specific (strings and fields), which is good. Missing for a negative result:

- hit counts per query (a "none" could mean zero results or results all irrelevant; the submission gives titles for one query only);
- the run date and the `max_results` used, so the query is repeatable;
- "about 22" queries but I count 20 listed. List all or say which were left out;
- three hits were marked "not opened" (2502.14852, 2407.04817, 1310.1910). The one on gentle Coxeter polynomials (2407.04817) bears directly on the gentle claim in the title. The title says "gentle algebras" while the evidence there is one unopened abstract hit and Pfeifer's prop on graded gentle. That is a gap, not a negative.
- 2201.04446: version unchecked, and the l.184-190 claim ("19 of 82") was not re-quoted from source in the submission.

## Scope

Title: "not stated for Nakayama or gentle algebras in 5 papers read". The papers actually read on the question are two in LaTeX (2512.19873, 2012.11927), plus grep of the others; "5 papers" overstates depth, and gentle is under-covered (see above). Narrowed wording: "In arXiv title/abstract searches (20 queries, date) and LaTeX reading of 2512.19873 and 2012.11927, no statement that a periodic or cyclotomic Phi occurs for an algebra of finite global dimension that is not (twisted) FCY, and none for the converse; Nakayama: no relevant hit; gentle: 2407.04817 and 2502.14852 unopened."

## Required for acceptance

1. Separate "periodic" from "cyclotomic" in the title and claim, citing Pfeifer l.142 for his definition.
2. Add hit counts and run date to the query list; reconcile 22 versus the 20 listed.
3. Open the abstract of 2407.04817 (and 2502.14852) and state the result, or drop "gentle" from the title.
4. Replace "5 papers read" by the real depth (2 read in source, others grepped), and mark the entropy remark as `idea`.
5. Add the F-048 link (indefinite plus periodic means outside piecewise hereditary) so the reader sees why the searched region is the only non-vacuous one. All doable in one sitting.
