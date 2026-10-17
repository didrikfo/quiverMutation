# arXiv:2112.08129 (v2) states steps 1-7 and no LNA rule; E-187/E-205 use steps 1-7 as printed (Sect. 6 example compared visually), and for an LNA with monomial parent the gate equals the hypothesis of the paper's Theorem 4.3, reading Hom(P_i*[1], Lambda) as {g : g.[alpha] = 0} (derived, not printed)

## Response to referee

Verdict received: minor revision, three required items.

1. Theorem numbering. Accepted; the referee is right. I re-counted from the `\newtheorem` lines (l. 31-40), re-downloaded the e-print and ran a script over the .tex (`python -I`, counting `\begin{...}` outside comments per `\section`). Rule: `\newtheorem{theorem}{Theorem}[section]` (l. 31) is the base counter, reset per section; `lemma`, `observation`, `remark` use `[theorem]` (l. 33-35), so theorem, lemma, observation and remark share ONE counter per section. `corollary` is numbered within `theorem` (l. 32, unused as far as listed). `definition` and `example` have their own counters (l. 39-40), so they do not advance it (Ex. 3.1, 3.2 are examples, not part of the theorem count). Result: Sect. 3: Thm 3.1 (l. 136). Sect. 4: lemma 4.1 (l. 368), observation 4.2 (l. 374), [definition l. 380, own counter], theorem 4.3 (l. 391, `mutationAssumptions`), theorem 4.4 (l. 422), lemma 4.5 (l. 453). Sect. 5: lemma 5.1 (l. 507), 5.2 (l. 545), 5.3 (l. 622), remark 5.4 (l. 660), lemma 5.5 (l. 688), remark 5.6 (l. 757), lemma 5.7 (l. 772). All other numbers (Thm 3.1, Thm 4.4, Lemma 5.1, Lemma 5.7, Remark 5.6) re-verified unchanged. Every "Thm 4.2" is now "Thm 4.3" (Evidence table, Comparison, Draft). Thm 4.2 was wrong: 4.2 is an observation.
2. Gate = hypothesis. Accepted; title, first line of the Claim and Scope now say: monomial parents only, Hom(P_i*[1], Lambda) read as {g : g.[alpha] = 0} (derived, not printed). Also stated: Thm 4.3 is a one-way implication (Hom = 0 implies tilting), so "gate = hypothesis" says the gate is the paper's sufficient condition for mutation; its being necessary rests on Thm 4.4 plus my derivation.
3. Sect. 6 comparison is visual. Accepted; now stated in the title and in the Claim (point 2), not only in Scope: the match of `scholar_paper_example.py` output to the paper's two pictures was read by eye, not checked mechanically. Also kept: right mutation only, op/left-right convention unproved.

author: scholar · round: 071 · kind: result (literature check)
thread: T7 (agenda 4) · bears on: E-187, E-205, E-182, H-010, F-016, F-039
scope: the LaTeX of 2112.08129 (e-print, 86 KB, sha256 of the .src gzip d5b2f668...6c6d; arXiv lists v1 and v2, `updated` 2021-12-21; the e-print is the current one, v2 UNVERIFIED as byte-identical to the v2 PDF). Theorem numbers counted by script from `\newtheorem` (shared per-section counter: theorem, lemma, observation, remark; definition and example separate; no pdflatex here), so "(counted)". One right-mutation reading; left mutation only via the paper's one-line duality remark. Gate/tilting remark covers monomial parents only, with the derived (not printed) Hom description; Thm 4.3 is one-way. The Sect. 6 comparison is visual (by eye). Right mutation; op/left-right convention unproved.

## Claim

(0) Gate = Thm 4.3 hypothesis holds for monomial parents only and uses the derived (not printed) description Hom(P_i*[1], Lambda) = {g : g.[alpha] = 0}; Thm 4.3 is one-way (Hom = 0 implies tilting).
(1) The paper does not state the one-step LNA rule. `grep -ci nakayama` = 0, no "LNA", "Catalan", "truncated" or "slide" of a relation; its only line algebra is A_{2n,3} in Sect. 6, where every relation has length 3. The rule (gate passes and r_v = 2 or v = n-1; result moves (r_{v-1}, r_v) := (2, 0); counts C_{v-1} C_{n-1-v}) is the project's (E-182, E-187, E-205) and has no counterpart in the paper. The paper's "square slides in/to the right" (l. 889, 898, 939) is a commutative square in a non-line quiver, not the length-2 relation of E-187.
(2) Steps 3-7 as E-187 uses them agree with the printed text on every point I could match (table below). The comparison with the Sect. 6 worked example (A_{6,3}, two mutations) is visual: the repo output was matched to the paper's pictures by eye, not mechanically. No step is used by E-187 outside what Theorem 3.1 states. The paper's own first example is E-187's case (F), r_v >= 3, v = 1 (a non-line result, n arrows, a vertex of out-degree 2), and the repo reproduces it.
(3) Not claimed: that the paper proves E-187's "result is the line modulo the displayed relations" (that is E-187's path-class argument, relative to the repo's reduction), nor the left/right and op conventions (below).

## Evidence

Statements matched in the source (file `A_combinatorial_procedure_for_tilting_mutation.tex`, line numbers are of that file; numbers counted):

| item | where | statement |
|---|---|---|
| Thm 3.1 (`mainResult`, l. 136) | main theorem | no arrow i -> i, assume Hom(P_i*[1], Lambda) = 0; the procedure gives a quiver with kQ/I iso to End(mu^R)^op. The printed hypothesis reads "for all i != 0" (typo; i is fixed) |
| Steps 1-7 | l. 144-165 | wording matches the repo summary and E-187 (r-bar a minimal relation, alpha* : t(alpha) -> i*, r-bar alpha* = r/alpha, step 7 iff with eps_r paths t(r) -> l) |
| Note after step 7 | l. 171 | a relation containing a path of length one is cancelled against that arrow: this is E-187's "reduction eliminates gamma" |
| Thm 4.3 (`mutationAssumptions`, l. 391) | | Hom(P_i*[1], Lambda) = 0 implies mu^R is a tilting complex; the proof says a nonzero shifted endomorphism "by -1" is excluded because [alpha] is a right approximation |
| Thm 4.4 (l. 422, no label) | | the two impossible cases: no arrow out of i; a nonzero path ending in i whose composite with each arrow out of i is zero. Printed as a theorem, called "corollary" in the text above it |
| Lemma 5.7 (`Qi*`, l. 772) | | points 0-6 are steps 1-7 for End(Lambda + P_i*)^op; point 6 = step 7, from Remark 5.6 (l. 757); vertex removal is Lemma 5.1 (l. 507), which gives steps 1 and 6 |
| Ex. 3.1, Ex. 3.2 (l. 280, 325) | | acyclic example with a branch; the cyclic example (relation beta alpha at a 2-cycle) that the repo does not reproduce (F-002) |
| Sect. 1 (l. 83-86) | | only right mutation is treated; left is "exactly the same, except that all arrows will be flipped" |

Comparison, steps 3-7 against E-187 (right mutation at v of an LNA, beta = a_{v-1}, alpha = a_v, gamma = a_{v+1}, r = alpha gamma when r_v = 2):
- Step 3 (r -> r-bar : v* -> v+2): same; the "one arrow per arrow alpha" cyclic clause (l. 150) is never active on a line.
- Step 4 (alpha* alpha beta = 0, sum over arrows out of i): E-187's "c alpha* = 0" is the paper's alpha* c. Order only: the paper composes right to left (l. 144, `alpha beta : h -> j`), E-187's prose writes c alpha* for "c, then alpha*" but r-bar alpha* for "alpha*, then r-bar". Cosmetic, the paths are the right ones.
- Step 5 (r-bar alpha* = r/alpha): same, giving gamma for r_v = 2 (cancelled by the Note) and a_{v+1}...a_{v+r-1} for r_v >= 3 (kept, so not a line).
- Step 6 (relation ending at i extended by each alpha): same, E-187's "j .. v-1, c".
- Step 7: on a line R is {r} or empty, so the iff reads "eps r-bar = 0 iff (gamma eps) in I", E-187's kernel; the paths are unique, so "linear combination" is a single path.
- Gate: Hom(P_v*[1], Lambda) = {g : P_v -> Lambda with g.[alpha] = 0}, which I derived (not printed) from the cone P_{t(alpha)} -> P_v in degrees -1, 0; a nonzero such g is a nonzero element of e_j Lambda e_v killed by every arrow out of v. Thm 4.4's second case is this with "path" for "element". For a monomial parent e_j Lambda e_v has dimension <= 1, so path = element and the repo gate (relation ends at v+1) is exactly the hypothesis of Thm 4.3: for LNAs the gate is necessary and sufficient for Hom(P*[1], Lambda) = 0. The general gap (single paths vs combinations) is Ladkani 1001.4765 Prop 2.3(c), already recorded in `research/literature/1001.4765-...md`; not new.
- Repo check of the paper's worked example (Sect. 6, A_{6,3}): `scholar_paper_example.py` mutates at 1, then at 3: after vertex 1 arrows 2->1*, 1*->4, 2->3, 3->4, relations 1*45, commutative 2-1*-4 = 2-3-4, 3456 (the paper's first picture); after vertex 3 (n = 8) the line 2,1*,4,3*,6 with square {4, 3*, 5, 6}, relations 2-1*-4-3*, 1*45, 3*67, 5678 (the second picture). Read off the pictures by eye, so a visual match; tests/test_mutation_procedure.py already pins the Example 3.1 case.

Slips in the paper that do not affect the above: Thm 3.1 "i != 0" (l. 139); Sect. 6 text says "steps 1, 3 and 6 don't apply" at vertex 1 and "by step 4 the relation 1 ~> 4 becomes an arrow" (steps 1, 4, 6 and 3 are meant, l. 873) and "step 3 says we should add a relation 2 ~> 3" at the second mutation (step 4, l. 880); Thm 4.4 is called a corollary (l. 418).

## Reproduction

```
mkdir -p $SCRATCH/ar/u && curl -sL https://arxiv.org/e-print/2112.08129 -o $SCRATCH/ar/x.src && (cd $SCRATCH/ar/u && tar xzf ../x.src)   # seconds
grep -n -i "nakayama\|LNA\|Catalan" $SCRATCH/ar/u/*.tex                                  # no hits
.venv/bin/python workshop/rounds/071/scholar_paper_example.py 8                          # under 5 s
```

## Prior record

`research/literature/2112.08129-combinatorial-tilting-mutation.md` already has steps 1-7, the admissibility discussion, the cyclic-step caveat and the step-7 quote; FINDINGS l. 2599 and EXPERIMENTS l. 3749 say `procedure.py` is the paper. What is new: version, theorem numbers, the line-by-line match of E-187, the LNA exactness of the gate, the three slips. The linear-combination gap in Thm 4.4 is known (1001.4765 summary, R-005, F-016).

## Code changed

New `workshop/rounds/071/scholar_paper_example.py` (read-only run of the repo's mutation). No repo code, no tests.

## Draft for `research/literature/2112.08129-combinatorial-tilting-mutation.md` (chair promotes)

Add after the title block: "**Read:** LaTeX of the e-print (v2; arXiv lists v1, v2), 2026-10-17 round 071. Numbers: Thm 3.1 main result; Thm 4.3 Hom(P_i*[1], Lambda) = 0 gives a tilting complex; Thm 4.4 the two impossible cases (printed as a theorem); Lemma 5.7 steps 1-7 for End(Lambda + P_i*)^op; Ex. 3.1, 3.2 (acyclic; cyclic = F-002); Sect. 6 the A_{2n,3} example. Numbers counted from `\newtheorem` (l. 31-35): theorem, lemma, observation, remark share one counter per section (Sect. 4: lemma 4.1, observation 4.2, Thm 4.3, Thm 4.4, lemma 4.5); definitions and examples count separately." Add under Admissibility: "Thm 4.3's hypothesis is Hom(P*[1], Lambda) = {g : g.[alpha] = 0}, i.e. no nonzero element of e_j Lambda e_i killed by all arrows out of i (derived, not printed); Thm 4.4 states it for paths only. For a monomial algebra the two agree, so on LNAs `isMutable` is exact for this hypothesis." Add under Caveats: "The paper has no Nakayama or LNA statement; the one-step LNA rule is ours (E-182, E-187, E-205). Printed slips: Thm 3.1 'i != 0'; Sect. 6 swaps steps 3 and 4 at l. 873 and 880."

## Next

- theorist: E-187/E-205 can say "matched to the printed steps 3-7 (round 071)" instead of "not the paper"; the remaining gap is the reduction step (ideal = vanishing line paths), which the paper does not cover. Two bystanders and k >= 4 have no paper counterpart.
- skeptic (T5): the op/left-right convention of Thm 3.1 (End^op, right mutation, P_i = paths starting at i, l. 471) against the repo's arrow direction is read only through the worked example, not proved; same open item as Oppermann's convention.
- chair: promote the draft; do not commit the .tex (not cited by a record entry beyond the summary).
