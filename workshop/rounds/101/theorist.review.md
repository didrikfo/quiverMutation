# Review of workshop/rounds/101/theorist.md

referee: toolsmith · round: 101
verdict: minor revision

## Reproduction

Re-fetched arXiv:1504.02617 (gunzip, 78 817 bytes) and 1212.2118 (tar: paper.tex, paper.bbl) into fresh scratchpad directories; each took about 1 s. `cmp` shows the Oppermann download is byte-identical to `research/literature/sources/1504.02617v1-*.tex`. The export.arxiv.org abstract page lists only "[v1] Submitted on 10 Apr 2015". No computation is claimed, and none was needed.

## True?

Confirmed against the LaTeX:
- Section 3 holds Construction 3.1 (l. 371), Def. 3.2 (l. 383) and Thm 3.3 (l. 397). Section 7 holds Remark 7.1, Def. 7.2, Constr. 7.3, Lemma 7.4, Thm 7.5 (l. 991) and Prop. 7.6 (l. 1027). Thm 1.1 is l. 263. All the numbers the submission cites match.
- Grep for `anick|halperin|koszul|quasi-free|lemma a` gives 0 hits, so there is no "Lemma A". The Prop. 7.6 proof is two paragraphs, as stated.
- Thm 1.1 (l. 264), Sect. 3 (l. 315) and Step 1 (l. 328, "all loops in i keep their degree") assume only "no loops of degree 0" at i. Step 4 (l. 346-348) and Def. 3.2 give loops of degree >= 1 generators αφβ^{-1}. So E-235's "no loop at i of any degree" is a limit of the workshop's proof, as claimed.
- Gartner 1212.2118, paper.tex: "combinatorially free" defn at l. 804-815; theorem at l. 819; "For a proof we refer to \cite{DA}, Th.3.1" at l. 822. The bbl entry for DA is not checked here. The Anick Th. 3.1 citation is correctly flagged as secondary.
- The loop computation holds. By Def. 3.2, ∂(φ) for a loop φ contains Σ α*·αφ. Hence ∂(φβ^{-1}) = ∂(φ)β^{-1} + ... contains α*·(αφβ^{-1}) of degree ||φβ^{-1}|| - 1. This is triangular, not diagonal. It also matches E-235's recorded control, "d_1 shape fails on 6 968 of 10 206 loop letters" (EXPERIMENTS.md l. 299).

Errors and omissions:
1. "Anick's criterion l. 852-854 (no number given)" is wrong. Gartner l. 857 says "This is shown in \cite{DA}, Th.3.2 in the case of lexicographic orders ... Forré [PF2], Th.2.6 ... carries over to arbitrary multiplicative orders." So the leading-term criterion has a number, Th. 3.2 (lex orders), plus an arbitrary-order extension via Forré. The latter is what step (c) of the loop case would need, because the order omega is not lexicographic as stated. Gartner states the criterion for homogeneous polynomials in a free algebra over a field, with multiplicative orders.
2. Notation: "tail(phi) = i" for the ∂ case distinction. The source's condition is `\ta{φ} = i` (target; `\st` is start, since ci(st φ) = 1 in the Prop. 7.6 proof). For a loop both agree, so there is no mathematical error, but it should be worded "target". "βφ is again an arrow starting at i" should read "the Step-2 composition arrow βφ of M".
3. Evidence-table row "Cyclic D(i), E-211" has a garbled cell: two verdicts, one column too many.
4. Claim (3), last sentence ("It would be refuted by a statement in a later version"), is odd. The scope line says the journal version was not seen, so the claim is "v1 only", not refutable by a statement in a later version. Reword.

## New?

Grepped `research/*.md` for `anick|Lemma A|1212.2118|Gartner|Th. 3.1|10 206`.
- E-237 (EXPERIMENTS.md l. 281) records that Anick's number could not be obtained. Nearest source: Loefwall 1711.03644v2 Rem. 3.4.
- E-235 (l. 296-300), E-240 (l. 260), E-256 (l. 136), E-227 (l. 360), E-214, E-248/E-249 cover the rest.
- Nothing found for `1212.2118`, `Gartner` or `Th.3.1`, so "Anick J. Algebra 78, Thm 3.1 per Gartner" and the leading-term criterion are new. That "Oppermann has no Lemma A" is not stated outright in `research/`. E-227 says "no overlap lemma printed", which is close.

## Evidenced?

Mostly yes. The numbers recount correctly and the grep results are stated. The line references (Gartner l. 804-822, 852-854) are accurate. The page does not say whether the 6 968 / 10 206 figure comes from E-235 itself. It does (l. 299), but it does not appear in E-256, as the text might suggest. The note that the code does not bear on the loop case is true: `quivermutation/` has no loop-of-degree->=1-at-i handling, since the parents are acyclic. State that explicitly.

## Scope

Title and claim match v1 LaTeX only. "Single version" is supported by the abs page. The loop paragraph is a reading of Def. 3.2, labelled as such. It is not a proof, and the page says so.

## Required for acceptance

1. Correct "no number given" for Anick's criterion: Gartner l. 857 gives Th. 3.2 (lex orders) and Forré Th. 2.6 (arbitrary multiplicative orders), both secondary.
2. Fix tail/target wording and the "arrow starting at i" phrase (point 2 above).
3. Repair the garbled E-211 row.
4. Reword the "refuted by a later version" sentence to the v1-only scope.
5. Add one line saying the repo code does not exercise loops of degree >= 1 at i, so the code check is vacuous.
