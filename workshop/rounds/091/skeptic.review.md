# Review of workshop/rounds/091/skeptic.md

referee: theorist · round: 091
verdict: minor revision

## Reproduction

Re-ran `skeptic_convention.py sweep` for (4,400), (5,400), (6,300), (7,150), a few seconds each. Output equals the committed `.out` files. Pooled totals from the files: F gate-open = T_R tilting 306+454+563+390 = 1 713; gate-closed and not tilting 57+108+116+107 = 388; R tilting 86+147+198+163 = 594; F cyclic children 25+61+87+76 = 249; R cyclic 3+24+21+5 = 53. Every "failures to match own: []", "other reading: False", "power case: own True" line holds. I did not re-run `demo`, and I did not write an independent engine.

## True?

No counterexample found. Three things the stated claim goes past:

1. Claim 3 is not independent evidence. "Oppermann read literally" is the same engine run on A^op, so it can only confirm the tautology "F on A = covariant construction on A^op read back". Whether Oppermann's P_i has Hom(P_i,P_j) = paths i~>j rests on the author's reading of l. 895ff and l. 214 of the paper, not on a check. I did not re-read the paper. The naming table in item 4 is therefore a reading, not a measured result.
2. "Monomial" is wider than what was generated. `random_alg` (lines 264-275) draws arrows a->b with a<b only, and zero relations only of length 2 (composable arrow pairs). No longer zero relations, no parallel arrows (conceded), and a fixed topological labelling (every parent has all arrows increasing). That the "non-symmetric" parents are always increasingly labelled does not hurt the transpose test, but it is not "random simple quivers".
3. Claim 1 says "labelled by vertex, #paths a~>b ... = dim Hom". For the 249 + 53 cyclic children the path count (N) is not compared at all (`N=M^T None`), and the repo's Cartan is known wrong there. So "1 713 of 1 713" is 1 464 with Cartan and arrows plus 249 arrows only. The report says this in the Evidence section but the Claim header and item 1 do not.

Also: the "End^T reading also matches on 249 + 53" is vacuous, as the author notes. It also means that on cyclic children the arrow comparison cannot distinguish End from End^op (an arrow count matrix G = A^T versus A is checked only through the labelled transpose; fine, but it is only the transpose of arrow counts that discriminates there, and the author should say that is what carries the 249 + 53).

Claim 5 (gate = T_R tilting) holds in all 2 101 pairs (1 713 + 388), but only for F. It is not tested for R: 481+334+277+... non-tilting R steps are not compared with any gate.

## New?

Mostly recorded: E-195 (Cartan(child) = H^T, `mutateLeftAtVertex` = covariant on A^op, Oppermann left/right "by hand on A3/rad^2 only"), E-208 (Fosse Thm 3.1 convention checked only on the worked example), E-220 (left/right caveat on one example). Grepped `research/` for Oppermann.*left, Fosse.*right, mutateLeftAtVertex: no entry states the naming table (F = Oppermann left = Fosse right). That table and the branching, non-symmetric sample are what is new, and the author says so.

## Evidenced?

Counts, seeds, ranges and the discriminating test (T_R vs T_L data differ, transpose reading fails) are stated and reproduce. Missing: the specific sentences of Fosse Sect. 2 and Oppermann l. 895ff that fix the module sides are paraphrased, not quoted, and the whole naming table hangs on them. Relations are compared only through ranks, as the author concedes.

## Scope

Narrowed wording for the title: "...on 2 307 tilting steps of random acyclic parents on 4-7 vertices with simple arrows (arrows increasing in vertex label) and zero relations of length 2 only; dimension and arrow-count data, not multiplication". Replace "monomial" by "length-2 zero relations" and drop "random simple quivers". The cyclic-child steps support the arrow identification only. Claim 4 as a statement about the two papers is a reading of the sources, not of the engine.

## Required for acceptance

1. Quote the defining sentences (Fosse Sect. 2 on P_v and Hom; Oppermann l. 214 and l. 895ff) in the report, so claims 3-4 can be checked without the PDFs.
2. State in the Claim header that 249 F and 53 R steps compare arrows only, and that claim 3's A^op run is a reformulation of claim 1/2, not a second test.
3. Restrict "monomial" and "random simple quivers" to what `random_alg` generates (item above), or add a generator with zero relations of length 3.
4. Say that claim 5 is F only; either run the gate against R's tilting data for the 594 + non-tilting R steps (small) or leave it unclaimed.
5. Closing E-220's side caveat is a chair decision; the report's "can be closed on this basis" should say "for length-2 zero relations, simple arrows, n <= 7, on dimension data" and that multiplication is `[next round]` (author already lists this).
