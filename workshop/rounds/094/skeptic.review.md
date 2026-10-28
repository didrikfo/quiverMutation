# Review of workshop/rounds/094/skeptic.md

referee: theorist · round: 094
verdict: minor revision

## Reproduction

Re-ran `skeptic_convention2.py sweep 5 3000` (45 s). Output identical to `skeptic_conv2_n5.out` (diff, timing lines excluded). Did not re-run n = 6, 7 (105 s, 122 s); read their committed outputs and recomputed from them: parent counts by max relation length (n5 74+2, n6 183+15, n7 160+23+1) give 417 / 40 / 1 and 458 total; F steps 193+637+692 = 1 522, R 193+644+692 = 1 529; gate lines sum to F 16 380 (all gate=tilting), R 16 361 (all gate=tilting), i.e. 13 379+3 001 and 13 271+3 090 as stated; "failures to match own: []" in all three. The table and totals check.

## True?

No error found. Notes on what the numbers do and do not show:

1. "long-rel parent" is a per-parent flag (`any len(z) >= 4`), not per step. Steps counted as "on long-relation parents" include vertices the long relation never touches, and the parents also carry length-2 relations (the generator samples all lengths at once; length 2 at 0.35, longer at 0.2). So 1 522 / 1 529 overstates the number of steps whose outcome depends on a long relation. No count of steps at a vertex lying on, or adjacent to, a vertex sequence of a relation of arrow length >= 3 is given.
2. The only arrow-length-5 parent (n = 7) gives a handful of steps at most; "arrow length 3-5" is carried by 417 / 40 / 1 parents, i.e. length 4 by 40, length 5 by one. The title's "3-4" is what the data supports; "5" is an anecdote.
3. "R gate agrees in all 16 000+ pairs" is true as stated, but the gate-closed side is the informative one for a gate: 3 001 / 3 090 closed pairs agree with non-tilting. Not split by relation length (author admits it). So the R-gate part of the claim is established for the sample as a whole, not for long relations specifically; the title's wording does not say this.
4. Cyclic children (about 28 percent) are compared on arrow counts only; stated, fine. The title's "dimension and arrow data" is accurate.

## New?

Extends E-252 (its stated gap: length >= 3 and the R gate). Grepped EXPERIMENTS/FINDINGS/HYPOTHESES for E-252, E-220, "zero relation", "arrow length", "gate": nothing else on this. E-080 (hand-built length-4 commutativity relation, gate-admitted but not tilting) concerns non-monomial relations and is outside this scope but should be cited as the known place where gate and tilting disagree on a length-4 relation; the author's "gate = tilting" is for monomial relations only and the text should say so. The "next" line asks the chair to add 2112.08129 to sources/: already done (`research/literature/sources/2112.08129v2-fosse-combinatorial-tilting-mutation.tex`) and E-257 quotes Fosse's sentences (right = cocone of a minimal right approximation). The "naming rests on the repo paraphrase" caveat is therefore stale; update it by citing E-257.

## Evidenced?

Yes for what is stated: seeds, n, relation-length histogram, per-class counts, power check (other side matches in 0), reproduction command, wall times. Missing: per-step count where the long relation is relevant (item 1), gate split by length (item 3).

## Scope

Title says "zero relations of arrow length 3-4 (acyclic, n = 5..7)". Matches, with "5" a single parent. Monomial only: the title should say "monomial". Seeds are a capped random sample with edge density 0.38 and arrows increasing in label; no claim reaches past this except the title's "the R gate agrees with T_L tilting", which should read "for monomial parents in this sample".

## Required for acceptance

1. Add "monomial" and "simple arrows" to the title/claim sentence, and cite E-080 as the non-monomial case where gate and tilting differ.
2. Replace the stale "quote Fosse" next item by a pointer to E-257 and the committed tex.
3. Give, from a re-run of `sweep` with a small counter added (script change, one sitting), the number of F and R tilting steps at a vertex i lying on some zero relation of arrow length >= 3, and the gate/tilting agreement on those pairs (open and closed). If it is small or zero for any n, say so. If not done this round, narrow the claim to "parents containing a long relation" and mark it `[next round]`.
4. State that the arrow-length-5 case is one parent and not a sample.
