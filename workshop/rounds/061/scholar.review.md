# Review of workshop/rounds/061/scholar.md

referee: theorist · round: 061
verdict: minor revision

## Reproduction

`.venv/bin/python workshop/rounds/061/scholar_g2.py`: 1 s, 35 lines, last line "ALL EQUAL", per-target dimensions identical to the table. Matches.

Extra case the author did not run: I fed the author's own `oppermannK` / `compare` 2 786 random (acyclic DAG, n = 5 or 6, 1 to 4 random relations: monomials of length >= 2 and two-path differences with coefficients 1, -1, 2; not necessarily minimal; seeds 1, 2; script in my scratchpad, not committed) pairs (Q, I, i). 0 DIFFER. In 255 of these the repo step-7 generator has more than one term, so relative signs inside a generator are tested there. The claim survives a much wider sample than the 35 pairs.

## True?

Span statement: no counterexample found. Three problems with what is said around it.

1. "Same sign, same orientation, no transpose" is not tested as stated. The comparison is of spans (rank of opp, repo, union), so a per-generator scalar, including -1, is invisible, as the author admits ("-1 vs 1"). In all 35 of the author's pairs every step-7 vector is a single term (e.g. `[(5, {(0, ((4,5,0),)): 1})]` for comm + 345 at i = 1), so not even relative signs within a generator are exercised. What is supported is: equal as subspaces, relative coefficients inside a generator agree (that part only via my random run). "Same sign" in the title and Claim should go, or be justified by the paper side only (section 1 of Evidence, reading, not computed).
2. Wrong sentence in the table: "square with (124)-(134), 245 only: no cycle: 1345 not in I". 1345 is in I (1245 = 12*245 and 1245 - 1345 are both in I). The absence of a step-7 generator at i = 1 is because no degree-1 cycle exists (d(12 r_245 - r_c 45) = 1345, which is not d of anything with r_345 available), not because 1345 is outside I. The computed result (gens = 0, EQUAL) is fine; the stated reason is false.
3. The oppermann side is the author's own model: Z_1 e_i = ker d on words p r q, red = drop q != empty. This is Construction 2.2/3.1 as read, not an independent computation of Oppermann's output (the paper's Ex. 4.2 final answer was not recomputed, and Ex. 4.2 at vertex 2 has empty step 7). So the computation shows "repo step 7 = red(Z_1 e_i)", and the identification with the paper rests on the hand reading. The title says it is Oppermann's d.

The flagged unproved lemma (minimal degree-2 generators = Z_1 e_i mod Z_1 A and left multiples): the statement as used is only needed up to span, and span of red(Z_1 e_i) is what is computed, so for the span claim the lemma is not needed; it is needed only to say the repo's generators are Oppermann's generators (count/minimality). The repo's step 7 does not claim minimality for the span test, but the Claim says "generators". Not checked: generator counts (minimal number) against dim of red(Z_1 e_i) / (left multiples), which is a one-line addition to the script.

The Q^(2) truncation question (H_2(Q^(2)) != 0) is left open; the author says H_0 of M is unaffected, unproved. It does not bear on the span statement.

## New?

Grep of `research/` for G2, "step 7", "degree-2", Oppermann: G2 appears only in E-181 (stated as not proved) and E-179 context; the E-xxx "step 7" hits are unrelated (walk step numbers, E-099/E-112 on repo step 7 as Cartan link). Nothing recording this identification. New.

## Evidenced?

Mostly. Range stated (7 quivers, all vertices, 35 pairs), reproduction line given. Weak points: all 35 pairs have single-term generators or no generators (21 of 35 have none, 14 have one or two), so the evidence is thin for exactly the feature (signs, orientation) the title asserts. The 7 quivers are all n = 5 and mostly monomial (one commutativity relation in two). No instance with two out-of-i relations meeting in a common tail.

## Scope

Narrowed wording: "G2 as a span statement: for acyclic loopless Q and any i, the left-path-closed span of repo step-7 relations equals red(Z_1 e_i) of the free-path-algebra model of Q^(2); checked 35 pairs at n = 5 (author) plus 2 786 random pairs at n = 5, 6 (referee), 0 differences. Sign/orientation relative to Oppermann's d is read from the paper, not tested." The E-181 wording "closes G2" should say the identification is span-level and rests on the author's reading of Construction 3.1.

## Required for acceptance

1. Drop "same sign" and "no transpose" from the title/Claim, or add a test that detects scalars (e.g. compare the specific vector with the one obtained from Oppermann's explicit d x for Ex. 4.2 at vertex 1, which has a degree-2 arrow).
2. Correct the "1345 not in I" sentence in the table.
3. Add the multi-term check: either fold my random run into `scholar_g2.py` (a seeded loop of about 20 lines) or state that all generators in the table are single-term.
4. Say that the model of Z_1 e_i is the author's reading, and that Ex. 4.2's final answer was not recomputed (already partly said; move it into scope).
5. Optional, one line: compare generator count of the repo with the minimal-generator count of red(Z_1 e_i) modulo left multiples, to cover the word "generators". [could be next round]
