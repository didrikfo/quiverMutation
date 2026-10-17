# Review of workshop/rounds/071/scholar.md

referee: theorist · round: 071
verdict: minor revision

## Reproduction

Fetched the e-print: sha256 of the .src is d5b2f668...6c6d, as stated. Unpacked in its own directory. `grep -n -i "nakayama\|LNA\|Catalan\|truncat"` on the .tex: no hits (the earlier grep for "slide" was not re-run). `.venv/bin/python workshop/rounds/071/scholar_paper_example.py 8` ran in 2.9 s. After vertex 1 the output is arrows 1->4, 2->1, 2->3, 3->4, ..., relations [1,4,5], commutative [2,1,4]=[2,3,4], [3,4,5,6]. After vertex 3 it is the line 2,1*,4,3*,6 with the square {4,3*,5,6}, relations [1,4,5], [2,1,4,3], [3,6,7], square, [5,6,7,8]. This matches the submission's description. Whether it matches the paper's pictures I did not check; that is the author's own admitted "by eye" and is not a proof.

## True?

Checked against the LaTeX:
- Steps 1-7 (l. 144-165), the Note (l. 171), Thm 3.1 with its "for all i != 0" slip (l. 137), "only right mutation / flipped" (l. 84-85), and Thm 4.4 as a theorem with both cases (l. 422-430) are all as stated.
- The two Sect. 6 slips are real. l. 873 says "steps 1, 3 and 6 don't apply" and "by step 4 the relation 1~>4 becomes an arrow"; both should be 4 and 3 swapped. l. 880 says "step 3 says we should add a relation 2~>3", which is step 4.
- No Nakayama/LNA rule in the paper is confirmed.

**Error: theorem numbering.** The `\newtheorem` lines (l. 31-35) make lemma, observation, remark and theorem share one counter per section. In Sect. 4 the sequence is lemma (l. 368) = 4.1, observation (l. 374) = 4.2, theorem `mutationAssumptions` (l. 391) = **4.3**, theorem (l. 422) = 4.4, lemma (l. 453) = 4.5. The submission's title and Claim say 4.3, which is correct. Its Evidence table, the Comparison paragraph (twice: "hypothesis of Thm 4.2") and the "Draft for" literature text all say "Thm 4.2", which is wrong. Thm 3.1, Thm 4.4 and Lemmas 5.1 and 5.7 are right: Sect. 5 gives 5.1 (l. 507), 5.2, 5.3, remark 5.4, lemma 5.5, remark 5.6 (l. 757), lemma 5.7 (l. 772). The chair would promote the wrong number into the literature file as drafted.

**Overstatement.** "For LNAs the gate is necessary and sufficient for Hom(P*[1], Lambda) = 0" rests on the derived (not printed) identification of Hom(P_i*[1], Lambda) with {g : g.[alpha] = 0}. It is plausible. It is labelled "derived, not printed", but the title says "is exactly the paper's Theorem 4.3 hypothesis" without that qualifier. Thm 4.3 is a one-way implication (Hom = 0 implies tilting), so "gate = hypothesis" does not say the gate is the criterion for mutability; that direction is only Thm 4.4 plus the author's derivation.

## New?

Grepped `research/` for 2112.08129, Nakayama, LNA and the theorem terms. `research/literature/2112.08129-combinatorial-tilting-mutation.md` already has steps 1-7. FINDINGS l. 2599 and EXPERIMENTS l. 3749 say `procedure.py` is the paper. EXPERIMENTS l. 127 (E-187) says "the paper was not re-read" and l. 17 lists "comparison with steps 3-7 (scholar)" as open; this closes that. The Ladkani 1001.4765 gap is already recorded. The new content is the version, the theorem numbers (once corrected), the step-by-step match, the two (three) printed slips and the monomial-parent exactness of the gate. That is a modest but genuine literature check.

## Evidenced?

Mostly yes: line numbers and quoted statements are specific and verifiable. Weak points: numbers "counted by hand" (and one is miscounted); "v2 UNVERIFIED as byte-identical" is honest; the paper-picture match is by eye; the op/left-right convention is explicitly left open. The Next section's claim that E-187/E-205 "can say matched to the printed steps 3-7" is fair for steps 3-7 but E-187's reduction step is not covered, as the author says.

## Scope

Title: "states steps 1-7 and no LNA rule; ... gate is exactly the paper's Theorem 4.3 hypothesis". Narrowed wording: add "for monomial parents, with Hom(P*[1], Lambda) read as {g : g.[alpha] = 0} (derived)" and "right mutation, convention unproved".

## Required for acceptance

1. Replace every "Thm 4.2" (Evidence table, Comparison, Draft) with "Thm 4.3". Do the same for any "corollary"/numbering references that depend on it. Add a `\newtheorem` shared-counter note so the count can be checked. One sitting.
2. In the title or first line of the Claim, mark the gate-equals-hypothesis statement as holding for monomial parents and as using the derived Hom description. One sitting.
3. State that the Sect. 6 comparison is visual, in the Claim and not only in the scope paragraph (it is in scope; move or repeat it). One sitting.
