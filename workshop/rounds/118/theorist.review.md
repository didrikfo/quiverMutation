# Review of workshop/rounds/118/theorist.md

referee: experimentalist · round: 118
verdict: minor revision

## Reproduction

Re-run (all matched the submission's numbers):
- `theorist_letters.py` (0.04 s): printed letters of M4, M4', M5, M5' as described.
- `theorist_m45.py` (0.05 s): 11 + 15 identities `[ok]`, "M4 True M5 True", real g raises on nothing, planted g raises on 3 of 3 in each.
- `theorist_bardzell.py` (0.6 s): A3 (1,0,0), A3/rad^2 (1,0,0), Kronecker (1,3,0), square ab=0 (1,1,1), square ab=cd=0 (1,1,0). Same.
- `theorist_nochildren.py` (10 s): 75 records, 36 / 18 / 2 split by cycle rank 1 / 2 / 3, NO children steps 57 and 61, both cycle rank 3. Same. NOTE: it loaded the author's existing `/tmp/skeptic_c1_13_ck1.pkl`. I did not rebuild the walk (the 9 + 2 min `069/skeptic_c1_13.py` runs); the 23 000 expansions / 75 J != 0 steps are therefore unverified by me.
- `115/theorist_steps.py` (0.3 s): 13 of 13 YES, corruptions 69 / 30 / 13. Same.

Beyond the author (item 13a): ran the Bardzell code on tree quivers with relations, which the submission cites only as "classical, source not checked": linear A_n/rad^2 for n = 5..8 and six branching trees with monomial relations (one to four relations, including rad^2 = 0 on a two-in two-out tree). All gave HH = (1,0,0,0,0). No counterexample to "tree => HH^m = 0, m >= 1" in this sample.

Beyond the author (item 9): by hand I re-expanded the M4 change of generators from the printed letters. d sg = sp v1 + r22 w - w fbar, d r22 = w u2, d sp = w u1p, d fbar = u1p v1 + u2 w all follow (I did the sums, including the c1, c2 cancellations in sp). Consistent with the script.

## True?

I found no error in the M4 argument as written. Points the author should state, none of which I could break:
1. The script checks the algebra of the substitution (that the new generators' d equals the claimed words). It does not check the three structural claims that carry the proof: (i) E = <t1,p1,s11,q11,s21,q21> is closed under d and kM = B * T(E) as a free product (needs that no generator of B has d landing in E, and that the substitution is invertible over the base, not only triangular); (ii) I = (E) + (Theta', d Theta') exactly (not merely contains); (iii) the g-filtered quotient B/(Theta', d Theta') is a filtered dg algebra with gr equal to the Lemma A quotient even though d fbar != 0 (the text says "D = 0 on the arrows" while fbar in Y has d fbar = u1p v1 + u2 w; true at gr level since g drops, but that is the point to say in the Lemma A match). Step 2 rests wholly on Lemma A with X1 = {w}, k^3 base; E-294 used k^2. Lemma A is stated for semisimple base in E-240 (as I read the cite in the submission), so this is probably fine, but "same use as E-294" should be checked against Lemma A's hypotheses for a base with three idempotents and an arrow i->k, k->j (M5), not asserted.
2. M5: the single step "the M4 shape verbatim" is the part least backed by the script: it checks 15 of 15 identities and acyclicity, which covers the substitution. The one non-obvious new item is that the extra arrow e acts on the left of bbar, so e is in no X / Y set; the text should say e is in B but outside Lemma A's X1, Y and why it does not break "one Theta for each composable pair (w, y), and no other" (words w e ... ? e starts at k, w ends at i, so w e is not composable; state it).
3. 13a: "none of the 56 refused steps has a tree child" holds (cycle rank >= 1 throughout). The inference "HH test never met the region where it must say YES" is fine. But "the 2 NO children are not an exception to the tree theorem" is trivially true once they have cycle rank 3; it is not evidence about HH^1 = 1 being a correct NO. What the data show is that parent and child, both with cycle rank 3, differ in HH^1 while the step was refused; whether the refusal is correct is separate and the submission does not claim it.
4. 13c is a bookkeeping fix. E-306 already carries "78 candidates, not 84" and "34 done in round 114, 44 now" (research/EXPERIMENTS.md E-306 Scope). So the only new content is the 69 / 30 / 13 label ("corruptions" row of E-305: numbers were flagged missing in E-305's Scope). That is reproduced.

## New?

grep of research/ for M4, M5, "loop", "cancel", "Lemma A", "pivot":
- E-284 (ten loop models, computation to weight 12, includes M4, M5), E-290 (M1 proved; text says "M4 shared pivot, M5 product pivot not covered"), E-294 (M2 proved by g-filtration), E-240 (Lemma A), E-279 (loops force a cycle). The M4 and M5 proofs are new; E-290's own caveat names exactly these two cases. HYPOTHESES.md H-023 (line 37-38) lists "loops of degree >= 1 at i" as an open joint (E-273) and mentions only M1/M2 as proved; this closes two more models, not the joint.
- 13a: E-170 / E-172 already record HH trivial on LNAs (the submission says so). The cycle-rank tabulation of the 75 steps and the 2 NO children is not in E-305. New.
- 13b: HH^2 = 1 for the ab=0 square is new as an independent method; E-305's Scope lists it as "unchecked".
- 13c: not new beyond the 69/30/13 reading (E-306 has 78).
- Proposed H-030: highest H in HYPOTHESES.md is H-027. Numbering H-030 skips H-028 and H-029; check whether those are taken by other rounds' proposals before the keeper assigns the number.

## Evidenced?

Item 9: specific and checkable, with identifiers of each pivot and substitution. Power line is honest (planted g control caught 3 of 3 per model), but that control only shows the g-check can fail; it does not test the three structural claims in "True?" 1. Item 13: counts and the cycle-rank table are given; the walk itself was not rebuilt by me. The 13b hand computation and the code agree, and the controls (ab=cd=0 gives (1,1,0)) show the code can return a different answer. The "Power" line for 13b is adequate. The statement "E-284 H(I) = 0 for M4, M5 to weight 12" is cited, not re-run; it is E-284's.

The proposed H-030 is broader than what was proved. It is stated for any number of out-arrows and paths r_p of any length. M4, M5 have 2 out-arrows and r_p of length 1 (and one with an extra arrow e). The "Status" sentence is honest (SUPPORTED by M1, M2, M4, M5 only), so this is acceptable as a hypothesis, but the pivot condition "one r_p ending in an arrow b with b-bar linear in d theta" is the thing that made M5 work only after a product-pivot fix (t_m' = t_m - e t_a). The general statement should say whether that fix always exists.

## Scope

Title says "(proved, ...)" for M4 and M5 "of H-023". The Claim line correctly narrows to the formal dg models and says H-023 is not touched. The title should carry that: "formal dg models M4, M5". "E-305's three open points are closed" is accurate only as 13a (data stated), 13b (second method), 13c (bookkeeping); 13a/13c do not bear on whether HH is a useful test. Narrowed wording: "kM -> N is a q.i. for the formal dg models M4 and M5 for all weights, given Lemma A for a three-idempotent base; HH^2 of the ab=0 square is 1 by two methods; the n = 7 class 1 refused steps have cycle rank 1-3 in every case."

## Required for acceptance

1. State, in the proof, the check that kM = B * T(E) and that I = (E) + (Theta', d Theta') (not only the d-formulas the script verifies); one paragraph, or extend `theorist_m45.py` to compare the dimension of kM and B + T(E) pieces to a small weight. Doable now.
2. Say explicitly that Lemma A's hypotheses (E-240) hold for the base k^3 of M5 (no new assumption on idempotents) and that fbar in Y with d fbar != 0 is handled by the g-filtration, with the composability of w and e addressed. Doable now.
3. Retitle/narrow per Scope; drop or reword "not an exception to the tree theorem" in 13a so it is not read as evidence about the NO. Doable now.
4. Check the H-number (H-028 / H-029) with the keeper, and state in H-030 whether the product-pivot substitution of M5 is claimed for all models of the form. Doable now.
5. The walk rebuild (`069/skeptic_c1_13.py`, two runs) was not re-run by me: `theorist_nochildren.py` read the author's /tmp checkpoint. A skeptic or the toolsmith should rebuild it. [next round]
