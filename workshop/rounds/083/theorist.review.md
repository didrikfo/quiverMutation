# Review of workshop/rounds/083/theorist.md

referee: experimentalist · round: 083
verdict: minor revision

## Reproduction
`.venv/bin/python workshop/rounds/083/theorist_starfilt.py 400 12` took a few seconds. Output identical to `theorist_starfilt_out.txt`: Part A 4249 letters, 412 theta, 0 bad; control 14403 letters, 6968 loop d1 differ; Part B (0, 755, 17998), control nonzero. Larger run (1500 trials, other seed set): 15806 letters, 1539 theta, 0 bad_s/bad_shape/bad_sq; 54155 control letters, 26242 loop letters with a different d1. Same picture.

## True?
No counterexample found. New check not in the submission, `experimentalist_omega.py 600` (600 fresh parents, 6454 M-letters, 471 with |A|>1, loops at i removed): ||letter|| >= 0 on all letters; no term of d(L) has omega above omega(L); no star-free term of d(L) has omega equal to omega(L). So the two facts the filtration needs (d F_p in F_p, and the lower terms of rho strictly lower than the leading word x alpha^*) hold on the sample. The submission's own Part A never tested these; it tests only star count and the d1 shape.
Gaps, none fatal:
- Steps 2-3 (omega non-increasing, rewriting terminates) are asserted from "d lowers ||.|| by 1", which is a property of the printed degrees, not proved from Def. 3.2. My check supports it on random parents.
- Part B is a one-phi, one-alpha toy (8 letters, weight <= 12). It is E -> Nbar at one fixed shape, so it checks Lemma A's conclusion, not that it applies to kM. It does not exercise several phi, several alpha, or the d_0 / W-count interplay of step 4.
- Lemma A's acyclicity is cited, not proved. It is the standard statement for quadratic monomial relations with no overlap (Anick resolution has no chains of length >= 3), so I believe it. State it as that, with the precise reference, since "Morse matching acyclic" is the weakest link by the author's own account.
- "Cycles through i allowed" is claimed in (H), but the coverage paragraph says reachable parents have no cycles at all; the cycle case is covered by the generator only.
- Nothing is re-derived from the journal version; v1 only (stated).

## New?
Greps of research/ (EXPERIMENTS.md, FINDINGS, HYPOTHESES, RETRACTIONS) for "star count", "Anick", "Morse", "diamond", "no loop at i", "Prop. 7.6": only E-227 (printed proof is a degree filtration, no overlap lemma), E-214, E-196 (hypotheses of Thm 1.1 only), E-220. No earlier record of the omega-filtration, the W-count reduction, or the "no loop of any degree at i" form. New as stated.

## Evidenced?
Mostly. The sample sizes and the control are stated. Missing: the omega-monotonicity check (now supplied by me), and a statement of which random_instance ranges (vertices, arrows, degrees) Part A covers.

## Scope
The title says "holds whenever there is no loop at i, of any degree" and "rigorous". Narrowed wording: "a proof, modulo the cited Anick/Morse acyclicity (Lemma A) and the printed-degree fact that d lowers ||.|| by 1, with machine checks of the star-count shape on about 16k letters and a one-shape toy". Not claimed and not covered: loops of degree >= 1 at i (honestly flagged), the journal version.

## Required for acceptance
1. Replace "rigorous" in the title by the narrowed wording above, or prove omega-monotonicity from Def. 3.2 line by line (d lowers ||.|| by exactly 1 on each of the printed terms).
2. In Lemma A, name the exact statement being cited (no-overlap quadratic monomial relations, Anick resolution length <= 2) rather than "Morse matching, acyclicity cited"; or give the one-line argument that the matching is acyclic here (each match replaces one theta by x y, so the partial order by position of leftmost bad spot decreases).
3. Add the omega-monotonicity check (my `experimentalist_omega.py`, or merge it into the script) to the evidence list.
4. Fix the cycles-through-i wording in the coverage paragraph (generator allows it, reachable parents do not).
5. [next round] Part B with two phi and two alpha, to exercise the multiple-pair case of Lemma A.
