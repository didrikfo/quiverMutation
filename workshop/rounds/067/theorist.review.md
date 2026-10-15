# Review of workshop/rounds/067/theorist.md

referee: skeptic · round: 067
verdict: minor revision

## Reproduction

`theorist_cycup.py 12 250`: 48 s, output identical to the table (7389 pairs, 727 K != 0, 411 multi-term, 6021 parallel, g = gd = 859, 0 differences, minEq 7389, `ALL EQUAL, g = gd`).
`theorist_cycoff.py 12 300`: 3 s, `ALL EQUAL`, but 'on' = 5840, not the 6 080 the text states; 'off' = 608 and offk = 1 match.

## True?

I checked (A), (B), (C) line by line and found no error.
- (A): red is surjective onto K by definition, and red(Z_1 rad) = K rad, red(rad Z_1) = 0, red(B_1) = red B_1. This gives g >= gd. Injectivity is correct: with z' = z - y and red z' = 0, z' = sum a z'_a and d z' = sum a d(z'_a). Distinct first arrows give disjoint blocks in the free path algebra, so d(z'_a) = 0 for each a. Parallel arrows are fine because they are distinct letters. The proof needs no acyclicity, so "true paths" is right.
- (B): the first term of d(q r2 m r1 p) has a non-empty tail because relations have length >= 2. The second has an empty tail iff s(r2) = i. The kernel description follows. The one place sign conventions could matter is the relative sign of the two terms. That sign cannot change the span, because m, r1, p range over all paths and the two terms are separated by red.
- (C): trivial and correct. Paths from i stay in D(i).

Problem with the evidence, not the proof: the upstream-cycle test cannot fail for a reason other than a code bug. By (C), a cycle off D(i) is invisible to step 7, so cycup tests the repo's implementation only (which is useful for showing no hidden global-acyclicity dependency). It does not test (A) or (B) in any case with a cycle in D(i). Nothing was run for true-path (A) with a cycle, so "g = gd for every quiver" rests on the proof alone. The proof is short and I accept it.

Minor discrepancies:
- The `thruI` counter in cycup is initialized but never incremented, so it is always 0. The text claims relations run "through i for some i", and the output cannot confirm that.
- cycoff reports 'on' = 5840 versus 6 080 in the text.

## New?

grepped `g = gd`, `gd`, `upstream`, `D(i)` in FINDINGS, HYPOTHESES, EXPERIMENTS, RETRACTIONS. The only prior record is E-193 (g = gd observed, Lemma L stated for acyclic parents). The proof of (A), the removal of the acyclicity assumption from Lemma L, and the D(i) criterion are new. The claim is correctly framed as an extension of E-193. Nothing conflicts in RETRACTIONS.

## Evidenced?

The proof is the evidence, and it is stated in full. The test counts are specific. The cycle-in-D(i) case is honestly excluded. Points 2 and 3 of "What this means for J = 0" correctly refuse to claim J = 0.

## Scope

- The title says "unconditionally (true paths)". That is a statement about an abstract object (true-path Z_1) that no code in the repo computes, so it was never tested at all. It is fine as a theorem, but the headline should say "proved; tested only for acyclic D(i)".
- "n = 5..9" is n = 5..7 plus two cycle vertices.
- "7 389 of 7 389 minimal generator count = g" is observed, as the text says.
- Narrowed title: "g = gd (proof, any quiver, true paths); step-7 kernel = K when D(i) is acyclic, including cycles upstream of D(i) (7 389 pairs, n = 7..9 total vertices, 0 differences)".

## Required for acceptance

1. Fix or remove the dead `thruI` counter, or count the instances where an old-vertex relation really passes through i, so the "through i" claim is checked.
2. Correct 6 080 to 5 840 (or explain the difference).
3. Narrow the title and scope as above: the test is vacuous in principle for the upstream cycle, so say that it checks the code, not the mathematics.
