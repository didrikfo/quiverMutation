# Review of workshop/rounds/069/scholar.md

referee: skeptic · round: 069
verdict: minor revision

## Reproduction

Both commands re-run (each under 1 s, not ~10 s). Output matched the tables exactly: sign script at 200 trials gives (0, 0, 7082) for the printed sign and (0, 788, 7082) for the control; lemma script gives s = -1: 610 agree, 0 disagree (337+182+75+16), s = +1: 330/280 (178+99+42+9 agree, 159+83+33+7 disagree).
Extra, not in the submission: 1000 trials give (0,0,36424) printed vs (0,4019,36424) control; 3000 trials give (0,0,109561) vs (0,12189,109561). Same seed 0 prefix, so these are more letters of the same stream rather than independent evidence, but no failure appeared.
I also re-read the .tex at the cited lines: Lemma l. 1055-1059 (with `\ta{\varphi}` typo), proof end l. 1104-1105 (underbrace `(-1)^{||p_2||} dp/alpha`, last line `+(-1)^{||p_2||}`), l. 1100-1105, Def. 3.2 l. 383-391, l. 1124 and l. 1140, Prop. l. 1109. These match the submission.

## True?

The claims survive. Checked by hand: Leibniz gives dp/alpha = (-1)^{|p_2|} p_2 (dp_1/alpha) with |p_2| the Q-degree, and ||p_2|| = |p_2| - 1 + c_i(t p_2) when s(p_2) = i, so the proof's last line differs from the statement by -(-1)^{c_i(t p)}. That agrees with the submission, and the "typo in proof, Lemma fine" reading is right.
d(alpha) = 0 by degree: alpha has degree 0, grading is non-negative, d has degree -1. Correct, and the paper's l. 1124/1140 do rely on it.

Defects found:
1. Wrong citation: the submission says l. 308 for "A the collection of all arrows of degree 0 which start in i". That sentence is at l. 321 (l. 308 is the Remark about R End of a silting complex). The grading citation l. 262 is fine (theorem statement, "concentrated in non-negative degree"). Also in claim 3(b) the "A not loops" fact is stated at the same l. 321, not l. 308.
2. Test self-consistency: the d² script and the lemma script are written by the same author who derived the signs, and both implement "Leibniz with ||.||-signs" via one shared M-degree function. A shared misreading of ||.|| (e.g. the c_i convention in `M.deg`, l. 44-50 of the script) would pass both. The control (+ sign gives 788 failures) shows the test has teeth against the sign, not against a wrong degree convention. Not a refutation, just a limit not stated.
3. Test class is narrow in a way the claim does not say: the random instance has fixed vertex shape (vertices 0..3, one fixed path 1->3->2 of degree-0 arrows, A drawn from {0->1, 0->2, 0->1 parallel}); degree-0 arrows into i come from 1,2,3 only; no degree-0 arrow 1->... that makes a degree-0 cycle through i, and d of degree-1 arrows is drawn from degree-0 paths of length <= 4. The submission's scope line says "4 vertices, 1-3 arrows in A" which is accurate; it does not say "fixed skeleton".

## New?

Grepped RETRACTIONS.md for 1504.02617 and Oppermann: nothing. E-196 (EXPERIMENTS.md l. 29-32) explicitly lists both observations (sign at l. 1105 vs l. 1058; d(alpha)=0 at l. 1124/1140) as "unresolved"; E-189 (l. 80) is the hypothesis list. So the new content is the resolution of those two plus the executable d² test, as the submission says. Not a duplicate. Note: E-196 also says the triangularity of the cancellation was checked only for one term; the submission rightly leaves Prop. 7.6 alone.

## Evidenced?

Mostly yes: counts, classes (end-at-i with start loop etc.), and controls are stated and reproduce. Missing: (a) which seed/RNG (seed 0, default argument) so the 200 is replicable; the submission says "200 random instances" without a seed, (b) the number of instances where A had 1, 2, or 3 arrows, (c) the claim "degrees <= 2" is stated honestly. The hand derivation in Claim 1 is the stronger evidence than the numerics and is complete enough to check from the text.

## Scope

Claim 3 ("all Sect. 8 uses is ...") is a reading claim over l. 1048-1144 and is stated as "used and visible", not "sufficient": fine. The title's "a d^2 test passes" is accurate for degree <= 2 arrows and random d; it should not be read as evidence for Thm 8.1 (the submission says so). Narrowed wording not needed beyond the two fixes below.

## Required for acceptance

1. Correct the citation l. 308 to l. 321 (both in Claim 2 and Claim 3(b), and in the Next bullet "checked against l. 262, 308").
2. State the seed (0) and the fixed skeleton of `random_instance` in the Scope line.
3. One sentence acknowledging that both scripts share the `M.deg` convention, so they test the signs given that reading of ||.||, not the reading itself. (Optionally: one hand-checked instance of ||.|| against l. 262-296 definitions.)
