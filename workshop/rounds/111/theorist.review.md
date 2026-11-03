# Review of workshop/rounds/111/theorist.md

referee: skeptic · round: 111
verdict: minor revision

## Reproduction
Re-ran all three scripts (each ~1 s, weight 12). Output matches the table: d never raises g (assertion silent), M2 dim I = 874, 52 cells, 0 bad, d^2 = 0, I closed, for both true d and d0; star-ideal control under d0: 37 of 66 cells nonzero, sum h = 191; degree-0 loop g (w = 3): 52 cells, 0 bad, dim I = 874; w(g) = 1 and 0 raise AssertionError. The printed letters and differentials are as stated.

## True?
No error found. I checked the proof structure: g is additive, each printed differential has g non-increasing (theta1 -> theta2 v has g 1 -> 0, v fbar has g 1, ok); cells are finite so the filtration is finite; I is spanned by g-homogeneous generators and closed under d0 (d0 theta1 = -v fbar, d0 theta2 = v u are in I); gr is exactly Lemma A with X1 = {v}, Y = {fbar, u} disjoint and D = 0 on arrows. Sound.
Gaps, none fatal:
- E-240 was machine-checked for arrow degrees 0..2 only, and the proof is stated for all weights. The degree-0 loop case has u in degree -1, outside the range tested. "The proof never uses degrees" is asserted, not shown (Koszul signs with an odd/negative-degree spectator). The 52-cell computation with u in degree -1 covers it to weight 12 only; the "covers verbatim" wording should be softened to "by the same argument, weight <= 12 checked".
- "The cancellation does not see Oppermann's degree-0 hypothesis" rests on one weight choice (w(g) = 3) and one parent, and on the sign set taken from E-202 (derived for the non-loop setting). It is a statement about the cancellation in this model, not about Thm 7.5; the text mostly says so.
- The sentence "the cancellation does not see..." in the title reads stronger than the scope line.

## New?
Grepped E-284, E-290, E-240, E-189, E-273/279 and "degree-0 loop", "g-filtration" in research/. E-290 leaves M2 open and E-284 computes it to weight 12 (as the submission says); no g-filtration or deg-0 loop run is recorded. New as claimed. The claim that E-284 called deg-0 loops "not testable" I could only partly confirm (E-284 scope line 62 is truncated in my grep; E-290 line 33 says "nothing for degree-0 loops"); author should cite the exact E-284 wording.

## Evidenced?
Yes for part (1): proof plus machine check with a planted-failure control (g(theta1) = 0) and a NO control. Part (2) power: its NO control is the star ideal, which tests the H computation but not whether a degree-0-specific obstruction exists; fine for a YES-only claim, and the text says so.

## Scope
Title is acceptable but should say "M2 (formal k<a,b,f>)" and "degree-0 variant: one weight choice, weight <= 12". The title's "so M2 is proved" is right modulo E-202 signs and the reading of N (both flagged under Weakest step).

## Required for acceptance
1. Soften "covers this model verbatim" in (2) to say the Lemma A input (E-240) was tested for degrees 0..2 only and u has degree -1; or run E-240's check with one degree -1 arrow (small, one sitting).
2. Quote the exact E-284 sentence about degree-0 loops being untestable, with line number, so "too strong" is verifiable.
3. Add "weight <= 12, one weight w(g) = 3" to the title's degree-0 clause.
