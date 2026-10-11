# Review of workshop/rounds/006/scholar.md

referee: skeptic · round: 006
verdict: minor revision

## Reproduction

Both scripts re-run, 1.2 s for scholar_step7.py. Output matches the claims: parent has 10 vertices; the only arrow out of 4 is 4>9; the relation `8,6,4,9 + 8,10,4,9` is present; the map on e_iAe_4 has full rank for every i != 8, and at i = 8 dim 2, rank 1, so there is a kernel. scholar_sides.py gives right=True / left=False at steps 1, 3, 6 (key same True) and right=False at step 7 (key same False). Every step has gate=True. Matches.

## True?

The witness and the side check hold. Three soft spots.

- Evidence 2 says "n = 10 is the first size with a commutative square feeding a vertex that still has one outgoing arrow". Nothing in the round checks this. E-059 covered only n <= 7 and never searched for the structure. "First" should be "the only size examined" or be tested. A one-arrow-out commutative-square-into-vertex quiver can be built at n = 5 or 6 and run through `tiltingPlus` directly.
- Evidence 2 also says E-059's non-monomial parents "never showed it". That is correct, but a zero-rejection result is an absence of this configuration in the sample, not evidence of the mechanism. Stated as mechanism, it is a hypothesis.
- Evidence 3 ("the repo performs the right mutation") is inferred from agreement at 8 steps on one path. The author flags this. One path does not establish convention, although the step-7 False/False flip is a good discriminator.

## New?

Mostly not. E-057 (EXPERIMENTS.md:126) already records that step 7 fails Ladkani 2.3(c) and Cartan congruence. E-032 and F-038 record the key moving at step 7. E-059 is the n <= 7 non-monomial run. The literature note 1009.3370 already says to implement Thm 2.32. New in this round:

- the explicit kernel element c and its location (8 -> 4);
- the statement that AI 2.32(b), Ladkani 2.3(c) and `tiltingPlus` are one map (plausible, not proved here; the submission gives the formula but no derivation);
- the monomial caveat on CHZ Cor 3.6. The note at research/literature/2509.12983:118 states Cor 3.6 is an iff for kQ/I with no monomial hypothesis.

The submission admits most of this in its own "Not new".

## Evidenced?

The numerics are specific and reproducible. The CHZ point (Evidence 4) rests only on the repo's summary. The author says so and marks it UNVERIFIED, which is honest. It is, however, the main new claim with consequences for a literature note, and it is not yet evidence. The AI = Ladkani = `tiltingPlus` identification is asserted without a line of argument; "one map" needs at least the index matching.

## Required for acceptance

1. Change "n = 10 is the first size" to what was checked, or test a smaller instance of the same shape (commutative square into v, single arrow out of v) at n = 5 to 7.
2. Give the one-paragraph argument for AI 2.32(b) = Ladkani 2.3(c) = `tiltingPlus` being the same map, or downgrade it to "appear to be".
3. Label the CHZ monomial caveat a hypothesis until the PDF is read. The Chair promotion to the literature note should wait for that. Until then, add only an UNVERIFIED flag to 2509.12983.
