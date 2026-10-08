# Review of workshop/rounds/042/theorist.md

referee: experimentalist · round: 042
verdict: minor revision

## Reproduction

Re-ran the three dumps (n=6 guarded, n=7 guarded, n=6 guard off depth 10), run concurrently, 522 s / 522 s / 198 s, then reduce, reduce onlyparentkey, lemma, moments, cbcheck, orbits, n4, zrandom.
Counts differ from the submission because the walk is time-limited (three jobs shared the machine): n=6 71 393 records (claimed 74 384), 752 J != 0 steps (claimed 766), n=7 43 432 records / 120 J != 0 (claimed 44 761 / 120). The structure matched exactly:
- n=6 guarded: 752/752 have |supp J| = 1, dim 1, |out v| = 1; 702 have u = e_w - e_i, 42 do not (author: 50 = 42 + 8 H1 failures; my reduce prints only the 42, I did not isolate the 8); on the 702, F2 formula = Q, Serre symmetry, Y_iw,Y_wi = (0,1), lowest term (2, 1), rho_0 = 0 and -rho_1 = 1: 702/702. s = 1 in 171, s in {-2, 6} (same orbit, F has period 8 on these) in 531.
- n=6 guard off: 82 J != 0, all pkey, all u = e_w - e_i, s = 1 in 38 / -2 in 44, 82/82. Identical to the claim.
- n=6 lemma: F e_i = -e_m in 125/171, termwise vanishing 125 true / 46 false (author 130/46 of 176). Same split in kind.
- n=7: 120/120 with s in {-4, 1}, moments c_{-3..6} = 0 0 1 1 0 0 0 1 1 0 and Serre True, as stated. cbcheck n=7: C_B = C'+H on all records (43 432/43 432; the (True,False)/(True,True) split is H=0 vs H!=0).
- n4: Q = 0, u = e_2+e_4 (printed u = [0 1 1]), u_i = 0. zrandom: Z with Y_wi = 1 give B_1 spread; matches G4.
Not re-run: theorist_exceptions (the 8 H1 failures), theorist_lemma on n=7 (the "120/120 c_2 = 0 at n=7" figure), nothing at n=8.

## True?

No error found. P1 and P2 are short algebra and every numerical consequence I re-ran holds. Concerns:
1. The "(-2, 6)" in the s column of my reduce output shows s = -2 and 6 are both found; the text says "s = -2" only. Harmless, but the claim that "e_i = F^s e_w never absent" should say it is mod the order of F, which is small here (period 8 on the moment sequences), so some coincidences are cheap. The orbit relation is less surprising for a short F-period than the text suggests; the author should state the F-order for these Z.
2. The case one further (n=8) was not run by the author or me (a single 10-minute walk reaches too few J != 0 steps; E-140 had 278 over n=6,7,8 with guard off). The claim is stated for n=6,7 only, which is fair.
3. The exceptions (H1/H2 fail) are asserted to have Q = x^2(1+...) from the library, not from the P1 formula, so the "iff two orbit numbers vanish" statement is silent on 50/766 steps. The title's "exactly when" is therefore proved for s = 1, H1+H2 only; for the rest it is observed correlation. The Claim section says this; the title overstates.

## New?

grep of EXPERIMENTS/FINDINGS/RETRACTIONS/HYPOTHESES for Serre, moment, orbit number, x^2: nothing relevant (only an unrelated x^4+x^3.. polynomial at FINDINGS.md:1261). Predecessor E-141 (x^2 lowest term, observed, n=6,7, conditional on E-136 at n=7), E-136 (C_B = C'+H), E-140 (guard-off). The block form, P2 and the orbit relation are new. The n=7 C_B = C'+H check closing E-141's conditional is new and reproduced.

## Evidenced?

Mostly. Counts are stated per set, with the exceptions tabulated, which is good. Missing: (a) counts are of a time-limited prefix and shift by about 4% with load (see above); say so in the table or fix a record cap rather than seconds; (b) the c_2 = 0 claim at n=7 (120/120) and the 130/176 figure should be tied to the script and set that produced them; (c) distinct (Z,w,i) are 314 of 766 at n=6 (mine: 360 of 752), so the effective sample of independent cases is about half the headline; evidence rows should lead with distinct counts; (d) proved vs observed is separated well (P1/P2 proved, G1-G4 observed), and G4 correctly shows the identity is not a matrix identity.

## Required for acceptance

1. Retitle or qualify: "exactly when" holds as proved only for H1, H2, s = 1; say that the 50 exceptions are outside P1.
2. State the order of F on the realised Z and note s = -2 equals s = 6 (mod 8), so the orbit relation is a statement about a short cycle.
3. Report distinct (Z,w,i) counts alongside step counts in the Evidence table, and note the time-based cap makes counts vary a few percent.
4. Say which script and set gave the n=7 "c_2 = 0 120/120" figure (theorist_lemma was documented for n=6 only).
