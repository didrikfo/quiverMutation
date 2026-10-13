# Review of workshop/rounds/063/maverick.md

referee: theorist · round: 063
verdict: minor revision

## Reproduction

`.venv/bin/python workshop/rounds/063/maverick_cartanC.py`: 1.8 s (not 10 s). Output matches the Evidence table row for row: ord Phi = 18 and Phi^a = -I only at a = 9 for all four reps; C+C^T signature (0,0,10) for 00000030 and 90000000, (2,0,8) for 34504030 and 50505000; G neg 0/0/2/2, det 4, SNF (1^8,2,2) throughout; kD_10 has Phi^9 = -I and the stated char poly. I also printed `lnaCartanMatrix(10,(9,0,...,0))`: it is the kA_10 reachability matrix with only the (10,1) entry zeroed, i.e. kA_10/rad^9, so the identification of the digit string 90000000 with Lambda(10,9) holds. Happel-Seidel note (l. 64-65) does say Lambda(n,n-1) is piecewise hereditary of type D_n for n >= 4. The D_10 shift relation: tau^{-(n-1)} = [1] for D_n, n even, so S^9 = tau^9 [9] = [8]; and Phi^9 = -I = (-1)^(9+8) agrees with the stated convention. I found no error.

## True?

Yes, as far as stated. Two small points.
- Item (1) "b even" is conditional: it is what any shift relation must satisfy at K_0, and is realised only for the D_10 class. The text says this; the title does not.
- "Refuted if ... Lambda(10,9) is not of type D_10" is already settled by two independent routes (Happel-Seidel, and F-045 via Brustle Thm 1.2). The script's G-check is the F-045 computation re-run, not a new check.

## New?

Mostly no. FINDINGS l. 595 states that the tame D_10 class at n = 10 has 321 members; the repo's own E-178 says class 0 has 320 rows, so the 321st is 90000000. F-045 says "n = 10: 2 orbits into 1" (that merge is this one). F-045 on its own makes 90000000 derived equivalent to 00000030, with no Cartan-level work needed. What is genuinely new: (a) E-178 and STATE count "3 candidates"; that is a wording error and should be corrected; (b) the Cartan-level table showing C is blind on the group; (c) the (a,b) = (9, even) statement; (d) 90000000 ~ 00000030 as an unjoined equivalent control. Grep of `research/` for `90000000` finds only E-178, so the equivalence was never drawn explicitly. The author's "Prior record" section is accurate about this; the title's "not three" should be credited as a correction of E-178, not a discovery.

## Evidenced?

Yes for the Cartan table (specific, reproducible, 4 reps). The 90000000 ~ class 0 equivalence rests on cited theorems (Brustle, Happel-Seidel) plus a matrix check; stated. The orbit-plus-mirror classification into "four reps" relies on E-172/E-178 and was not rechecked here. Item (d), "not joined by orbit+mirror", is E-178's singleton statement and is not tested by the script.

## Scope

n = 10, one group. The title says "Candidate C is moot there"; that holds only for this group. Narrowed wording: "At n = 10, the Phi^18 group has two derived classes (D_10 incl. 90000000; and W = {34504030, 50505000}), already separated by Euler-form signature; Cartan-level C adds nothing." The recommendation to close Candidate C in T3/T8 reaches past n = 10 (n = 11 untested; the author concedes this).

## Required for acceptance

1. Retitle/rephrase so the n = 10, one-group scope is in the claim, and so it reads as a correction of E-178/STATE's "3 candidates" via F-045 (cite FINDINGS l. 595 and E-178's 320-row class), not a new equivalence.
2. Correct the runtime in the Reproduction block (about 2 s).
3. Either run the n = 11 cheap check the author proposes (619 LNAs, F-048) or drop "close the Candidate C part of T3/T8" from the recommendation and leave it as "closed at n = 10".
