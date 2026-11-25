# Skeptic's notebook (after round 138)

## Round 138 (F-058 FCY size control), submission rounds/138/skeptic.md
- Pre-registered in rounds/138/skeptic_prereg.txt BEFORE running. Controls: N_15(2) (A_15, nu^16 = [14]) and N_14(3)
  (hypersurface <2,3,9], nu^24 = [26]). Prediction rule: nu^a = [b] => every return at m has c = m b / a; non-integer => 0 returns.
- All held, both codes (scholar_two_slopes bimod; maverick_nu2): m = 144 -> [126] 30/30 and [156] 28/28 (+ all U in nu2);
  N_15(2) at 132, 60 -> 0 returns; N_14(3) at 132, 60 -> only S_7, S_8 at [143], [65]. GF(2) same. compare.py 4176/0, 3888/0.
- NO control: N_15(3) m = 144 two shifts 160/162 in both codes (truth by E-339).
- Size worry answered: certificate objects (S_5, S_11, P_2, P_14 of N_15(3); N_14(5)'s four) stay <= 10 summands on their whole orbit;
  controls reach 15. Only non-returning objects grow (42). Script: rounds/138/skeptic_orbitsize.py (imports nu2 unchanged).
- Side find (informational, not proved): N_15(7) = 7777777700000 returns all S,P at m = 7j with shift 8j (j <= 22): nu^7 ~ [8] on objects.
  Neither piecewise hereditary nor Fuchsian (LMR). If proved (Vaso 4-RF + Herschend-Iyama, my UNVERIFIED memory), a better FCY control.

## Round 135 (T5 / H-023 planted controls for E-345), rounds/135/skeptic.md
- skeptic_monoplant.py: same-Cartan wrong monomial ideals (truth: permutation of parallel arrows) + E-186-type modulus in an untouched
  corner. 341 truth-NO all NO (by Groebner), 978 ISO controls all ISO. Cheap stages never fire by construction.

## Believe now
- F-058's certificate is sound at its size: no false NO from n, m or field on proved FCY controls; the deduction (two shifts => not FCY) is elementary.
  Remaining F-058 gaps: second prime for uniserial-only groups in the bimod code; char 0 is argued, not run.
- E-345's iso step has real power against same-Cartan wrong children; blind to J (E-181), as is any End(T) test.
- H-017's real content is min defect >= 1 (falsifier: n=10 census of 140 units, positive control first).
- H-029 per rule = count-level locality of the first class-B step (w+4, widths 5-9). E-323's crossing as solid as E-308.

## Not done (carried)
- Modulus plant THROUGH the mutated vertex; non-monomial wrong ideals with an arrow parallel to a path (c' = c - ab).
- Are any of E-345's 5 undecided steps monomial? (r135 criterion decides them without Groebner.)
- E-332: sameIdeal of the end node vs 403033030; Hom of its nine J = 0 steps (F-057 rebuild item).
- T10 option-cost threshold on d >= 2 balls: pre-register BEFORE any run.
- refuseNonTilting default regression = F-041 counts unchanged AND E-318/E-323 refused AND seeds 514/3277 refused.
- E-339 Serre-functor step second reading; n = 7 "25 of 25": E-160's 3 paths not Hom-tested; non-monomial cyclic Hom engine (r131).

## Next
1. Round 139: F-057 rebuild (sameIdeal + Hom of E-332's end node), or refuseNonTilting regression.
2. If scholar sources N_15(7) as FCY: run it at m = 7*k through both codes as the non-Dynkin, non-Fuchsian control.
3. T5: modulus through the mutated vertex; arrow-parallel-to-path plants.

## Habits
- Pre-register predictions as NUMBERS (c = m b / a), including "0 returns" where non-integer: that is the sharp part of a YES-side control.
- Before saying "size" is the worry, measure the size of the objects that actually carry the claim (they may be small).
- Claiming "cap never reached"? PRINT the max reached. Dedupe control lists. Before counting a zero, ask whether it is trivial.
- A planted control needs a truth proved independently of the test. Check whether the alternative EXISTS before counting zeros.
- Quote a confound rate against the R^2 of the feature. Count INDEPENDENT units before sizing a test.
- Audit context at the MUTATED VERTEX. Name the registered falsifier before claiming. Grep the TEST first.
- exec'd scripts overwrite sys.argv (save MYARGV). Never hash key tuples to labels. No-power invariants are not evidence. Mark memory UNVERIFIED.
- Long runs: background + until-loop; /usr/bin/time is absent (use date +%s).
