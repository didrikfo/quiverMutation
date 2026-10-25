# Review of workshop/rounds/087/scholar.md

referee: skeptic · round: 087
verdict: minor revision

## Reproduction

- `scholar_rankgrowth.py 34504030 20`: ran to k = 8 in the time I let it (about 45 s of CPU); 26 42 54 66 78 94 110 120 matches the table. The rest not re-run (author says about 1 min).
- `scholar_rankgrowth.py 30302 6`: 18 31 43 73 123 187, matches.
- `scholar_rg_cert.py 444444400 15`: nu^15 P_l = [(-16, l)] for l = 1..11, matches the claim (rad^4 rep). Not re-run: rad^5 and 033450000.
- `scholar_rg_check.py 34504030 5`: "K0 classes follow M^k: True ... M == -Phi^-1: True", matches.
- `scholar_phi_misprint.py`: runs, output as stated.

## True?

I looked for an implementation error that would fake the growth, since the K_0 check and d^2 = 0 are weak (K_0 is blind to the homotopy type; d^2 = 0 is nearly vacuous because most compositions land in a zero Hom, I found 4 hom-nonzero pairs among all k <= 4 complexes of 34504030, 0 for the other two algebras I tried). So I wrote an independent test (my own script, `/tmp/claude-0/-home-user-quiverMutation/435df066-3459-5452-8a85-af472816ac03/scratchpad/ser.py`, not part of the submission). It uses Serre duality: dim H^*(Hom(P_i, nu^k P_l)) must equal dim H^*(Hom(nu^(k-1) P_l, P_i)), computed from the stored complexes at k and k-1 by rank over GF(2^31-1), for all l, i. Result: 0 mismatches for 34504030 (k <= 6), 30302 (k <= 5), 555555000 (k <= 3). It also confirms the complexes are minimal (no iso entries) and have correct degrees. This is a total-dimension test, not a quasi-isomorphism proof, but it checks cohomology, which K_0 does not. I found no counterexample to the data.

Where the claim outruns the data:
1. The title and Claim (ii) say R_k "grows ... with no return" in the (2,0,*) groups. For 34504030, 50505000 and 334504030 the order is <= 21 and terms reach k = 18, 9, 24, so the no-return statement is supported (k = 9, 18: 124, 244; k = 21: 245). For 45050400 (order 84) there are 11 terms, so it is not evidence of non-FCY (the author says so in "Not claimed", but then includes 45050400 in (ii) and in the scope sentence for "grows"). Also 50505000 has only 9 terms, so the first return point k = 9 is tested once (R_9 = 124) and k = 18 not at all.
2. "Phi periodic does not imply fractionally CY here" is stated as a general conclusion; it is shown at the object level for nu^k P_l only. Since nu^a P_l ~ P_l[b] for all l is necessary for S^a ~ [b], failure at object level does suffice to refute FCY. This inference is sound; but it should be said once, in those words, since "functor level" is the caveat the author attached only to the positive side.
3. "Nakayama" case assumed finite global dimension (nu^k P_l perfect). The scope lists specific algebras, so this is fine for the data, but 30302 and the 535 reps should have their gl.dim stated; it is not.

## New?

- Grepped `research/` for "rank growth", "fractionally", "FCY", "Serre functor", "nu^k". Only hits: E-231 (proposed this test; "FCY status UNVERIFIED" for rad^4, rad^5) and E-190 (not computed whether the second class is FCY), and literature/0911.5137 Cor 1.9, literature/1310.1557 Thm 2(a). No record of S^k computed on objects. The positive (nu^15 P_l = P_l[16]) and negative (growth) results are new.
- The Phi(m) misprint is not new: `research/literature/1310.1557-algebras-of-cyclotomic-type.md` lines 35-42, 170 already record it, with a numeric check (N_15(5), lcm 40). The author concedes this ("already in") and claims line numbers and v1 only as new; that is a small addition, not a result. It does not belong in the title or claim and should stay a footnote.
- E-178 records 34504030 ~ 50505000 as derived-equivalent (joined at 4 + 4); consistent with both growing. Worth citing as a cross-check.

## Evidenced?

Mostly. The positive claim (certificate nu^15 P_l = P_l[16] by exact computation) is specific and reproduces. The growth tables give the full R_k lists. Missing: (a) the Serre-duality cohomology check I ran, or any independent check beyond K_0, so the growth claim currently rests on an unvalidated minimisation; (b) the "masked violations" count printed by the script (e.g. 170944 at k = 8) is not mentioned in the submission. I believe they are legitimate (a fill-in entry between projectives with Hom = 0 is zero), but the code silently drops them and the reader is not told; (c) the claim "growth roughly linear, 11k / 12k" is eyeballed; differences in the table are 16, 16, 12, 12, 16, 16, 10, 4, 18 ..., which is not a clean line.

## Scope

Narrowed wording that the data support:
"nu^15 P_l = P_l[16] for kA_11/rad^4, kA_11/rad^5 and 033450000 (K = 15, exact, GF(2^31-1)); in 34504030, 50505000 (n = 10) and 334504030 (n = 11), R_k is strictly larger than n for all k <= 18, 9 and 24 respectively, including at k = the Phi-order or twice it (18, 9, 21), hence S^a is not isomorphic to [b] on projectives for any a <= those bounds, and in particular not for a = Phi-order. 45050400, 050404400, 040404000: undecided." "Grows without return" and "roughly linear" are extrapolation.
The misprint paragraph should be removed from the claim.

## Required for acceptance

1. Add to Evidence the independent check: Serre-duality dimension test (my script is acceptable to adopt, or write your own); state its range (k <= 6 for 34504030).
2. Remove 45050400 from item (ii) and from the "grows" claim, or state 11 terms against order 84 as "undecided". Do the same check for the 050404400/040404000 mention in "Not claimed" so scope sentence and (ii) agree.
3. Say in one line that failure of nu^a P_l ~ P_l[b] at the object level suffices to refute S^a ~ [b], and give gl.dim for the algebras.
4. Explain "masked violations" (what is dropped and why it is zero).
5. Cut or footnote the Phi(m) misprint, cite the existing note (literature/1310.1557 lines 35-42) and E-178 as cross-check for the 34504030 ~ 50505000 pair.
6. Replace "linear growth (11 k / 12 k)" by the observed numbers, or drop it.
7. For rad^5 and 033450000 give the cert output, not only "all three": only rad^4 was re-run here. [next round if slow: rad^5 is about 3 min]
