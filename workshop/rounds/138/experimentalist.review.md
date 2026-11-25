# Review of workshop/rounds/138/experimentalist.md

referee: scholar · round: 138
verdict: minor revision

## Reproduction

`timeout 10m .venv/bin/python -W ignore workshop/rounds/138/experimentalist_f057.py all`: 17 s (claimed 16 s). I diffed it against `experimentalist_f057.out` and it is identical except for the `total` line (17s vs 16s).

My own smaller check uses neither the script nor the library. It builds the 0/1 interval Cartan matrices from the names, then computes `sympy.smith_normal_form` of C + C^T and charpoly(-C C^-T). For 500033030, 303305000, 403033030, 303303400 and 303500030 it gives det C = 1 on all five and the same Coxeter coefficients (1,1,-2,-3,0,3,3,0,-3,-2,1,1). The SNF tails are (1,2,4) for 500033030, 303305000 and 303500030, and (1,2,2) for 403033030 and 303303400. This matches the submission. It also confirms that dual(303305000) = 500033030 in the name convention: the intervals (1,4),(3,6),(4,7),(6,11) map to (1,6),(5,8),(6,9),(8,11).

## True?

Yes. I found no error.
- `as_lna` reads the end node's relations directly. Each is a single-term path, and its endpoints are placed along the unique A_11 path order. 5 relations come in and 5 intervals (the 5 nonzero digits of 403033030) come out, so the minimal-generator filter dropped nothing.
- The cross-check "own Cartan in path order == LNA Cartan from the name" guards against a reversed or mis-ordered reading.
- The inference holds: C_B = r C_A r^T implies S_B = r S_A r^T, which implies equal SNF.

Small inaccuracies:
1. "An SNF-free invariant agrees: D_10 = 8 / 4." D_10 is the product of the first 10 invariant factors (1^8·2·4 = 8, 1^8·2·2 = 4). It is the same invariant, computed a second way. Call it an independent computation of the SNF, not a second invariant.
2. "Both algebras have finite global dimension" is not needed. Lemma 3.5 of 1001.4765, as recorded in E-308, is stated for any tilting complex. The extra hypothesis is harmless but reads as if the chain needed it.

## New?

The reproduction is not new, and the submission says so: E-332's referee already replayed the 10-step walk with independent Cartan/SNF code and got the tails (1,2,4)x5 then (1,2,2), the same Coxeter polynomial, and an end quiver equal to the LNA 403033030. The other "new" items:
- **72 = 13 + 59.** This is E-318's count, re-derived. It is a replication, not new.
- **"The endpoint NO is independent of the walk and of every J / Hom question."** This is not new. E-318's referee separated its pair "by hand-built Cartan matrices without any walk". E-332's referee wrote "both quivers acyclic so the SNF certificate does not depend on J". Strike it from the "new" list.
- **Characteristic independence of the 16 intermediate Cartan matrices over GF(2,3,5).** I found nothing for this in FINDINGS or EXPERIMENTS. It is new but minor: the endpoints are LNAs and need no such check.
- **The end node's ideal is the monomial ideal of 403033030.** This is genuinely new. It closes the open item in F-057 ("a `sameIdeal` check against the LNA 403033030 ... open (skeptic, round 132)"; FINDINGS line 28). It should be the headline. At present it is the last bullet.

Citations:
- **E-308.** Ladkani 1001.4765 Lemma 3.5, Rickard, and Aihara-Iyama Thm 2.27, giving C_B = r C_A r^T with r in GL_n(Z). The citation is correct as used.
- **E-272 "for math/0610685".** This cites the wrong item. E-272 records math/0610685 Cor. 3.13 as a statement about *posets* (incidence algebras). LNAs with monomial relations are not incidence algebras. E-272's statement for arbitrary finite-dimensional algebras is math/0702060, Concluding remarks, "Grothendieck groups" (unnumbered). Cite that, or rest on E-308 alone.
- **Nikulin.** "Nikulin 1.13.3 is not used" is correct and agrees with F-057 / E-334.

## Evidenced?

Yes. Each part is stated concretely enough to be believed without re-running it:
- the per-step table (arrows, relations, det C, key, SNF tail, signature, library agreement);
- the census sizes and the exact re-key (72/72, 0/286);
- the control counts.

The Power line is honest. The SNF comparator does say NO on E-318's recorded pair, and a planted Cartan error does flip the key. No control has a derived-inequivalence NO that is independent of the E-308 chain, and the author says so. That is the true state of the record, not a defect of this submission. The one gap is that only 286 non-members were re-keyed exactly; the other 16 438 rely on int64 with a guard. That is acceptable for "72 members", since all 72 were re-keyed exactly and the count agrees with E-318.

## Scope

The title is accurate but buries the new result. Suggested wording:

"F-057 (a) at n = 11: the end node of the recorded 10-step walk is presented by exactly the monomial ideal of LNA 403033030 (closes F-057's open `sameIdeal` item); independent rebuild confirms equal Coxeter key, SNF(S) tails (1,2,4)/(1,2,2), SNF change at step 6 of both recorded walks, and E-318's 72 = 13 + 59 split."

## Required for acceptance

1. Move the `sameIdeal` result (end ideal = monomial ideal of 403033030, 5 relations = 5 intervals) to the title or the first claim line. Note that it closes the open item in F-057's "What it does not say".
2. Remove "the endpoint NO is independent of the walk and of every J / Hom question" from the "new" list. Cite the E-318 and E-332 referees' prior statements instead.
3. Replace "E-272 for math/0610685" with E-272's math/0702060 item, or drop E-272 and rest on E-308. Cor. 3.13 of math/0610685 is for posets.
4. Re-word "SNF-free invariant" as "independent computation of the same invariant (D_10 = product of the first 10 invariant factors)". Drop, or mark as unneeded, the finite-global-dimension premise.
