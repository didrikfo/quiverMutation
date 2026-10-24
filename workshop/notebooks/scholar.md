# Scholar's notebook (after round 085)

## Believe
- Candidate C = S^a(A) ~ A[b] (fractionally CY, exponent a derived invariant). 1310.1557 Thm 2(a) FCY => Phi periodic; no theorem in 0911.5137 or 1310.1557 gives periodic => FCY. Sharp Phi^a = (-1)^(a+b) I is ours (E-190).
- 0911.5137 (LaTeX checked R082): Thm A/B/C l. 272/363/526; Cor 1.2 l. 439, 1.3 l. 488, 1.9 l. 671 (silent at prime n = 11), 1.13 l. 967, Remark 1.14 l. 972.
- 1310.1557: Thm 1 l. 182, Thm 2 l. 217, Lemma 2.9 l. 407, Thm 2.10 l. 425 (its (1)=>(3) proof is wrong; true period lcm{m}), Thm 6.1 l. 1103. n = 11 groups 535/64/24/15; group 535 (kA_11/rad^4, rad^5) FCY status UNVERIFIED.
- arXiv e-print works from the cloud (curl -sL https://arxiv.org/e-print/ID; single gzip or tar.gz). Export API search works (--data-urlencode search_query).
- Older: AI Thm 2.31, 2.32(b); Oppermann 1504.02617 v1 Thm 1.1 / Prop 7.6; T5 gap = repo child = H_0(kM) (E-217).
- R085 (T5): E-235's "Lemma A" (Anick, monomial quadratic no-overlap => strongly free) has no theorem number reachable: Anick J. Algebra 78 (1982) 120-140 and Trans AMS 296 (1986) are not on arXiv. Nearest statement: Loefwall arXiv:1711.03644v2 Sect. 3, Dfn 3.2, Rem 3.4 (DGA: strongly free iff A<S;dS=w> acyclic), unnumbered para l. 730 (monomial criterion, wording "common beginning and ending" must mean overlap); field T(V) only, not path algebras over semisimple R. Halperin-Lemaire Thm 3.3 (via math/0601421) is Lie-algebra inertness, not Lemma A.
- R085: F-047's 2400230 / 3304000: 3304000 is NOT almost separate (Fosse's theorem does not apply to it); named through lnaMoves mates 3004000, 3024000 (almost separate, P^(1,2)_(1,1,2)); 2400230 ~ those only by Fosse thm:QuipuToAn "if". classPreservingOrbit closures do not join them. So "naming only" = lnaMoves link (computed) + Fosse (cited).

## Did
- R001..R082 audits and summaries (see lit summaries). R085: Lemma A source hunt (scholar_quipu.py in rounds/085), F-047 basis.

## Next
- If someone has Anick 1982/86 in hand: the theorem number of the monomial/no-overlap case.
- Check quipuForAlmostSeparateLNA against Fosse's theorem statement (only docstring read in R085); implement Fosse Cor ops fully so 2400230 and 3004000 join (positive control for the 3060000 vs 3304000 NO).
- Source for FCY of kA_n/rad^r (Herschend-Iyama, Lenzing-de la Pena).
- Lessons: grep summaries first; counters: remark/prop/cor share the lemma counter; a paper's remark ("Observe that...") is not a theorem, do not number it.
