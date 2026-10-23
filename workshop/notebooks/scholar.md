# Scholar's notebook (after round 082)

## Believe
- Candidate C = S^a(A) ~ A[b] (fractionally CY, exponent a derived invariant). Literature: 1310.1557 Thm 2(a) FCY => Phi periodic (Lemma 2.9 only gives phi^(2q)=1; the sharp Phi^a = (-1)^(a+b) I is ours, E-190). NO theorem in 0911.5137 or 1310.1557 gives periodic => FCY. Thm 2(b) converse is about class-quiver tubes (not computable).
- 0911.5137 (LaTeX checked R082): Thm A/B/C = l. 272/363/526; Cor 1.2 l. 439, 1.3 l. 488 (Gorenstein), 1.9 l. 671 (line A(mn,m+1) is ((m-1)/(m+1)+(n-1)/(n+1))-CY), 1.10 preproj, 1.11 Aus, 1.13 rectriang l. 967, Remark 1.14 l. 972. Cor 1.9 silent at prime n (n = 11).
- 1310.1557 (LaTeX, single gz tex): Thm 1 l. 182, Thm 2 l. 217, Lemma 2.9 l. 407, Thm 2.10 l. 425, Thm 6.1 l. 1103, Sect. 6.2 l. 1138. Misprint phi(m) vs m in Sect. 2.2, Thm 2(c)(d), 2.10(3) is real, and its proof of (1)=>(3) is wrong; true period lcm{m}. n = 11 groups: 535 Phi_2Phi_6Phi_30 (ord 30), 64 Phi_2Phi_8Phi_18 (72), 24 Phi_2Phi_4Phi_8Phi_12 (24), 15 Phi_2^3Phi_6Phi_14 (42). Group 535 contains kA_11/rad^4 and rad^5 (only anchor with a known family); their FCY status is UNVERIFIED.
- arXiv e-print works from the cloud (curl -L https://arxiv.org/e-print/ID); 1310.1557 is a single gzip, not a tar; rate-limit once, retry after 20 s.
- Older: AI Thm 2.31, 2.32(b); Oppermann 1504.02617 v1 Thm 1.1 / Prop 7.6; T5 gap = repo child = H_0(kM) (E-217). Prop 7.6 and "leading words" notes in R079.

- R082 reply: prime rank 11 forces a trivial factor in tensor/T_n constructions, so only lines/Dynkin kA_11 from 0911.5137 at n = 11; windows are minimal a; misprint already recorded (lit summary l. 37-42).

## Did
- R001..R079 audits and summaries. R082: LaTeX theorem-number audit of both papers (no summary errors), scholar_cyclo.py (n = 11 cyclotomic table), named an object-level test (minimal complexes of nu^k P_i, rank growth as the NO signal).

## Next
- If maverick/toolsmith run the nu^k test: check its controls (A(10,3) fraction 1, D_10 a = 9 b = 8). Look for a source stating FCY (or not) for kA_n/rad^r beyond mn composite (Herschend-Iyama 0908.? "n-representation-finite", Lenzing-de la Peña): needed for the 535-group anchor.
- Lessons: grep summaries first; the summaries here were already right, so a literature round on an existing summary is a short audit; counters: remark/prop/cor share the lemma counter.
