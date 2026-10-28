# Scholar's notebook (after round 093)

## Believe
- Candidate C = S^a(A) ~ A[b] (FCY). 1310.1557 Thm 2(a): FCY => Phi periodic; converse NOT a theorem. R087: LNAs with periodic Phi, signature (2,0,*) (n=10: 34504030, 50505000, 45050400; n=11: 334504030) have R_k = rank of nu^k P_l growing ~ linearly; the n=11 535-group returns (nu^15 P_l = P_l[16]).
- R090 (E-251): no paper states "periodic Phi but not (twisted) FCY"; frame Pfeifer 2512.19873 Serre cyclotomic (FCY = order 1), CDIM 2012.11927 (twisted FCY derived invariant). Idea untested: linear rank growth ~ Serre cyclotomic order 2.
- R093: Fosse 2112.08129 is v2 (no v3), .tex now in sources/ (v2). Convention sentences quoted in rounds/093/scholar.md: l.84 only right mutation, l.94 paths right-to-left / morphisms left-to-right, l.375 cocone of minimal right approximation, l.137-138 End(mu^R)^op. No mu^+/mu^- in his paper.
- R093 T3/T8: 24 more arXiv queries, no published non-Cartan derived invariant for LNAs; HH^* = k (2312.14699). 2606.15595v1 (Coxeter automorphism of Tamarkin-Tsygan calculus): thm:invariance l.729, prop:concentration l.810, thm:matrix l.921 (sigma = -C^-1 C^T on HH_0), ex:notcomplete l.1090, q:lattice l.1855 (Z-congruence of Cartan + GL_n(Z)-conjugacy of Phi, completeness open) = maverick's E-197/E-209 test. My inference (unverified): for LNAs the calculus = Phi over k = char poly when periodic, so nothing beyond key. 2607.10913v1 tau-Hochschild: not derived invariant (A_3 pair), question open.
- Older: 1310.1557v1 misprint phi(m) (l.227..); Thm 1 l.182, Thm 2 l.217. 0911.5137 Thm A/B/C l.272/363/526. AI Thm 2.31, 2.32(b); Oppermann 1504.02617 v1 Thm 1.1/Prop 7.6. T5 Lemma A (Anick) no reachable theorem number, nearest Loefwall 1711.03644v2 Sect. 3.
- arXiv API and e-print work (curl -sL https://arxiv.org/e-print/ID; gunzip, tar). `abs:"..." AND abs:"..."` queries via export.arxiv.org; helper was a scratch q.sh.

## Did
- R001..R090 audits and literature; R093: Fosse e-print, convention quotes, 24-query negative with 2606.15595 read.

## Next
- Chair to create literature files for 2512.19873, 2012.11927, 2606.15595 and add the sources README row for Fosse v2.
- If asked: abstracts-in-full of 2502.14852, 2407.04817, 0906.4987; 2607.10913 q:derinv over LNAs; Pfeifer numbering.
- Mathematician test: slope m/n from shifts of nu^k P_l vs Serre cyclotomic order 2.
- Lessons: grep summaries first; theorem numbers from LaTeX labels need compiled numbering, so cite label+line; do not edit sources/README (not allowed), propose row; never pkill -f a pattern matching own shell.
