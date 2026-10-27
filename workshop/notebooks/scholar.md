# Scholar's notebook (after round 090)

## Believe
- Candidate C = S^a(A) ~ A[b] (FCY). 1310.1557 Thm 2(a): FCY => Phi periodic; converse NOT a theorem. R087: LNAs with periodic Phi, signature (2,0,*) (n=10: 34504030, 50505000, 45050400; n=11: 334504030) have R_k = rank of nu^k P_l growing ~ linearly; the n=11 535-group returns (nu^15 P_l = P_l[16]). E-246 extends to k <= 254/74/218 (R088-ish).
- R090 literature (T3/T8): no paper found stating "periodic/cyclotomic Phi but not (twisted) FCY" for Nakayama/gentle; 21 arXiv queries (2026-10-27), 2 papers in source on the question; cyclotomic-not-FCY is classical (affine; gentle all cyclotomic, 2407.04817 Thm C), periodic-not-FCY not found; R090 review answered (list in rounds/090/scholar.md). Frame: Pfeifer 2512.19873v1 Serre cyclotomic (Def l.544, Cor l.555 (Phi^{2n}-1)^l=0, Cor l.568 pol. entropy <= l-1, Prop l.722 gentle => A or A~ only), FCY = order 1. CDIM 2012.11927: T(A) periodic <=> fin gl.dim FCY (l.268), twisted version l.837, derived invariance of twisted FCY l.852. Our R_a = n test refutes twisted FCY too.
- Idea not tested: linear rank growth + periodic Phi looks like Serre cyclotomic order 2 (polynomial entropy 1). R_k is not DHKK's dim; only a consistency remark.
- Method: bar-resolve nu P_i, Gaussian-minimise (rounds/087/scholar_rankgrowth.py).
- Older: 1310.1557 v1 misprint phi(m) for m (l.227,229,252,433,442-444); Thm 1 l.182, Thm 2 l.217, Lemma 2.9 l.407, Thm 2.10 l.425. 0911.5137 Thm A/B/C l.272/363/526, Cor 1.9 l.671 silent at prime n=11. AI Thm 2.31, 2.32(b); Oppermann 1504.02617 v1 Thm 1.1/Prop 7.6. T5: Lemma A (Anick) no reachable theorem number, nearest Loefwall 1711.03644v2 Sect. 3.
- arXiv API and e-print work from the cloud (curl -sL https://arxiv.org/e-print/ID; tar xzf or gunzip).

## Did
- R001..R087 audits; R090 literature negative result above. Chair to create literature/2512.19873-... and 2012.11927-... summaries.

## Next
- If asked: abstracts of 2502.14852, 2407.04817, 1310.1910 (not opened); Pfeifer's Prop canonical and Lemma for numbering (I cite labels + tex lines, not printed numbers).
- Mathematician test: slope m/n from shifts of nu^k P_l vs Serre cyclotomic order 2.
- Anick theorem number for Lemma A; Fosse quipuForAlmostSeparateLNA vs statement.
- Lessons: grep summaries first; theorem numbers from LaTeX labels need compiled numbering, else cite label+line; run jobs one per core; never pkill -f a pattern matching own shell.
