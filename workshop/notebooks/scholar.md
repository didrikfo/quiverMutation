# Scholar's notebook (after round 055)

## Believe
- The two papers are now read from LaTeX (sources/). arXiv is still unreachable; no longer needed for these two.
- AI Thm 2.31 (mutation of silting is silting) and Thm 2.32(b) (tilting iff Hom(g,D) injective, M tilting, D contravariantly finite) are Theorems, as cited. Thm 2.50(c) is the vertex form: Hom(eA/eA(1-e)A,(1-e)A)=0; for loopless v this is Hom(S_v, sum_{j!=v}P_j)=0, i.e. all J_j=0.
- J=0 is a whole-T condition: 2.32's proof splits Hom(T,T[!=0]) into four terms; Hom(N,N[<0])=0 follows from Hom(N,D[<0])=0. Generation and finiteness are automatic in K^b(proj A). "Silting not tilting iff some J_i!=0" = 2.31 + 2.32(b), not a sentence in the paper.
- AI mu^-_{P_v}(A) = N + D with N = (D'->P_v), P_v in degree 1; its [1]-shift is the K^{[-1,0]} object of the torsion pair (filt S_v, S_v^perp). At S={v}: AI 2.50(c) = CHZ Prop 3.5(2) (v not in supp soc P_j, j!=v). Remaining unproved: End(T) = repo's rewrite (Oppermann).
- CHZ Cor 3.6 is for kQ/I (not kA_n/I); the paper has no "monomial"; Ex 3.4 silently needs it. 5-vertex counterexample: a:1>2,b:2>4,c:1>3,d:3>4,e:4>5, I=<(ab+cd)e>, S={4}: Prop 3.5 fails, Cor 3.6 holds (scholar_chz_nonmonomial.py). Iff for monomial I (all LNAs).
- Open wording issue: E-364/E-419 say "left approximation" for the out-arrow cover P_tb -> P_v, which is a right approximation (AI Def 2.30).

## Did
- R055 response: referee minor revision answered (rerun ok; 'iff for monomial I' marked as own derivation; row 7 = wording slip right/left; Lemma L1 AI 2.50(c) = CHZ 3.5(2) at S={v} numbered; both summaries updated).
- R001..R053 as before; R055: full citation audit (table in scholar.md, 20 rows), edited both literature summaries (provenance lines, corrections marked r055), wrote the counterexample script.

## Next
- If the chair applies the table, check the edited research/ lines read as proposed.
- Read AI section 4 only if a question needs it (t-structures; not needed so far).
- Settle E-364/E-419 convention (left/right) with the skeptic.
- Test the Cartan-defect statement on E-131 rows; ask toolsmith for a socle-based (all-paths) gate on E-066 parent.
- Lessons: a theorem number from memory is often right but the hypotheses are what to read; LaTeX numbering is one counter per section (use awk on environments); a paper's "in other words" lines hide hypotheses; do not copy "unread/UNVERIFIED" flags forward once a source is in hand.
