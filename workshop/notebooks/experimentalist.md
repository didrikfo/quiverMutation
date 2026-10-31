# Experimentalist notebook (rewritten round 103)

## What I now believe
- 103 literature (T7, S-1): arXiv search finds NO source for the Catalan factorisation, the K-threshold law K0 = 3, 4, 5, or |orbit B_n| = n 2^(n-5): "not addressed" (arXiv title/abstract search, 5 sources opened; OEIS unreachable, HTTPS 403). Adjacent: LNAs <-> Dyck paths (2602.04571 Prop. 2.1, number by my count), Catalan counts (1910.01494).
- Gate part of the single-mutation rule is independently published for gentle LNAs: Deng-Geng-Liu 2609.39536 Thm 4.2(3) / 7.1 (numbers from LaTeX counters, PDF unchecked). Repo gate = paper criterion on all 3 168 gentle cells n = 3..9 (paper = repo LEFT gate; mirror = RIGHT gate; crossed controls fail). Nothing published for r_v >= 3.
- Fosse 2305.06642 Cor. 4.6 (r = 0, n_0 = n-3, l = 3) gives A_n ~ B_n for all n >= 5, so E-218/E-247 "equivalent known only at n = 9, 10" is too weak; orbit(B_n) is a proper subset of the derived class (A_n missing). Orbit at n = 6..8 hits only table class D_n (2 of its 3 shapes). Chair to correct wording.
- Older: 099 pair (8:4)(10:4) has no lowering control (T7 uninformative there); 097 n = 15 key scan (S-1 item 3 key level); row convention: n-vertex LNA is an n - 2 entry row, entry i = relation of that many arrows starting at vertex i+1.

## What I tried
- 103: `rounds/103/experimentalist_arxivq.py` (14 queries), `_arxivabs.py` (abstracts and id lookups), `_gentle_gate.py`, `_fosse_check.py`. arXiv API rate-limits (429): space calls 6 s, retry.
- Earlier rounds: 099 `_pairs/_single/_single4`; 097 `scan15`; 089 `orb13`; see those rounds.

## What I would do next
1. Scholar: PDF-verify the three "my count" theorem numbers; OEIS lookup of 12, 28, 64, 144, 320, 704 (blocked here).
2. Extend the gate check to rows with a length-3 relation using Ladkani 1001.4765 Prop 2.3(c) (library) as the source.
3. T7: lowering control for start-gap-2 pairs (stepped chain 8, 10, 12) before more two-bystander cells; S-1: is A_n reached from B_n by any tilting mutation at n = 6 (outside the move set)?
