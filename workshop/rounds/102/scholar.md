# Among the sources read, mixed-digit derived-equivalent LNA pairs occur only in Fosse 2305.06642 (Cor. 4.3, 4.5, 4.6, almost separate); its Cor. 4.6 has a misprint in 2 of 8 algebras; "rank >= l_p+2" is Nikulin Thm 1.14.2 in 3 restatements, Cor. 1.13.3 UNVERIFIED

author: scholar · round: 102 · kind: literature
thread: T3/T8 (H-024) · bears on: H-024, E-265, E-267, E-269, E-274
scope: n = 9..12; conclusions hold for the sources read only; Nakayama/gentle/monomial literature named in the request; e-prints read as LaTeX (latest version at fetch, 2026-10-31). Not read: Happel-Seidel (not on arXiv), Bobinski-Geiss-Skowronski (not on arXiv; API search 0 hits), Ladkani math/0610685 Sec. 4 beyond the existing project summary, Nikulin 1979 itself. Computation only checks a printed statement against the repo.
power: n/a (literature). The Cartan-level checks in Evidence are positive controls (YES expected) and cannot say NO.
backlog: n/a


## Response to referee

Verdict minor revision; items 1-4 done.

1. Title and claim 1 narrowed to the sources read (listed in Scope; Happel-Seidel, BGS, Ladkani math/0610685 Sec. 4 and Nikulin 1979 are unread); "exist in print only" dropped. Cor. 4.6 needs all l_i >= 3 and almost separate relations, so mixed pairs outside that class are not addressed. The title and claim 3 no longer assert "not Cor. 1.13.3": what is established is only that the rank/length statement is Thm 1.14.2 in the three restatements read; the content of 1.13.3 is UNVERIFIED.
2. E-267 (research/EXPERIMENTS.md l. 62-66) is now cited. It calls 1.13.3 "Prop.", E-265 and this submission "Cor."; I cannot reconcile them, since no source read quotes 1.13.3, so the type is unverified and I write "1.13.3" bare. Cor. 1.13.4: re-read 2202.08030 l. 863 this round; it is cited (not quoted) for uniqueness of the even lattice of signature (2,10) with a given discriminant form. It is a third candidate for "unique in its genus", statement UNVERIFIED, and not the rank >= l_p+2 statement; nothing here depends on it.
3. Moved into the claim line: (a)+(b) of Thm 1.14.2 is weaker than a single "rank >= l_2+2" at p = 2 ((b) allows rank = l_2 with q = u_1+q' or v_1+q'), and is not the single bound E-265/E-267 state.
4. Misprint arithmetic, checkable by hand: vertices are 1..n and a relation of length L starting at s ends at s+L. Printed entry: start n-n_0-l_0+1, length n_0+l_0, end n+1 > n, off the line. With n_0+2: end n-l_0+3 <= n since l_0 >= 3. The scholar_mixedpairs.py count (1222 of 4888 fail) is the same arithmetic. The PDF check stays with the Skeptic.

## Claim

(1) Among the sources named in the request, the only published source of derived-equivalent LNA pairs with MIXED digits (>= 2 distinct nonzero relation lengths) is Fosse arXiv:2305.06642 (within the class l_i >= 3, almost separate; Cor. 4.3 = `corollary:2Rels`, Thm 4.5 = `thm:QuipuToAn`, Cor. 4.6 = `cor:EquivNakayamaAlgebras`). Happel-Seidel (via LMR 2112.15587 Prop. 4.1, Thm 4.4), LMR, Ueda 2302.02880 (Thm 1.1-1.5, Cor. 1.3, 1.6) and Ladkani 0911.5137 (Cor. 1.2) are all radical powers kA_n/rad^l, digit string l^k 0^m, so none gives a mixed pair. Dong-Lin-Ruan 2203.15735 (Prop. 4.5) likewise. BGS gentle one-cycle algebras are not linear Nakayama. (2) Cor. 4.6 as printed has a misprint: in the third and fourth right-hand entries the last relation length is printed `n_0 + l_0` and must be `n_0 + 2`. (3) The statement "indefinite even, (a) rank >= l_p + 2 for p != 2, (b) if rank = l_2 then q = u_1+q' or v_1+q' => class = genus" is Nikulin Theorem 1.14.2 in the three restatements read (2003.08914, 1704.00950, 1508.05251); (a)+(b) is WEAKER than a single "rank >= l_2+2" at p = 2. Whether Cor./Prop. 1.13.3 (E-267 says Prop.) or Cor. 1.13.4 says the same is UNVERIFIED; refuted if Nikulin's own text places it elsewhere. See E-267.

## Evidence

Mixed-digit pairs (one per n, from Cor. 4.6; all four members of each orbit have equal `coxeterKey`; digits = relation lengths, Kupisch series from `LinearNakayamaAlgebra.kupischSeries`):

| n | pair (digits) | Kupisch A | Kupisch B | theorem |
|---|---|---|---|---|
| 9 | 3030400 ~ 4003030 | 343544321 | 454343321 | Cor. 4.6 (1st ~ 2nd column, last-rel. swap + dual) |
| 10 | 30305000 ~ 50003030 | 3,4,3,6,5,5,4,3,2,1 | 5,6,5,4,3,4,3,3,2,1 | Cor. 4.6 |
| 11 | 303030400 ~ 400303030 | 3,4,3,4,3,5,4,4,3,2,1 | 4,5,4,3,4,3,4,3,3,2,1 | Cor. 4.6 |
| 11 | 304005000 ~ 003040030 | 3,5,4,7,6,5,5,4,3,2,1 | 5,4,3,5,4,5,4,3,3,2,1 | Cor. 4.6 |
| 12 | 3030305000 ~ 5000303030 | 3,4,3,4,3,6,5,5,4,3,2,1 | 5,6,5,4,3,4,3,4,3,3,2,1 | Cor. 4.6 |
| 12 | 3040050000 ~ 0040040030 | 3,5,4,7,6,5,6,5,4,3,2,1 | 6,5,4,6,5,4,5,4,3,3,2,1 | Cor. 4.6 |

Cor. 4.3 (a relation of length 2 can be added or deleted inside the hypotheses) gives further mixed pairs: any source with a free admissible start vertex. Theorem numbers: counted from `\newtheorem`/section order in the LaTeX (section 4 "Quipus and Nakayama algebras": 4.1 lemma, 4.2 `thm:AnToQuipu`, 4.3 `corollary:2Rels`, 4.4 lemma, 4.5 `thm:QuipuToAn`, 4.6 `cor:EquivNakayamaAlgebras`); no engine run, so the numbers rest on the counter order. Hypotheses: almost separate relations (n_{i+1} >= n_i + l_i - 1), all l_i >= 3 in Cor. 4.6, field not restricted in the statement (not checked in the proof).

Check of Cor. 4.6 against the repo (`scholar_mixedpairs.py`, all 611 LNAs at n = 9..12 with >= 2 relations, all l_i >= 3, almost separate, >= 2 distinct digits; the 8 algebras re-implemented from the printed formulas, not from `classPreservingOrbit`):
| check | result |
|---|---|
| members as PRINTED that are admissible LNAs | 1222 of 4888 fail, exactly the 3rd/4th right-hand entry for every source: start n-n_0-l_0+1 plus length n_0+l_0 ends at vertex n+1, off the line |
| same with last length `n_0+2` | 0 failures; all 8 are admissible and have the source's `coxeterKey` |
| set of the 8 vs repo `classPreservingOrbit` | equal for all 611 |
| Cartan-level, 48 pairs (orbit members and length-2 additions, sample of 24 sources) | Phi charpoly equal, S = C+C^T signature equal, det S equal, discriminant forms isomorphic (`disc_iso`): 48 of 48; Z[Phi]-module search not run (cypari2 not installed in .venv) |

So the published statement is SUPPORTED with the corrected entry, CONTRADICTED as printed (for 2 of 8 algebras), and the repo already implements the corrected form (`swapFirstRelation` etc.). This strengthens E-274: mixed-digit YES controls for the H-024 tests exist, but the 48 pairs check only the necessary conditions, not the E-269 search itself (maverick's job; not rerun).

Claim-by-claim table:
| claim | verdict | where |
|---|---|---|
| mixed-digit derived-equivalent LNA pairs in Happel-Seidel / LMR | not addressed (radical powers only, LMR Prop. 4.1, Thm 4.4 via existing summary; HS unread) | 2112.15587 |
| ... in Ueda | not addressed (Thm 1.1 needs L(s,t,u), not Nakayama; Cor. 1.3, 1.6 radical powers) | 2302.02880 (grep "Kupisch": 0 hits) |
| ... in Ladkani 0911.5137 | not addressed (Cor. 1.2 radical powers) | grep "Nakayama" only in Cor. 1.2 context |
| ... in BGS | not reachable | no e-print |
| Fosse Cor. 4.6, 8 algebras | supported after correcting `n_0+l_0` -> `n_0+2` | 2305.06642 sec. 4 |
| Nikulin class = genus, rank >= l_p+2 | supported by (primary restatement) as **Thm 1.14.2**: 2003.08914 `thm:Nikulin-uniqueness` (l. 1554): "(a) rank T >= l_p(q(T)) + 2 for all p != 2, (b) if rank T = l_2(q(T)) then q(T) = u_1 + q' or v_1 + q'" => one class in the genus. Same theorem cited as 1.14.2 by 1704.00950 (l. 1589, with the equality-of-length case treated separately) and 1508.05251 (l. 366-380, commented-out paraphrase: "rank >= l_2+2 or U(2^k)/V(2^k) summand") | |
| Nikulin Cor. 1.13.3 statement | not matched. Cited (without the rank condition) as "T_X is determined by signature and discriminant form" in math/0602015 l. 1732, l. 1795 (van Geemen-Sarti), and for "U + <-2k> is unique in its genus" in 2008.12077 l. 356. It is not quoted with a statement in any source read | |
| E-265's "rank >= l+2 and indefinite = Cor. 1.13.3" (via 1404.5460 l. 336 in round 101) | CONTRADICTED by the three sources above as the number of that statement; 1404.5460 not re-read this round | |

Consequence for H-024: the hypothesis is "rank >= l_p + 2 for all odd p, plus the p = 2 condition (b)"; E-265/E-267 state a single bound "rank >= l(A_q)+2" which is stronger than (a)+(b) at p = 2, so using it is safe but not what Nikulin states. 1.13.3 is not needed for any NO (round 101 conclusion stands).

E-267 (research/EXPERIMENTS.md) already records that 1.13.3 is not found in 28 e-prints and that 1.14.2 and Cor. 1.13.4 are what others cite.

Searched: arXiv API (ti:"discrete derived categories": 4 papers by BPP, none BGS; au:Bobinski/Geiss/Skowronski: 0 hits; ti:Nakayama AND ti:derived equivalen: 0 hits; API index thin, 0 is not absence); web search for Nikulin 1.13.3 (statement not returned). E-prints read: 2305.06642, 2302.02880, 2112.15587, 0911.5137 (grep only), 2003.08914, 2008.12077, 1704.00950, 1508.05251, math/0602015.

## Reproduction

```
curl -sL https://arxiv.org/e-print/<id> -o src; tar xzf src || gunzip -c src > main.tex    # own scratchpad dirs, ~1 min
timeout 10m .venv/bin/python workshop/rounds/102/scholar_mixedpairs.py 9 12     # ~4 min; prints the counts above
timeout 10m .venv/bin/python workshop/rounds/102/scholar_table.py               # seconds; the pair table
```

## Prior record

`research/literature/2305.06642-quipu-quivers-nakayama.md` records Cor. "cor:EquivNakayamaAlgebras" and Cor. 2Rels; the repo implements them (`nakayama.py` swapFirstRelation/swapLastRelation/relationDual, `tests/test_quipu_symmetry.py`) and lists "radical powers only" as the caveat of every other summary. The printed misprint is not recorded (grep "n_0+l_0", "misprint" in research/: none). E-274 already says no mixed-digit control existed; this supplies the source and 6 explicit pairs, but they are the repo's own classification theorem, so they are positive controls of the Cartan tests, not independent of the project's authorship (Fosse is the project author).

## Code changed

None. Scripts: `workshop/rounds/102/scholar_mixedpairs.py`, `scholar_table.py`.

## Next

- Chair: promote a note on 2305.06642 (misprint `n_0+l_0` -> `n_0+2` in Cor. 4.6 right column entries 3 and 4; numbering 4.3/4.5/4.6) and change the Nikulin citation to Thm 1.14.2 with hypotheses (a),(b) (primary-text restatement 2003.08914), keeping 1.13.3 UNVERIFIED.
- Maverick: run the E-269 search on the 6 pairs above (mixed digits, includes n = 11, 12).
- Skeptic: confirm the misprint against the PDF of v-latest (I read LaTeX only).
- Open: genuinely independent mixed pairs would come from Ladkani math/0610685 Cor. 4.6 / Prop. 4.7 (A_Y for closed subsets Y of a poset) when A_Y is a line with relations of different lengths; not searched.
