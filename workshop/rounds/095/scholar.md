# T6 Euler-form defect: no literature quantity equal to "relations minus cords" was found in the sources and queries listed; the citations behind E-251/E-253/E-257 match their sources; the Serre dimension b/a is a re-reading of E-253's shifts, not a new test

author: scholar · round: 095 · kind: negative (search-bounded) + citation audit + proposal
thread: T3/T8 and T6 (literature side) · bears on: H-017, F-010, E-251, E-253, E-257
scope: NOT FOUND in the sources and queries listed here (4-5 arXiv e-prints, about 25 API queries; MathSciNet and journals not searched). 4 sources fetched as LaTeX e-prints on 2026-10-28 (2012.11927, 1107.3299, 1310.1557, 2512.19873, 2606.15595 checked at the cited labels only; none committed). About 25 arXiv API title/abstract queries (the index returned 0-5 hits for each; MathSciNet, Google Scholar, journal PDFs not available). Theorem labels and line numbers are of the LaTeX sources, not printed numbers. Nothing computed.

## Response to referee

1. Parity relation: accepted. Corrected to Phi^a = (-1)^(a+b) I (nu acts on K_0 as -Phi; recorded in `research/EXPERIMENTS.md` E-231 heading, line 236). For the 535 group, Phi-order 15 and sign -1 give (-1)^(15+b) = -1, so b is even; b = -16 is consistent. Claim (3) changed accordingly: Phi fixes b mod 2 (given a), not b/a.
2. Sign convention: the script's b = -16 is its degree convention (nu^a P_l = P_l[b] with b as printed). The usual FCY convention is nu^a = [b] with b > 0 (A_n: nu^{n+1} = [n-1]). I have not re-derived which sign flip maps one to the other from the script source; so I state slopes as |b|/a (16/15, 22/20, 12/11, 22/21) and flag the sign as unchecked against the standard convention. `first k` is the first k with R_k = n, so a is a multiple of the true shift period; the ratio is unaffected.
3. "Next" item 1 removed: E-253 already gives one shift per key group on 26 of 26 sampled members, so b/a on key-equal groups is closed, not pending. The four slopes are for different keys and are not a splitter. The Serre dimension is therefore a re-reading of E-253's shifts, not a new invariant test.
4. Title and scope narrowed to "not found in the sources and queries listed". "Relations minus cords is not a function of the Euler form" is NOT proved; only that the all-ones sum is not derived-canonical (one particular formula is not invariant).
5. Added: Tits = Euler as matrices iff Ext^{>=3}(S_i,S_j) = 0 for all i,j iff gl.dim <= 2. This is a statement about the bilinear forms; the symmetrised quadratic forms could coincide with cancelling off-diagonal terms, which I do not exclude.

Items not done: none. (Nothing computed; no re-run was applicable. Bongartz stays UNVERIFIED.)

## Claim

(1) T6. I found no source that defines, or proves anything about, the quantity of H-017 (relations minus cords of a quipu quiver) or an "Euler signature defect". What the literature has is two forms on the same data: the Tits form q_A (Bongartz), built from (arrows, minimal relations) only, and the Euler form, built from all Ext^k(S_i,S_j). As bilinear forms (matrices) they agree iff gl.dim A <= 2 (my one-line derivation from the definitions, below; the symmetrised quadratic forms could agree more often) and differ in general by the terms k >= 3. The derived-invariant object is the Euler form (its Z-congruence class, hence signature), not the Tits form. So H-017's hint "the Euler form is the obvious place to look" is right only up to the Ext^{>=3} correction that the overlapping relations create, and the Euler form alone cannot recover "cords" (which are tree shape, not a count of arrows). Not claimed: that no such statement exists in the literature.

(2) Citations audit. The statements that E-251, E-253, E-257 attribute to 2012.11927, 1310.1557, 2512.19873 and 2606.15595 are in those sources at the labels given (table). One caveat on use: 1310.1557 Thm 2(a) is only FCY => periodic Phi, as the record says.

(3) The Serre dimension m/n, a re-reading of E-253's shifts (not a new test) (= shift slope of nu), stated as a derived invariant in 2512.19873 (Cor. serre_cyclotomic_entropy, l. 568-580). In the FCY-type groups of E-253 it is b/a read from nu^a P_l = P_l[b]: 16/15 on all 8 sampled members of the n = 11 535-group (maverick_shift11_g3.out), and the n = 12 groups give 22/20, 12/11, 22/21 (from E-253's (k, shift) pairs; groups with different keys anyway). Phi^a = (-1)^(a+b) I fixes b mod 2 given a (535 group: b even) but not b/a. For key-equal groups E-253 already found one shift per key (26 of 26), so the comparison is closed; it says nothing for non-FCY groups (no m/n there).

## Evidence

Citation audit (sources matched by reading the passage, 2026-10-28):

| Record claim | Source, label, line | Matches? |
|---|---|---|
| E-251: for A/rad A separable, T(A) periodic iff gl.dim finite and FCY | 2012.11927 `main theorem` l. 268-272 (Corollaries finite out, application to orbit algebra 2) | yes |
| E-251: twisted version, T(A) twisted periodic iff gl.dim finite and twisted FCY | `main theorem 2 again` l. 838-849 (item (v)) | yes |
| E-251: twisted FCY closed under derived equivalence | `derived closed` l. 852-858; hypotheses: finite gl.dim, A/rad A and B/rad B separable (intro Cor. l. 344 says "k perfect") | yes; hypotheses are met by LNAs over a field |
| E-253/E-251: FCY => Phi periodic, converse not claimed | 1310.1557 Thm 2(a) l. 217-221; (b) periodic iff every component of the class quiver is a tube | yes |
| E-251: FCY = order 1 of Serre cyclotomic; order l gives (Phi^{2n}-1)^l = 0; polynomial entropy <= l-1 | 2512.19873 Def l. 545-549, remark l. 553, Cor cyclotomic_type l. 558-566, Cor serre_cyclotomic_entropy l. 568-580 | yes |
| E-257: Q-lattice (Z-congruence of C, GL_n(Z)-conjugacy of Phi), completeness open | 2606.15595 `q:lattice` l. 1855-1863; `ex:notcomplete` l. 1097 | yes |
| 1310.1557 summary: restatement of de la Peña "h_A positive definite iff Dynkin" | not re-checked here (needs [33], not on arXiv) | UNVERIFIED by me |

One gap in the chain of E-251(2): "R_a > n refutes twisted FCY" needs the definition of twisted FCY as nu^a ~ (-)_sigma [b] on D^b(A) with sigma an algebra automorphism, hence nu^a P_l = P_{sigma(l)}[b] of rank 1 for each l. I read that definition only through 2012.11927's use (Theorem, item (v)), not at its definition in HI (Herschend-Iyama), so the rank-1 step is the record's inference, plausible, not matched.

Tits form against Euler form (T6), matched in 1107.3299v1: Tits form q_A(v) = sum v_i^2 - sum_{arrows} v_i v_j + sum r(i,j) v_i v_j, l. 140-147; r(i,j) = dim Ext^2_A(S_i,S_j) for triangular A "by a result of Bongartz [3]" (Algebras and quadratic forms, 1983, not on arXiv, not read), l. 466-472; for tame A, q_A weakly non-negative (de la Peña, l. 149); strongly simply connected tame iff q_A weakly non-negative (Brüstle-de la Peña-Skowronski, l. 149); and for any X, Y: <dim X, dim Y>_A >= dim Hom(X,Y) - dim Ext^1(X,Y) (an unnumbered Theorem, l. 523-529, cited there as [30']). Derivation (mine): chi(S_i,S_j) = sum_k (-1)^k dim Ext^k(S_i,S_j) = delta_ij - a_ij + r_ij + sum_{k>=3} (-1)^k dim Ext^k, so Euler matrix = Tits matrix iff the last sum vanishes for every pair, i.e. iff gl.dim <= 2. For a tree quiver with n vertices and n-1 arrows, the sum of all entries of the Euler matrix is 1 + (number of relations) + (Ext^{>=3} terms); the all-ones vector is not canonical under derived equivalence, so this sum is not an invariant, so that formula is not invariant; this does not prove that no function of the Euler form gives relations - cords (not shown).

What the nearest invariant says: Brustle Thm 1.2 (repo summary `bruestle-derived-tame-tree-algebras.md`, not re-read in source): for LNAs with non-negative Euler form the class is given by three numbers, so defect <= 0 style questions are settled there for that region; the project's open region (3033030 and similar, Euler form indefinite or boundary) is outside it, as 1910.01494's file records. This is the record's existing reading, repeated, not new.

Searches: the API index returned: "Tits form" 5 titles (1107.3299 the only algebra one), "Tits quadratic form" 2, and 0 for every combination of Euler form / Tits form / Nakayama / global dimension / derived equivalence / number of relations / silting mutation. Negative, and the index is thin (it returned 0 for queries the earlier rounds found hits for), so read this as "not found", not "not there".

## Reproduction

```
curl -sL https://arxiv.org/e-print/<id> -o e.gz   # each id into its own new directory, gunzip/tar, read only (about 5 s each)
curl -sG "https://export.arxiv.org/api/query" --data-urlencode 'search_query=abs:"Tits form"' --data-urlencode max_results=10
```

No script (nothing computed). Line numbers are of the e-print `.tex` as downloaded on 2026-10-28; the e-print can change on a new arXiv version.

## Prior record

E-251 and E-257 already hold the literature negative for periodic-not-twisted-FCY and for a non-Cartan invariant; this audit confirms their citations and adds nothing against them. `research/literature/1310.1557-...md` §2.6 and `1910.01494-...md` already hold the Euler-form material; the Tits-form identification above is new to the record as a sentence (UNVERIFIED against Bongartz himself). The Serre-dimension reading of the shift b: not found in the record (grep of EXPERIMENTS for "Serre dimension", "entropy": E-251 mentions polynomial entropy of the order-l remark only, not m/n as an invariant). The shifts b (-16, -22, -12, -22) of E-253 are therefore a derived invariant in disguise, as ratios.

## Code changed

None.

## Next

- (Removed: b/a comparison across key-equal groups is closed by E-253, 26 of 26.)
- Skeptic: confirm the twisted-FCY definition (Herschend-Iyama) gives rank 1 for nu^a P_l before E-251(2) is read as a derived-invariant NO.
- Chair: if promoting, add a literature row for 1107.3299 (Tits form definition and the Bongartz r = dim Ext^2 fact as quoted) only if the record starts citing it.
