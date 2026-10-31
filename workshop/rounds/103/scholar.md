# Anick Th. 3.1 / 3.2 and Forre Th. 2.6: the same numbers appear in a second secondary source (Efrat 2604.25789v1, joint tag 3.1+3.2, primaries unread; Efrat cites Gartner elsewhere); for J != 0 steps none of the 3 papers read (AI 1009.3370, CKL, Dugas) addresses derived equivalence of End(T)

author: scholar · round: 103 · kind: literature
thread: T5/T10 · bears on: H-023, H-015, E-273, E-277
scope: arXiv LaTeX of 6 papers (list below), statements only; the second Anick source is SECONDARY, not independent (Efrat cites Gartner15 15 times, incl. l.906, l.912, but NOT at the Anick sentence l.929-930, which cites Anick82 and Forre11 directly); Efrat's tag is joint 3.1+3.2. Part (2) covers only AI, CKL, Dugas. Anick 1982 and Forre 2011 primaries NOT read (not on arXiv). The sources state the criterion for the free algebra F<X> (tensor algebra over a field), not for path algebras over a semisimple ring; that transfer is not checked here. The J != 0 part is a negative search (API index is thin: 0 hits is not absence).
power: n/a (literature)
backlog: n/a

## Response to referee

Re-ran this sitting: Efrat fetch (`grep -c Gartner15 mild_paper.tex` = 15; lines listed below); Dugas, CKL, AI e-prints re-fetched and grepped (strings below). Referee verdict minor revision; all 4 required items done, item 5 deferred.

1. Retitle/rescope: DONE. Title, Scope line and Claim now say "second secondary source", "primaries unread", "joint tag 3.1+3.2", and "none of the 3 papers read". "Independent" is withdrawn. Fact for the record: the 15 Gartner15 hits are at l.240, 246, 317, 395, 906 (Golod-Shafarevich inequality, "Lemma 2.5"), 912 (strongly free, "Prop. 2.9"), 1066, 1461, 1467, 1481, 1500, 1775, 1943 (bib), 2345, 2387; the Anick sentence l.929-930 and Proposition l.935-943 cite only Anick82 and Forre11. That is weak evidence of independence at that spot, not proof; I do not claim it.
2. Tag recommendation: KEEP UNVERIFIED, annotate "2 secondary sources (Gartner 1212.2118 splits 3.1/3.2; Efrat joint)". Do NOT strike. Striking needs Anick 1982 / Forre 2011 primaries. Claim text changed to say this (the earlier "un-tagged ... strike if chair accepts" wording is withdrawn).
3. Efrat's tau hypothesis (Prop. l.935-943): tau: X -> Z>=0 (weights 0 allowed on letters), gamma_j tau-homogeneous of POSITIVE tau-degree, order a total order making X* an ordered monoid, field F. Against E-235 (non-negatively graded parent, letters of degree >= 0, no loop at i; E-227 order omega = degree + star count): weight 0 on letters is allowed, but positivity of each relation's tau-degree is NOT checked here: a relation built only from weight-0 letters (omega = 0) would fall outside the Proposition. Open for the theorist; the Anick/Forre numbers are also moot until the path-algebra transfer (k^n base) is made.
4. Part (2) grep strings, all on the LaTeX main.tex fetched from arxiv.org/e-print (case-insensitive, `grep -ci`): counts for AI 1009.3370 / CKL 1305.2576 / Dugas 1110.1679: "silting" 203 / 26 / 3; "tilting" 112 / 76 / 28; "derived equiv" 5 / 30 / 25; "Anick" 0 / 0 / 0; "strongly free" 0 / 0 / 0. In AI all 5 "derived equiv" hits are in the commented-out or bibliography block (l.2597 onward, titles of cited papers), so AI states no derived equivalence of End(T) for non-tilting silting T in the body; I did not read the 203 silting hits line by line, so "not addressed" for AI rests on the "derived equiv" grep plus Thm 2.32 as already summarised. CKL: 30 "derived equiv" hits read by context only at l.223-249, 429-474, 638 (lifting of stable equivalences, Rickard's Nakayama permutation); not read exhaustively, mark "grepped, context-read at those lines". Dugas: the l.158 text (inside the proof of the Lemma mu^-(mu^+(A)) = A) reads "...G(A) is concentrated ..." with G induced by T_U; its "we know that G(A) is a tilting complex" claim is the line quoted in Claim (2) (grep "we know that" -> l.158). Dugas never states "End(T) derived equivalent" for a non-tilting T: tilting is assumed (Aihara's tilting mutation). Unchanged conclusion; string sets now stated. Search queries remain API-index-thin (not re-run).
5. Reframe (J != 0 derived-equivalence counterexample via Cartan data on silting-not-tilting pairs): next-round item for toolsmith/experimentalist, not done this sitting (needs code and a pair-generation design, outside a literature response). Already listed under Next.

## Claim

(1) Efrat, "Mild Pro-p Groups and Ordered Monoids", arXiv 2604.25789v1 (math.GR, 2026-04-28), states in its own words that Anick [J. Algebra 78 (1982) 120-140] proved the combinatorial-freeness criterion in "Th. 3.1 and Th. 3.2" for the length-lexicographic order, and that Forre [Crelle 658 (2011) 173-192] "Th. 2.6" extends it to arbitrary total orders making X* an ordered monoid. This agrees with the numbers E-273 took from Gartner 1212.2118 (Anick Th. 3.1 = monomials: strongly free iff combinatorially free; Th. 3.2 = leading-term criterion, lex; Forre Th. 2.6 = arbitrary multiplicative orders): SUPPORTED as a citation, but still secondary: both papers cite the primaries and do not reprove them. Efrat's tag is a joint "Th. 3.1 and Th. 3.2", so it does not separate the two; Gartner does. Verdict for the tags: the numbers move from "Gartner only" to "two secondary sources agree" (not shown independent: Efrat cites Gartner15 elsewhere). Recommendation: KEEP UNVERIFIED with the note "2 secondary sources"; do not strike (STEERING says "after matching them in the source": the source here is not Anick).

(2) Derived equivalence of the result of a J != 0 step: NOT ADDRESSED by any source read. What is stated: Aihara-Iyama 1009.3370 Thm 2.32 (already in `research/literature/1009.3370-silting-mutation.md`, matched in R055) gives tilting iff Hom(g,D) injective; a J != 0 step is silting-not-tilting, so Keller's equivalence (Prop. 2.3) is not available for End(T). No source I found says End(T) is, or is not, derived equivalent to A by some other functor. Chan-Koenig-Liu 1305.2576v1: not addressed (self-injective, simple-minded systems; derived equivalence only enters as lifting of stable equivalences). Dugas 1110.1679v2: not addressed (weakly symmetric algebras, mutation at a loopless vertex always uses a tilting complex; the proof of the unnumbered lemma after l.151 (l.158) says "we know that G(A) is a tilting complex"); tilting is automatic there, so the J != 0 phenomenon cannot occur in the hypotheses. This does not contradict H-015 either way.

## Evidence

Queries (export.arxiv.org API unless noted) and results:
- `all:"strongly free" AND all:Anick` -> 1711.03644v2 (Cyclic homology, gl.dim <= 2). Its l.730 says only: strongly free monomial sets = no mutual containment and no common beginning/ending; leading-monomial criterion "if there is an order of the variables". No theorem numbers: "not addressed" for numbers, consistent in content. Cites Anick as [anick] J. Algebra 78, 120-140.
- `all:"mild pro-p groups"` -> 1212.2118v1 (Gartner, = E-273 source, l.796-857 re-read: Thm "combfreetheo" l.818-822 "see [DA] Th.3.1"; "Anick's criterion" l.852 "DA Th.3.2 lexicographic; Forre PF2 Th.2.6 arbitrary multiplicative orders"), 2106.03675v2 (Mild pro-p groups and Koszulity conjectures: l.931 "Anick's criterion [anick]*{Thm. 3.2}", uses lex order, matches Th. 3.2 = leading-term criterion), 2604.25789v1 (Efrat: l.929-945 Proposition "Anick criterion", as in Claim (1)). math/0602189 and 2112.05027: no "Anick" in source.
- Forre 2011 itself: title/author/ti queries return nothing; not on arXiv index (or index gap).
- Dugas: `au:Dugas AND ti:tilting` -> 1110.1679v2. CKL: `all:Chan AND Koenig AND Liu AND simple-minded` -> 1305.2576v1; grep for silting/tilting shows only Thm D (strongly tilting connected, RFS) and lifting of stable equivalences.
- The three papers (1212.2118, 2106.03675, 2604.25789) agree on all three numbers where each gives them: Th. 3.1 (monomials), Th. 3.2 (lex leading terms), Forre Th. 2.6 (ordered monoid).

Hypothesis check against E-273 / E-235: all three state the criterion for homogeneous elements of F<X> with a grading tau, "strongly free" via the Hilbert-series equality H(z)(1 - sum z^tau(x) + sum z^deg rho_j) = 1 (Efrat l.905-911). The cancellation use in E-235/Lemma A is for dg path algebras over a semisimple base k^n: Anick's original 1982 paper works over a field with free algebra on X; the passage to a path algebra (X = arrows, with idempotent bookkeeping) is routine but NOT in these sources. That is the one hypothesis mismatch to watch (blind spot).

## Reproduction

```
mkdir -p $SCR/e && cd $SCR/e && curl -sL https://arxiv.org/e-print/2604.25789 -o s.src && tar xzf s.src && grep -n "Th\.\\\\ 3\.1\|Forr" mild_paper.tex
```
(<10 s; any scratch dir; sed -n 917,945p mild_paper.tex shows the statement.) Same for 2106.03675 (grep -n -i anick a.tex) and 1212.2118 (sed -n 796,860p paper.tex). Versions read: 2604.25789v1, 2106.03675v2 (id list returned v2), 1212.2118v1, 1711.03644v2, 1110.1679v2, 1305.2576v1. No repo code involved.

## Prior record

E-273 (Anick numbers via Gartner only, UNVERIFIED), STATE "UNVERIFIED citations" line, 1009.3370 summary (Thm 2.32, Thm 2.50(c); J_i = Hom(S_v, e_i A)). The Anick part is a second-source confirmation, not new mathematics. The J != 0 result is a null: nothing in AI/CKL/Dugas decides it; this does not move H-015.

## Proposed summary for the chair to promote (short)

`research/literature/2604.25789-efrat-mild-ordered-monoids.md`: "Anick criterion (Efrat Prop. at l.935 of v1): leading terms combinatorially free for an ordered-monoid order => strongly free; credits Anick Th. 3.1+3.2 (length-lex) and Forre Th. 2.6 (arbitrary). Definition of strongly free by Hilbert series equality; combinatorially free = no w_i a middle factor of w_j (i != j), no proper left factor of w_i a proper right factor of w_j. Not a source for path algebras."

## Code changed

None.

## Next

- Human/chair: if the bar for striking UNVERIFIED is the primary, Anick 1982 and Forre 2011 PDFs are needed (not on arXiv).
- theorist: write the path-algebra version of the criterion for Lemma A (arrows as generators over k^n; idempotent-graded combinatorial freeness) or cite a source for it; this, not the theorem numbers, is the real gap. Also loops of degree >= 1 (unitriangular, E-273).
- J != 0: only an own proof or an explicit pair can decide; sources do not. A pair A, End(T) with J != 0 and derived-inequivalent Cartan data is the cheapest test (toolsmith/experimentalist), and a derived-equivalent one would show the theorem is only sufficient.
