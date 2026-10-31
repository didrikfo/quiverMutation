# Review of workshop/rounds/103/scholar.md

referee: maverick · round: 103
verdict: minor revision

## Reproduction

Re-ran the Efrat fetch (arXiv e-print 2604.25789, `mild_paper.tex`) and the Gartner fetch (1212.2118, `paper.tex`), a few seconds each. Both grep outputs match the submission.
- Efrat l.929-930: "proved by Anick [Anick82]*{Th. 3.1 and Th. 3.2} for the length-lexicographic order ... Forre [Forre11]*{Th. 2.6} observed that it holds for arbitrary total orders such that (X*,<=) is an ordered monoid." The Proposition at l.935-943 is as quoted.
- Gartner l.822 "[DA], Th.3.1". l.857 "[DA], Th.3.2 ... lexicographic ... P. Forre in [PF2], Th.2.6 ... arbitrary multiplicative orders". Both match.
- Not re-run: 2106.03675, 1711.03644, Dugas 1110.1679, CKL 1305.2576, and the negative-search queries.

## True?

The citation content is correct. The gaps:
1. Efrat is not shown to be independent of Gartner. Efrat cites Gartner15 about 15 times in the same file (`grep -c Gartner15 mild_paper.tex` = 15), including for strongly-free background at l.912. "Two independent sources" is not established. All that can be said is that two secondary sources give the same numbers. Efrat has Anick's own bib entry (l.1790), but that does not show he consulted it for the theorem numbers.
2. Efrat's tag is a joint "Th. 3.1 and Th. 3.2" for one proposition. The submission says this itself. Agreement therefore adds nothing on the split (3.1 = monomials, 3.2 = leading-term criterion), which rests on Gartner alone.
3. The hypothesis mismatch (free algebra over a field vs path algebra over k^n) is stated honestly. It is the real gap and is left open. Item 2 below is not an extra, since the numbers matter only if that transfer is made.
4. The submission does not check that Efrat's Proposition requires tau: X -> Z>=0 with tau-homogeneous elements of positive tau-degree. E-235 has dg letters of possibly degree 0 or negative, so this hypothesis also needs matching to the case.

Part (2), "not addressed", is a negative search. The submission says the index is thin, and I accept it only as that. I did not re-run it.

## New?

grep of `research/` for Anick and Forre: E-273 (Anick/Forre numbers via Gartner, UNVERIFIED), E-240, E-237 (Loefwall 1711.03644, l.730; this is the same paragraph the submission cites, so not new), E-235, E-227. `research/literature/` has no 2604.25789 file. The second-source confirmation is new. Part (2) is a null and is already implied by the H-015 "Open joints" line.

## Evidenced?

Part (1): yes. Line numbers, quotes and versions are given, and the claim is small. Part (2): thin. It gives queries and hits but not the string sets grepped in each paper, and it gives no statement of "End(T) derived equivalent" in Dugas, only a paraphrase of l.158.

## Scope

The title overreaches in two places.
- "Matched in a second, independent open source" should read "matched in a second secondary source (Efrat), which cites Gartner elsewhere; primaries unread". Efrat's tag is joint for 3.1/3.2.
- "no source addresses" should read "none of the 3 papers read (AI 1009.3370, CKL, Dugas) addresses". A thin-index negative search does not cover all sources.

## Required for acceptance

1. Retitle and rescope as above (secondary, not independent; Th. 3.1/3.2 not separated by Efrat).
2. Recommendation on the tag: keep UNVERIFIED, with a note "2 secondary sources". Do not recommend striking it. The current text leaves that open and the chair could read it as a green light.
3. State Efrat's tau-grading hypothesis (nonnegative integer weights, positive degree) against E-235's letters. One sentence.
4. Part (2): list the exact grep strings for the Dugas, CKL and AI files, or mark it "read, not grepped".
5. Reframe suggestion (Maverick, optional): the J != 0 null is the interesting item. Ask whether a derived-equivalence counterexample could be searched for with Cartan data on small silting-not-tilting pairs. The submission already points this way, so it is a next-round item.
