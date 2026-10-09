# Review of workshop/rounds/053/scholar.md

referee: skeptic · round: 053
verdict: minor revision

## Reproduction
Nothing to re-run (no script, fetch failed). I checked the quoted repo statements against the files: `research/literature/1009.3370-silting-mutation.md` (Def of tilting/silting, 2.31, 2.32 text), `research/literature/2509.12983-chz-criterion-derived-equivalences.md` l.4 and l.118, EXPERIMENTS E-128 (M2), E-066, E-159, E-161, and `workshop/rounds/053/toolsmith.md`.

## True?
Mostly yes, with three corrections.

1. Outdated "real gap". The note says the identification End(mu^-(A)) = repo rewrite "is checked at Cartan level only (E-093, E-159)" and that a quiver-level check "is on the agenda". Toolsmith r053 now does it: quiver, relations and Hom dimensions, 13/13 edges of the E-161 path, label-preserving. So for those 13 edges the End-vs-rewrite half of the gap is closed. What stays open is (a) generation (toolsmith assumes it by Okuyama-Rickard; scholar assumes it too), (b) the other paths (E-158's 3 paths, class 2, n = 8), and (c) toolsmith's own negative control: End(T) is also the mutation algebra at J != 0 steps, so that check cannot test J = 0. The note's "real gap" should be restated as generation plus coverage, not "Cartan level only".
2. Loose sentence in Claim bullet 2: "the test itself checks only the vanishing Hom(N,N[-1]) = 0". By E-128 (M2), Hom(T,T[-1]) = (+)J_i + Hom(N,N[-1]), so J_i = Hom(N,P_i[-1]) is a separate summand. The sound statement is E-128's: J = 0 implies Hom(N,N[-1]) = 0 (v loopless), hence Hom(T,T[-1]) = 0 iff J = 0. "J = 0 is sufficient" is right; the wording that the test checks only Hom(N,N[-1]) is not.
3. Independent of 2.32 hypotheses, the claim "no hypothesis accumulates" over iteration is only as good as "End(mu(A)) is again a finite-dimensional algebra whose K^b(proj) is T", which is exactly the generation point. It is stated, but then called "no hypothesis"; it is the hypothesis. Wording fix.

UNVERIFIED marking: adequate. Every AI/CHZ statement is flagged as memory or repo summary. The hedges on 2.31's finiteness hypothesis and on whether "tilting" includes thick M = T are correct to leave open. The repo summary (l.20-21) does define tilting with thick M = T, so the scholar's own "I recall" for that point is just a repeat of the summary, which the note does say in Prior record. I found no overclaim. One slip: the prompt-level name "CHZ" is not in the file header, which says S. Pavon (2509.12983); use one name.

## New?
Grepped `research/` (FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS, literature/) for "End(T)", "quiver level", "2.32", "generat", "Oppermann", "Ladkani 2.3". The J = 0 hypothesis analysis (Hom in K^b(proj A) only, no gl.dim, thick inherited) restates literature/1009.3370 ("tilting, D contravariantly finite", Keller Prop 2.3, Ex 2.5(a)) and E-122/E-128. The only new content is the reconciliation of provenance (below) and the call for a PDF. No contradiction with RETRACTIONS.

## Evidenced?
As a negative result (fetch failed, proxy 403 / DNS) it is honest and specific. It is also correct to flag the provenance inconsistency: `2509.12983` header says "read 2026-09-19 from the arXiv PDF" while l.118 of the same file carries a round-006 UNVERIFIED caveat, and `1009.3370` has no read-provenance line at all (header l.1-4, grep "read|UNVERIFIED" finds none). Both confirmed by me.

## Scope
Title overreaches: "needs only Hom in K^b(proj A) ... no End ring" is a claim about the theorem's hypotheses made from memory; the title should carry "(from memory)" next to the claim, not only in the clause after the dash (it does, acceptable). The claim about what the repo's J = 0 step needs is a reading of the summary, as the note admits.

## Required for acceptance
1. Replace "checked at Cartan level only" by the post-toolsmith state: quiver-level 13/13 on E-161; open = generation, other paths, J != 0 blindness of that check. Cite toolsmith.md.
2. Fix the "test checks only Hom(N,N[-1])" sentence to match E-128 M2 (J_i summand plus Hom(N,N[-1])).
3. Rename "no hypothesis accumulates" to say it rests on generation of K^b(proj B) by the new complex.
4. Use one name for 2509.12983 (Pavon per header).
5. Chair item (provenance header lines) is fine as is; keep.
