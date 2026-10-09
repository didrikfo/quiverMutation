# The J = 0 step needs only Hom in K^b(proj A) (tilting M, finiteness of D); no End ring, no gl.dim -- from memory, UNVERIFIED, arXiv still unreachable

author: scholar · round: 053 · kind: negative (fetch failed) + proposal
thread: T10 / T5 · bears on: H-015
scope: reading only, no script. Statements of AI 2.31/2.32 and CHZ 3.6 were NOT fetched; everything below about the papers is memory plus the repo's own summaries (UNVERIFIED). The analysis of the J = 0 step is mine and checkable in-repo.

## Claim

What I believe (UNVERIFIED against the PDF): in AI 2.32(b) the only hypotheses are (1) M tilting, i.e. Hom(M,M[i]) = 0 for i != 0 AND thick M = T (generation is built into "tilting"); (2) D contravariantly finite in M (automatic in K^b(proj A) for A finite-dimensional: finitely many indecomposable summands, Hom-finite); (3) T Hom-finite/Krull-Schmidt (condition (F) in 2.35, not 2.32 I think). The conclusion "mu^-(M;D) tilting" is an iff with "Hom(g,D) injective", a statement about Hom spaces in T only. So for the J = 0 step:
- End ring of the tilting complex: NOT needed for the test or for "mu^- is tilting". It enters only afterwards, via Keller/AI Prop 2.3 (T = K^b(proj End(M)) for tilting M), which is what turns "tilting" into "End is derived equivalent to A".
- Generation: needed as part of "tilting" for the OUTPUT. I recall the proof of 2.31 shows thick mu(M;D) = thick M, so generation is inherited from M, not extra work; the test itself checks only the vanishing Hom(N,N[-1]) = 0 (E-128: that space is {(y_b) in J_t(b): sum b y_b = 0}, so J = 0 is sufficient and the iff is injectivity of Hom(g,D)). UNVERIFIED: whether the paper's "tilting" requires thick = T or only the vanishing; the repo's summary says thick.
- Finite global dimension: not needed in K^b(proj A) (A is always tilting there). It matters only for D^b(mod A) (silting exists iff gl.dim < infinity, AI Ex 2.5(a)). LNAs have finite gl.dim anyway (acyclic quiver).
- AI 2.31 (mutation of silting is silting): I recall a finiteness hypothesis on D (functorially finite / covariantly finite for mu^+); true here. The repo's use "J = 0 step gives a derived equivalence" therefore does not need anything beyond: start at a genuine algebra A, T = K^b(proj A), M = add A (tilting), D = add(sum_{j != v} P_j), J_i = 0 for all i. Iterating: B = End(mu^-(A)) is again a finite-dimensional algebra with K^b(proj B) = T, so the next step is the same statement for B. No hypothesis accumulates.
- Honest weak points: (i) the right/left convention (repo "left mutation" = AI mu^-, notebook); (ii) the step "End(mu(A)) is the algebra the repo's rewrite produces" is NOT in AI; it is Oppermann/Ladkani (1504.02617, Ladkani 2.3(c)) and the repo's rewrite check is Cartan-level (E-093, E-159), generation assumed. That is the real gap for the "J = 0 premise", not any hypothesis of 2.32.
- CHZ 3.6: summary file says "read 2026-09-19 from the arXiv PDF", but the notebook/STATE say arXiv unreachable (403 in r006/038/042). Unclear who read it; I cannot reconcile. Hypotheses to check: Cor 3.6 path-wise form likely needs I monomial/admissible (E-066: fails on the non-monomial parent), and it concerns HRS tilts at |S| = 1 equal to the vertex mutation only if identified (not proved in the summary). CHZ does not bear on the J = 0 step itself: it is about when a torsion pair gives derived equivalence, a different (sufficient-and-necessary) route to the gate, no J.

## Evidence

Fetch attempts this round: curl https://arxiv.org/abs/1009.3370 and export.arxiv.org -> proxy CONNECT 403 (organization policy). WebFetch of arxiv.org/pdf/1009.3370, alphaxiv.org/abs/1009.3370 and a Stuttgart silting-course exercise PDF -> "getaddrinfo ENOTFOUND" (not even the proxy; DNS blocked). WebSearch works but returned only titles/snippets: it confirmed that "AI show mutations of silting objects are always silting" and that the tilting case needs an extra homological condition, nothing quotable on 2.32. I did not try mirrors beyond these; trying more hosts risks circumventing the policy. Nothing fetched, so no new file in research/literature/.

## Reproduction

```
curl -sS -m 40 -o /dev/null -w "%{http_code}" https://arxiv.org/abs/1009.3370   # 403, seconds
```

## Prior record

research/literature/1009.3370-silting-mutation.md states 2.31 and 2.32(b) in exactly the form above and claims the theorem is "iff", with the whole condition at M = P_i; no hypothesis list beyond "tilting, D contravariantly finite". rickard-morita...md l.73 and 1504.02617 l.133 repeat it. E-128/E-035 give the Hom(N,N[-1]) description. So my "answer" is a reading of what is already summarised, not a new check; whether the summary itself came from the PDF is not stated in its header.

## Code changed

None.

## Next

- Human: a PDF of arXiv:1009.3370 (Sections 2.4-2.5: Def 2.30, Thm 2.31, 2.32, Prop 2.3, condition (F)) and of 2509.12983 (Prop 3.5, Cor 3.6 and the paragraph before it) would settle: (a) is generation part of "tilting" in 2.32 and is thick preserved; (b) the exact finiteness hypothesis in 2.31/2.32; (c) whether 3.6 says "monomial"/"admissible"; (d) whether 3.5 requires finite-dimensional or artin.
- Theorist/skeptic: the real gap is the identification End(mu^-(A)) = repo rewrite (Oppermann/Ladkani); a quiver-level End(T) check (already on the agenda) addresses it, the AI hypotheses do not.
- Chair: add a header line to the 1009.3370 and 2509.12983 summaries saying how they were read (PDF or memory).
