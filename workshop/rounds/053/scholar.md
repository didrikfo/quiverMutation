# The J = 0 step needs only Hom in K^b(proj A) (tilting M, finiteness of D); no End ring, no gl.dim -- but the iteration rests on generation of K^b(proj B) by the new complex; from memory, UNVERIFIED, arXiv still unreachable

author: scholar · round: 053 · kind: negative (fetch failed) + proposal (revised after referee)
thread: T10 / T5 · bears on: H-015
scope: reading only, no script. Statements of AI 2.31/2.32 and Pavon 2509.12983 (Prop 3.5, Cor 3.6) were NOT fetched; everything below about the papers is memory plus the repo's own summaries (UNVERIFIED). The analysis of the J = 0 step is mine and checkable in-repo. The End-vs-rewrite identification is taken from toolsmith r053 (label-preserving, 13/13 edges of the E-161 path only; class 2, the 3 E-158 paths and n = 8 not covered).

## Response to referee

Verdict was minor revision; all four required items are done, item 5 kept.

1. **Outdated "real gap" (required 1): accepted.** The old text said the End(mu^-(A)) = repo-rewrite identification is "checked at Cartan level only (E-093, E-159)" with a quiver-level check on the agenda. Replaced (Claim, last bullet, and Next) by the state after `workshop/rounds/053/toolsmith.md`: quiver, relations and Hom dimensions, 13 of 13 edges of the E-161 path, label-preserving (no vertex permutation tried), over the algebraic closure. Open: (a) generation of K^b(proj B) by the new complex (toolsmith assumes it by Okuyama-Rickard, as do I); (b) coverage: class 2, the 3 E-158 paths, n = 8, the 25 E-155 paths; (c) toolsmith's own negative control: End(T) is also the mutation algebra at J != 0 steps (8 of the 16 failing steps decided, iso in all 8; 8 undecided because of parallel arrows), so this check cannot test the J = 0 premise. The gap is restated as generation plus coverage, not "Cartan level only".
2. **Sentence on what the test checks (required 2): accepted.** I wrote "the test itself checks only the vanishing Hom(N,N[-1]) = 0". Wrong as worded. By E-128 (M2), Hom(T,T[-1]) = (+)_i J_i + Hom(N,N[-1]) with J_i = Hom(N,P_i[-1]) a separate summand. Sound statement, now in the Claim: J = 0 implies Hom(N,N[-1]) = 0 (v loopless), hence Hom(T,T[-1]) = 0 iff J = 0; "J = 0 is sufficient" stands.
3. **"No hypothesis accumulates" (required 3): accepted.** Reworded: the iteration is as good as "End(mu^-(A)) is a finite-dimensional algebra B with K^b(proj B) = T", i.e. generation of T by the new complex. That is a hypothesis (inherited from thick mu(M;D) = thick M in AI 2.31, which I recall but did not read), not the absence of one.
4. **One name for 2509.12983 (required 4): done.** "CHZ" replaced by "Pavon (2509.12983)" per the file header. (The same replacement in my own Prior record; the repo file name still contains "chz", I did not touch it.)
5. Chair item (provenance header lines for 1009.3370 and 2509.12983): kept as is.

The referee's other remarks (UNVERIFIED marking adequate; no new literature file since nothing was fetched) need no change.

## Claim

What I believe (UNVERIFIED against the PDF): in AI 2.32(b) the only hypotheses are (1) M tilting, i.e. Hom(M,M[i]) = 0 for i != 0 AND thick M = T (generation is built into "tilting"); (2) D contravariantly finite in M (automatic in K^b(proj A) for A finite-dimensional: finitely many indecomposable summands, Hom-finite); (3) T Hom-finite/Krull-Schmidt (condition (F) in 2.35, not 2.32 I think). The conclusion "mu^-(M;D) tilting" is an iff with "Hom(g,D) injective", a statement about Hom spaces in T only. So for the J = 0 step:
- End ring of the tilting complex: NOT needed for the test or for "mu^- is tilting". It enters only afterwards, via Keller/AI Prop 2.3 (T = K^b(proj End(M)) for tilting M), which is what turns "tilting" into "End is derived equivalent to A".
- Generation: needed as part of "tilting" for the OUTPUT. I recall the proof of 2.31 shows thick mu(M;D) = thick M, so generation is inherited from M, not extra work. What the test decides (E-128 M2): Hom(T,T[-1]) = (+)_i J_i + Hom(N,N[-1]) with J_i = Hom(N,P_i[-1]); J = 0 implies Hom(N,N[-1]) = 0 (v loopless), so Hom(T,T[-1]) = 0 iff J = 0, and the iff of the theorem is injectivity of Hom(g,D). The test is therefore not "only Hom(N,N[-1])"; J is a separate summand. UNVERIFIED: whether the paper's "tilting" requires thick = T or only the vanishing; the repo's summary (literature/1009.3370, l.20-21) says thick M = T.
- Finite global dimension: not needed in K^b(proj A) (A is always tilting there). It matters only for D^b(mod A) (silting exists iff gl.dim < infinity, AI Ex 2.5(a)). LNAs have finite gl.dim anyway (acyclic quiver).
- AI 2.31 (mutation of silting is silting): I recall a finiteness hypothesis on D (functorially finite / covariantly finite for mu^+); true here. The repo's use "J = 0 step gives a derived equivalence" needs: start at a genuine algebra A, T = K^b(proj A), M = add A (tilting), D = add(sum_{j != v} P_j), J_i = 0 for all i. Iterating: B = End(mu^-(A)) is a finite-dimensional algebra, and the next step is the same statement for B provided K^b(proj B) = T, which holds iff mu^-(A) generates T. That generation is the one hypothesis carried along the iteration (inherited from M by thick mu(M;D) = thick M, if my recollection of 2.31 is right); it is not free.
- Open points, restated after toolsmith r053: (i) the right/left convention (repo "left mutation" = AI mu^-, notebook); (ii) the step "End(mu(A)) is the algebra the repo's rewrite produces" is NOT in AI; it is Oppermann/Ladkani (1504.02617, Ladkani 2.3(c)). It is now checked at quiver level: toolsmith.md reports End_K(T) (quiver, relations, Hom dimensions) equal to the rewrite on 13 of 13 edges of the E-161 path, label-preserving. This closes the End-vs-rewrite half for those 13 edges only. What stays open for the "J = 0 premise" is generation (assumed by toolsmith by Okuyama-Rickard, and by me), coverage of the other paths (class 2, E-158's 3 paths, n = 8), and the fact that End(T) is the mutation algebra also at the J != 0 steps (8 of 16 decided), so that comparison cannot discriminate J = 0; the discriminating object is Hom(T,T[-1]) (E-126/E-159), not End(T).
- Pavon 2509.12983 (Cor 3.6): summary file says "read 2026-09-19 from the arXiv PDF", but the notebook/STATE say arXiv unreachable (403 in r006/038/042). Unclear who read it; I cannot reconcile. Hypotheses to check: Cor 3.6 path-wise form likely needs I monomial/admissible (E-066: fails on the non-monomial parent), and it concerns HRS tilts at |S| = 1 equal to the vertex mutation only if identified (not proved in the summary). Pavon does not bear on the J = 0 step itself: it is about when a torsion pair gives derived equivalence, a different (sufficient-and-necessary) route to the gate, no J.

## Evidence

Fetch attempts this round: curl https://arxiv.org/abs/1009.3370 and export.arxiv.org -> proxy CONNECT 403 (organization policy). WebFetch of arxiv.org/pdf/1009.3370, alphaxiv.org/abs/1009.3370 and a Stuttgart silting-course exercise PDF -> "getaddrinfo ENOTFOUND" (not even the proxy; DNS blocked). WebSearch works but returned only titles/snippets: it confirmed that "AI show mutations of silting objects are always silting" and that the tilting case needs an extra homological condition, nothing quotable on 2.32. I did not try mirrors beyond these; trying more hosts risks circumventing the policy. Nothing fetched, so no new file in research/literature/.

## Reproduction

```
curl -sS -m 40 -o /dev/null -w "%{http_code}" https://arxiv.org/abs/1009.3370   # 403, seconds
```

## Prior record

research/literature/1009.3370-silting-mutation.md states 2.31 and 2.32(b) in exactly the form above and claims the theorem is "iff", with the whole condition at M = P_i; no hypothesis list beyond "tilting, D contravariantly finite". rickard-morita...md l.73 and 1504.02617 l.133 repeat it. E-128/E-035 give the Hom(N,N[-1]) description. So my "answer" is a reading of what is already summarised, not a new check; whether the summary itself came from the PDF is not stated in its header. The quiver-level End(T) result is toolsmith r053 (13/13, label-preserving); before it only E-093/E-159 (Cartan / Hom dimensions).

## Code changed

None.

## Next

- Human: a PDF of arXiv:1009.3370 (Sections 2.4-2.5: Def 2.30, Thm 2.31, 2.32, Prop 2.3, condition (F)) and of 2509.12983 (Prop 3.5, Cor 3.6 and the paragraph before it) would settle: (a) is generation part of "tilting" in 2.32 and is thick preserved; (b) the exact finiteness hypothesis in 2.31/2.32; (c) whether 3.6 says "monomial"/"admissible"; (d) whether 3.5 requires finite-dimensional or artin.
- Theorist/skeptic: the open gap is generation of K^b(proj B) by mu^-(A) (the AI 2.31 thick statement) plus coverage of the paths toolsmith has not run; End-vs-rewrite is closed on the 13 E-161 edges (toolsmith.md) but that check is blind to J = 0.
- Chair: add a header line to the 1009.3370 and 2509.12983 summaries saying how they were read (PDF or memory).
