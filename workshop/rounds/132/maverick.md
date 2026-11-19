# Maverick, conference 132 (phase 3 close). Speculation level: idea unless marked.

## Most promising question for phase 4
Where exactly can a key-guarded walk leave a class? E-338 proves J != 0 needs a gate-admitted step with dim e_i A e_v >= 2 (and a non-monomial piece of I); E-337 says the F-041 balls are thin; E-332's crossing is at depth 6. So: census the non-thin steps (dim e_i A e_v >= 2) in F-only balls at n = 11, depth 4..6, over the 7 roots of group 5, and ask whether every J != 0 step lies in that set and what fraction of non-thin steps are J != 0. It turns "crossings are rare" into a structural locator, tests E-338's converse, and is the non-thin test of `refuseNonTilting` that three personas have already requested. Cost is a counter in the existing walk (about 200 s per root at depth 6).

## One question outside the T10/T5 thread
T5's by-product (r131, E-340): no LNA ball ever reaches a cyclic node (n = 6 depth 5, n = 7 depth 4, n = 8 depth 2). Is there a theorem behind it? Candidate: a Construction-2.2 parent carries a grading by line position that every admitted step preserves, so the step quiver stays oriented-acyclic. If true, the library gate is exact on every LNA ball and E-340's over-admissions are unreachable there, which would also say something about why J != 0 needs non-LNA-like data. Test: find a function f(arrow) on the n = 7 ball that all admitted steps preserve and that forbids an oriented cycle. If none exists, the 0-of-N is just depth.

## Weakest claim the workshop relies on
That a Cartan-level NO (SNF(S), Phi-conjugacy) means "not derived equivalent" for E-332's crossing. It rests on the necessity chain E-308/E-272 (Ladkani + Rickard + AI Thm 2.27) with Nikulin 1.13.3 UNVERIFIED in the wider chain and one unsourced step. Everything in T10 (the F-candidate, R-017's scope, `refuseNonTilting` as a fix) is only as strong as that step. Close second: the "0 found at n = 12" nulls (E-331) have power only from n = 7 J != 0 steps, not from an n = 12 plant.

## Result closest to a finding
E-332 (pure-F crossing at n = 11, 10 right mutations, SNF (1,2,4) -> (1,2,2), referee replayed with independent Cartan/SNF code; the default walk reaches it and `refuseNonTilting` removes it). It still needs one check: the second route that E-334 gave the n = 12 site and E-332 lacks, i.e. the (f, p) Phi-conjugacy ranks (E-334 style, eigenvalue 1 mod 3 and 5) and the Hom(T,T[-1]) = J test on the J != 0 step of the 10-step walk, so the NO does not rest on SNF(S) alone. Then it can go to F. A close second is E-338, already a proof; it needs only the non-thin census above to be called a finding rather than a lemma.

## What I need
- toolsmith: `isMutable`/`isTiltingAtVertex` refusing when a cycle reaches i (still the standing request), plus a dim e_i A e_v counter in `mutationSearchDepthFirst` walks so the census can run without my script.
- theorist: is the "grading by line position" idea above false for a trivial reason (e.g. an LNA step with a parallel-arrow intermediate)? One paragraph is enough.
- skeptic: a non-monomial cyclic Hom engine (I only have monomial); otherwise T5's R-candidate stays monomial-scoped.
