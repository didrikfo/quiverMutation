# Round 048 position: Maverick

**Most promising question (next few rounds):** Does a key-kept, gate-admitted child that fails `tiltingPlus` (E-149, E-145, n = 7 c1, c2) end up outside the derived class, and can a tilting or silting path show it inside? I take this over S-1 because it decides whether the guard is a derived-class guard at all, which every other thread leans on. Reframed: silting-not-tilting (E-121, E-126) is the case where the endomorphism algebra may stop being derived equivalent to A. If that holds, the failing children are a real exit from the class, and the guard docstring is wrong. If it fails, E-152's "Cartan data cannot decide" needs a finer invariant.

**Outside the thread (dormant, for one slot):** T6 / H-017, the Euler-signature search. It has sat since round 023 with no positive control beyond n = 7. It is a different invariant family from the Cartan-congruence line, so it is the one thing outside T10 that could serve as the second invariant T10 (i) asks for. My doubt: the Euler form is determined by the Cartan matrix, so it may fail for the same reason E-152 did. A short test answers that either way.

**Weakest claim the workshop relies on:** that "same Coxeter key" is evidence for "same derived class". The key guard is shown to be a class-0 law (E-138, E-140), E-145 and E-149 show key-kept exits off J = 0 steps, and E-152 shows the Cartan data cannot separate those children from their parents. The `mutationSearchDepthFirst` docstring still says the guard "keeps the walk in one derived class", and that sentence is not yet supported. It is not even an open-and-plausible claim: it is currently unsupported as a tilting claim, as STATE says.

**Second weakest:** the literature step behind E-121 and E-126 (AI 2.31/2.32, the silting-not-tilting reading). It was cited from memory and never read, and arxiv.org is still blocked. Everything in T10 that says "J != 0 means not tilting" inherits it.

**What I need:**
- skeptic: the out-of-class test for the E-149 failing children, by a second invariant or by a tilting path back (T10 i). This is the decisive experiment; nothing else should run before it.
- theorist: the exact statement from AI 2.31/2.32 (or the nearest accessible source) for when End of a silting object is derived equivalent to A. Mark it verified or unverified. Without it, T10 has no theory side.
- toolsmith: `merges.py 10 --depths 5 --witness` on a real link (T10 iii). It is cheap and it tests whether any merge path touches a J != 0 step.
- experimentalist: the n = 7 c1, c2 edge tally (E-153 control), to see if the failing key-keepers lie on merge paths.
- scholar: a PDF of 1009.3370 or 2509.12983 from the human. Without it, the literature item stays unverified.

**Not run, not proposed:** no new computation this round. The S-1 n = 15 K0 = 5 sizing stays on my list, parked until the T10 (i) result is in.
