# Scholar, round 048 (conference)

Read: personas/scholar.md, notebooks/scholar.md, STATE.md, DIGEST.md (rounds 044-047), STEERING.md. Nothing run.

**Most promising question.** Do the E-151 / E-147 key-keeping, `tiltingPlus`-failing children (n = 7 c1, c2) lie in the derived class of an LNA? Test it constructively: search for a sequence of tilting mutations (gate plus `tiltingPlus`, key guard off) from each child back to any LNA of that class. A tilting path is a derived equivalence, so a hit settles membership positively and is cheap to check. A miss is only a bounded miss, and then a finer invariant is needed. This is the question STATE already names as "needs a finer invariant or a tilting path back"; the path search is the one that can decide it now.

**Outside the thread.** T6 / H-017 (Euler signature, dormant since round 023). Before anyone spends a run on it: the signature of the symmetrised Euler form is a function of the Cartan congruence class (Sylvester), and E-154 already shows the failing children are congruent to their keys. So if H-017 is that signature, it cannot separate them and is only a control. Check the H-017 header against this before it is revived; if it is a different invariant, say which.

**Weakest claim.** That the Coxeter key is a proxy for derived class, in the guard docstring of `mutationSearchDepthFirst` ("keeps the walk in one derived class"). The support is: class 0 only (E-140, E-142), E-147 class-1 counts from one script, and AI 2.31/2.32 and CHZ 3.6 cited from memory (arXiv still 403). The key equality is a necessary condition for nothing the workshop has verified. The docstring should stay unreworded until the path test in the first question reports.

**What I need.**
- toolsmith: BFS from the 25 E-154 children over tilting-only edges (gate + `tiltingPlus`, key guard off) to any LNA, with a positive control (a known tilting-connected pair at n = 7) and closure flags.
- skeptic: independent hand check of two or three children, as the tilting path is found.
- theorist: is the Coxeter polynomial a derived invariant for these algebras (finite global dimension check for the children), and does that rule the children out, or does it agree with the key?
- experimentalist: the E-155 edge tally on n = 7 c1, c2 (already requested).
- human: PDFs of 1009.3370 and 2509.12983 remain the only way to replace the UNVERIFIED flags.
