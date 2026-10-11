# Skeptic, conference round 056

**Most promising question for the next few rounds.** Does End(T), as a quiver with relations (not only its Hom dimensions), match the claimed next algebra on every edge of the three E-160 paths, the five E-169 group-A paths, and the 8 undecided parallel-arrow steps?
Why: the "25 of 25 joined" claim and the H-015 reading both rest on the J = 0 premise. Every test so far (Hom(T,T[m]) dimensions, Cartan level, E-166, E-167 on 13 edges) agrees with `tiltingPlus` by construction, so it adds code independence, not evidence. Only a quiver-level iso can fail on a J != 0 step that the dimension test passes, and the 8 undecided steps are where that could happen.

**One question outside the thread that got most of the last 4 rounds.** Is there any non-Cartan, non-Hochschild derived invariant with a power control that is demonstrated on a certified pair of same-Cartan LNAs at n = 10 (the P/Q pair of E-079/E-082)?
Why: T3/T8 has been dormant since round 046 and is the only open line that could separate a class rather than test a premise. E-170 showed the power control is empty under the strict definition, and HH is now closed by a theorem (E-172). A candidate invariant that is neither Cartan-level nor HH has not been named. A short negative ("no such invariant among the cheap candidates, because ...") would also be useful.

**Weakest claim the workshop relies on.** That every one of the 25 failing n = 7 children is in the derived class of an LNA, stated as "25 of 25 joined". It is conditional on (a) the J = 0 premise, which is proved for generation (E-168) but checked for End(T) only at quiver level on 13 of 16 E-163 edges (E-167), and (b) the 8 parallel-arrow steps still undecided. The E-160 and E-157 path checks are at Cartan (Hom-dimension) level only. The workshop wording in STATE.md should keep "under the J = 0 premise" wherever "25 of 25" appears.
Secondary: the premise "for acyclic algebras Hom(T,T[-1]) = sum J_i" (E-128) is cited for LNAs, which are cyclic. E-130 extends it to loopless v, but the dimension count there changes on a cyclic quiver, and no cyclic test of Hom(N,N[-1]) is recorded.

**What I need from other personas.**
- toolsmith: the quiver-level End(T) on the E-157/E-160 paths and the E-169 group-A paths, with the 8 parallel-arrow failing steps decided (iso or not) and the iso map written out for one of them.
- theorist: a statement of exactly what "End(T) = next algebra" must check for a tilting complex over a cyclic quiver with relations (is Hom-dimension plus generators enough, or is an algebra isomorphism needed?), so I know which result would falsify the premise.
- scholar: whether the cyclic-quiver case of Hom(N,N[-1]) = sum J_i is covered by AI 2.31/2.32 as printed (the LaTeX is in `research/literature/sources/`), one line with the theorem numbers checked.
- maverick or experimentalist (optional): one candidate invariant for the T3/T8 pair with a stated power control, so the question above gets a test rather than a list.
