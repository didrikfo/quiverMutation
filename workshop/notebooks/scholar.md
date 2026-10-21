# Scholar's notebook (after round 079)

## Believe
- AI Thm 2.31, 2.32(b) are Theorems as cited (AI v3 .tex l. 945, 982); Thm 2.50(c) (l. ~1651-1660) is the vertex form = CHZ Prop 3.5(2) at S={v} (Lemma L1, ours). "Silting not tilting iff some J_i != 0" = 2.31 + 2.32(b), not printed.
- T5, "J = 0 => repo child derived equivalent to parent" has three links: (a) J=0 <=> T tilting (AI 2.32(b)); (b) tilting => End(T) derived equiv (Rickard 1989 Thm 6.4, per repo note; not re-read); (c) repo child = End(T) = H_0(kM): the named gap = Oppermann Thm 1.1/7.5 + Gamma_red=N (E-217, hand) + E-199 (acyclic D(i) only) + E-189 generator count + left/right convention. Prop 7.6 is on none of these for J=0 steps; it only reads H_1 and up off kM.
- Oppermann 1504.02617 v1 ONLY (sha 7cdf284f). Thm 1.1 l. 263, Setup l. 283, Def 3.2 l. 383, Constr 7.3 l. 958, Lemma 7.4 l. 963, Thm 7.5 l. 991, Prop 7.6 l. 1027 (proof 1034-1041, two sentences, degree filtration "strictly smaller degree"). Hypotheses: non-negative, free degree-wise finite, no degree-0 loop at i; nothing on cycles.
- R079: the paper NEVER says "leading words" or "overlap" (grep). That phrase is ours. Shape argument: leading words are (non-star)(star), no star letter begins one, so no combinatorial overlaps. Inertness of partial(phi alpha^-1) (Anick/Dicks, UNVERIFIED memory) is the unproved part; off T5's critical path.
- E-211 (cyclic D(i)): J = 0 yet child fails = repo step 7, not J, not the paper.
- Fosse 2112.08129 (v2): Thm 3.1, 4.2, 4.4, Lemma 5.7 = steps 1-7; slips Thm 3.1 "i != 0", Sect. 6 steps 3/4 swapped.
- Side: AI mu^+ = left (2.32(a)), mu^- = right (2.32(b), 2.50(a)); Oppermann Thm 1.1 = left = mu^+ by name; repo J is the mu^- count; matches Oppermann's A_3/rad^2 example only (one example, not a convention proof).

## Did
- R001..R077 (audits, summaries, litfixes, Oppermann Thm 1.1/Sect. 8, 2112.08129, cancellation, J check). R079: third write-up of Prop 7.6 (E-214, E-220 already had it); new: the paper has no overlap claim, shape argument.

## Next
- Theorist/skeptic: corrected E-211 child = H_0(kM)? machine check of Theta/H_0 on random parents (E-217 open).
- Journal version of Oppermann: does 7.6 gain a proof? Repo J vs AI Hom(S_v,-) convention on one non-monomial case. Rickard 6.4 re-read.
- Lessons: grep EXPERIMENTS for the exact assignment FIRST (7.6 was done in E-214, again in 077, again now); grep the source for the assigned phrase before assuming the paper contains it; a title must not say "not used" unless every compared object avoids it; find an independent example when a name and a formula disagree; count theorem numbers from \newtheorem.
