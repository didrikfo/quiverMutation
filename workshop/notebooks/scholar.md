# Scholar's notebook (after round 077)

## Believe
- AI Thm 2.31, 2.32(b) are Theorems as cited (AI v3 .tex l. 945, 982); Thm 2.50(c) (l. ~1651-1660) is the vertex form = CHZ Prop 3.5(2) at S={v} (Lemma L1, ours). "Silting not tilting iff some J_i != 0" = 2.31 + 2.32(b), not printed.
- T5, "J = 0 => repo child derived equivalent to parent" has three links: (a) J=0 <=> T tilting (AI 2.32(b); citation, abstract T); (b) tilting => End(T) derived equiv (Rickard 1989 Thm 6.4, note says so; not re-read); (c) repo child = End(T) = H_0(kM): THE named gap, = Oppermann Thm 1.1/7.5 + Gamma_red=N (E-217, hand) + E-199 (acyclic D(i) only) + E-189 Q^(<=2)/generator count + left/right convention. Prop 7.6 is on none of these for J=0 steps; it only reads H_1 off kM.
- Oppermann 1504.02617 v1 ONLY (sha 7cdf284f). Per-section numbers: Thm 1.1 l. 263, Setup l. 283, Def 3.2 l. 383, Constr 7.3 l. 958, Lemma 7.4 l. 963, Thm 7.5 l. 991, Prop 7.6 l. 1027 (proof 1034-1042, two sentences). Hypotheses: non-negative, semi-free degree-wise finite, no degree-0 loop at i; nothing on cycles. They enter only in (c): generators of M non-negative, H_0(kM) needs no cancellation (J_0 in d kM_1).
- E-211 (cyclic D(i), source vertex 1): J = 0 (R077 script, repo convention x in e_jAe_v killed by all arrows out of v), so T tilting and End(T) ~ parent; failure = repo step 7 (vertex-simple tails), not J, not the paper. J cannot serve as a guard for it.
- R074: triangularity of every term of d(phi alpha^-1); leading words never overlap; dg inertness step (Anick, UNVERIFIED) not proved.
- Fosse 2112.08129 (v2 e-print): Thm 3.1, 4.2, 4.4, Lemma 5.7 = steps 1-7; no LNA rule; slips Thm 3.1 "i != 0", Sect. 6 steps 3/4 swapped.
- Left/right and op convention of Thm 3.1 vs repo arrows: matched only through examples.

- Side: AI mu^+ = left (left approx, 2.32(a)), mu^- = right (2.32(b), Okuyama-Rickard 2.50(a)); Oppermann Thm 1.1 = left = mu^+ by name. Repo J = Hom(S_v,e_jA) is the mu^- (2.50(c)) count; it agrees with Oppermann's printed A_3/rad^2 example (tilting at 1, not at 2); the dual count (paths out of v killed by arrows in) says not tilting at both, so it is not Oppermann's side. A match by one example, NOT a convention proof (R077 response). E-211: J^dual=(1,0,0).

## Did
- R001..R074 (audits, summaries, litfixes, Oppermann Thm 1.1/Sect. 8, 2112.08129, cancellation). R077: 7.6 written out again (already E-214!), factored the J=0 implication; script rounds/077/scholar_jcheck.py (J on E-211 and A_3/rad^2 controls).

## Next
- Theorist/skeptic: corrected E-211 child (rbar.b.c = 0) = H_0(kM)? Machine check of Theta/H_0 on random parents (E-217 open); the dg inertness step.
- Ask whether the journal version of Oppermann corrects 7.6; check the repo's J vs AI's Hom(S_v,-) convention against the repo's gate on one non-monomial case.
- Trace module/composition conventions of AI vs Oppermann line by line; re-read Rickard 6.4.
- Lessons: a title must not say "not used" unless every object compared avoids it; when a name (left/right) and a formula disagree, find the independent example before claiming a match; grep EXPERIMENTS for the exact assignment first (E-214 had done 7.6 already; I re-derived); check which object each theorem computes (Gamma, N, kM, End(T)) before saying where a proposition is used; write a model before arguing signs; count theorem numbers from \newtheorem and \numberwithin.
