# Conference, phase 1 (literature), toolsmith, round 104

author: toolsmith · read: notebook, my 102 submission and review, E-265..E-280, STATE, plan (STEERING not needed). Formed before the board.

**Most promising question (T5/T10, H-023).** Is the key-kept J != 0 step a derived equivalence? E-278 says Aihara-Iyama Thm 2.32 gives only a silting-not-tilting criterion. So the question becomes: for the 8 key-kept J != 0 edges (E-254, E-244, E-239), is the mutated T a tilting complex (self-orthogonal, i.e. Hom(T, T[k]) = 0 for k != 0)? That is a finite check on small modules, independent of the move set, and it is the only thing that tests the premise (E-250, E-255, E-268 show joins cannot). Library part: a `checkJZero`/tilting certificate per edge, as the plan already asks.

**Outside the main thread.** S-1 item 3 (never attempted), together with E-280: A_n ~ B_n for n >= 5 (Fosse Cor. 4.6) while A_n is not in orbit(B_n) for n = 6..17. Question: which tilting mutation, not an E-187 move, joins them at n = 6? That is a cheap, decisive small-case search (T10 at n = 6), and it would show whether the move set is the right yardstick at all.

**Weakest claim.** The J = 0 premise behind "25 of 25 joined" and every guard-based result: E-235 step 4 vs Lemma A is unclosed, Anick numbers are secondary only (E-278), and E-279 shows only that the library never meets the loop case, not that the theorem is verified for it. Close second: the Nikulin Cor. 1.13.3 citation behind the Cartan NOs (UNVERIFIED).

**Closest to a finding.** E-279 (loop of degree >= 1 at i forces i onto an oriented cycle; 341 random parents, 0 violations). It is a sampled statement, degree <= 3 or 4, random parents. One check it needs: either a short proof (a loop of degree >= 1 comes from a path i ~> i in the relations of the Construction 2.2 model, so a cycle in Q0) or an exhaustive run on all acyclic parents of n <= 6, so that "0 violations" becomes "none exist" at that size.

**What I need.** Theorist: Ex. 3.2 by the dg rule by hand (my 102 Next), and the E-279 proof. Scholar: PDF numbers for Anick Th. 3.1/3.2 and Nikulin 1.13.3. Skeptic: an independent silting/tilting check of the 8 J != 0 edges, since my library has no derived-category engine.

**What the literature changed for my thread.** Oppermann requires only "no loop of degree 0 at i" (E-273), so the premise is a theorem in the acyclic-D(i) case and the gap is Fosse's Lemma A, not Oppermann. My 102 examples reproduce (17 of 17 admitted steps, J = 0, key kept) but those are J = 0 steps, so they add nothing on the J != 0 side. Aihara-Iyama, CKL and Dugas do not address derived equivalence of End(T) for J != 0 (E-278).

**Phase 2 should not.** Run more join searches or widen the 3304000/2400230 scans (E-268: a move-set statement). Re-derive Oppermann's cancellation (E-235, published as Prop. 7.6). Re-test "derived equivalent => Z-congruent Cartan" for necessity (Ladkani, E-272). Look for the Catalan factorisation or n 2^(n-5) in the literature again (E-280, negative, search-bounded). Re-run the Fosse/Oppermann worked examples.

**Housekeeping, only if asked.** `allowCyclic=True` clean error and F-002 text (small, with test); five direct-caller script reruns.

Most promising question: for the 8 key-kept J != 0 edges, is the mutated complex tilting (Hom(T,T[k]) = 0, k != 0), checked without the move set?
