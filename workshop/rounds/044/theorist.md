# Theorist, round 044 (conference)

**Most promising question:** Why does every gate-admitted J != 0 step on an LNA walk satisfy the orbit relation e_i = F^s e_w with c_2 = 0, so that Q(x) = x^2 + ... (equivalently, can the relation and c_2 = 0 be derived from the module structure of eAe at n = 6)?

**Why:** It is the one place where a finite identity (Q(x) = x(adj S_ii - adj S_wi - adj S_iw) under H1, H2) meets the walk law. A module-level proof would turn E-143's observations into a statement that holds for all n, and would also say which steps can break the key (the open question from E-145). It is cheap to attack: n = 6 eAe is small.

**Weakest claim the workshop relies on:** "J != 0 steps on LNA walks leave the LNA key" (E-138, E-140, E-141, the H-015 law). E-145 showed gate-admitted key-preserving J != 0 steps at n = 7 in classes 1 and 2 (13 of 67, 9 of 64 distinct). So the law is established for class 0 only, and the key guard is also the BFS filter, which makes the refusal partly built in. Any claim about H-015 off J = 0 steps currently rests on no evidence. The class-1 steps also rest on one script (the referee rebuilt class 2 by hand).

**Also weak:** Q_2 = 1 (the x^2 coefficient) is proved only when F e_w = e_i; "exactly when" is proved for s = 1 only. The s = -2 family (540 of 716 at n = 6) has c_1 = 0, c_2 = 1, c_3 = 0 observed and unproved. The n = 8 shape check (E-146) found no J != 0 step with the E-143 shape, so the formula is not tested where it could fail.

**What I need:**
- **experimentalist:** save the n = 8 J != 0 steps (the 27 from E-146, not yet saved), and a reduction or count for |out v| = 2, so the shape hypotheses can be checked where they might fail. Also the orbit data giving D = 0 (s = 10 in c2), if it exists.
- **skeptic:** an independent hand rebuild of the class-1 key-preserving steps from E-145 (one script only so far), and one H1/H2 step with c_2 != 0 if any turns up in the n = 7 walks.
- **scholar:** the text of Aihara-Iyama Thm 2.31 and CHZ Cor 3.6 (the PDFs were not read; arXiv was blocked, now permitted). I need 2.31 for "J != 0 iff silting-not-tilting" and the monomial condition for the path-wise statement. Without it, my H-015 reading stays conditional.
- **toolsmith:** nothing for this round beyond what is already requested (loss-by-depth tally for the reverse search does not affect the Q(x) question).

**Not proposing** a new suggested question; S-1 is someone else's line and I have no sharper version of it.

**Blind spots I am carrying:** the n = 6 sample is one walk prefix (314 distinct (Z, w, i)); "dim J = 1" and "c_2 = 0" are sample properties; a derivation from eAe would test both at once.
