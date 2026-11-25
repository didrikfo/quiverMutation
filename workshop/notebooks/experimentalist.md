# Experimentalist notebook (rewritten round 138)

## What I now believe

F-057 (a) rebuild, round 138 (`rounds/138/experimentalist_f057.py`, 16 s):
- The n = 11 crossing reproduces with my own Cartan / SNF / key / LNA code. Only the mutation comes from the library.
- Endpoints: LNA 500033030 = dual(303305000) and LNA 403033030.
  - Same key; SNF(C + C^T) tails (1,2,4) / (1,2,2).
  - D_10 (gcd of 10x10 minors) is 8 / 4.
  - The SNF changes once, at step 6.
- The endpoint NO is about two named LNAs. It needs no walk, no J and no Hom; it rests only on the E-308 / E-272 chain.
- The end presentation IS the monomial ideal of 403033030. This settles F-057's open sameIdeal item for the end node.
- Key group census: 16 796 LNAs at n = 11 (no length-1 relations), 79 keys. The group is 72 = 13 + 59, matching E-318.
- All 16 intermediate Cartan matrices are identical over Q, GF(2), GF(3) and GF(5).
- Library `integerCartanMatrix` is the TRANSPOSE of dim e_i A e_j (documented). Remember this when comparing.
- Not checked: that the J = 0 steps are derived equivalences (open in F-057), E-333's Hom computation, and the library search itself.
- No control exists where derived-inequivalence NO is known apart from the Cartan chain. Say so in Power lines.

S-1 / H-026, round 135 (details in rounds/135):
- Key changes happen only through K = 0, on eight classes (n = 13, 15). This is near-definitional; the tail-deletion key control gives thousands of hits.
- Boundary changes are `rules` moves, directed K = 0 -> K = 1, onto rows starting 03.
- At h >= 4 every boundary count is 0.
- Bridge depth: the majority-key claim fails on two small mut bridges. Print every bridge's map.

T7 / H-010, round 133:
- No clean interior (8,4)(10,3) cell lowers at k = 4: 15 decided, 9 capped at 560 s.

Practical:
- About 2 effective cores.
- Fraction charpoly costs about 30 ms per 11x11. Use integer Faddeev in numpy int64 with an overflow guard, and re-key the hits exactly.
- Watch for name clashes: `p` was both a path and a prime in round 138.
- Write the submission early. Print vacuous counts. Power controls belong on the same side of the question.

## What I tried
- 138: F-057 (a) rebuild (walks, census, controls).
- 135: bridge and boundary censuses.
- 133: band / pick / reach / cox.
- 130: screen.
- 127: gap2fix.
- 125: transit.
- 123: classlevel.

## What I would do next
1. F-057 (b) rebuild at n = 12: the group of 47 (14 + 33) and the F3 step at 0040033030 · F3 R9. The script needs N = 12 and R-step support (dualPathAlgebra around the F step).
2. Hom check of the nine J = 0 steps of walk b (with skeptic). This makes "leaves its class at step 6" not depend on the mutation code.
3. S-1 at n = 17 for h >= 4 (overnight, or with a faster key: the int64 Faddeev now exists in the round-138 script).
4. T7: the 9 capped cells at a 1500 s cap.
