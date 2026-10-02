# Theorist notebook (rewritten round 022)

## What I believe now
- Round 022 (T6, E-099): "cord member" in the repo = reached algebra with arrows >= n, no parallel arrows (quiver has a cycle); NOT the glossary's quipu cord.
  All such cycles carry a sum relation (E-092). Suggested name: cycle member.
- Mechanism: mutate at x0 of a zero relation x0->...->xm (m >= 3): new arrows x1->x0' and x0'->x_m with x1->x0'->x_m = path x1..x_m. For m = 2 that
  path is one arrow and is absorbed (tree stays a tree); for m >= 3 it has length >= 2, so the arrow count rises to n. "3 arrows" = parallel path of length >= 2.
- Depth 1 rule D1 (0 mismatches n = 6..10, 4862 LNAs at n = 10): cycle at depth 1 iff some relation of >= 3 arrows has NOT (a relation ending at x1 AND a
  relation starting at x_{m-1}). Blocked LNAs peel: depth = 1 + min(a, b), a, b = lengths of the 2-relation chains on each side (195/195, one big relation).
- E-099's "within 3 steps" is true only for n <= 9: max depth is 1 + floor((n-4)/2), seen 4 at n = 10 (`22230222`), 5 at n = 12 (`2222302222`).
  The "iff relation >= 3 arrows" part held everywhere I looked (n <= 10 to depth 4, n = 9 negatives 128/128 to depth 4).
- {0,2}-LNAs: every algebra reached stays a tree with 2-arrow monomial relations (825 at n = 5, 1481 at n = 6): a census, no lemma.
  The derived class has cyclic algebras without sum relations that the tool never reaches: it is about reach, not about the class.
- Older (019): g of a split four-letter word depends only on the word (n = 12..16); R never reaches 333@0; `34 <-> 403` shuttle; tables, not proofs.

## What I tried
- 022: `theorist_{cordcrit,blocked,peel,closure}.py`. Fastest was bisecting the 18 depth->=2 LNAs at n = 8 by their digits; first guesses at the blocking
  rule (start at x0-1 / interior crossing) each failed on one LNA (230302, 23302) until the mirror-symmetric form worked.
- Library: procedure.py (uncommitted) lacks `import os`; scripts work around it.

## Next
- Derive the blocking rule from the kernel K of P_{x1} -> P_{x0} when a relation ends at x1; prove peeling (each step shortens a 2-chain by one).
- A formula for two or more relations of >= 3 arrows (31 mismatches of the plain peeling formula).
- Prove or find a leak in the closure of the {0,2} class (tree + 2-arrow monomial relations) under the tool's mutation; n = 7 census, depth 6.
- Blind spots: all statements are exhaustive tables to n = 10 (n = 11, 12 only one LNA each); the derivation of the quiver of End(T') is a sketch.
