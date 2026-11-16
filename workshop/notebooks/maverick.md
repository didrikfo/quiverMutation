# Maverick's notebook

## What I now believe (after round 129)
- T10, pure F (r129): **the library's F-only key-guarded walk crosses at n = 11.** From LNA 303500030 (·F1 F1 F4 F3 F1, then F8) and from dual 303305000 (·F6 F11 F7 F11 F7, then F7) five J = 0 steps reach a parent whose key-kept J != 0 non-tilting step changes SNF(S) (1,2,4) -> (1,2,2). `mutationSearchDepthFirst` with its defaults reaches the same nodes at depth 6 (maverick_libwalk.py). `refuseNonTilting=True` removes them. These 2 parents (up to iso) are new, not among the 13 F/R parents of r125. **Pure-F LNA-to-LNA join across classes:** dual(303305000) · F6 F11 F7 F11 F7 F7* F3 F3 F4 F7 = LNA 403033030 (59-class), checked with `walk`. The other child (P_a's) reaches no LNA within F depth 5.
- Round 127's "every crossing needs an R move" was a depth artefact. My r127 guess was "pure-F crossing would be deeper than 10 min at n = 12". At n = 11 it is cheap (depth 6, about 4 min).
- E-331 (n = 12, F-only depth 5, 0 J != 0) stands. The 14-class's 7 orbit reps at depth 5 also give 0.
- **R inverts F**, but only up to the ideal: the referee's 3 "non-invertible" edges are the same algebra with other generators. (5,6)(6,7)(7,8) vs (5,6)(6,3)(3,8) are equal mod the shared commutativity relation. Over 3 135 J = 0 edges in 4 balls, 0 lack an inverse and 55 are ideal-only. Struct iso and canonicalKey compare GENERATORS (F-049). Use `sameIdeal` (maverick_inv.py) when an iso miss matters. E-149's n = 6 "lost in reverse" (10.3 %) is probably the same thing (untested).
- The round-127 backward R-only balls are only a lower bound (presentation dedup, LNA recognition by monomial generators).
- Reframing (idea): J = 0 F-walks are descending chains in the silting poset of L. The new sites are non-tilting irreducible steps 5 tilting steps below an LNA, the cleanest test case for "J != 0 <=> silting-not-tilting".
- E-318 pair (r125): X_a = coker(Phi - a) derived invariant; X_+-3 separate; fails as a function of (key, SNF) for general forms. S-1 / H-026: K0 = min(h, n-4-h) at odd n.

## What I tried
- r129: workshop/rounds/129/maverick_inv.py (ctrl7 | inv | cls | diag | jx | after | idealctrl) and maverick_libwalk.py; outputs in maverick_inv_out.txt. Costs: cls n = 11 depth 6 is 75-273 s per root; n = 12 depth 5 is 61-186 s; inv n = 12 depth 4 is 389 s; libwalk n = 11 depth 6 is about 200 s.
- Power controls: n = 7 class 1 depth 8 gives E-204's 4 J != 0 steps exactly (multi-seed, 274 s; `fwd 7 8 30400` / `40030` give 2 + 2). The FR depth-3 balls find 12 of 12 parents. idealctrl: 6 of 6.
- r127 pureF fwd/back/ctrl; r125 X_a invariants; r121 census47.

## Watch for
- Venv loses python-flint/cypari: `uv pip install --python .venv/bin/python cypari python-flint`. Never use sympy SNF/det at n >= 11 (INV in toolsmith_jxsearch is sympy).
- Exec-importing scripts: split on "\nif __name__ == '__main__':". A str.replace on that marker hits the exec line too (cost me a broken file).
- Never start a job whose output file a queued job also writes. Check `ps` first.
- `procedure.relationsFrom` ignores a hand-set arrowRels unless `rels` agrees. Build test algebras with `procedure.toPathAlgebra`.
- A NO needs a control that says NO. An iso NO needs the ideal test.

## Next
- Overnight: n = 11 13-class F-only depth 7 (all 13 roots); n = 12 14-class depth 6. Count the pure-F crossing sites; do any reach a 59-class LNA (pure-F LNA-to-LNA join)?
- Theorist: Hom-level silting check at P_a·F8, P_b·F7.
- Rerun E-149's lost edges with sameIdeal; if all are recovered, E-149's 10.3 % is an iso artefact.
- Lemma (LNA-only): is coker(Phi - 3) determined by (key, SNF(S))? n = 13.
