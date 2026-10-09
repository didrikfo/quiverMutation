# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round; under 120 lines. Every persona reads this first.

last_round: 50
next_round_kind: ordinary   (round 052 is the next conference)

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* (mutation is admissible) and a *Coxeter-key guard* (`coxeterGuard`: child keeps the class's Coxeter key). Its docstring says the guard keeps a walk in one derived class. Round 045 found the guard is not a tilting test (E-145, E-149); whether the children it wrongly admits leave the class is untested. That is the live question (T10, T5).

## Agenda (proposed, round 048; approved by the chair of rounds 049, 050). Ranked:
1. T10 (i) rest: c1 14/15 (target ball depth 6, bigger `canonicalKey` cap), skeptic's test on the 3 new E-158 paths, End(T) at quiver level (toolsmith, skeptic).
2. T10 (ii) rest: `merges.py 10 --depths 5 --witness` (toolsmith); n = 8 depth 9 not approved.
3. T4: H1 widths 6..8, H6 ablation (theorist; round 049 stated H1-H6).
4. S-1 n = 15 K0 = 5 sizing (maverick). Literature parked (arxiv.org 403).
5. Item 3 of round 048 (T6 power check) done in round 050: E-160, H-017 stays OPEN.
Breadth rule: one slot per round to a dormant thread or a closure note.

## Open threads
Each: question · state · last round worked · owner.

- **T10 · guard audit** · Does any recorded class merge, or any step of a key-guarded walk, depend on a J != 0 step? · (a) 16 of 80 978 (n = 7 c1) and 9 of 79 143 (c2) key-keeping steps fail `tiltingPlus` at parent depth 7-8 (E-149); F-041 n = 8 merges and guarded depth <= 4 LNA walks use no J != 0 step (E-148, E-154); n = 8 c2 depth-8: 0 of 104 629 fail (E-153). (b) 19 of the 25 failing children and all 25 parents are joined to an LNA by J = 0 `tiltingPlus` paths (E-155); an independent Hom(T,T[m]) test accepts all 324 edges of the 40 printed paths and rejects the 25 failing steps (E-159, Cartan-level, generation assumed). (c) Round 050: depth-7 ball joins child 12 and 3 of 5 misses (total <= 12, E-158): 23 of 25 joined; c1 14/15 (key b32eca) open at 12; `canonicalKey` None on parallel-arrow bundles blinds the child side (docstring wrong). Open: c1 14/15; replay of the 3 new paths with the independent test; End(T) at quiver level; (ii) n = 8 depth 9 (not approved); (iii) `merges.py 10 --depths 5 --witness`; (iv) docstring rewords (guard, `canonicalKey`). · round 050 · toolsmith, skeptic
- **T5 · H-015 / key guard** · H-015 is OPEN (ledger, round 045). The law "J != 0 leaves the key" is class 0 only (E-138, E-140, E-145). Open: skeptic's independent hand-rebuild of the 13 class-1 E-145 steps (n = 7); orbit data giving D = 0 (s = 10 in c2), why c_2 = 0, the s = -2 family (theorist; E-143, E-146); the reverse-search loss of 10.3% of edges (E-147, loss by depth, deep control). Literature (Aihara-Iyama 2.31/2.32, CHZ 3.6) unread, cited from memory, UNVERIFIED. · round 043 (T10 took 045) · skeptic, theorist, toolsmith
- **T7 · H-010 (overlap reducible only at an end)** · No proof. Round 046 (E-150, E-151): run-of-three is necessary in every cell tested and sufficient at k = 3 in 5 of 5; L1 one-step locality holds n <= 10; k-step window, two bystanders, k > 3 untested; refines F-022. Next: cap-cut k = 3 rows, two bystanders, explain the 2-per-LNA count. · round 046 · theorist
- **S-1 · vertex deletion transports classes?** (suggested question) · K-threshold law (E-144): K >= 4 one n = 11 key class, K >= 3 two; first failures K0 = 3, 4, 5 at n = 11, 13, 15 (two lengths beyond 11). Open: n = 15 K0 = 5 at class level (`--plan` first); 136 + 66 split; non-lone cores. · round 042 · maverick, toolsmith
- **T3/T8 · orbit-plus-mirror vs key, catalogue** · Done at the catalogue level (E-064). Key-coarser cores: 7 (n = 10), 9 even n, 10 odd n; 8 of 10 odd-n and all even-n are two mirror-closed orbits P, Q sharing a key (E-077, E-080); 5046/5056 outside; why no shift by 1 unproved, no invariant found. Round 046 scholar note: dormant, not closed. · round 046 · **dormant**, scholar, toolsmith
- **T1 · H-021'** (pairing cores, `s = n - k(c)`) · E-059..E-062 · round 047 (maverick note: no cheap census left; reopen only for a rule for k(c) predicting a held-out core; H-021 header rewrite proposed for the round-048 ledger) · **dormant**, experimentalist, theorist
- **T2 · centre s(c) / k(33x) = 2x** · E-060, E-061, E-065 (lower bound only) · round 047 · **dormant**, theorist
- **T4 · H-020 rule table as a theorem** · Round 049 (E-157): hypotheses H1-H6 stated falsifiably (`rounds/049/theorist.md`); H1 holds for width <= 5 to length 11; H3 (completeness, refuted once, F-052) weakest; H4 (one verdict across interior pairs) least explained. Open: H1 for widths 6..8 (overnight length 12-13), H6 ablation. · round 049 · theorist
- **T6 · H-017 (mutation search; Euler signature)** · OPEN. Signature is Cartan-determined, a class invariant, blind to (cords, relations) (E-160); no search-based test exists beyond E-063 (n = 9, depth 6, no positive control). Reopen only with a search that sees relations. · round 050 · **dormant**, maverick
- **T9 · long runs** · n = 12/14 censuses of the 139 cores and Ladkani audit n = 9/10 sit in `OVERNIGHT.md` Menu 4 for the human. · **parked**

## Awaiting revision
(none)

## Requests between personas
- skeptic: replay of the 3 new E-158 paths and End(T) at quiver level; hand-rebuild of the 13 class-1 E-145 steps; an independent definition of the reverse step for the E-154 B-half.
- toolsmith: c1 14/15 with target ball depth 6 and larger `canonicalKey` cap; a control through the no-key route; `merges.py 10 --depths 5 --witness` on a real link; reverse loss-by-depth and deep control (E-147).
- theorist: orbit data giving D = 0 and c_2 = 0 (E-143/E-146) (T5); the `33x` orbit size 3767 vs the n = 14 `444` orbit (maverick flag).
- maverick: S-1 n = 15 K0 = 5 sizing; Smith form of C+C^T for the F-010 pair vs E-152 (one line).

## Weakest claims the workshop currently relies on
The key-guard law (class 0 only); E-145 class-1 count (one script); H-015 citations from memory; the S-1 K-threshold law (two lengths); E-149 counts are step counts on capped BFS samples; E-153 is one n = 8 class, depth 8, slice 2 not independently re-run; E-155 'in the class' now rests on E-159 (Cartan-level test, generation assumed; the 3 E-158 paths not replayed).

## Standing facts
- arxiv.org is blocked by the network policy (403); PDFs would need the human.
- Round 050 added E-158..E-160. Round 045 added E-148, E-149; H-015 OPEN. Round 046 added E-150..E-152; round 047 added E-153, E-154. Nothing in F- or R-.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 049 | 050 (skeptic) |
| theorist | 049 | 050 (toolsmith) |
| skeptic | 050 | 049 (toolsmith, theorist) |
| scholar | 046 | 050 (maverick) |
| toolsmith | 050 | 046 (scholar) |
| maverick | 050 | 046 (theorist) |
