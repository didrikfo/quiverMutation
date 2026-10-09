# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round; under 120 lines. Every persona reads this first.

last_round: 52
next_round_kind: ordinary   (round 053)

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* (mutation is admissible) and a *Coxeter-key guard* (`coxeterGuard`: child keeps the class's Coxeter key). Its docstring says the guard keeps a walk in one derived class. Round 045 found the guard is not a tilting test (E-147, E-151). Under the J = 0 premise all 25 failing n = 7 children are now joined to an LNA (E-157, E-160, E-163); the premise itself is the live question (T10, T5).

## Agenda (proposed, round 052; ordinary rounds work from it until the human changes it). Ranked:
1. T10 premise, non-Cartan test: Hom replay of the 3 E-160 paths and quiver-level End(T) on the 13 E-163 edges (skeptic, toolsmith); which hypothesis of AI 2.31/2.32 the J = 0 step needs, generation (theorist).
2. T10 (ii): group-A witness path 05040330 -> 33460000 by a tilting-only, relabelling-aware search, then `tiltingPlus` replay (experimentalist, toolsmith). No overnight.
3. T3/T8 (outside T10): finer non-Cartan derived invariant separating the P and Q orbits at n = 10; power control first (equal-Cartan, different-class pairs at n <= 9) (maverick, theorist, toolsmith).
4. T4/H-020: cores with >= 3 relations at n = 13, 14 where rules = [] changes a verdict, else restate H-020 for the move set; H6 in a rules-only walk (theorist).
5. Breadth T7/H-010: run-of-three bystander at k = 4, two arrows away, with positive control (toolsmith/theorist). S-1 n = 15 sizing after item 1.
Breadth rule: one slot per round to a dormant thread or a closure note. Toolsmith round: `canonicalKey` cap and docstring rewords (guard, `canonicalKey`).

## Open threads
Each: question · state · last round worked · owner.

- **T10 · guard audit** · Does any recorded class merge, or any step of a key-guarded walk, depend on a J != 0 step? · (a) 16 of 80 978 (n = 7 c1) and 9 of 79 143 (c2) key-keeping steps fail `tiltingPlus` at parent depth 7-8 (E-151); F-041 n = 8 merges and guarded depth <= 4 LNA walks use no J != 0 step (E-150, E-156); n = 8 c2 depth-8: 0 of 104 629 fail (E-155). (b) 19 of the 25 failing children and all 25 parents are joined to an LNA by J = 0 `tiltingPlus` paths (E-157); an independent Hom(T,T[m]) test accepts all 324 edges of the 40 printed paths and rejects the 25 failing steps (E-161, Cartan-level, generation assumed). (c) Round 050: depth-7 ball joins child 12 and 3 of 5 misses (total <= 12, E-160); `canonicalKey` None on parallel-arrow bundles blinds the child side (docstring wrong). Round 051: c1 14/15 joined at total 13 (7 + 6, E-163; skeptic replayed with Hom test), so 25 of 25 joined under the premise; Hom-tested only for E-157 paths and this one. (d) Round 051 (E-164): `merges.py 10 --depths 3 4 5` found no link; F-037 (one of two n = 10 merges) 5 paths use no J != 0 step; group-A merge (05040330 -> 33460000) unreplayed, no path recorded; `meetingPoints` uses labelled `quiverKey`. Open: Hom replay of the 3 E-160 paths; End(T) at quiver level; group-A witness path; n = 8 depth 9 (not approved); docstring rewords (guard, `canonicalKey`). · round 051 · toolsmith, skeptic, scholar
- **T5 · H-015 / key guard** · H-015 is OPEN (ledger, round 045). The law "J != 0 leaves the key" is class 0 only (E-140, E-142, E-147). Open: skeptic's independent hand-rebuild of the 13 class-1 E-147 steps (n = 7); orbit data giving D = 0 (s = 10 in c2), why c_2 = 0, the s = -2 family (theorist; E-145, E-148); the reverse-search loss of 10.3% of edges (E-149, loss by depth, deep control). Literature (Aihara-Iyama 2.31/2.32, CHZ 3.6) unread, cited from memory, UNVERIFIED. · round 043 (T10 took 045) · skeptic, theorist, toolsmith
- **T7 · H-010 (overlap reducible only at an end)** · No proof. Round 046 (E-152, E-153): run-of-three is necessary in every cell tested and sufficient at k = 3 in 5 of 5; L1 one-step locality holds n <= 10; k-step window, two bystanders, k > 3 untested; refines F-022. Next: cap-cut k = 3 rows, two bystanders, explain the 2-per-LNA count. · round 046 · theorist
- **S-1 · vertex deletion transports classes?** (suggested question) · K-threshold law (E-146): K >= 4 one n = 11 key class, K >= 3 two; first failures K0 = 3, 4, 5 at n = 11, 13, 15 (two lengths beyond 11). Open: n = 15 K0 = 5 at class level (`--plan` first); 136 + 66 split; non-lone cores. · round 042 · maverick, toolsmith
- **T3/T8 · orbit-plus-mirror vs key, catalogue** · Done at the catalogue level (E-066). Key-coarser cores: 7 (n = 10), 9 even n, 10 odd n; 8 of 10 odd-n and all even-n are two mirror-closed orbits P, Q sharing a key (E-079, E-082); 5046/5056 outside; why no shift by 1 unproved, no invariant found. Round 046 scholar note: dormant, not closed. · round 046 · **dormant** (proposed agenda item 3, round 052), scholar, toolsmith
- **T1 · H-021'** (pairing cores, `s = n - k(c)`) · E-061..E-064 · round 047 (maverick note: no cheap census left; reopen only for a rule for k(c) predicting a held-out core; H-021 header rewrite proposed for the round-048 ledger) · **dormant**, experimentalist, theorist
- **T2 · centre s(c) / k(33x) = 2x** · E-062, E-063, E-067 (lower bound only) · round 047 · **dormant**, theorist
- **T4 · H-020 rule table as a theorem** · Round 049 (E-159): hypotheses H1-H6 stated falsifiably; H1 holds for width <= 5 to length 11. Round 051 (E-165): H1 spot check at one length each for widths 6..8 (0 failures in 43 116), gaps at lengths 10-12 and widths 9..11; H6: anchored rules redundant only in the reduced walk, matter in a rules-only walk. H3 (completeness, F-052) weakest; H4 least explained. · round 051 · theorist
- **T6 · H-017 (mutation search; Euler signature)** · OPEN. Signature is Cartan-determined, a class invariant, blind to (cords, relations) (E-162); no search-based test exists beyond E-065 (n = 9, depth 6, no positive control). Reopen only with a search that sees relations. · round 050 · **dormant**, maverick
- **T9 · long runs** · n = 12/14 censuses of the 139 cores and Ladkani audit n = 9/10 sit in `OVERNIGHT.md` Menu 4 for the human. · **parked**

## Awaiting revision
(none)

## Requests between personas
- skeptic: Hom replay of the 3 E-160 paths and End(T) at quiver level; hand-rebuild of the 13 class-1 E-147 steps; an independent definition of the reverse step for the E-156 B-half.
- toolsmith: a control with minimum exactly 13; an orbit-based canonical form for the 12 high-cost keyless nodes; a relabelling-aware `meetingPoints` (E-164); a control through the no-key route; `merges.py 10 --depths 5 --witness` on a real link; reverse loss-by-depth and deep control (E-149).
- theorist: orbit data giving D = 0 and c_2 = 0 (E-145/E-148) (T5); the `33x` orbit size 3767 vs the n = 14 `444` orbit (maverick flag).
- maverick: S-1 n = 15 K0 = 5 sizing; Smith form of C+C^T for the F-010 pair vs E-154 (one line).

## Weakest claims the workshop currently relies on
The key-guard law (class 0 only); E-147 class-1 count (one script); H-015 citations from memory; the S-1 K-threshold law (two lengths); E-151 counts are step counts on capped BFS samples; E-155 is one n = 8 class, depth 8, slice 2 not independently re-run; 'in the class' for the 25 children rests on the J = 0 premise, Hom-tested (Cartan level, generation assumed) only for the E-157 and E-163 paths; the E-165 length check is one length per rule.

## Standing facts
- arxiv.org is blocked by the network policy (403); PDFs would need the human.
- Round 052 was a conference: no change to any H-, F-, R- status (ledger in proceedings). Round 051 added E-163..E-165. Round 050 added E-160..E-162. Round 045 added E-150, E-151; H-015 OPEN. Round 046 added E-152..E-154; round 047 added E-155, E-156. Nothing in F- or R-.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 051 | 050 (skeptic) |
| theorist | 051 | 050 (toolsmith) |
| skeptic | 050 | 051 (toolsmith, theorist) |
| scholar | 046 | 051 (experimentalist) |
| toolsmith | 051 | 046 (scholar) |
| maverick | 050 | 046 (theorist) |
