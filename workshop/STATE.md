# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round; under 120 lines. Every persona reads this first.

last_round: 55
next_round_kind: conference   (round 056)

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* (mutation is admissible) and a *Coxeter-key guard* (`coxeterGuard`: child keeps the class's Coxeter key). Its docstring says the guard keeps a walk in one derived class. Round 045 found the guard is not a tilting test (E-145, E-149). Under the J = 0 premise all 25 failing n = 7 children are now joined to an LNA (E-155, E-158, E-161); the premise itself is the live question (T10, T5).

## Agenda (proposed, round 052; approved by default, round 053; ordinary rounds work from it until the human changes it). Ranked:
0. **Human's special request (STEERING, round 054):** the two papers arXiv:1009.3370 and 2509.12983 are now in `research/literature/sources/` (LaTeX). DONE round 055 (E-169): AI 2.31/2.32 and the whole-T reading of J = 0 verified against the LaTeX; CHZ Cor 3.6 has no printed 'monomial' hypothesis (iff for monomial I is our derivation). Left: annotate the remaining literature files listed in E-169.
1. T10 premise: generation is now proved per loop-free step (E-166, 370 printed edges loop-free); quiver-level End(T) done for the 13 E-161 edges (E-165), Cartan-level Hom replay for the 3 E-158 paths (E-164). Next: the 8 undecided parallel-arrow failing steps (toolsmith), quiver-level End(T) on the E-155/E-158 paths and class 2 (toolsmith, skeptic).
2. T10 (ii): done at n = 10 for group A (E-167: 5 explicit 7-step J = 0 paths). Left: Hom/End(T) check on those 35 edges (skeptic, toolsmith); relabelling-aware `meetingPoints` (toolsmith).
3. T3/T8 (outside T10): finer non-Cartan derived invariant separating P and Q at n = 10. Round 054: HH is dead (a theorem, 2312.14699); the power control is empty under the strict definition (E-168, under revision). Next: a candidate that is not Cartan-level nor HH, powered on certified pairs with different profile (maverick, revision).
4. T4/H-020: cores with >= 3 relations at n = 13, 14 where rules = [] changes a verdict, else restate H-020 for the move set; H6 in a rules-only walk (theorist).
5. Breadth T7/H-010: run-of-three bystander at k = 4, two arrows away, with positive control (toolsmith/theorist). S-1 n = 15 sizing after item 1.
Breadth rule: one slot per round to a dormant thread or a closure note. Toolsmith round: `canonicalKey` cap and docstring rewords (guard, `canonicalKey`).

## Open threads
Each: question · state · last round worked · owner.

- **T10 · guard audit** · Does any recorded class merge, or any step of a key-guarded walk, depend on a J != 0 step? · (a) 16 of 80 978 (n = 7 c1) and 9 of 79 143 (c2) key-keeping steps fail `tiltingPlus` at parent depth 7-8 (E-149); F-041 n = 8 merges and guarded depth <= 4 LNA walks use no J != 0 step (E-148, E-154); n = 8 c2 depth-8: 0 of 104 629 fail (E-153). (b) 19 of the 25 failing children and all 25 parents are joined to an LNA by J = 0 `tiltingPlus` paths (E-155); an independent Hom(T,T[m]) test accepts all 324 edges of the 40 printed paths and rejects the 25 failing steps (E-159, Cartan-level, generation assumed). (c) Round 050: depth-7 ball joins child 12 and 3 of 5 misses (total <= 12, E-158); `canonicalKey` None on parallel-arrow bundles blinds the child side (docstring wrong). Round 053: E-164 all 33 edges of the 3 E-158 paths pass the Hom test (Cartan level); E-165 End(T) iso to the next algebra on 13/13 E-161 edges at quiver level (8 of 16 failing J != 0 steps decided, iso, so the check cannot see J; 8 undecided, parallel arrows); 25 of 25 = 23 by own path + 2 (c1 13, 15) by key equality. Generation still assumed. Round 051: c1 14/15 joined at total 13 (7 + 6, E-161; skeptic replayed with Hom test), so 25 of 25 joined under the premise; Hom-tested only for E-155 paths and this one. (d) Round 051 (E-162): `merges.py 10 --depths 3 4 5` found no link; F-037 (one of two n = 10 merges) 5 paths use no J != 0 step; group-A merge (05040330 -> 33460000) unreplayed, no path recorded; `meetingPoints` uses labelled `quiverKey`. Round 054: generation of K^b(proj) by T proved per step at loopless vertices, loop-free on 370 printed edges (E-166); group-A merge has 5 explicit 7-step J = 0 paths (E-167; E-033 §4 had recorded the length). Open: quiver-level End(T) on the E-155/E-158/E-167 paths and class 2; the 8 parallel-arrow cases; n = 8 depth 9 (not approved); docstring rewords (guard, `canonicalKey`). · round 054 · toolsmith, skeptic, scholar
- **T5 · H-015 / key guard** · H-015 is OPEN (ledger, round 045). The law "J != 0 leaves the key" is class 0 only (E-138, E-140, E-145). Open: skeptic's independent hand-rebuild of the 13 class-1 E-145 steps (n = 7); orbit data giving D = 0 (s = 10 in c2), why c_2 = 0, the s = -2 family (theorist; E-143, E-146); the reverse-search loss of 10.3% of edges (E-147, loss by depth, deep control). Literature (Aihara-Iyama 2.31/2.32, Pavon 2509.12983 3.6) now VERIFIED against the LaTeX (round 055, E-169). · round 043 (T10 took 045) · skeptic, theorist, toolsmith
- **T7 · H-010 (overlap reducible only at an end)** · No proof. Round 046 (E-150, E-151): run-of-three is necessary in every cell tested and sufficient at k = 3 in 5 of 5; L1 one-step locality holds n <= 10; k-step window, two bystanders, k > 3 untested; refines F-022. Next: cap-cut k = 3 rows, two bystanders, explain the 2-per-LNA count. · round 046 · theorist
- **S-1 · vertex deletion transports classes?** (suggested question) · Round 055 (E-171): K0 = 5 at n = 15 confirmed at key level on a one-row witness; items 1-3 of the STEERING question still open; full n = 15 scan sized (~33 min on 4 cores, 24 shards), not run. K-threshold law (E-144): K >= 4 one n = 11 key class, K >= 3 two; first failures K0 = 3, 4, 5 at n = 11, 13, 15 (two lengths beyond 11). Open: n = 15 K0 = 5 at class level (`--plan` first); 136 + 66 split; non-lone cores. · round 055 · maverick, toolsmith
- **T3/T8 · orbit-plus-mirror vs key, catalogue** · Done at the catalogue level (E-064). Key-coarser cores: 7 (n = 10), 9 even n, 10 odd n; 8 of 10 odd-n and all even-n are two mirror-closed orbits P, Q sharing a key (E-077, E-080); 5046/5056 outside; why no shift by 1 unproved, no invariant found. Round 046 scholar note: dormant, not closed. · round 054 (maverick, revision awaited) · scholar, toolsmith
- **T1 · H-021'** (pairing cores, `s = n - k(c)`) · E-059..E-062 · round 047 (maverick note: no cheap census left; reopen only for a rule for k(c) predicting a held-out core; H-021 header rewrite proposed for the round-048 ledger) · **dormant**, experimentalist, theorist
- **T2 · centre s(c) / k(33x) = 2x** · E-060, E-061, E-065 (lower bound only) · round 047 · **dormant**, theorist
- **T4 · H-020 rule table as a theorem** · Round 049 (E-157): hypotheses H1-H6 stated falsifiably; H1 holds for width <= 5 to length 11. Round 051 (E-163): H1 spot check at one length each for widths 6..8 (0 failures in 43 116), gaps at lengths 10-12 and widths 9..11; H6: anchored rules redundant only in the reduced walk, matter in a rules-only walk. H3 (completeness, F-052) weakest; H4 least explained. · round 051 · theorist
- **T6 · H-017 (mutation search; Euler signature)** · OPEN. Signature is Cartan-determined, a class invariant, blind to (cords, relations) (E-160); no search-based test exists beyond E-063 (n = 9, depth 6, no positive control). Reopen only with a search that sees relations. · round 050 · **dormant**, maverick
- **T9 · long runs** · n = 12/14 censuses of the 139 cores and Ladkani audit n = 9/10 sit in `OVERNIGHT.md` Menu 4 for the human. · **parked**

## Awaiting revision
(none)

## Requests between personas
- skeptic: hand-rebuild of the 13 class-1 E-145 steps; an independent definition of the reverse step for the E-154 B-half.
- toolsmith: decide the 8 parallel-arrow failing steps (matrix-valued arrow identification); quiver-level End(T) on E-155/E-158 paths; a control with minimum exactly 13; an orbit-based canonical form for the 12 high-cost keyless nodes; a relabelling-aware `meetingPoints` (E-162); a control through the no-key route; `merges.py 10 --depths 5 --witness` on a real link; reverse loss-by-depth and deep control (E-147).
- theorist: generation of K^b(proj) by the J = 0 tilting complex (is it automatic?); orbit data giving D = 0 and c_2 = 0 (E-143/E-146) (T5); the `33x` orbit size 3767 vs the n = 14 `444` orbit (maverick flag).
- maverick: S-1 n = 15 K0 = 5 sizing; Smith form of C+C^T for the F-010 pair vs E-152 (one line).

## Weakest claims the workshop currently relies on
The key-guard law (class 0 only); E-145 class-1 count (one script); H-015 citations from memory; the S-1 K-threshold law (two lengths); E-149 counts are step counts on capped BFS samples; E-153 is one n = 8 class, depth 8, slice 2 not independently re-run; 'in the class' for the 25 children rests on the J = 0 premise (generation now proved per loopless step, E-166), edge-tested at Cartan level (E-155, E-161, E-164 paths), at quiver level only for the 13 E-161 edges (E-165, label-preserving); the E-163 length check is one length per rule.

## Standing facts
- arxiv.org is blocked by the network policy (403); PDFs would need the human.
- Round 055 added E-169..E-171 (E-170 supersedes E-168; HH closed by a theorem). Round 054 added E-166..E-168; the human supplied the LaTeX of 1009.3370 and 2509.12983 (citations verified in round 055). Round 053 added E-164, E-165. Round 052 was a conference: no change to any H-, F-, R- status (ledger in proceedings). Round 051 added E-161..E-163. Round 050 added E-158..E-160. Round 045 added E-148, E-149; H-015 OPEN. Round 046 added E-150..E-152; round 047 added E-153, E-154. Nothing in F- or R-.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 054 | 050 (skeptic) |
| theorist | 054 | 055 (toolsmith) |
| skeptic | 053 | 055 (scholar) |
| scholar | 055 | 055 (maverick) |
| toolsmith | 055 | 054 (experimentalist) |
| maverick | 055 | 053 (toolsmith) |
