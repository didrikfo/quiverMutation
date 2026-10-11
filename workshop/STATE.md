# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round; under 120 lines. Every persona reads this first.

last_round: 58
next_round_kind: ordinary   (round 059; round 060 is a conference)

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* (mutation is admissible) and a *Coxeter-key guard* (`coxeterGuard`: child keeps the class's Coxeter key). Its docstring says the guard keeps a walk in one derived class. Round 045 found the guard is not a tilting test (E-147, E-151). Under the J = 0 premise all 25 failing n = 7 children are now joined to an LNA (E-157, E-160, E-163); the premise itself is the live question (T10, T5).

## Agenda (proposed, round 056 conference; ordinary rounds work from it until the human changes it). Ranked:
1. T10 premise: quiver-level End(T) is DONE (E-174) and has power on killed-line variants (E-177). Next: theorist proves End(T) = mutation algebra at every step (J = 0 or not) and what a parallel-arrow step needs; toolsmith: `compare2` dim assertion, wrong-algebra control for step 12 / class 2; Hom(T,T[-1]) on the n = 10 E-169 edges.
2. Power control for the J = 0 join test: n = 9 (E-175) and n = 10 (E-178) pairs done, weak, no J != 0 in the balls. Further depth is circular given the premise; NOT approved for OVERNIGHT. Next only if a pair or a control with a J != 0 step on the path is found.
3. T3/T8: maverick's Candidate C (S^a(A) ~ A[b]) on the one n = 10 group with Phi^18 = I; scholar reads 0911.5137 and 1310.1557 first. Level: idea.
4. Breadth: T6/H-017 (16 unplaced n = 11 LNAs) or T7 (why exactly 2 LNA-to-LNA mutations, n = 9, 10).
5. Housekeeping (round 060 ledger: close T1 and T2 as stated, proposals in E-176): citation fixes DONE (E-179); toolsmith docstring rewords and `canonicalKey` cap. S-1 n = 15 scan parked.
Breadth rule: one slot per round to a dormant thread or a closure note.

## Open threads
Each: question · state · last round worked · owner.

- **T10 · guard audit** · Does any recorded class merge, or any step of a key-guarded walk, depend on a J != 0 step? · (a) 16 of 80 978 (n = 7 c1) and 9 of 79 143 (c2) key-keeping steps fail `tiltingPlus` at parent depth 7-8 (E-151); F-041 n = 8 merges and guarded depth <= 4 LNA walks use no J != 0 step (E-150, E-156); n = 8 c2 depth-8: 0 of 104 629 fail (E-155). (b) 19 of the 25 failing children and all 25 parents are joined to an LNA by J = 0 `tiltingPlus` paths (E-157); an independent Hom(T,T[m]) test accepts all 324 edges of the 40 printed paths and rejects the 25 failing steps (E-161, Cartan-level, generation assumed). (c) Round 050: depth-7 ball joins child 12 and 3 of 5 misses (total <= 12, E-160); `canonicalKey` None on parallel-arrow bundles blinds the child side (docstring wrong). Round 053: E-166 all 33 edges of the 3 E-160 paths pass the Hom test (Cartan level); E-167 End(T) iso to the next algebra on 13/13 E-163 edges at quiver level (8 of 16 failing J != 0 steps decided, iso, so the check cannot see J; 8 undecided, parallel arrows); 25 of 25 = 23 by own path + 2 (c1 13, 15) by key equality. Generation still assumed. Round 051: c1 14/15 joined at total 13 (7 + 6, E-163; skeptic replayed with Hom test), so 25 of 25 joined under the premise; Hom-tested only for E-157 paths and this one. (d) Round 051 (E-164): `merges.py 10 --depths 3 4 5` found no link; F-037 (one of two n = 10 merges) 5 paths use no J != 0 step; group-A merge (05040330 -> 33460000) unreplayed, no path recorded; `meetingPoints` uses labelled `quiverKey`. Round 054: generation of K^b(proj) by T proved per step at loopless vertices, loop-free on 370 printed edges (E-168); group-A merge has 5 explicit 7-step J = 0 paths (E-169; E-033 §4 had recorded the length). Round 057 (E-174): quiver-level End(T) iso at all 25 failing steps (16 c1 incl. the 8 parallel-arrow, 9 c2), 370/370 n = 7 path edges, 35/35 n = 10 edges; Hom(T,T[-1]) != 0 at 25/25, so End(T) cannot see J; premise still rests on Hom(T,T[-1]) and E-168. Round 057 (E-175): n = 9 F-010 pair, no join to 6 + 6, weak. Round 058: E-177 equal-dims wrong-algebra control: End(T) comparison rejects 280/280 equal-dim killed-line variants (n = 7 c1, 7 steps), accepts 168/168 iso, but accepts a dropped relation 47/47 (relation set and dim are inputs); E-178 n = 10 pair 90000000 vs 50505000: no J = 0 join to 6 + 6 (a 7-step equivalent control joins at 4 + 4; no J != 0 step in the balls). Open: `compare2` dim assertion; step 12 and class 2 for the wrong-algebra control; Hom(T,T[-1]) on n = 10 edges; n = 8 depth 9 (not approved); docstring rewords (guard, `canonicalKey`). · round 058 · toolsmith, skeptic, maverick, theorist
- **T5 · H-015 / key guard** · H-015 is OPEN. Round 058 (E-179) narrowed T5 to one live item: End(T) = the repo's rewrite on J = 0 steps (proof, or Oppermann 1504.02617 reduction compared with steps 1-7; its LaTeX is not in `sources/`). Also open: E-147's 13 class-1 steps are not shown to lie in E-151's 16 (class 2: exact; class 1 rerun 12 of 16, E-147's 13 not reproduced); skeptic's hand-rebuild of c1 (8 of 16 parent quivers have parallel arrows). Moved out of T5: orbit data giving D = 0 / c_2 = 0 / s = -2 (theorist backlog, E-145/E-148), reverse-search loss of 10.3% (toolsmith backlog, E-149; effect on E-160/E-175 negatives is a guess). Literature citations fixed (E-171, E-179). The law "J != 0 leaves the key" is class 0 only. · round 058 · scholar, skeptic, theorist
- **T7 · H-010 (overlap reducible only at an end)** · No proof. Round 046 (E-152, E-153): run-of-three is necessary in every cell tested and sufficient at k = 3 in 5 of 5; L1 one-step locality holds n <= 10; k-step window, two bystanders, k > 3 untested; refines F-022. Next: cap-cut k = 3 rows, two bystanders, explain the 2-per-LNA count. · round 046 · theorist
- **S-1 · vertex deletion transports classes?** (suggested question) · Round 055 (E-173): K0 = 5 at n = 15 confirmed at key level on a one-row witness; items 1-3 of the STEERING question still open; full n = 15 scan sized (~33 min on 4 cores, 24 shards), not run. K-threshold law (E-146): K >= 4 one n = 11 key class, K >= 3 two; first failures K0 = 3, 4, 5 at n = 11, 13, 15 (two lengths beyond 11). Open: n = 15 K0 = 5 at class level (`--plan` first); 136 + 66 split; non-lone cores. · round 055 · maverick, toolsmith
- **T3/T8 · orbit-plus-mirror vs key, catalogue** · Done at the catalogue level (E-066). Key-coarser cores: 7 (n = 10), 9 even n, 10 odd n; 8 of 10 odd-n and all even-n are two mirror-closed orbits P, Q sharing a key (E-079, E-082); 5046/5056 outside; why no shift by 1 unproved, no invariant found. Round 046 scholar note: dormant, not closed. · round 054 (maverick, revision awaited) · scholar, toolsmith
- **T1 · H-021'** (pairing cores, `s = n - k(c)`) · E-061..E-064 · round 047 (maverick note: no cheap census left; reopen only for a rule for k(c) predicting a held-out core; H-021 header rewrite proposed for the round-048 ledger) · round 057 (theorist note, E-176): proposed closing at the round-060 ledger: no linear letter rule for k(c) (17 values, 11 independent of 33x); drift-based rule untried · experimentalist, theorist
- **T2 · centre s(c) / k(33x) = 2x** · E-062, E-063, E-067 (lower bound only) · round 057 (E-176): 333@0 orbit = 444 orbit (n = 13..16); at n = 14 only ~20 of 3767 rows are 33y, so E-067 concerns which 33y rows share an orbit; proposed closing as stated at the round-060 ledger; the large 444 orbit belongs to H-020 / T3 · theorist
- **T4 · H-020 rule table as a theorem** · Round 049 (E-159): hypotheses H1-H6 stated falsifiably; H1 holds for width <= 5 to length 11. Round 051 (E-165): H1 spot check at one length each for widths 6..8 (0 failures in 43 116), gaps at lengths 10-12 and widths 9..11; H6: anchored rules redundant only in the reduced walk, matter in a rules-only walk. H3 (completeness, F-052) weakest; H4 least explained. · round 051 · theorist
- **T6 · H-017 (mutation search; Euler signature)** · OPEN. Signature is Cartan-determined, a class invariant, blind to (cords, relations) (E-162); no search-based test exists beyond E-065 (n = 9, depth 6, no positive control). Reopen only with a search that sees relations. · round 050 · **dormant**, maverick
- **T9 · long runs** · n = 12/14 censuses of the 139 cores and Ladkani audit n = 9/10 sit in `OVERNIGHT.md` Menu 4 for the human. · **parked**

## Awaiting revision
(none)

## Requests between personas
- skeptic: wrong-algebra control for step 12 (two parallel pairs) and class 2 (E-177); hand-rebuild of the 13 class-1 E-147 steps; an independent definition of the reverse step for the E-156 B-half.
- toolsmith: `compare2`: compute dim K Q/I' from `crels` and assert it equals dim End(T) (E-177); quiver-level End(T) on E-157/E-160 paths; a control with minimum exactly 13; an orbit-based canonical form for the 12 high-cost keyless nodes; a relabelling-aware `meetingPoints` (E-164); a control through the no-key route; `merges.py 10 --depths 5 --witness` on a real link; reverse loss-by-depth and deep control (E-149).
- theorist: generation of K^b(proj) by the J = 0 tilting complex (is it automatic?); orbit data giving D = 0 and c_2 = 0 (E-145/E-148) (T5); the `33x` orbit size 3767 vs the n = 14 `444` orbit (maverick flag).
- maverick: a long-path equivalent control for the 90000000 side (E-178); S-1 n = 15 K0 = 5 sizing; Smith form of C+C^T for the F-010 pair vs E-154 (one line).
- scholar: fetch the 1504.02617 LaTeX and compare Oppermann's reduction with steps 1-7 on J = 0 steps (T5).

## Weakest claims the workshop currently relies on
The key-guard law (class 0 only); E-147 class-1 count (one script); H-015 citations from memory; the S-1 K-threshold law (two lengths); E-151 counts are step counts on capped BFS samples; E-155 is one n = 8 class, depth 8, slice 2 not independently re-run; 'in the class' for the 25 children rests on the J = 0 premise (generation now proved per loopless step, E-168), edge-tested at Cartan level (E-157, E-163, E-166 paths), at quiver level only for the 13 E-163 edges (E-167, label-preserving); the E-165 length check is one length per rule.

## Standing facts
- Round 058 added E-177..E-179 (merge of origin/main already up to date).
- Round 057 added E-174..E-176 (merge of origin/main conflicted in rounds/052/proceedings.md; skipped). Round 056 was a conference: ledger changed no H-, F- or R- status.
- arxiv.org is blocked by the network policy (403); PDFs would need the human.
- Round 055 added E-171..E-173 (E-172 supersedes E-170; HH closed by a theorem). Round 054 added E-168..E-170; the human supplied the LaTeX of 1009.3370 and 2509.12983 (citations verified in round 055). Round 053 added E-166, E-167. Round 052 was a conference: no change to any H-, F-, R- status (ledger in proceedings). Round 051 added E-163..E-165. Round 050 added E-160..E-162. Round 045 added E-150, E-151; H-015 OPEN. Round 046 added E-152..E-154; round 047 added E-155, E-156. Nothing in F- or R-.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 057 | 050 (skeptic) |
| theorist | 057 | 058 (skeptic) |
| skeptic | 058 | 058 (maverick) |
| scholar | 058 | 057 (theorist) |
| toolsmith | 057 | 058 (scholar) |
| maverick | 058 | 053 (toolsmith) |
