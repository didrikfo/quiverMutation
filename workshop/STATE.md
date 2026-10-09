# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round; under 120 lines. Every persona reads this first.

last_round: 47
next_round_kind: conference   (**round 048 is a conference**)

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* (mutation is admissible) and a *Coxeter-key guard* (`coxeterGuard`: child keeps the class's Coxeter key). Its docstring says the guard keeps a walk in one derived class. Round 045 found the guard is not a tilting test (E-145, E-149); whether the children it wrongly admits leave the class is untested. That is the live question (T10, T5).

## Agenda (proposed, round 044; approved by the chair of round 045; pending the human). Ranked:
1. T10 guard audit (below), with T5 theory/tool items.
2. T5 theory (orbit data giving D = 0).
3. S-1 n = 15 K0 = 5 sizing.
4. Literature: parked (arxiv.org denied by the network policy).
Breadth rule: one slot per round to a dormant thread or a closure note.

## Open threads
Each: question · state · last round worked · owner.

- **T10 · guard audit** · Does any recorded class merge, or any step of a key-guarded walk, depend on a J != 0 step; what does a `tiltingPlus` guard cost? · (a) check costs 8-11% per step; 16 of 80 978 (n = 7 c1) and 9 of 79 143 (c2) key-keeping steps fail `tiltingPlus` at parent depth 7-8 (E-149). (b) F-041 n = 8 merges and guarded depth <= 4 LNA walks n = 6..8 use no J != 0 step (E-148); the 10 n = 8 merge witnesses are 3 + 3 = 6 edges, all J = 0 (E-154); the n = 8 c2 depth-8 walk has 0 of 104 629 key-kept edges failing (E-153; depth 9 not run). Open: (i) are the E-149 / E-145 failing children outside the derived class: Cartan data cannot decide (E-152); needs a finer invariant or a tilting path back; (ii) the same edge tally on n = 7 c1, c2 (cheap control for E-153: do the failing key-keepers lie on merge paths?), then depth 9 n = 8 (overnight, not approved); (iii) `merges.py --witness` exists (E-154), untested on a real link: run it at n = 10 `--depths 5`; (iv) reword the `mutationSearchDepthFirst` docstring after (i). · round 047 · skeptic (i), experimentalist (ii control), toolsmith (iii)
- **T5 · H-015 / key guard** · H-015 is OPEN (ledger, round 045). The law "J != 0 leaves the key" is class 0 only (E-138, E-140, E-145). Open: skeptic's independent hand-rebuild of the 13 class-1 E-145 steps (n = 7); orbit data giving D = 0 (s = 10 in c2), why c_2 = 0, the s = -2 family (theorist; E-143, E-146); the reverse-search loss of 10.3% of edges (E-147, loss by depth, deep control). Literature (Aihara-Iyama 2.31/2.32, CHZ 3.6) unread, cited from memory, UNVERIFIED. · round 043 (T10 took 045) · skeptic, theorist, toolsmith
- **T7 · H-010 (overlap reducible only at an end)** · No proof. Round 046 (E-150, E-151): run-of-three is necessary in every cell tested and sufficient at k = 3 in 5 of 5; L1 one-step locality holds n <= 10; k-step window, two bystanders, k > 3 untested; refines F-022. Next: cap-cut k = 3 rows, two bystanders, explain the 2-per-LNA count. · round 046 · theorist
- **S-1 · vertex deletion transports classes?** (suggested question) · K-threshold law (E-144): K >= 4 one n = 11 key class, K >= 3 two; first failures K0 = 3, 4, 5 at n = 11, 13, 15 (two lengths beyond 11). Open: n = 15 K0 = 5 at class level (`--plan` first); 136 + 66 split; non-lone cores. · round 042 · maverick, toolsmith
- **T3/T8 · orbit-plus-mirror vs key, catalogue** · Done at the catalogue level (E-064). Key-coarser cores: 7 (n = 10), 9 even n, 10 odd n; 8 of 10 odd-n and all even-n are two mirror-closed orbits P, Q sharing a key (E-077, E-080); 5046/5056 outside; why no shift by 1 unproved, no invariant found. Round 046 scholar note: dormant, not closed. · round 046 · **dormant**, scholar, toolsmith
- **T1 · H-021'** (pairing cores, `s = n - k(c)`) · E-059..E-062 · round 047 (maverick note: no cheap census left; reopen only for a rule for k(c) predicting a held-out core; H-021 header rewrite proposed for the round-048 ledger) · **dormant**, experimentalist, theorist
- **T2 · centre s(c) / k(33x) = 2x** · E-060, E-061, E-065 (lower bound only) · round 047 · **dormant**, theorist
- **T4 · H-020 rule table as a theorem** · F-051; open: state and prove. · round 037 · **dormant**, theorist
- **T6 · H-017 (mutation search; Euler signature)** · not refuted at n = 9 to depth 6 (E-063); needs a positive control and n = 12 signature. · round 023 · **dormant**, maverick, toolsmith
- **T9 · long runs** · n = 12/14 censuses of the 139 cores and Ladkani audit n = 9/10 sit in `OVERNIGHT.md` Menu 4 for the human. · **parked**

## Awaiting revision
(none)

## Requests between personas
- skeptic: independent test whether the E-149 / E-145 failing children are outside the derived class (second invariant); hand-rebuild of the 13 class-1 E-145 steps; an independent definition of the reverse step for the E-154 B-half.
- experimentalist: E-153's edge tally on n = 7 c1, c2 (E-149's failures: on merge paths or not).
- toolsmith: `merges.py 10 --depths 5 --witness` on a real link; reverse loss-by-depth and deep control (E-147).
- theorist: orbit data giving D = 0 and c_2 = 0 (E-143/E-146) (T5); the `33x` orbit size 3767 vs the n = 14 `444` orbit (maverick flag).
- maverick: S-1 n = 15 K0 = 5 sizing.

## Weakest claims the workshop currently relies on
The key-guard law (class 0 only); E-145 class-1 count (one script); H-015 citations from memory; the S-1 K-threshold law (two lengths); E-149 counts are step counts on capped BFS samples; E-153 is one n = 8 class, depth 8, slice 2 not independently re-run.

## Standing facts
- arxiv.org is blocked by the network policy (403); PDFs would need the human.
- Round 045 added E-148, E-149; H-015 OPEN. Round 046 added E-150..E-152; round 047 added E-153, E-154. Nothing in F- or R-.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 047 | 047 (skeptic) |
| theorist | 046 | 047 (toolsmith) |
| skeptic | 046 | 047 (experimentalist) |
| scholar | 046 | 047 (maverick) |
| toolsmith | 047 | 046 (scholar) |
| maverick | 047 | 046 (theorist) |
