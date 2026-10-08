# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round; under 120 lines. Every persona reads this first.

last_round: 45
next_round_kind: ordinary   (round 046 ordinary, 047 ordinary, **048 is the next conference**)

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

- **T10 · guard audit** · Does any recorded class merge, or any step of a key-guarded walk, depend on a J != 0 step; what does a `tiltingPlus` guard cost? · (b) the 10 F-041 n = 8 merges and guarded depth <= 4 LNA walks at n = 6..8 use no J != 0 step (E-148); (a) the added check costs 8-11% per step, and 16 of 80 978 (n = 7 c1) and 9 of 79 143 (c2) key-keeping steps fail `tiltingPlus` at parent depth 7-8 (E-149). Open: (i) do those children lie outside the derived class (independent invariant, not the Cartan congruence of the step's own map); (ii) deep merges (distance >= 5: E-094 n = 8 c2 depth-8 walk, pipeline depth 5-6) with per-edge `tiltingPlus` tally (overnight candidate, not in OVERNIGHT.md); (iii) `merges.py` stores no witness paths: store them; (iv) reword the `mutationSearchDepthFirst` docstring after (i). · round 045 · skeptic (i), toolsmith (iii), experimentalist (ii)
- **T5 · H-015 / key guard** · H-015 is OPEN (ledger, round 045). The law "J != 0 leaves the key" is class 0 only (E-138, E-140, E-145). Open: skeptic's independent hand-rebuild of the 13 class-1 E-145 steps (n = 7); orbit data giving D = 0 (s = 10 in c2), why c_2 = 0, the s = -2 family (theorist; E-143, E-146); the reverse-search loss of 10.3% of edges (E-147, loss by depth, deep control). Literature (Aihara-Iyama 2.31/2.32, CHZ 3.6) unread, cited from memory, UNVERIFIED. · round 043 (T10 took 045) · skeptic, theorist, toolsmith
- **T7 · H-010 (overlap reducible only at an end)** · *Awaiting revision.* No proof found; "run-of-three" criterion (bystander shares >= 2 arrows) tested to k = 2, not k = 3; mostly restates F-022. Closure not agreed. · round 045 (major revision) · theorist
- **S-1 · vertex deletion transports classes?** (suggested question) · K-threshold law (E-144): K >= 4 one n = 11 key class, K >= 3 two; first failures K0 = 3, 4, 5 at n = 11, 13, 15 (two lengths beyond 11). Open: n = 15 K0 = 5 at class level (`--plan` first); 136 + 66 split; non-lone cores. · round 042 · maverick, toolsmith
- **T3/T8 · orbit-plus-mirror vs key, catalogue** · Done at the catalogue level (E-064): orbit-plus-mirror refines the key for all 139 cores at n = 10, 12..16. Open: why the 9/10 key-coarser cores are parity classes; `--max-word 5` at n = 14. · round 006 · **dormant**, toolsmith
- **T1 · H-021'** (pairing cores, `s = n - k(c)`) · E-059..E-062 · round 021 · **dormant**, experimentalist, theorist
- **T2 · centre s(c) / k(33x) = 2x** · E-060, E-061, E-065 · round 021 · **dormant**, theorist
- **T4 · H-020 rule table as a theorem** · F-051; open: state and prove. · round 037 · **dormant**, theorist
- **T6 · H-017 (mutation search; Euler signature)** · not refuted at n = 9 to depth 6 (E-063); needs a positive control and n = 12 signature. · round 023 · **dormant**, maverick, toolsmith
- **T9 · long runs** · n = 12/14 censuses of the 139 cores and Ladkani audit n = 9/10 sit in `OVERNIGHT.md` Menu 4 for the human. · **parked**

## Awaiting revision
- theorist: `rounds/045/theorist.md`, answer `rounds/045/theorist.review.md` (items 1-7: 21 of 21 not 22; retitle without "Lemma"/"holds at 3 mutations"; cite F-022; run `(8:3)(9:3)(10:4)` at k = 3; mark locality as conjecture or run L1 for n <= 9; include the left bystander `(7:3)(8:3)(9:3)`; status recommendation matching the evidence).

## Requests between personas
- skeptic: independent test whether the E-149 / E-145 key-keeping failing children are outside the derived class (a second invariant); hand-rebuild of the 13 class-1 E-145 steps; check E-084 class-index agreement.
- toolsmith: store witness paths in `merges.py`; reverse loss-by-depth and deep control; promote `tiltingPlus` only if the human agrees (Questions, round 045).
- experimentalist: deep replay of E-094's n = 8 c2 depth-8 walk with per-edge `tiltingPlus` tally (size with `--plan`; overnight if over 10 min).
- theorist: orbit data giving D = 0 and c_2 = 0 (E-143/E-146), after the T7 revision.
- maverick: S-1 n = 15 K0 = 5 sizing.

## Weakest claims the workshop currently relies on
The key-guard law (class 0 only); E-145 class-1 count (one script); H-015 citations from memory; the S-1 K-threshold law (two lengths); E-149 counts are step counts on capped BFS samples.

## Standing facts
- arxiv.org is blocked by the network policy (403); PDFs would need the human.
- Round 045 added E-148, E-149; H-015 status changed to OPEN. Nothing in F- or R-.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 045 | 043 (skeptic) |
| theorist | 045 | 043 (skeptic) |
| skeptic | 043 | 045 (experimentalist, toolsmith) |
| scholar | 042 | 041 (theorist) |
| toolsmith | 045 | 043 (maverick) |
| maverick | 042 | 045 (theorist) |
