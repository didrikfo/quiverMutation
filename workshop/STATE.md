# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 44
next_round_kind: ordinary

## Open threads

Round 006 (ordinary: toolsmith, theorist, scholar) recorded E-064..E-066; threads T3/T8, T2/T4, T5 are updated below. Round 005 was a conference; its ranked agenda (below, `rounds/005/proceedings.md`) is **proposed, round 005**, pending the human's approval in STEERING.md; ordinary rounds work from it. Priority order: T3 (with T8) first, then T2/T4, then T6's control, then T5 reading; T1 and T9 wait on overnight. Round 004 recorded E-061..E-063. Each line: id · question · suited to · status.

- **T1** · H-021: restated without the mirror clause (H-021': a pairing core has `s = n - k(c)`, `k`, `d` independent of n; covers only cores that pair). The 7 cores of E-059 pair at even n, fail at odd; the 12 cores of E-060 keep `k`, `d` at 15 and 16 and pair at 13, 14, 15; at 16 `46 3355 3445` lose the fit through an unmerged equal-size middle pair (size 20300, as `4056` in E-058), so parity is not the whole story (E-062). Open: is 20300 one shared orbit; the 3 at n = 17, 18; the other 127 cores (census at 12, 14 in Menu 4, overnight); why the outside block folds. · experimentalist, theorist · open
- **T2** · (proposed rank 2, with T4: theorist derives `k(33x) = 2x` from the rule table; experimentalist gives `34x`, `44x`, `45x`) Centre `s(c) = n - k(c)`. Interior blocks: `s` = first + last outside offset (13/13, E-060, E-061); end-touching failures `4045 3556 4506 4556` pick the smaller consistent centre for three and the larger for `4506` (no min/max rule); all-outside `4046 5046 5056` maybe a period-2 translation, not a reflection (unverified). `k(33x) = 2x`, `d = x - 3` for x = 3..7 (E-061). Open: why `k = 2x`; `k(c)` for other prefixes (`34x`, `44x`, `45x`); a null test for `|R| <= 4` fits; T4 may explain the top `x - 3` singleton offsets. · theorist, skeptic · open
- **T3** · (proposed rank 1, with T8: toolsmith commits the whole-catalogue orbit task, experimentalist runs 14..16, skeptic referees; first check: is the 20300 pair of `4056` the same orbit as that of the E-060 cores at 16) `3346` never pairs (proved by key, E-054). `4056` at n = 16: `{1,2}` are two mirror-image orbits, key pairs them (E-058), so "orbit = key class" is false there; compare orbit-plus-mirror against key over the `--max-word 4` catalogue at 14..16. Open: other cores that split; the other four H-020 failures (ledgers not in repo). · experimentalist, skeptic · open
- **T4** · H-020 as a theorem: the rule table acts the same at every interior position (F-051). State and prove. · theorist · open
- **T5** · (proposed rank 4: scholar reads Aihara-Iyama Thm 2.32, arXiv:1009.3370, and CHZ, arXiv:2509.12983, against E-032 step 7) H-015: Ladkani 2.3(c) agrees with the gate everywhere tested, incl. non-monomial parents at n = 5..7 (E-055, E-057); only gate-admitted rejection is E-032 step 7. `isTilting` not promoted (STEERING q3). Open: the audit at n = 9/10 (overnight, Menu 4). · scholar, toolsmith · waiting on overnight
- **T6** · (proposed rank 3: toolsmith builds the positive control first; skeptic: null test for `|R| <= 4` fits) H-017: not refuted at n = 9 to depth 6; the Coxeter polynomial and Euler form do not predict (cords, relations); the Euler signature `pos <= n-2` separates outside-every-quipu-class for n = 8..11 (E-063). Open: n = 12 signature (about 15 min in shards); an LNA in the n = 13 class of `P^(1,0,3,0,1)_(1,1,1,1)`; a positive control for the mutation search; a tree eigenvalue argument for `relations >= cords + 1`; depth 7 at n = 9 (overnight, pending question). · maverick, toolsmith, theorist · open
- **T7** · H-010, overlap reducible only at an end: a proof from step 7 of arXiv:2112.08129. Do not run `probe.py --steps 7`. · theorist · open
- **T8** · Tooling: `batch.py orbits N` exists (round 003; `45`, `344` pinned). Open: run it over the whole catalogue (`--jobs 4`), commit the fit/slide step as a task, orbit-plus-mirror prefilter (not key). · toolsmith · open
- **T9** · H-019, H-013: long runs. Overnight waiting on the human: n = 12 and n = 14 censuses of the 139 cores (Menu 4; `workshop/rounds/002/experimentalist_census.py` and `experimentalist_fit.py`); Ladkani audit n = 9/10. · any · parked

Next round (007) is ordinary; 008 is the next conference. Overnight run added in round 005: H-017 depth 7 at n = 9 (Menu 4). Still waiting on the human: n = 12/14 censuses, n = 14 census, Ladkani audit n = 9/10.

## Round 044 updates (supersede everything below where they differ)

Round 044 was a conference (all six personas; `rounds/044/`). **Round 045 is ordinary; 048 is the next conference.** Ranked agenda **(proposed, round 044)**, replacing the round-040 agenda, pending the human; ordinary rounds work from it:
1. T5 key guard: skeptic hand-rebuilds the 13 class-1 E-145 steps (n = 7); experimentalist guard-off n = 7 census (`--plan` first) with per-step `tiltingPlus` + Cartan congruence, saves the n = 8 J != 0 steps.
2. T5 theory: theorist, orbit data giving D = 0 (s = 10 in c2), why c_2 = 0, the s = -2 family; which property the meet needs (tilting, silting, key).
3. T5 tools: toolsmith, can a key-preserving walk step fail `tiltingPlus` (tally); loss-by-depth for the reverse search, deep control. No overnight.
4. Literature: scholar retries arXiv 1009.3370 / 2509.12983 (standing permission); else park.
5. S-1 n = 15 K0 = 5 sizing (maverick, toolsmith); min(h, K) vs core key.
Weakest claims named: the key-guard law (class 0 only; guard is also the BFS filter); E-145 class-1 count (one script); H-015 citations from memory; S-1 law (two lengths beyond n = 11).
Requests added: skeptic <- experimentalist, theorist, scholar, toolsmith, maverick: independent rebuild of class-1 steps. experimentalist <- toolsmith, skeptic: guard-off n = 7 census, saved n = 8 steps. theorist <- scholar, experimentalist: orbit data for D = 0.

## Round 043 updates (supersede everything below where they differ)

Round 043 (ordinary: skeptic, experimentalist, toolsmith; referees theorist, skeptic, maverick; all minor revision, accepted with qualifications) recorded E-145..E-147. **Round 044 is a conference** (conference_every = 4). The round-040 agenda stands, pending the human.
- **T5 / key guard (E-145):** at n = 7 classes 1 and 2 there are gate-admitted J != 0 steps that keep the class key (13 of 67, 9 of 64 distinct; walk descendants, tiltingPlus False, Cartan incongruent), so E-140's "none keeps the key" and the E-138/E-141 law hold for class 0 only; the key guard is no evidence for H-015 off J = 0 steps. Orbit relation e_i = F^s e_w fails on most H1/H2 steps in c1, c2. Open: which orbit data give D = 0 (s = 10 in c2); rebuild the class-1 steps by hand; E-140's command at 500 s; whether D = 0 children are derived-equivalent (none claimed).
- **T5 / n = 8 (E-146):** no n = 8 J != 0 step has the E-143 shape (|out v| = 2 or |supp J| = 2, so s, c_2 undefined); Q's lowest term is x^3 for some class-1 steps; E-143's "orbit relation never absent" holds only at n = 6, 7 c0. Shallow sample (cap, depth 7-9). Open: reduction for |out v| = 2; save the steps; guard-off depth 6 (>10 min, overnight candidate); n = 9 shape.
- **T5 / reverse meet (E-147):** reverse positive control at depth 2-4 passes (12/12); 154 of 1500 edges (10.3%) are lost in reverse, not to a filter (same-vertex opposite step lands on a different same-key algebra). Control is shallow; loss by depth not tallied. Open: loss-by-depth tally, deep control, reason for the loss; the 3 h reverse job stays out of OVERNIGHT.md.
- Agenda item 3 (literature) parked: no PDFs. S-1 n = 15 K0 = 5 sizing not worked.
- Requests added: theorist: orbit data giving D = 0 (B = 0 = 1 - c_s - c_{-s}); E-143 reduction for |out v| = 2. skeptic: explain the 154 lost edges (module condition). experimentalist: E-140 command at 500 s; n = 8 steps saved; toolsmith: loss by depth, deep reverse control.

## Round 042 updates (supersede everything below where they differ)

Round 042 (ordinary: scholar, theorist, maverick; referees skeptic, experimentalist, toolsmith; all minor revision) recorded E-143, E-144. **Round 043 is ordinary; 044 is the next conference.** The round-040 agenda stands, pending the human.
- **T5 / Q(x) (E-143):** under shape hypotheses H1, H2 on the v-row, Q(x) = x(adj S_ii - adj S_wi - adj S_iw) (proved); with F e_w = e_i the x^2 coefficient is 1 iff c_2 = 0 (proved). Observed only: H1/H2 from the gate; e_i = F^s e_w (s = 1 or -2) in every step with H1, H2; c_2 = 0; the 50 off-shape steps. Open: module-level derivation (eAe at n = 6); an H1/H2 step with c_2 != 0 (would be a key-preserving candidate); n = 8 tally.
- **T5 / literature (item 3):** arXiv is blocked by the proxy (403) in rounds 006 and 042; AI 2.31/2.32 and CHZ 3.6 stay UNVERIFIED. By hand, "monomial" is needed for path-wise Cor 3.6 (already flagged in literature/2509.12983 and E-066/E-122). Parked unless the PDFs are supplied.
- **S-1 (E-144):** n = 12 confirms the K-threshold law (K >= 4 one n = 11 key class, K >= 3 two); single-3 key table: first failures K0 = 3, 4, 5 at n = 11, 13, 15. Open: n = 15 K0 = 5 at class level (enumerate from the lone-3 orbit, size first); E-115-type label to confirm the 136 + 66 split is two derived classes; non-lone cores; why the split is by min(h, K).
- Requests added: theorist/skeptic: derive H2 and the orbit relation; skeptic: an H1/H2 step with c_2 != 0. experimentalist: n = 8 tally for E-143. toolsmith/experimentalist: n = 15 K0 = 5 sizing.

## Round 041 updates (supersede everything below where they differ)

Round 041 (ordinary: experimentalist, theorist, toolsmith; referees skeptic, scholar, maverick; all minor revision, accepted with qualifications) recorded E-140..E-142. **Round 042 is ordinary; 044 is the next conference.** The round-040 agenda (below) stands, pending the human.
- **T5 / key guard (E-140):** 278 gate-admitted J != 0 steps at n = 6, 7, 8 (classes 0-2; 5 of 9 cells non-vacuous): none keeps the class key. Random non-LNA parents do (55 of 1264). No LNA-parent analogue found. Open: per-cell tally of child keys; n = 8 c1/c2 deeper; why c2 has no J != 0 steps.
- **T5 / theory (E-141):** "J != 0 implies the key moves" is false for general parents (n = 4 example, non-LNA key); on LNA walks Q(x) = det(xC_B+C_B^T) - det(xC'+C'^T) has lowest term x^2 (coefficient 1) at n = 6, 7. No proof; the x = 0, infinity and trace routes are closed. Open: why Q_2 = 1; n = 7 C_B = C' + H check; is the n = 4 pair tilting-equivalent.
- **T5 / meet (E-142):** `toolsmith_n6meet.py` (round 041 copy) has `--tilting-only`, `--control`, `--revcontrol`, `--reverse`, closure flags. 16/16 hits share 0 keys, no BFS closes: bounded miss. Reverse search loses about 12% of edges and has no positive control. Open: a reverse control at depth >= 2; why edges are lost; the 3 h job (hits 0, 4, 13) not yet in OVERNIGHT.md.
- **Agenda items 3 (scholar, literature) and 4 (S-1 n = 12) were not worked this round.** Next round: scholar on item 3, maverick/experimentalist on item 4, theorist or skeptic on Q_2.
- Requests added: skeptic: independent build of the n = 4 child and search for an LNA-key J != 0 child at n <= 6 using Q(x). toolsmith: reverse positive control. experimentalist: per-cell tally of child LNA keys.

## Round 040 updates (supersede everything below where they differ)

Round 040 was a conference (all six personas; `rounds/040/`). **Round 041 is ordinary; 044 is the next conference.** Ranked agenda **(proposed, round 040)**, replacing the round-036 agenda, pending the human's approval in STEERING.md; ordinary rounds work from it:
1. T5: does any gate-admitted J_i != 0 step keep the LNA key (key guard off; n = 6, 7, 8, classes 0, 1, 2; positive control; confirm child key is computed on the child)? experimentalist (`--plan` first), skeptic, theorist (R(x) = 1 + t for j = e_i).
2. T5: `--tilting-only` meet in `toolsmith_n6meet.py` with positive control and closure flag; the other 13 hits; reverse search from the hits (J = 0 vertices). toolsmith, skeptic.
3. T5: AI Thm 2.31/2.32 and Ladkani 2.3(c) against `tiltingPlus`; the 2.31 citation is from memory until read. scholar, skeptic.
4. S-1: n = 12 key class after `--plan`; fix `maverick_single.py` label block; null from a core-length statistic. experimentalist, maverick, toolsmith, skeptic.
Weakest claims named: E-134 membership (unsupported); E-138 refusal pattern (guard is the BFS filter); d_i = 2 (capped, E-131 counterexamples); H-015 as the tilting test; the S-1 K-threshold law (two lengths).
Rounds 036-038 update sections dropped as superseded (see their proceedings).

## Round 039 updates (superseded by 040 where they differ)

Round 039 (ordinary: skeptic, theorist, maverick; referees experimentalist, experimentalist, skeptic; all minor revision, accepted with qualifications) recorded E-137..E-139. **Round 040 is a conference** (conference_every = 4). The round-036 agenda stands (proposed, pending the human).

- **T5 (E-137):** E-134's meeting needs a gate-admitted, key-preserving, non-tilting step out of the hit (kernel dim 2); tilting-only BFS (LNA side 27 518 in 150 s, not closed; hits 323-501) shares nothing. The meeting is unsupported, not refuted. Open: positive control and closure flag for the tilting-only run; replay of the other 13 hits; whether the Coxeter guard as such admits the step (`coxkey-pres = True`); reverse search from the hits using J = 0 vertices only.
- **T5 (E-138):** all 192 J != 0 steps at n = 8 c0 (and 123 at n = 6, 7) leave the class key; out(i) = 3 is a sample property, the gate allows out(i) = 2 (hand algebra). Open: c1, c2, n = 9; a gate-admitted J != 0 step whose child keeps the key (would break the law); prove R(x) != 1 + t when j = e_i; machine-check row 7798's kernel element.
- **S-1 (E-139):** n = 13 failure confirmed at the orbit level (4349 + 674). Open: n = 12 (K >= 4 holds, K >= 3 fails; `--plan` first) and n = 15; which backward moves join (4,5) and (5,4); the {h, K} dependence from H-020.
- Requests added: skeptic: a J != 0 step whose child passes the key guard (n <= 7). toolsmith: `--tilting-only` in `toolsmith_n6meet.py` with a positive control; reverse-direction search from the hits. experimentalist: n = 8 c1, c2 for the J != 0 / key-guard column; n = 12 key class for S-1. theorist: R(x) = 1 + t for j = e_i.

## Round 006 updates (supersede the thread text above where they differ)

- **T3/T8 done at the catalogue level (E-064):** orbit-plus-mirror refines the key in all 139 cores at n = 10, 12..16; equals it in 129-132; the rest are key-coarser cores where the key merges the two parity classes (9 at even n, 10 at odd n). The 20300 pairs of `4056` and the 3 cores of E-062 are one orbit and its mirror, so the "unmerged middle pair" is resolved by the mirror join. Open: why those 9/10 cores are parity classes, and whether they are the 7 of E-059; `--max-word 5` at n = 14 (needs `--plan`); a prefilter by key difference (true in 139 x 6 cases, not proved); `toolsmith_orbitclass.py` has no `--budget-hours`.
- **T2/T4 (E-065):** `k(33x) = 2x`, `d = x - 3` from the drift of the double mutation plus seed `333` (lower bound derived; end link and upper bound computed only). Fails for `44x` (one orbit with `333@0`), `34x` only `345` explained, `45x` no drift. Open: the full-orbit size comparison for 33x; the end link as a proof; n = 17 for 34x/45x; n = 18 rule-table claim has no output file; x = 9 excluded by a script cap.
- **T5 (E-066):** step 7 is rejected at a commutativity element; AI 2.32(b), Ladkani 2.3(c), `tiltingPlus` appear to be one map (no derivation); CHZ Cor 3.6 path-wise wording may need "monomial" (UNVERIFIED, PDF unread, arXiv blocked by the proxy this round). Open: derive the one-map identity; build a smaller instance (commutative square into a vertex with one outgoing arrow) at n = 5..7 and run `tiltingPlus`; compute End of the two-term complex directly for an independent check.
- **T7 (H-010):** step 7 is a concrete non-monomial case to test the "overlap reducible only at an end" proof on.

Requests added: toolsmith to theorist: why the key-coarser cores are the parity classes. scholar to theorist: the n = 5..7 instance for T5. scholar/toolsmith: read the CHZ PDF where the network allows.

## Awaiting revision

(none)

## Requests between personas

- theorist (to experimentalist): `k(c)` for `34x`, `44x`, `45x` at n = 14..17; the 4046-type cores at 14 (period-2 translation?).
- maverick (to toolsmith): a positive control for the H-017 mutation search.
- maverick (to experimentalist): n = 12 signature census in shards.
- experimentalist (to theorist): what separates unmerged `{4}{5}` of `3355` at 16 from merged `{3,4}` at 14.
- experimentalist: orbit-plus-mirror vs key classes over the catalogue at 14-16 (skeptic).
- theorist: at n = 16 `4056` has a self-dual pair `{0,3}` and a mirror pair `{1,2}`: is the split a function of parity or of position relative to the middle? (skeptic)
- theorist: the centre `s(c)` formula (round 001 request, still open).
- toolsmith: committed census task; key prefilter with the E-058 caveat.

## Rota

<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 043 | 042 (theorist) |
<!-- 044 was a conference: all six wrote position statements -->
<!-- 036 was a conference: all six wrote position statements -->
<!-- 012 was a conference: all six wrote position statements -->
<!-- 008 was a conference: all six wrote position statements -->
<!-- 005 was a conference: all six wrote position statements; nobody refereed -->
| theorist | 042 | 043 (skeptic) |
| skeptic | 043 | 043 (experimentalist) |
| scholar | 042 | 041 (theorist) |
| toolsmith | 043 | 042 (maverick) |
| maverick | 042 | 043 (toolsmith) |
