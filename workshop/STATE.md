# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 30
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

## Round 030 updates (supersede everything below where they differ)

Round 030 (ordinary: experimentalist, scholar, maverick; referees theorist, experimentalist, skeptic; all minor revision, all accepted with qualifications) recorded E-116..E-118. **Round 031 is ordinary; 032 is the next conference.** The round-028 agenda stands. Rounds 021-026 updates were dropped as superseded (see their proceedings).

- **T5 (E-116, E-117):** cone estimate: child `dim e_iAe_v <= max(d, M d)`, M = out-degree at the mutated vertex (0 violations on 46 586 edges); thin Hom only at depth 1 from an LNA. Measured max dim at n = 8 c0/c1 is 1 to depth 3, 2 at depth 4, then grows (8 at depth 8): it does NOT bound E-113's circuits. Caveats: depth = BFS depth through acyclic algebras only; dim count unvalidated for parallel arrows (a doubled relation gives dim 2). Open: restrict to non-parallel pairs; the mechanism behind "components <= 2 edges" must come from the two-term kernel structure (g_i), not a dimension bound; the 4 coefficient-2 pairs; match the 22 to the 28; n = 9.
- **S-1 (E-118):** with E-115 labels the free-end K >= 3 rule holds at n = 8, 9, 10 (3/3, 6/6, 12/12) and fails at n = 11 (20/21; one class of 1305 LNAs splits inside one orbit); K >= 4 holds on the 10 resolved n = 11 classes. Open: per-end table of the K = 3 failure; K >= 4 among the 418 unresolved sources; derive the rule (F-028/H-020); is the threshold K >= 4 at n = 12; core-with-room question 3.
- Requests added: experimentalist: dim count at non-parallel pairs; Gamma_i shapes of the dim-3 out-2 rows; n = 9 c0 prefix. maverick/toolsmith: per-(orbit, K) image table for the n = 11 failure; resolve more n = 11 sources by another route. theorist: derive the K threshold; two-term kernel route to "no circuit >= 3".

## Round 029 updates (supersede everything below where they differ)

Round 029 (ordinary: theorist, skeptic, toolsmith; referees skeptic, theorist, experimentalist; all minor revision, all promoted with the referees' qualifications) recorded E-113..E-115. **Round 030 is ordinary; 032 is the next conference.** The round-028 agenda stands.

- **T5 (E-113):** at n = 8 c0/c1 every `Gamma_i` component has <= 2 edges (no nn 2-cycle, no circuit >= 3; scalar 1). The skeptic's 22 rows are half-W (20) or two loose pendants (6), not length-2 ground paths. Open: the mechanism (is `dim e_iAe_v` <= 1 forced by the tilting complex, so that a bound follows?); the 4 coefficient-2 pairs; match the 22 to the 28; a labelled example of each shape; n = 9.
- **T5 (E-114):** the 61 D' rejects are admitted because the gate tests single paths (near-tautology). Loose D' shape: 84 rows at n = 8 c0, 61 rejected; at n = 6, 7 c0 (5, 20 rows) all accepted; no out-degree 2 reject elsewhere in capped walks. Open: why the 25 loose-shape W-false rows at n = 6, 7 fail W.
- **S-1 (E-115):** the F-047 profile places the 16 n = 9 and 176 n = 10 unresolved LNAs (n = 9 is F-047's own table); it does not separate at n = 11 (24 of 442). The E-112 labels at n = 9, 10 are now independently supported (necessary evidence); n = 11 labels need another route.
- Requests added: skeptic: the 25 loose-shape W-false rows at n = 6, 7. theorist: bound `dim e_iAe_v` (Hom between thin modules / mutation depth); the coefficient-2 pairs; match 22 to 28. experimentalist: max `dim e_iAe_v` by depth at n = 8, 9; S-1 K >= 3 rule with the corrected labels (rerun with E-115 classes) and n = 11 by another route. toolsmith: a hand-built nn algebra whose key is an LNA key (proof candidate); an independent n = 11 resolver.

## Round 028 (conference) -- proposed agenda (proposed, round 028; supersedes the round-024 agenda; ordinary rounds work from it until the human changes it in STEERING.md)

Details: `rounds/028/proceedings.md`. Round 029 is ordinary; next conference 032. T1-T9 and the round updates below stay current.
1. **No long circuits / converse of W** (T5): theorist + scholar prove or refute "no nn 2-cycle, no circuit >= 3 on LNA-derived algebras"; can D, G, H occur inside a walk? First: the skeptic's 22 n = 8 c0 rows (both J nonzero, empty intersection).
2. **Controls for W** (T5): experimentalist + toolsmith hand-build a parallel-arrow positive control and a "cancels" example; skeptic checks gate admission of the 61 D' rejects.
3. **Resolve the 16 n = 9 and 176 n = 10 cospectral LNAs** (S-1 validation): toolsmith + experimentalist; size with `--plan` first.
4. **S-1: derive K >= 3** (maverick + theorist): why K = 2 fails; relation to H-020; then n = 11, S-1 question 3.
5. Carried, lower rank: mirror chain (T6), `k(33x) = 2x` (T2/T4), n = 15..17 (T1/T2).

## Round 027 updates (supersede everything below where they differ)

Round 027 (ordinary: scholar, experimentalist, maverick) recorded E-110..E-112, all accepted with caveats. **Round 028 is a conference** (multiple of 4). The round-024 agenda stands.

- **T5 (E-110):** "reject => W" is not derivable from step 7 (J is defined by the parent alone). For monomial + two-term relations (scalar 1) J != 0 iff the circuit graph Gamma_i has a circuit; W is its length-2 ground-path case. Hand-built D (cancels / nn), G (out-degree 1, no long square), H (3-term) are gate-admitted rejects W misses. On walks J != 0 is always a length-2 ground path. Open: why no nn or longer circuit on LNA-derived algebras; balanced circuits for scalars != 1; >= 3-term relations.
- **T5 (E-111):** W has 0 mismatches on 32 132 fresh out-degree 2 rows (n = 8 c0, 2, 3; n = 9 c0 prefix; 1 135 parallel) but only 61 positives (n = 8 c0). Out-degree >= 3: 0 rejects. Open: parallel-arrow positive control; classes 4..10 at n = 8; n = 9 c0 beyond the prefix (no overnight adopted).
- **S-1 (E-112):** no positional rule for deleting a vertex (same-class pairs kept 0.43 vs 0.19 chance at n = 9); deleting a free end vertex with K >= 3 free vertices there gives an image class depending only on the source class at n = 8, 9, 10 (3, 5, 10 resolved classes). Open: derive it (F-028/H-020?); '?' ends at n = 10; independent validation of the class labels; n = 11; the 16 unresolved n = 9 LNAs; question 3 of S-1 (core with room to move).
- Requests added: theorist: why no nn/longer circuit on LNA-derived algebras; the balanced-circuit lemma; why K >= 3. toolsmith/experimentalist: resolve the 16 n = 9 and 176 n = 10 cospectral LNAs; hand-built parallel and "cancels" controls for W. experimentalist: n = 11 for the K >= 3 rule.

## Round 024 (conference) -- proposed agenda (proposed, round 024; supersedes the round-020 agenda; ordinary rounds work from it until the human changes it in STEERING.md)

Details: `rounds/024/proceedings.md`. Round 025 is ordinary; next conference 028. T1-T9 and the round updates below stay current.
1. **The 42 n = 8 class-0 rejects** (T5): scholar classifies them (revision due); skeptic: reachable by key-preserving walks?; toolsmith: Cartan test through the rewrite; experimentalist: capped walks at n = 8 classes 1, 3 and n = 9 class 0.
2. **Mirror chain** (T6): maverick sorts the 24 failing LNAs by big blocker; experimentalist n = 11/12 deeper shapes; theorist a derivation.
3. **Split words / shuttle / `k(33x) = 2x`** (T1/T2/T4): experimentalist n = 15..17; theorist shuttle row.
4. **Cord = commutativity cycle** (T6): theorist, toolsmith.

## Round 020 (conference) -- proposed agenda (proposed, round 020; supersedes the round-016 agenda; ordinary rounds work from it until the human changes it in STEERING.md)

Details: `rounds/020/proceedings.md`. Round 021 is ordinary; the next conference is 024. T1-T9 and the round updates below stay current.
1. **A5-shape check and kernel/cokernel identity** (T5): experimentalist checks all rejecting parents are A5-shaped; toolsmith replays the 10 E-084 parents and adds the Cartan assertion; scholar/theorist derive "(k,i) = dim coker" from step 7.
2. **Split words, gap, shuttle** (T1/T2): theorist, shuttle as a rule-table row and why the gap is word-only; experimentalist, R-terminal and right gap for `3334 2455 3335` and 4-letter words at n = 12..17; skeptic, selectivity null.
3. **Cords** (T6): toolsmith, non-MONO `--plan` over the 429 n = 8 LNAs; theorist, cord = commutativity cycle.
4. **`k(33x) = 2x`** (T2/T4): theorist, a mechanism.



(Round 023 updates dropped as superseded; see `rounds/023/proceedings.md`. Round 015-017 updates and the round-016 agenda are dropped from this file as superseded; see `rounds/015..017/proceedings.md` and E-085..E-090.)

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
| experimentalist | 030 | 030 (scholar) |
<!-- 012 was a conference: all six wrote position statements -->
<!-- 008 was a conference: all six wrote position statements -->
<!-- 005 was a conference: all six wrote position statements; nobody refereed -->
| theorist | 029 | 030 (experimentalist) |
| skeptic | 029 | 030 (maverick) |
| scholar | 030 | 022 (toolsmith) |
| toolsmith | 029 | 027 (maverick) |
| maverick | 030 | 015 (toolsmith) |
