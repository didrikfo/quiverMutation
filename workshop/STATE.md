# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 7
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

## Round 007 updates (supersede everything below and above where they differ)

Round 007 (ordinary: experimentalist, skeptic, maverick) recorded E-067..E-069. Next round (008) is a **conference**.

- **T2 (E-068):** `34x` pairs `o <-> hi - o` (`k = x + 3`, not `2x`) at n = 14..17 for x = 4, 5, 7, 8, 9; `346` is one orbit; `45x` has no reflection (parity splits for `455`, `458`, `459`); `4046` is a reflection `k = 11` (size-paired), the parity translation is `5046/5056` at odd n. E-060's `4046@13` = `{0,2},{1,3}` did not reproduce (`{0,2},{1},{3}`). Open: mirror-join check of the equal-size singleton pairs; `44x`; x >= 10; n >= 18.
- **T1/T2 null (E-067):** 39 of 109 n = 13 fits are one-orbit (vacuous); informative fits mostly chance-level; the interior-core centre formula and the 12 E-060 cores at n = 15, 16 survive. Open: a neighbour-aware null (none committed); a random (not pre-selected) sample of cores at 14-15; null for `k(33x) = 2x`.
- **T6 (E-069):** the H-017 search finds a class iff a member is within its depth (round trips 273/273 at n = 7, 91 of 132 LNAs; 0 at depth L-1). The depth-4 negative for the 16 below-diagonal candidates is therefore weak. Open: L = 5 control at n = 6/7 (and check that recorded paths are shortest); depth 5-6 on the 16 candidates at n = 9 (size first; overnight only if unavoidable).

Requests added: theorist to experimentalist/toolsmith: mirror-join on `344`, `348`, `349` at n = 15..17 and `4046` at 14..16; maverick: L = 5 control; skeptic: null for `k(33x)`.

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
| experimentalist | 007 | 007 (maverick) |
<!-- 005 was a conference: all six wrote position statements; nobody refereed -->
| theorist | 006 | 007 (skeptic) |
| skeptic | 007 | 007 (experimentalist) |
| scholar | 006 | 004 |
| toolsmith | 006 | - |
| maverick | 007 | - |
