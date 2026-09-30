# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 4
next_round_kind: conference

## Open threads

Round 004 records E-061, E-062, E-063. (Round 005 is a conference: round 003's STATE wrongly said 004 was ordinary.) Each line: id · question · suited to · status.

- **T1** · H-021: restated without the mirror clause (H-021': a pairing core has `s = n - k(c)`, `k`, `d` independent of n; covers only cores that pair). The 7 cores of E-059 pair at even n, fail at odd; the 12 cores of E-060 keep `k`, `d` at 15 and 16 and pair at 13, 14, 15; at 16 `46 3355 3445` lose the fit through an unmerged equal-size middle pair (size 20300, as `4056` in E-058), so parity is not the whole story (E-062). Open: is 20300 one shared orbit; the 3 at n = 17, 18; the other 127 cores (census at 12, 14 in Menu 4, overnight); why the outside block folds. · experimentalist, theorist · open
- **T2** · Centre `s(c) = n - k(c)`. Interior blocks: `s` = first + last outside offset (13/13, E-060, E-061); end-touching failures `4045 3556 4506 4556` pick the smaller consistent centre for three and the larger for `4506` (no min/max rule); all-outside `4046 5046 5056` maybe a period-2 translation, not a reflection (unverified). `k(33x) = 2x`, `d = x - 3` for x = 3..7 (E-061). Open: why `k = 2x`; `k(c)` for other prefixes (`34x`, `44x`, `45x`); a null test for `|R| <= 4` fits; T4 may explain the top `x - 3` singleton offsets. · theorist, skeptic · open
- **T3** · `3346` never pairs (proved by key, E-054). `4056` at n = 16: `{1,2}` are two mirror-image orbits, key pairs them (E-058), so "orbit = key class" is false there; compare orbit-plus-mirror against key over the `--max-word 4` catalogue at 14..16. Open: other cores that split; the other four H-020 failures (ledgers not in repo). · experimentalist, skeptic · open
- **T4** · H-020 as a theorem: the rule table acts the same at every interior position (F-051). State and prove. · theorist · open
- **T5** · H-015: Ladkani 2.3(c) agrees with the gate everywhere tested, incl. non-monomial parents at n = 5..7 (E-055, E-057); only gate-admitted rejection is E-032 step 7. `isTilting` not promoted (STEERING q3). Open: the audit at n = 9/10 (overnight, Menu 4). · scholar, toolsmith · waiting on overnight
- **T6** · H-017: not refuted at n = 9 to depth 6; the Coxeter polynomial and Euler form do not predict (cords, relations); the Euler signature `pos <= n-2` separates outside-every-quipu-class for n = 8..11 (E-063). Open: n = 12 signature (about 15 min in shards); an LNA in the n = 13 class of `P^(1,0,3,0,1)_(1,1,1,1)`; a positive control for the mutation search; a tree eigenvalue argument for `relations >= cords + 1`; depth 7 at n = 9 (overnight, pending question). · maverick, toolsmith, theorist · open
- **T7** · H-010, overlap reducible only at an end: a proof from step 7 of arXiv:2112.08129. Do not run `probe.py --steps 7`. · theorist · open
- **T8** · Tooling: `batch.py orbits N` exists (round 003; `45`, `344` pinned). Open: run it over the whole catalogue (`--jobs 4`), commit the fit/slide step as a task, orbit-plus-mirror prefilter (not key). · toolsmith · open
- **T9** · H-019, H-013: long runs. Overnight waiting on the human: n = 12 and n = 14 censuses of the 139 cores (Menu 4; `workshop/rounds/002/experimentalist_census.py` and `experimentalist_fit.py`); Ladkani audit n = 9/10. · any · parked

Round 005 is a conference (every persona writes a position statement; the chair proposes an agenda). Then ordinary rounds work from the agenda.

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
| experimentalist | 004 | 003 (theorist) |
| theorist | 004 | 003 |
| skeptic | 002 | 004 (x2) |
| scholar | 002 | 004 |
| toolsmith | 003 | - |
| maverick | 004 | - |
