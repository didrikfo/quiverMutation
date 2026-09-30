# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 3
next_round_kind: ordinary

## Open threads

Round 002 revised T1, T3, T5 (records E-056, E-057, E-058). Each line: id · question · suited to · status.

- **T1** · H-021: restated without the mirror clause (round 003, H-021': a pairing core has `s = n - k(c)`, `k` and `d` independent of n; covers only cores that pair). The 7 mirror-without-reflection cores pair at even n 12..18 and fail at odd n 13..17 (E-059): an odd-n effect for them. Open: does parity govern the other 132 cores (census at n = 12 and 14, then 15); why the outside block folds; the n-independence of `s`, `d` rests on 12 cores at one step. · experimentalist, theorist · open
- **T2** · Centre `s(c) = n - k(c)`: for interior outside blocks `d = t - h`, `s` = first + last outside offset (13/13 at 13, E-060). Open: the 62 all-inside cores (`d` invisible in the slide; `33x` gives `|d| = x - 3`), the 4 end-touching failures `4045 3556 4506 4556`, the 3 all-outside failures `4046 5046 5056`, and the formula for `k(c)`. Commit the interior/end-touch split as a column of `theorist_shortfall.py`. Test interior blocks with `m >= 5` at n >= 15. · theorist, skeptic · open
- **T3** · `3346` never pairs (proved by key, E-054). `4056` at n = 16: `{1,2}` are two mirror-image orbits, key pairs them (E-058), so "orbit = key class" is false there; compare orbit-plus-mirror against key over the `--max-word 4` catalogue at 14..16. Open: other cores that split; the other four H-020 failures (ledgers not in repo). · experimentalist, skeptic · open
- **T4** · H-020 as a theorem: the rule table acts the same at every interior position (F-051). State and prove. · theorist · open
- **T5** · H-015: Ladkani 2.3(c) agrees with the gate everywhere tested, incl. non-monomial parents at n = 5..7 (E-055, E-057); only gate-admitted rejection is E-032 step 7. `isTilting` not promoted (STEERING q3). Open: the audit at n = 9/10 (overnight, Menu 4). · scholar, toolsmith · waiting on overnight
- **T6** · H-017, a quipu carries more relations than cords: the Euler-form count; test on recorded classes at small `n`. · maverick, theorist · open
- **T7** · H-010, overlap reducible only at an end: a proof from step 7 of arXiv:2112.08129. Do not run `probe.py --steps 7`. · theorist · open
- **T8** · Tooling: `batch.py orbits N` exists (round 003; `45`, `344` pinned). Open: run it over the whole catalogue (`--jobs 4`), commit the fit/slide step as a task, orbit-plus-mirror prefilter (not key). · toolsmith · open
- **T9** · H-019, H-013: long runs. Overnight waiting on the human: n = 14 census of the 139 cores (about 90 min, 4 procs; use `workshop/rounds/002/experimentalist_census.py` and `experimentalist_fit.py`); Ladkani audit n = 9/10. · any · parked

Suggested next round (the chair may change it): theorist on T2 (the column, the 4 + 3 failures, `33x`); experimentalist on the 7 cores plus the 12 of E-060 at n = 15, 16; skeptic on interior blocks with `m >= 5`; maverick on T6 (never yet worked); toolsmith on the whole-catalogue orbit run if time allows.

## Awaiting revision

(none)

## Requests between personas

- theorist (to experimentalist): n = 15, 16 for the 12 cores of E-060 to test the shift of `s`, `d`.
- experimentalist: orbit-plus-mirror vs key classes over the catalogue at 14-16 (skeptic).
- theorist: at n = 16 `4056` has a self-dual pair `{0,3}` and a mirror pair `{1,2}`: is the split a function of parity or of position relative to the middle? (skeptic)
- theorist: the centre `s(c)` formula (round 001 request, still open).
- toolsmith: committed census task; key prefilter with the E-058 caveat.

## Rota

<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 003 | 003 (theorist) |
| theorist | 003 | 003 |
| skeptic | 002 | 003 (x2) |
| scholar | 002 | - |
| toolsmith | 003 | - |
| maverick | - | - |
