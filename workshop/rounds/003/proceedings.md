# Round 003 -- proceedings

Step 0.5: no open question was unanswered (round 002 q1, q2 answered by the human; both applied this round: theorist restated H-021 without the mirror clause, experimentalist ran the 7 survivors). `main` merge: no-op.

## experimentalist -- T1
* **Claim.** At n = 14 the 7 survivors all pair by reflection with no strict mirror; over n = 8..18 they pair at every even n 12..18 and fail at 13, 15, 17.
* **Referee (theorist): minor revision**; reproduced (table exact; `344`, `4405` re-run at 14, 15, 16). Wanted: title qualified to these 7, `d` per core, reconcile 7 vs 8, mark "would pair" as interpretation.
* **Decision: accept.** All four points are wording or a table I built from the committed fit script (`d` = 0/1/2 by core, constant across even n; 8 strict-true = 7 + `406`). Written into E-061.
* **Promoted.** E-061.

## theorist -- T1, T2
* **Claim.** H-021 restated without the mirror clause; for interior outside blocks the shortfall is `t - h`, `s` = first + last outside offset (13/13 at n = 13, 3/3 at n = 14); end-touching 17/21.
* **Referee (skeptic): minor revision**; reproduced. Errors: "m is 3 to 5" (2 to 4 at 13), "exactly" a block at an end (one direction only). Partly known: F-053 has it for `45`. Wanted the interior/end split as committed code.
* **Decision: accept with caveats.** Wording fixed in E-062; the missing column is recorded there as not done and goes to next round (theorist). The result is near-tautological by (a), and the entry says so.
* **Promoted.** E-062; H-021 status line updated (still OPEN).

## toolsmith -- T8
* **Claim.** `batch.py orbits N` reproduces E-052's report for `45` at 13 in 8 s; test pins it.
* **Referee (skeptic): minor revision**; reproduced, `3346` also matches E-056 data. Wanted a `344` cross-mirror pin, a partial-ledger resume test, a docstring caveat.
* **Decision: accept, after I added all three** (`tests/test_orbits_task.py`: 5 pass; docstring in `batch.py`). No E entry: tool, no new result. Not done: the whole-catalogue run (`--jobs 4`, not run this round); no Coxeter prefilter (key is not orbit, E-060).

Glossary: shortfall/overhang, odd-n effect.

## Decisions taken for the steering committee
None.

## Questions for the steering committee
1. The experimentalist's finding suggests the full n = 12 and n = 14 censuses of the 139 cores are more informative than n = 14 alone (does the 30-core "no fit" set change with parity?). Add both to `OVERNIGHT.md` Menu 4 (each about 90 min on 4 procs, via `batch.py orbits`)? **Recommend: yes** (not done this round, since the human runs them). If unanswered the next chair adds them to the menu.
