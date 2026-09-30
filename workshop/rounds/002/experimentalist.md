# H-021's mirror clause fails under every reading tried: loose (20 counterexamples), strict (108), and weak-strict (7)

author: experimentalist · round: 002 · kind: negative (revision of round 001)
thread: T1 · bears on: H-021, F-053, E-052, H-020

## Response to referee

Referee (skeptic) verdict was minor revision; three items.

1. **State which reading is refuted; run the strict reading over the 139 cores.** Done, at n = 13, all 139 cores re-walked with the committed script (395 orbit walks, 0 caps, all closed). Three readings, defined in `experimentalist_fit.py`: *loose* (some orbit holds the mirror of some placement, as round 001), *strict* (the mirror of `c@p` lies in an orbit whose held offsets do not include `p`: the mirror of `c@p` in the orbit of `c@q`, `q != p`, different orbits; STEERING's reading), *strict2* (an orbit holds the mirror of `c@p` and also some `c@q`, `q != p`; this also counts a pair `{p,q}` holding its own mirrors). Results below: the loose reading is refuted in one direction (20 cores); strict is refuted in the other (108 pairing cores have no strict mirror); strict2 is the only one where "pairing => mirror" survives, and "mirror => pairing" fails for 7 cores. So no reading gives an "exactly when". The referee expected `344 4044 4403` to survive strict: they do, with `366 4404 4405 4605`, 7 of the 20. The other 13 (`3033 3034 3044 3303 3346 3456 3466 3566 3606 4056 4406 506 6005`) do not.
2. **Overhang-fit code.** Now committed: `workshop/rounds/002/experimentalist_fit.py` (rule in its docstring), input `workshop/rounds/002/experimentalist_census_n13.jsonl` (the raw census, 139 lines), output `experimentalist_fit_n13.txt` (one line per core). It reproduces round 001's 62/28/17/2/30 exactly. Caveat kept: `d` up to 6 on ranges of 4 to 7 offsets is generous; only "no fit" is strong.
3. **`3346` sentence.** Reworded: under the loose reading `3346` and `4056` hold their own placements' mirrors (each orbit its own offset), so the loose test cannot separate them from `45`. Under the strict reading they hold no mirror and agree with "no pairing", so H-021's text about `3346` is *correct* under the strict reading and misstated only under the loose one. Round 001 claimed the misstatement without giving the reading.

## Claim

At `n = 13`, over the 139 placeable single-cluster cores of `--max-word 4` (all orbits closed), let P = "orbits pair by some reflection `o <-> s - o`, overhang `d <= 6`" (109 cores). Then: loose mirror holds for 129 cores; P => loose (109/109) but loose => P fails for 20. Strict mirror holds for 8 cores; strict => P fails for 7 (`344 366 4044 4403 4404 4405 4605`), P => strict fails for 108. Strict2 holds for 116; P => strict2 (109/109); strict2 => P fails for the same 7. H-021 survives on no reading as an iff. Not claimed: anything at n >= 14, or that the 7 are not "pairing with a defect" (`344` pairs `{0,6}`,`{1,5}` and swaps `2`,`4` by mirror only); nor two-cluster words.

## Evidence

| P (reflection fit) | loose | strict | strict2 | cores |
|---|---|---|---|---|
| yes | yes | no | yes | 108 |
| yes | yes | yes | yes | 1 (`406`) |
| no | no | no | no | 10 (`30xy`, `330x`) |
| no | yes | no | no | 13 (`3033 3034 3044 3303 3346 3456 3466 3566 3606 4056 4406 506 6005`) |
| no | yes | yes | yes | 7 (`344 366 4044 4403 4404 4405 4605`) |

Agreement of P with each reading: loose 119/139, strict 24/139, strict2 132/139. The 7 survivors have small orbits (sizes 18 to 50) and cross-offset mirrors, e.g. `344`: `{0,6} {1,5}` own mirrors, `{2}` holds the mirror of `4` and `{4}` that of `2`. `45` at 13 again `{0,5} 2386 {1,4} 1127 {2,3} 4217 {6} 447`, each orbit holding exactly its own offsets' mirrors (why strict is false for every clean pairing: the mirror acts inside the orbit).

## Reproduction

```
# 4 shards, each timeout 10m; wall ~5 min at n = 13 (139 cores, 395 walks)
for k in 0 1 2 3; do timeout 10m .venv/bin/python workshop/rounds/002/experimentalist_census.py 13 --shard $k/4 --out logs/t1-census-n13-shard$k.jsonl & done; wait
.venv/bin/python workshop/rounds/002/experimentalist_fit.py logs/t1-census-n13-shard*.jsonl   # seconds
.venv/bin/python workshop/rounds/002/experimentalist_fit.py workshop/rounds/002/experimentalist_census_n13.jsonl   # same, from the saved census
```

## Prior record

E-052/F-053 (pairing for nine cores), H-021 ("exactly when"), round 001's E-053 (loose reading fails). Nothing in RETRACTIONS on H-021. The strict-reading table and 7-core list are new; the census script is now in the repository (`experimentalist_census.py` is round 001's, unchanged apart from location).

## Code changed

Added `workshop/rounds/002/experimentalist_census.py` (copy of round 001's), `experimentalist_fit.py`, data files above. No library changes, no tests touched.

## Next

* Chair/theorist: restate H-021 without the mirror clause, or as a statement about `d(c)`; the mirror does not discriminate. Open: is `d(c)` = `|head - tail|` (H-020)? Not tested.
* Skeptic: the 7 strict cores are a family (`344`-type, `44xy`): a reflection with a defect or a different mechanism?
* Overnight (already in Menu 4): n = 14 census with this same script, then `experimentalist_fit.py`.
