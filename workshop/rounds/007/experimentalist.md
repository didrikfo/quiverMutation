# 34x: the offsets pair as o <-> hi - o at n = 14..17 (x = 4,5,7,8,9, not 6), equivalently k = x + 3 (not 2x); 45x has no reflection; 4046 is a size-paired reflection, the period-2 translation belongs to 5046/5056 at odd n

author: experimentalist · round: 007 · kind: result (data; no mechanism)
thread: T2 · bears on: H-021, H-020, E-060, E-061, E-065 (theorist request: `34x`, `45x` at n = 14..17; one 4046 core at 14)

## Claim

**34x.** For x = 4, 5, 7, 8, 9 at every n = 14, 15, 16, 17 (20 words, all orbits closed) the offsets of `34x` pair as `o <-> s - o` with `s = hi = n - x - 3`,
so **`k = n - s = x + 3`**, which is equivalent to `s = hi` (the pairing `o <-> hi - o` is the content; the n-independence of k adds nothing beyond 0 and hi sharing an orbit). **`d = 0`** (no offset above `s`) is true by `s = hi`; it is not the raw singleton count `d` of the table, which is 0..5, nor `d_eff`. The pairing of the equal-size middle singletons is **by size only**, not verified as same-orbit (no mirror-join run). `34x` is `k = x + 3`, **not** `2x`: the gap `2x - k = x - 3` is
1, 2, 4, 5, 6 at x = 4, 5, 7, 8, 9. `346` never fits: one orbit of all offsets at n = 14..17 (the `333@0` orbit of E-065, unverified here). So the answer to "which x fit `s = n - k`
with n-independent k": all x except 6; none fit with `k = 2x`.
**45x.** No reflection at any x, n = 14..17. `456`, `457`: one orbit of all offsets (as `346`). `455`, `458`, `459` split by **parity of the offset** at some n and are one orbit at others
(table below). The `k` the table script printed for `458@15`, `459@14, 16` is an artefact (a pair sum taken from a class that also has three members); it is not a fit.
**4046.** Not a period-2 translation. `4046` at n = 12..16 is a reflection `s = n - 11`, `k = 11`, `d = 1`. The parity-class translation is seen in `5046` and `5056`, and only at odd n.
Not claimed: any mechanism; x >= 10; n >= 18; mirror-joined orbits (I use orbit only and treat equal-size singletons `o`, `s - o` as a mirror pair, see caveat).

## Evidence

Script `experimentalist_kd.py` (held offsets of the same word, `freeMoves.REDUCED` walk, limit 300000, as `theorist_link.py` but prints `s`, `k`, `d`). "d_eff" pairs singleton
orbits `o`, `s - o` of equal size (the mirror pair of E-064) and drops the centre `2o = s`. `d` raw counts those singletons.

| word | n=14 | n=15 | n=16 | n=17 | k | 2x |
|---|---|---|---|---|---|---|
| 344 (x=4) | s7 d0 | s8 d3 | s9 d0 | s10 d5 | 7 | 8 |
| 345 | s6 d1 | s7 d0 | s8 d1 | s9 d0 | 8 | 10 |
| 346 | one orbit | one orbit | one orbit | one orbit | none | 12 |
| 347 | s4 d1 | s5 d0 | s6 d1 | s7 d0 | 10 | 14 |
| 348 | s3 d0 | s4 d1 | s5 d2 | s6 d3 | 11 | 16 |
| 349 | s2 d1 | s3 d0 | s4 d3 | s5 d2 | 12 | 18 |

(`sN dM`: pair sum N, raw singleton count M; `d_eff` = 0 in all 20 cells, every raw singleton is the centre or one of an equal-size mirror pair, e.g. `348@17`: offsets 2, 4 both 9622, centre 3 alone.)
Power (pair counts only; the referee notes a wrong k is not independently constrained per n): pairs with distinct orbit sizes per cell range from 1 (`349@14`, `348@14`: 2) to 5 (`344@16`, `345@17`); Low-power cells: `349@14` (one pair), `348@14`, `349@15` (2).

45x (offset classes, orbit sizes in the `.txt` files):

| word | n=14 | n=15 | n=16 | n=17 |
|---|---|---|---|---|
| 455 | {0,2,4,6}{1,3,5} | one | {0,2,4,6,8}{1,3,5,7} | one |
| 456 | one | one | one | one |
| 457 | one | one | one | one |
| 458 | one | {0,2,4}{1,3} | one | {0,2,4,6}{1,3,5} |
| 459 | {0,2}{1} | one | {0,2,4}{1,3} | one |

Parity split at (455: even n), (458: odd n), (459: even n). 456, 457 never split. No fixed `n + x` rule covers 456/457, so I do not offer one.

4046-type cores (`experimentalist_core4046.py`; classes of held offsets, orbit size):

| word | n | classes (size) |
|---|---|---|
| 4046 | 12 | {0,1}(261) {2}(1766) |
| 4046 | 13 | {0,2}(63) {1}(104) {3}(2116) |
| 4046 | 14 | {0,3}(383) {1}(388) {2}(388) {4}(11820) |
| 4046 | 15 | {0,4}(81) {1}(149) {3}(149) {2}(85) {5}(8794) |
| 4046 | 16 | {0,5}(521) {1}(530) {4}(530) {2}(534) {3}(534) {6}(77735) |
| 5046 | 11 | {0}(454) {1}(218) |
| 5046 | 12 / 14 / 16 | one orbit (1766 / 11820 / 77735) |
| 5046 | 13 | {0,2}(4217) {1,3}(2116) |
| 5046 | 15 | {0,2,4}(18157) {1,3,5}(8794) |
| 5056 | 13 / 14 | as 5046 at 13 / one orbit (11820) |

`4046`: hi = n - 10; pairs `0 <-> n-11`, `1 <-> n-12`, ... with the middle ones as equal-size singletons (388,388; 149,149; 530,530; 534,534) = mirror pairs; the last offset `hi` is a
lone large orbit (size equal to that of `5046`'s one orbit at n = 14, 16). So `s = n - 11`, `d = 1` at 12..16, with the middle singletons paired by equal size only. `5046`/`5056` at 13 and 15: two orbits by parity of offset, at 12, 14, 16 one orbit:
period-2 translation (steps of 2 stay in the class, steps of 1 leave it) at odd n, consistent with E-059's parity split. Notice `4046@13 offset 3` (2116) and `5046@13 {1,3}` (2116) have the same size (probably one orbit, not checked).
**Discrepancy to flag:** E-060 and STATE say `4046 5046 5056` all give `{0,2},{1,3}` at n = 13. I get that for `5046`, `5056`; for `4046` at 13 I get `{0,2},{1},{3}`. Either the old note listed the family loosely or it used another orbit notion; I did not find the ledger.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/007/experimentalist_kd.py 16 34 3-9      # 34x at n = 16: 8 min; n = 14 1 min; n = 15 2 min; n = 17 needs x = 8, 9 separately (3 min each), 347 6 min
timeout 10m .venv/bin/python workshop/rounds/007/experimentalist_kd.py 17 45 3-9      # 45x n = 17, 7 min
.venv/bin/python workshop/rounds/007/experimentalist_table.py                        # the table, reads experimentalist_kd_*.txt (1 s)
timeout 10m .venv/bin/python workshop/rounds/007/experimentalist_core4046.py 14 4046  # 4046 at 14 (11 s); n = 16 99 s
```
Outputs kept: `experimentalist_kd_{34x,45x}_n{14..17}.txt` (`34x_n17.txt` stops after x = 7; x = 8, 9 in `34x_n17_8.txt`, `_9.txt`), `experimentalist_core_*.txt`.
Four processes in parallel, each under 10 min; `34x_n17` alone hit the 10 min budget at x = 8 and was rerun per x.

## Prior record

E-061 (`k(33x) = 2x`, `d = x - 3`), E-065/theorist round 006 (`34x` only `345` explained; `k = x + 3` observed at x = 5, 7, 8 with low power; `346` joins the `333@0` orbit (UNVERIFIED here; taken from round 006); `45x` parity at 455 n = 16, no reflection).
New here: `k = x + 3` at x = 4, 5, 7, 8, 9 over n = 14..17, including the powered cells the theorist asked for (n = 16, 17) and x = 9; `d_eff = 0` throughout; `458`, `459` parity splits at
n = 14..17; `4046` as a reflection (k = 11, d = 1) and `5046/5056` as the odd-n translation. Grepped `research/` for `4046`: only the E-060 and E-056 lines above; nothing on `k(34x)`. `RETRACTIONS.md` has nothing on these.

## Code changed

None in the library. New: `experimentalist_kd.py`, `experimentalist_table.py`, `experimentalist_core4046.py` (all `workshop/rounds/007/`). No tests run (library untouched).

## Next

* theorist: `k(34x) = x + 3` means `hi = s`, i.e. the reflection covers the whole placement range (`d = 0`): this is a different kind of pairing from `33x` (overhang `x - 3`). Is it the mirror of the LNA (reversal sends `34x` to `x43`, not itself), or a chain as in E-065 with a seed? `346` joining the `333@0` orbit is the exception to explain; so is `4046` having `hi` as a lone big orbit.
* skeptic: the pairing is by orbit + equal size for the middle singletons; a mirror-join run (E-064 `toolsmith_orbitclass.py`) on `344`, `348`, `349` at n = 15..17 would convert "equal size" into "same orbit". Also a null: random `s` with the same orbit sizes?
* experimentalist (me): `34x` at n = 18 (347 at 17 took 6 min, 18 is over budget, OVERNIGHT proposal: `347..349` at n = 18, 19, about 40 min each); `44x` was not asked and not run.
