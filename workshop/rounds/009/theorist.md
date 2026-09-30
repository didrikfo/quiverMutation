# The drift family `aax` closes (computed criterion) when its seed `aaa` does not collapse onto the word `34`; `33x`, `55x`, `66x` close, `44x` does not

author: theorist · round: 009 · kind: result (criterion, computed; not proved)
thread: T2/T4 · bears on: H-020, H-021, F-032, F-051, E-061, E-065, E-068

## Claim

Drift families are `aax` (a = 3..6): the double mutation sends `aax@o -> aa(x-1)@(o+1)` in the interior (verified for `33x` in E-065; for `44x`, `55x`, `66x`
see the neighbour table below), so `c = x+o` is conserved along a chain that ends at the seed `aaa`. The chain closes into a pair `{c, n-c}` (an orbit of
only chain rows, `k = 2x + 3 - a`) unless the chain reaches the big orbit `O*` (the orbit of `333@0`, size 3767 at n = 14). **Criterion (interior
collapse):** the seed `aaa` has a double-mutation neighbour of span 2, its collapse `(a-1)a`. For a = 3, 5, 6 this is `23`, `45`, `56`; the orbit of `aaa@o`
(o interior) then holds no span-<=2 word at all its offsets, and the family is rigid. For a = 4 the collapse is `34`, and the orbit of `444@o` holds `34`
at every offset (`34` is joined to the sliding word `44` by a computed 7-step path); the whole family `44x` is one orbit with `333@0`. Tested: n = 14,
prediction made before running for `55x` (x = 5..8) and `66x` (x = 6..8): both predicted rigid, both rigid, all orbits closed. Not claimed: a proof
of why `34` (and only `34`) reaches the slider; n other than 14 for `55x`, `66x`; `34x`, `45x` (no drift: `45x` is not a drift family, see below).

## Evidence

Local facts at n = 20, offset 6 (`theorist_nbrs.py`, `theorist_translate.py`):

| seed | double-mutation neighbours | collapse (span 2) |
|---|---|---|
| 333 | 23@6, 302@7, 334@5, 4033@6 | 23 |
| 444 | 34@6, 403@7, 445@5, 5044@6 | 34 |
| 555 | 45@6, 504@7, 556@5, 6055@6 | 45 |
| 666 | 56@6, 605@7, 667@5, 7066@6 | 56 |

Each seed also has the drift neighbour `aa(a+1)@(o-1)` and two rows with a zero letter (`302@7`, `4033@6`, ...); every `aa` slides (`aa@o -> aa@(o+-1)`, all a = 2..9),
so sliding alone does not separate the cases. `aa -> aaa` expansions: `23 -> 333`, `34 -> 444`, `45 -> 555`, `56 -> 666` (one neighbour each).

Orbit test at n = 14 (`theorist_closure.py`: for each offset the reduced-walk orbit, closed in every case; "held" = offsets of the word in one orbit):

| family | interior collapse reaches full translator | prediction | outcome (held sets) |
|---|---|---|---|
| 33x, x=3..7 | no (`23`); only `333@0,8` reach `34` (edge) | rigid | pairs `o <-> n-2x-o` (E-065 rechecked; sizes 320, 886, 364, 491) |
| 44x, x=4..7 | yes (`34`, all 8 offsets of `444`) | merged | one orbit, size 3767, all offsets |
| 55x, x=5..8 | no (`45`) | rigid | s = 6, 4, 2, 0 (`[0,6],[1,5],[2,4],[3]` for 555), k = 2x-2 |
| 66x, x=6..8 | no (`56`) | rigid | s = 5, 3, 1 (`[0,5],[1,4],[2,3]` for 666), k = 2x-3 |

The rigid values fit the E-065 formula `k = 2x + w0 - x0` with w0 = 3 and x0 = a: `k = 2x + 3 - a`. The end link is uniform: in `33x`, `55x`, `66x` the rows in
the size-3767 orbit `O*` are exactly the seed at offset 0 (`333@0`, `555@0`, `666@0`) and the right-touching rows `aax@hi`; every other orbit is a pure chain pair. So the
unpaired offsets (o > s) are the right-touching ones and the seed at the source end: that is the observed form of the end link, for three values of a
instead of one. Sizes are shared across families (11820 = `55x` offset 1..2 = `66x` offset 1 = `45x` orbits), i.e. chain orbits are not family-specific.

A wider check (`theorist_translator.py`): "the orbit holds a span-<=3 word at every offset" does NOT predict merging (10 of 24 cells agree for `33x`), because the
orbits of `334@1` etc. hold `36`, `66` at every offset yet stay rigid. Only a translator reached from the seed by a collapse works. The shortest path
`34@5 <- 403@6 ... 34@7 <- 3333@7 <- 3403@7 <- 44@7 <- 44@6 <- 44@5` (`theorist_path34.py`, n = 14) shows which moves join `34` to `44`; the step I cannot explain is why
`3333 <-> 3403` (a free move) has no analogue for `23`, `45`, `56`. That step is the weakest link, and it is computed, not derived.

`45x` is not a counterexample: its seed `455` collapses to `35`, and its drift is absent (neighbours `5555`, `5666`), so the `45x` orbits are parity classes
(`455`: `{0,2,4,6},{1,3,5}` sizes 320, 272) and `456..458` merge into `O*`/11820 at n = 14 (the 456 row is the same 3767 orbit). The criterion is silent there.

## Reproduction

```
.venv/bin/python workshop/rounds/009/theorist_nbrs.py 44 4-6            # seconds: local neighbour table (also 33, 55, 66, 34, 35)
.venv/bin/python workshop/rounds/009/theorist_translate.py              # seconds: aa self-slides; xy -> xyz expansions
timeout 10m .venv/bin/python workshop/rounds/009/theorist_closure.py 14 55 5-8   # ~4 min; outputs theorist_closure_{33,44,45,55,66}_n14.txt
timeout 10m .venv/bin/python workshop/rounds/009/theorist_path34.py     # seconds: 34 to 44 path
timeout 10m .venv/bin/python workshop/rounds/009/theorist_translator.py 14 33 4-8   # the failed broader criterion (10/24)
```

## Prior record

E-065 states the drift and `33x` closure and says `44x` fails without a reason; F-051/H-020 give the interior-uniform rule; E-068 records `34x`, `45x`. `grep` of
`research/` for "collapse", "translator", "55x", "66x" finds nothing: the `55x`/`66x` closure and the collapse-to-`34` reason are new. They are not in
`RETRACTIONS.md`. The E-065 statement "the argument does not extend to 44x" becomes: the argument extends to every `aax` with a != 4, and fails at a = 4
because of `444 -> 34`.

## Code changed

None in the library. New scripts only, in `workshop/rounds/009/` (`theorist_*.py`); no tests touched.

## Next

- Explain `34`: why `34@o ~ 34@(o+2)` through `403` and `3333 -> 3403 -> 44`, while `23`, `45`, `56` have no such route (hand computation with F-032's double mutation; this is the proof gap).
- Run `55x` at n = 15 and 16 (the n = 14 run took ~4 min; try `--plan`-sized first) and `77x` (needs n >= 15) to confirm `k = 2x + 3 - a` and "no collapse to 34".
- experimentalist: does `44x` also fail at n = 17+ (it is one orbit at 14..16, E-065).
- skeptic: a neighbour-aware null for "a = 4 is special" (four values of a tested, one exception).


## Chair note (round 009, after referee)
Referee's extra runs agree: `55x` at n = 12, 13; `66x` at n = 13; `44x` at n = 12 (merged, 1410); `77x` at n = 15 (rigid, pairs [0,5],[1,4],[2,3]). The drift for x > a is evidenced only through the pair structure, not shown separately. The criterion fits one positive datum (a = 4); `PREDICT` lines are printed by `theorist_closure.py` before the run. Computed, not proved.
