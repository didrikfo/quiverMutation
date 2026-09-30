# The interior/end-touch split is committed and explains none of the 7 failures; but k(33x) = 2x exactly, n = 13..17, so d = x - 3 is a fact about that family

author: theorist · round: 004 · kind: result (column + partial negative + family formula)
thread: T2 · bears on: H-021, H-020, F-053, E-060

## Claim

(1) *The column (round 003 referee's item 1).* `theorist_shortfall.py` now prints `cls` = allI | allO | int | end0 | endhi | endboth
(outside block touches neither / offset 0 / offset hi / both), `pred` = first + last outside offset, and `cons` = all
slide-consistent centres. n = 13, 109 cores with a fit: **int 13/13** (pred = fitted `s`), **end0 7/10, endhi 10/11**,
**allO 10/13**, allI 45/62. Reproduces E-060's 13/13 and 17/21 (end0 + endhi = 21). Definition used: block = the set of offsets
with slide `o`; "touches" = contains offset 0 or hi.
(2) *What the split does not explain.* All 7 failures are cases where the slide is consistent with the fitted `s` (fitted `s` is in `cons` for
**109/109** cores) -- the slide only rules centres out; it never selects. In the 4 end-touching failures the orbit picks the *smaller*
consistent centre for `4045 3556 4556` (`cons` = {1,2}(+7); `s = 1`: only `{0,1}` merges, pred 2 would merge `{0,2}`) and the *larger* for
`4506` (`cons` = {4,5}, `s = 5`; pred 4). So no rule "min/max of `cons`" either. The 3 all-outside failures `4046 5046 5056` have
a silent slide (every centre 1..5 is consistent) and orbits `{0,2},{1},{3}` / `{0,2},{1,3}`: the second is a *parity* pattern (translation
by 2), not a reflection, and both `s = 2` and `s = 4` fit it (`5046`, `5056`, `5006`, `5066`: `allfit = [2,4]`). I read these 3 as cores where
the "pairing" is not a reflection at all. All 7 have `|R| = 4` (offsets 0..3) or 5, where a fit of one adjacent pair has almost no power; the
same holds for many passes (`3446 3666 3445`: a single 2-orbit). I have no test that separates "block folds" from "one coincidental merge" at `|R| <= 4`.
(3) *`33x`, the all-inside family.* Fitted centre and shortfall, n = 13..17 (`33` is one orbit for n >= 14, vacuous):

| core | hi | s | k = n - s | d = hi - s | unpaired offsets |
|---|---|---|---|---|---|
| 333 | n-6 | n-6 | 6 | 0 | the middle only |
| 334 | n-7 | n-8 | 8 | 1 | `hi` |
| 335 | n-8 | n-10 | 10 | 2 | `hi-1, hi` |
| 336 | n-9 | n-12 | 12 | 3 | `hi-2 .. hi` |

**`k(33x) = 2x`, `hi = n - x - 3`, `d = x - 3`, with the `d` top offsets `s+1..hi` singleton orbits**; identical at n = 13, 14, 15, 16, 17 (20
fits, 0 deviations; orbits are exact reflection pairs `{o, s-o}`). Not tested: `x >= 7` (not tried),
other families, and n >= 18. What it does not say: why `k = 2x`. It is a description; the mechanism is not found.

## Evidence

`theorist_shortfall_n13.txt` (per-core: slide, cls, s, hi, d, pred, ok, in_cons, cons; last lines the class table),
`theorist_k_n13.txt` (per core: digit sum, length, hi, `w = n - hi`, s, k, d, every centre that closes the orbits; the fitted centre is unique for
most cores, so d is not a fit artefact except `36 405 5004 5006 5046 5056 5066`, several centres), `theorist_33x_n{14..17}.jsonl`,
`theorist_33x.txt`. Side data from `theorist_k_n13.txt` (not a result): in the four-digit families `3344 3355 3366` and `3404 3405 3406` `d = -2`
constant, `3444 3445 3446 3555 3556 3666` `d = 2`, i.e. increasing the last two digits together moves `hi` and `s` in step; in `33x` alone `s` moves twice as fast as `hi`.
Weakness: slides are n = 13 only (`theorist_slides_n13.jsonl`); the 33x fits at 14..17 use orbits only, which is all (3) needs.

## Reproduction

```
.venv/bin/python workshop/rounds/004/theorist_shortfall.py > workshop/rounds/004/theorist_shortfall_n13.txt   # 1 s
.venv/bin/python workshop/rounds/004/theorist_k.py > workshop/rounds/004/theorist_k_n13.txt                    # 1 s
for n in 14 15 16 17; do timeout 9m .venv/bin/python workshop/rounds/002/experimentalist_census.py $n --cores 33,333,334,335,336 --out workshop/rounds/004/theorist_33x_n$n.jsonl & done; wait   # ~2 min in parallel
.venv/bin/python workshop/rounds/004/theorist_33x.py > workshop/rounds/004/theorist_33x.txt                    # 1 s
```

## Prior record

E-060 (rules of H-020/`d = t - h`, 13/13 and 17/21), F-053 (`45`), E-056 (109/139). The column and its counts reproduce E-060 and close
the referee's item. The `33x` law `k = 2x` is not in `research/` (grepped `k(c)`, `shortfall`, `33x`); `HYPOTHESES.md` H-021 leaves `k(c)` open. No overlap with RETRACTIONS.

## Code changed

New only: `theorist_shortfall.py`, `theorist_k.py`, `theorist_33x.py` (+ output files). No library change, no tests run.

## Next

- experimentalist: `k(c)` for other three-letter prefixes at n = 14..17 (`34x`, `44x`, `45x`, `4033`-type) to see whether `k` is linear in a last digit generally; and the 4-core parity cases `5046 5056` at 14 (is `{0,2},{1,3}` a translation of period 2 for all n? If so H-021 does not apply).
- skeptic: is the `|R| <= 4` fit powerless? A null test would randomise the merge structure.
- theorist: find why `k = 2x` from the rule table (T4): the top `x - 3` offsets are exactly the singletons, so they might be placements whose reduced row has no partner (edge effect of the sink end).
