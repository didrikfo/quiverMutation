# H-021 restated without the mirror clause: the reflection centre is the sum of the first and last outside offset, d(c) = tail - head, when the outside block is interior

author: theorist · round: 003 · kind: result (restatement + partial mechanism)
thread: T1, T2 · bears on: H-021, H-020, F-053, E-052, E-058

## Claim

**Restatement (T1).** For a single-cluster core `c` and length `n`, let `R(c,n)` be the offsets with a row, `lo = 0`,
`hi = max R`. Say `c` *pairs at n* if there is an integer `s` (centre) such that the reduced-walk orbit partition of
`{c@o : o in R}` is closed under `o -> s - o` wherever `s - o in R`, with at least one orbit holding some `o != s - o`
together with `s - o`. Its *shortfall* is `d(c) = hi - s` (signed; `s` is the fit with least `|d|`, `|d| <= 6`).
**H-021':** for every core that pairs, `s(c,n) = n - k(c)` with `k(c)` independent of `n`, and `d(c)` is independent of `n`.
No clause about a mirror in the orbit. What H-021' says about non-pairing cores: nothing (30 of 139 at `n = 13`, 20 with a
loose mirror; they are not "exceptions", they are outside the statement). At `n = 13`: pairing 109 of 139 (E-058, unchanged).
Not tested by me: `n`-independence beyond 13 -> 14 for 12 cores (all 12 fit again, `s` up by 1, `d` unchanged; table below).

**Mechanism for `d` (T2), partial.** Write the slide (H-020) as `i^h o^m i^t` (`h` = head, `t` = tail), `m >= 1`.
(a) *Proved, trivial:* the verdict is constant on an orbit (E-052). So a centre `s` that pairs is slide-consistent, and if the
reflection maps the outside block `[h, h+m-1]` onto itself then `s = 2h + m - 1` = (first outside offset) + (last outside
offset) `= hi + h - t`, i.e. **`d(c) = t - h`** (signed; `d > 0` means the reflection overhangs the sink end, `d < 0` the source end).
(b) *Data:* the fitted `s` equals first-outside + last-outside for **13 of 13** cores whose outside block touches neither end
of the slide (`n = 13`), and 3 of 3 at `n = 14`; for cores whose outside block touches offset 0 or `hi` it holds for 17 of 21 at
13 and 4 of 7 at 14. The 4 (+3) exceptions `4045 3556 4506 4556` (at 14 the same except `3556`, not run) are cores where the block touches
an end, so the slide-consistent centres are not unique (e.g. `4045 = oooii`: `s = 1` and `s = 2` are both slide-consistent) and the
orbits choose the other one, by 1. So the H-020 head/tail difference **is** the shortfall when the outside block is interior,
and is not forced when it touches an end.
(c) *What is not explained:* the 62 cores whose slide is all inside carry no verdict information, so `d` is not visible in
head/tail there; `|d|` is 0 for 45 of them, 1 for 13, 2 for 3, 3 for 1 (`336`; `33x` gives `|d| = x - 3`). The 30 cores
with no fit. And why the o-block folds onto itself at all (that is F-053 again, not derived). Nothing here says why `k(c)` has the values it has; only that it is `w(c) + t(c) - h(c)`
with `w = n - hi` constant, which trades T2's unknown for H-020's two numbers.

## Evidence

Slides were computed here (the census had orbits only): for each recorded orbit the verdict of its least offset by
`batch._verdictFor(n, c, o, free=REDUCED)`, 395 verdicts at `n = 13` (about 5 min, 4 procs), 12 cores at 14. No undecided, no
core with a mixed-verdict orbit. Table (`n = 13`, cores with a fit, slide class):

| slide class | cores | `s = first o + last o` (or `hi` for no `o`) |
|---|---|---|
| mixed, outside block interior | 13 | 13 |
| mixed, block touches an end | 21 | 17 (fails: `4045 3556 4506 4556`) |
| all outside | 13 | 10 (fails: `4046 5046 5056`, `d = 1`, orbits `{0,2}`) |
| all inside | 62 | not testable (slide silent) |

`n = 14`, 12 cores: `45 46 504 3344 3355 3445 3444 3345` pass; `4556 4045 4506` fail as at 13; `334` all inside, `d = 1`
as at 13. Data: `theorist_shortfall_n13.txt` (per-core), `theorist_rule_n14.txt`. Weakness: the 14 sample is 12 chosen
cores (mixed slides and the 13 failures), not a census; `d <= 6` fit is generous (only failures of the rule are strong).

## Reproduction

```
for k in 0 1 2 3; do timeout 10m .venv/bin/python workshop/rounds/003/theorist_slides.py $k/4 logs/th-s$k.jsonl & done; wait  # ~5 min; cat -> theorist_slides_n13.jsonl
.venv/bin/python workshop/rounds/003/theorist_shortfall.py > workshop/rounds/003/theorist_shortfall_n13.txt   # seconds
.venv/bin/python workshop/rounds/003/theorist_rule.py workshop/rounds/003/theorist_census_n14_sample.jsonl workshop/rounds/003/theorist_slides_n14_sample.jsonl   # seconds
```
The `n = 14` census sample: `experimentalist_census.py 14 --cores 45,504,46,3344,3444,3355,3445,4506,4556,4045,334,3345` (about 20 s per core).
The "first + last o" split by touching an end is a 12-line snippet, reproduced by `theorist_shortfall.py` output columns `s hi h t interior`.

## Prior record

F-053 says the slide of `45` is a palindrome on `0 .. n-8` and the offset past it is the difference: this is the same statement for one
core (`d = t - h = 1`); the general form (`d = t - h` whenever the outside block is interior, with 13/13 and 3/3) is not in `research/`
(grepped `shortfall`, `overhang`, `first outside`). H-021 asks "is `d(c)` the head/tail difference" (untested until now). The
`s = first o + last o` reading is close to a tautology (a); the content is that the block is folded onto itself in 13/13 cases and
that the exceptions are exactly a block at an end. No overlap with RETRACTIONS.

## Code changed

None to the library. New: `theorist_slides.py`, `theorist_shortfall.py`, `theorist_rule.py`, data `theorist_slides_n13.jsonl`,
`theorist_slides_n14_sample.jsonl`, `theorist_census_n14_sample.jsonl`, `theorist_shortfall_n13.txt`, `theorist_rule_n14.txt`. No tests.

## Next

- experimentalist: the same 12 (and the 4 failures) at `n = 15` to test the `n`-shift of `s` and `d`; slides for the n = 14 census when it runs.
- skeptic: is the 13/13 interior-block result an artefact of small blocks (`m` is 3 to 5)? Find an interior-block core with `m >= 6` at `n >= 15`.
- theorist next: why an all-inside core has `d > 0` (`33x`); the orbit-side reason the block folds (T4).
