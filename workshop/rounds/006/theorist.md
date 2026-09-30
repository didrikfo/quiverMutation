# k(33x) = 2x is a conservation law of the double mutation plus one self-dual seed; the argument does not extend to 44x, and explains 34x only at x = 5

author: theorist · round: 006 · kind: result (mechanism, one step unproved) + negative for the extension
thread: T2/T4 · bears on: H-021, H-020, F-051, F-032, E-061

## Claim

**(1) Not the rule table.** For `x >= 5` the floating half of the rule table does nothing to `33x` in the interior (orbit of size 1 at
n = 18, offsets 4-6, x = 5, 6). For `x = 4` the table has the two rules `334 <-> 333` (size-4 orbit). What moves `33x` is the double
mutation of arXiv:2310.08346 (F-032), so H-020/F-051's "the rule table acts the same at every interior position" is true for `33x` only
through `doubleMutation`, which is interior-uniform by construction.

**(2) The derivation.** Write `33x@o` for the intervals `A=[a,a+3], B=[a+1,a+4], C=[a+2,a+2+x]`, `a = o+1`.
(i) *Drift.* `L_t` on `B` (a relation starts at `a`; none starts at `t-1 = a+3`) sends `A` to `[a,a+4]` (contains `B`, dropped),
`C` to `[a+3,a+2+x]`, adds `[a+2,a+5]`: result `B, N, C'` = **`33(x-1)@(o+1)`** for every `x >= 4`, every `n`, every interior `o`.
The dual `R` gives `33(x+1)@(o-1)`. So **`c = x + o` is conserved** (the right end of the long relation is fixed) and the chain
`S_c = { 33y@(c-y) : 3 <= y <= c }` lies in one orbit. (Checked against `doubleMutation.rewritesOf` for all 33 placements at n = 14: 0 missing.)
(ii) *Seed.* `333` is self-dual, and `R` on `333@m` gives `4033@m`, the mirror of `334@(n-7-m)`. So the chain `S_c` (`m = c-3`) is joined to the
mirror of the chain `S_{n-c}`: the conserved label is identified under `c <-> n - c`.
(iii) *Consequence.* `33x@o` and `33x@o'` share an orbit iff `o + x` and `o' + x` are `{c, n-c}`, i.e. `o' = o` or `o' = n - 2x - o`.
That is the reflection with **`s = n - 2x`, `k = 2x`**. It needs `o' >= 0`, so offsets `o > s` have no partner: their number is
`hi - s = (n-x-3) - (n-2x) = x - 3 = d`. **`d = x - 3` is the count of offsets where the second chain no longer has room for `33x`.**
Both E-061 numbers follow; `hi = n - x - 3` is the footprint `x + 3`.
(iv) *Generalisation of the formula.* For any drift family `F_x` (footprint `x + w0`) with a self-dual member `x0`:
`k = 2x + w0 - x0`. For `33x`: `w0 = x0 = 3`.

**(3) Status of the proof.** Lower bound "orbit contains `S_c` and the mirror of `S_{n-c}`" = steps (i), (ii), by hand, plus the code check.
**Weakest step, not derived:** (a) that the *mirror* chain (rows `D(y)@m`) is joined to the 33-rows `33y@(n-c-y)`, i.e. the end link
(at n = 13, x = 4 the walk does `D(9)@1 -> 0(10)000000030 -> 3(10)0.. -> 339@0`, anchored/edge moves, `theorist_path.py`);
(b) the upper bound "nothing else is in the orbit" (only computed). Both are verified, not proved, below.

**(4) Extensions, stated so they can fail.**
* **`44x`** has the same drift (`445@6 -> 446@5, 444@7`, `theorist_nbrs.py`) and the self-dual seed `444`, so the argument *predicts*
  `k = 2x - 1`, `d = x - 4`. **This is false:** at n = 14, 15, 16 every `44x`, x = 4..8 (all offsets), lies in one orbit, the one of `333@0`
  (sizes 3767, 5648, 8134). The chain is only a lower bound; `33x` is exact because nothing else attaches, `44x` attaches to the `c = 3` orbit.
* **`34x`**: `L` on the middle relation sends it to `444(x-1)@o`, not to `34(x-1)`: no drift inside the family. The data: `k = x + 3 = w` (footprint),
  `d = 0`, i.e. `34x@o ~ 34x@(hi-o)`, for x = 5 (n = 14, 15, 16, all offsets), x = 7 (n = 14, 15 offsets 0-2 of 5; 16 offset 0) and 8 (n = 14), and `344` at
  even n (14, 16). `346` never pairs (one orbit with `333@0`, n = 14..16). **Only `345` is explained**: `334@o -> 4444@o <- 345@o`
  (both one double mutation, n = 20 o = 6), so `345@o` carries `c = o + 4` and pairs with `o' = n - 8 - o = hi - o`: `k = 8` is the `33x` law at `x = 4`
  shifted by the footprint. `347, 348, 344` carry no `33y` row (empty c-set): their reflection is observed, not derived.
* **`45x`**: no drift neighbour (`455 -> 5504, 35@68, 5055`). Observed: `456` and `457` merge (into the `333@0` orbit and the `347@0` orbit), `455` merges at 15 and
  at 16 is a **parity** translation `{0,2,4,6,8}, {1,3,5,7}` (the `4046` pattern of E-061), `458` at 15 is `{0,2,4}`. Predicted by the argument: no reflection. Consistent.

## Evidence

Exact-membership test (`theorist_chain.py`): for each `33x@o`, x = 3..8 (x = 9 excluded: the word string "3310" is misread, the script caps at 9), every offset, the
walk (reduced, whole table, closed in every case) is compared with the **prediction**
"the 33y@p rows in the orbit are exactly those with `y + p in {x + o, n - x - o}`":

| n | placements | `33y@p` held == predicted | closed | orbit sizes |
|---|---|---|---|---|
| 14 | 39 | 39 | 39 | 320 .. 3767, equal for equal `{c, n-c}` |
| 15 | 45 | 45 | 45 | |
| 16 | 51 | 51 | 51 | |

135/135, no deviation, for `x <= 8`; includes the singletons `o > s`. This is a statement about the forward orbit of the reduced walk (closed), not derived equivalence.
Table for the extensions (`theorist_link.py`, held offsets of the same word; c-set = `c = y + p` of the `33y` rows in the orbit; `M` = one orbit of all offsets):

| word | n = 14 | n = 15 | n = 16 |
|---|---|---|---|
| 345 | {0,6}{1,5}{2,4}{3}; c {o+4, 10-o} | {0,7}{1,6}{2,5}{3,4} | {0,8}{1,7}{2,6}{3,5}{4} |
| 346 | M (3767) | M (5648) | M (8134) |
| 347 | {0,4}{1,3}{2} | {0,5}{1,4}{2,3} | {0,6}, rest unfinished |
| 344 | {0,7}{1,6}{2,5}{3,4} | {0,8}{1,7}{3,5}, 2,4,6 single | {0,9}..{4,5} |
| 44x, x=4..8 | M (3767) | M (5648) | M (8134), x = 4..6 done |
| 455/456/457 | not run | M / M / M (2290) | {0,2,4,6,8}{1,3,5,7} / M / M |

What this does not cover: `x >= 9`; `n >= 17`; `n = 13` fits were not redone; the `34x`, `45x` rows at n = 16 for 347, 348, 457, 458 are incomplete
(one job still running when I stopped). Which n informs: the pairing tests have power only when `s = n - k` leaves at least 2-3 pairs, so **n = 16 and 17**
(not 14: `336`, `337` have 1-2 pairs); even and odd `n` both matter because `344` and `455` break parity. A test that would kill the `34x` law
`k = x + 3`: any `34x`, x = 5, 7, 8, at n = 17 whose orbits are not `{o, hi - o}`; for `x = 9`: `349` at n = 16, 17 (`hi = 4, 5`).

## Reproduction

```
.venv/bin/python workshop/rounds/006/theorist_rule.py                        # rules with LHS 33x in the table (1 s)
.venv/bin/python workshop/rounds/006/theorist_step.py 18 5 6                 # one-step double / table neighbours of 335@6 (1 s)
.venv/bin/python workshop/rounds/006/theorist_local.py 18 5 6                # floating-only orbit size 1; doubles alone give 2482 (20 s)
.venv/bin/python workshop/rounds/006/theorist_path.py 13 4 0                 # shortest path 334@0 -> 334@5, rule by rule (1 s)
timeout 10m .venv/bin/python workshop/rounds/006/theorist_chain.py 14        # the 135-placement test (52 s at 14; ~4 and ~6 min at 15, 16 in parallel)
.venv/bin/python workshop/rounds/006/theorist_nbrs.py 44 4-8                 # drift neighbours of 44x (1 s)
timeout 10m .venv/bin/python workshop/rounds/006/theorist_link.py 14 34 4-8  # held offsets + 33y c-set (1-4 min per n)
```
Outputs kept: `theorist_chain_n14.txt`, `_n15.txt`, `_n16.txt`, `theorist_link_{34,44,45}x_n{15,16}.txt` (n = 14 34x/44x in this file's table).

## Prior record

E-061 states `k(33x) = 2x`, `d = x - 3` as a description (round 004, mine); "why" was open. F-032 has the double mutation but not this use of it; F-053/E-052 have
reflection pairs for `45`; E-059 has parity for 7 cores (the `455` parity at 16 is new data, in the `4046` family, unexplained). Grepped `RETRACTIONS.md`, `conserved`, `drift`:
nothing. The claim "it is the rule table" from the assignment is corrected in (1). `S_c` / `c <-> n - c` is new.

## Code changed

None in the library. New scripts in `workshop/rounds/006/`: `theorist_{rule,step,local,path,chain,nbrs,words,link}.py`. No tests run (no library file touched).

## Next

* theorist: prove the end link (a) (a chain-top `33c@0` joined to its mirror `D(c)@(n-c-3)`) and the upper bound by listing every move out of the orbit rows for one `(n, c)` (orbit size 310 at n = 16).
* experimentalist: `33x` at n = 17, x = 7, 8 (exact-membership test of `theorist_chain.py`, edit its cap 9); `34x`, `45x`, x = 5..9, n = 16, 17 against the predictions above (`k = x + 3`, `d = 0` for `34x`; no reflection for `45x`); the n = 15 `344` odd-n singletons (E-059 parity).
* skeptic: `44x` merges into the `333@0` orbit: is that orbit an end-touch artifact (the `c = 3` chain has `33c@0` at the source) that also swallows the `34x`/`45x` members at x = 6? `346`, `456`, `446` all land there.
