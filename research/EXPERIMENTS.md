# Experiments

Runs made, newest first, with parameters and outcome — including runs that found
nothing, which are recorded precisely so they are not repeated. See
[`README.md`](README.md).

---

## E-099 — At n = 8, every one of the 365 LNAs with a relation of at least 3 arrows has a cord member within 3 mutation steps and none of the 64 with only 2-arrow relations (or none) has one at depth <= 5; no monomial cord at depth 3 for any of the 429
*2026-10-01* · **`workshop/rounds/021/maverick_predict.py 8 3 LO HI` over the four shards (indices 0-428; `maverick_criteria.py` on `maverick_predict_L3_s*.txt`: 365 predicted and found, 64 predicted not and found not, 0 mismatches); the 64 are exactly the sequences over {0,2} (2^6); L = 5 for those 64: 0 cord members (`maverick_predict_L5_zero_s*.txt`); `MONO=1` at L = 3: 0 of 429. Minimum cord depth over the 365: 1 for 347, 2 for 17, 3 for 1 (LNA 205, `223022`).** · *workshop round 021, maverick, refereed by theorist*

**Limits.** An observed fit at n = 8, L = 3 (positives) and L <= 5 (negatives), not a theorem; LNA 205 reaches its first cord at depth 3, so L = 3 is the edge and no L = 4 check was made. The 64 are cordless only to depth 5: E-087 finds n = 8 sum-cord controls at depth 6, so "cordless" is not established. `MONO` has no positive control at n = 8 (E-092), so 0 of 429 is not a measured rate. The heuristic (a relation of >= 3 arrows becomes a sum relation) is unproved and does not explain why 18 LNAs need more than one step; the proposed proof route for the rad^2 = 0 direction is unchecked. Refines E-087's "a 3 early in the sequence" to relation length. Referee re-ran 20 LNAs (match) and the criterion table on the committed shards; did not re-run the L = 5 negatives. The proposed non-MONO L = 5 overnight walk of the 365 (about 4.5 CPU-hours) is not added to `OVERNIGHT.md`: the criterion already says which 365 carry cords.

**Reproduction.** `for i in 0 1 2 3; do timeout 10m .venv/bin/python workshop/rounds/021/maverick_predict.py 8 3 $((i*108)) $((i*108+107)) > /tmp/p$i.txt; done` (about 3 min per shard); `.venv/bin/python workshop/rounds/021/maverick_criteria.py workshop/rounds/021/maverick_predict_L3_s*.txt`.

---

## E-098 — Over all nondecreasing words with 4-6 letters in 2..9 at n = 12..14, "exactly one placement in the `444` orbit" is not enriched over a pooled-rate binomial null and occurs for words without a 4 as often or more often; the right gap g <= 1 of that placement beats a uniform-position null
*2026-10-01* · **`workshop/rounds/021/skeptic_null.py N` (N = 12, 13, 14, seconds each) and `skeptic_null_gap.py N`. Exactly-one counts for 4-letter words: with a 4 6/9/12, without a 4 6/13/22 at n = 12/13/14 (pools 20/35/56 and 15/35/70); binomial expectation at the stratum rate 7.4/11.8/20.3 and 6.0/14.1/27.1. Exactly-one and g <= 1: n = 14, 11 of 12 (uniform null 4.5) with a 4, 19 of 22 (9.1) without; n = 13, 8 of 9 (3.6); 6-letter words 7/10 (p = .09) with a 4, 7/8 (p = .02) without.** · *workshop round 021, skeptic, refereed by experimentalist*

**Limits.** The words are exhaustive, not sampled; p-values are indicative (words share placements; n = 12..14 are not independent). E-091's 6/9/12 reproduce exactly, so its claim is untouched; what is shown is that it does not separate "has a 4" from "four-letter word" (no-4 counts are larger, partly because the pool is larger). "Not enriched" is relative to the pooled-rate binomial null only. The submission's "same counts" was wrong and is corrected here. `|S|` is 3767 at n = 14; orbit closure (versus the 300000 cap) was not stated by the author or checked by the referee. n >= 15 not run, although a run takes seconds. The g <= 1 concentration describes the orbit's right-end shapes, not 444-flavoured words (the null ignores that placements near the right end are likelier to be in S for any reason). Reads with E-079 (4 versus collapse to `34` inseparable).

**Reproduction.** `for n in 12 13 14; do timeout 10m .venv/bin/python workshop/rounds/021/skeptic_null.py $n; .venv/bin/python workshop/rounds/021/skeptic_null_gap.py $n; done` from the repository root.

---

## E-097 — The rewrite's (k,i) entry is dim coker g_i (from steps 4 and 6) and, given that step 7 returns generators of the whole ideal out of k*, its (i,k) entry is dim ker psi_i; the Cartan congruence then fails exactly when some dim ker g_i != 0; rejecting parents are not all strict-A5
*2026-10-01* · **`workshop/rounds/021/scholar_step7_entries.py`: `--e078` (75 tilting / 6 non-tilting), `5 --class 0` (16 620 steps, 0 non-tilting), `6 --class 0 --budget-sec 500` (99 267 + 1 123 distinct rejecting parents), `7 --class 0 --budget-sec 500` (56 365 + 156); both entries compared on every step, 0 violations. Shapes of the rejecting parents: strict A5 square 767 of 1 123 at n = 6 (156 of 156 at n = 7, strict test only), every one of the 1 123 has a relation of two paths from one vertex ending x, v, e with v having one out-arrow.** · *workshop round 021, scholar, refereed by skeptic*

**Limits.** A derivation sketch for (a) and (b) checked on the ranges above, not a proved theorem. (a) uses e_iBe_{t alpha} = e_iAe_{t alpha}, which the text does not justify (it is tested by the script). (b) is conditional on step-7 completeness, which no script tests separately from the final dimension. Nothing is said about the diagonal or rows/columns i, j != k (E-093 observes no difference), nor about whether B is the true End(T). The shape check is one-sided: no control of how many tilting parents have the long square. The `hasA5` test needs no relation and `hasLongSquare` is loose, so "neither: 0" is weak. Counts at n = 6, 7 depend on the cap (E-095: 907 / 143). "A5-shaped" in E-084 and E-095 is to be read as this long-sided square (356 of 1 123 at n = 6 fail the strict test). Referee re-ran `--e078`, `5 --class 0` (identical) and `6 --class 0 --budget-sec 120` (116 rejecting parents, 0 violations).

**Reproduction.** `.venv/bin/python workshop/rounds/021/scholar_step7_entries.py --e078` (2 s); `timeout 10m .venv/bin/python workshop/rounds/021/scholar_step7_entries.py 5 --class 0` (60 s); `... 6 --class 0 --budget-sec 500`.

---

## E-096 — For the 18 split four-letter words of E-091 the one placement in the `444` orbit has a fixed right gap g (0 or 1) at n = 12..17 for 17 words, `3344` is the exception; lemma R alone reaches `333@0` from none of 45 word-n cases
*2026-10-01* · **`workshop/rounds/019/theorist_gaps.py n` (n = 12..16; referee also n = 17): g = 0 for `2224 2334 4556 4667 4778 4889`, g = 1 for `224x` (x >= 5) and `344x` (x >= 4), g independent of n; `3344` does not fit (right gap 7 at n = 16, 8 at n = 17; in-orbit placement at left gap 1). `theorist_rchain.py n` (n = 12, 13, 14, 16; 6+9+12+18 = 45 cases): lemma R, with every step filtered against `doubleMutation.rewritesOf`, stops at boundary shapes `4`, `4y`, `334`, `3 b (d+1)` (`357 368 379 38(10)`) that are themselves in the orbit, never at `333@0`; the BFS paths printed at n = 12, 14 (`theorist_split.py`) go through the `34@k <-> 403@(k-1)` shuttle (7..13 steps). R is valid only when the shortened interval does not swallow a relation on its left (the `3344` case; E-088's R test covered isolated runs).** · *workshop round 019, theorist, refereed by skeptic*

**Limits.** A relabelling of E-091's "last or second-to-last for 17 of 18": the 18 words are E-091's list, so the evidence is not independent of it; the new content is the n-independence of g, the R-terminals and "R alone does not reach `333@0`". `3344` is 17 of 18 plus an exception, not a rule at every n (its left-gap reading is a relabelling, not derived). Why only one placement per word (short shapes are in the orbit only at the boundary) is read from tables for n = 12..16, not proved. R-chains: n = 15 not run; the shuttle claim rests on n = 12, 14 BFS output. `3334`, `2455` not checked.

**Reproduction.** `for n in 12 13 14 15 16; do timeout 10m .venv/bin/python workshop/rounds/019/theorist_gaps.py $n; done` (about 15 s each); `for n in 12 13 14 16; do timeout 10m .venv/bin/python workshop/rounds/019/theorist_rchain.py $n; done`; `theorist_split.py 12`, `theorist_shapes.py`, `theorist_jlabel.py 12` likewise.

---

## E-095 — On the guarded walks at n = 5..7 and the E-078 family, dim ker is 1 at exactly one vertex of each of 1 050 rejecting parents (907 at n = 6, 143 at n = 7, all distinct), never above 1 per vertex, and no non-tilting step has dim ker 0
*2026-10-01* · **`workshop/rounds/019/experimentalist_kerhist.py`: `--e078` (75 tilting / 6 non-tilting; totals 3/2/1 occur only here), `5 --all` (30 300 steps, 0 non-tilting), `6 --class 0 --budget-sec 480` (83 591 tilting + 907 non-tilting), `7 --class 0 --budget-sec 480` (53 502 + 143); 0 disagreements between `tiltingPlus` and the Cartan congruence; per-vertex dim ker <= 1 throughout.** · *workshop round 019, experimentalist, refereed by theorist*

**Limits.** Fills the gap named in E-093's Limits (histogram of dim ker, distinct parents). "Tilting <=> dim ker 0" is E-093's identity read on more steps, not new content. "Per-vertex dim ker <= 1" is about the walks' parents (all A5-shaped, E-084; no shape check made), not the gate; "one bad vertex per parent" and "dim ker per vertex" are different quantities. n = 6, 7 counts depend on the 480 s cap (n = 6 907 here vs 696 in E-093, attributed to load, unchecked); n = 6 classes 1-3, n = 7 classes 1+ not run.

**Reproduction.** `.venv/bin/python workshop/rounds/019/experimentalist_kerhist.py --e078` (3 s); `... 5 --all` (82 s); `timeout 10m ... 6 --class 0 --budget-sec 480`; `timeout 10m ... 7 --class 0 --budget-sec 480`.

---

## E-094 — With checkpoint/resume the n = 8 class 2 depth-8 guarded walk completes under the fixed library: 24 316 expansions, 63 221 distinct algebras, 2 rejections, 0 key-moved steps; the 10 key-moved steps of E-084 are equal in count to the steps that now keep the key
*2026-10-01* · **`workshop/rounds/019/toolsmith_walk.py 8 --class 2 --depth 8 --budget-sec 480 --ckpt FILE` (two slices, 493 s + 208 s; checkpoint 38 MB kept outside the repository); over E-084's first 20 899 expansions: guard-tilt 89 189 now vs 89 179 + 10 key-moved before, noguard-NOTtilt and the rejection line identical, distinct 54 333 vs 54 326; the second rejection (path (17, 8, 5, 6, 8, 8, 2, 5), vertex 5, gate True, `tiltingPlus` False, guard-refused) not inspected. Resume equals uninterrupted `scholar_walk.py` byte-for-byte at n = 7 classes 0, 1 (depth 6) and, by the referee, class 2 (depth 5).** · *workshop round 019, toolsmith, refereed by scholar*

**Limits.** Closes E-090's loose end only for the n = 8 class 2 guarded walk, full depth 8 (the frontier is not closed; depth 9 not run), under the E-089 fix; it says nothing for other classes. The 10 key-moved steps are matched to the new steps by count (89 179 + 10 = 89 189), not by replaying the 10 E-084 parents; the +7 distinct algebras is unexplained. n = 8 resume equivalence rests on count agreement with E-084 (depth 1-7 lines identical), not on an n = 8 slice-vs-uninterrupted run. Referee did not re-run the depth-8 walk (about 11 min).

**Reproduction.** See E-090 for the uninterrupted walker; `timeout 10m .venv/bin/python workshop/rounds/019/toolsmith_walk.py 8 --class 2 --depth 8 --budget-sec 480 --ckpt /tmp/tw_n8c2.ckpt` (run twice); resume check: `workshop/rounds/019/toolsmith_walk_resume_check.txt`, recipe in `workshop/rounds/019/toolsmith.md`.

---

## E-093 — On every gate-admitted step tested (guarded BFS from the LNAs, n = 5..7) the Cartan congruence `R C R^T = Cartan(child)` fails exactly where `tiltingPlus` fails, and then the discrepancy is row k, off the diagonal, equal to minus the kernel dimension of the map `p -> (p beta)_beta`; no step disagrees
*2026-10-01* · **`workshop/rounds/018/scholar_cartan_vs_tilt.py` and `scholar_e078_diff.py`: 18 E-078-family algebras (75 tilting+congruent steps, 6 non-tilting and non-congruent), n = 5 both key classes closed (11 700 algebras, 30 300 steps, 0 non-tilting), n = 6 class 0 (65 914 / 696, stopped at 420 s), n = 7 class 0 (39 901 / 111, stopped at 420 s); 0 steps with tilting and congruence disagreeing; in all 807 non-tilting steps X - Y is row k off-diagonal and equals minus dim ker g_i (E-078 at n = 5, vertex d: a single entry -1 at (d, a)). All non-tilting steps are guard-refused.** · *workshop round 018, scholar, refereed by experimentalist*

**Limits.** An observation on guarded-walk steps, not a theorem: the walks are not independent (all rejecting parents of the A5 shape, E-084), and the distribution of dim ker g_i and the number of distinct parents were not reported, so "equals minus dim ker" may be shown only for dimension 1. The n = 6, 7 counts depend on the 420 s wall-clock cap (referee's n = 7 rerun: 24 187 algebras, 44 241 / 120, same pattern, 0 disagreements). The `illegal` and `dim` buckets of the script were 0. The derivation (`chi(P_i, T^+_k) = dim coker g_i - dim ker g_i`) is a sketch and its step "the rewrite's (k, i) entry is dim coker g_i" is read off the data, not from step 7 of `procedure.mutateAtVertex`; the diagonal and column-k entries are not covered; no independent End(T). The n = 8 class 2 steps of E-085 are not re-run (E-090). Consequence, as far as it survives: E-085's congruence is not a second, independent piece of evidence against non-tilting steps; on tilting steps it checks the rewrite, which `tiltingPlus` cannot.

**Reproduction.** `.venv/bin/python workshop/rounds/018/scholar_cartan_vs_tilt.py --e078` (3 s); `... 5 --all` (86 s); `... 6 --class 0 --budget-sec 420`; `... 7 --class 0 --budget-sec 420`; `.venv/bin/python workshop/rounds/018/scholar_e078_diff.py`.

---

## E-092 — At n = 8 (six LNAs, path length <= 5) all 2376 cord members carry a sum (commutativity) relation and no member is monomial; `MONO=1` finds no monomial cord member for 6 of 14 LNAs with cords at n = 5 (depth 7) and 1 of 5 at n = 4 (depth 8); the filter is positive on a hand-built monomial cord seed
*2026-10-01* · **`workshop/rounds/018/maverick_producer.py 8 5 4 9 10 11 12 13` (2376 members: 2298 with one cycle, 78 with two; sum-relation shapes 2+2 1642, 3+2 358, 2+3 190, 4+2 114, ...; 0 monomial-only); `MONO=1 workshop/rounds/015/toolsmith_cords.py` at n = 4 L = 8 and n = 5 L = 7 (0 members); `maverick_filtercheck.py 5` (4 members from a hand-built monomial seed); `maverick_monocord.py 5`.** · *workshop round 018, maverick, refereed by skeptic*

**Limits.** The n <= 5 negatives are informative only for the LNAs that have cords at all (6 of 14 at n = 5, 1 of 5 at n = 4); the rest have none under any filter. n = 8 is depth 5, shallower than the depth 6 at which E-087 first found its n = 8 sum-cord controls, so "no monomial cord near the LNAs" is not excluded at greater depth. "Cycle covered by a sum relation" is a union-of-supports test and does not show the cycle is one relation's commutativity cycle (78 two-cycle members). The heuristic "a cord needs a sum relation" is still unproved and still has no positive MONO control at n = 8; the Coxeter polynomial cannot test it (15 of 42 and 190 of 736 monomial unicyclic algebras at n = 4, 5 share a polynomial with an LNA). The filter control seed is not LNA-like ("ILLEGAL RELATION" on several mutations). The E-076 candidates are monomial quipus with cords, the control members have sum cords: the two differ in kind (speculation).

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/018/maverick_producer.py 8 5 4 9 10 11 12 13` (2.5-5 min); `MONO=1 timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 5 7 0 1 0 100 --plan`; `MONO=1 timeout 5m .venv/bin/python workshop/rounds/018/maverick_filtercheck.py 5`.

---

## E-091 — At n = 12..17, among 4-letter words (letters <= 9, a 4, at least 4 placements) 6, 9, 12, 15, 18, 18 are split across the `444` orbit, each with exactly one placement in it; `3334` and `2455` have no placement in it; E-086's "0 partial" is vacuous for merged words
*2026-10-01* · **`workshop/rounds/018/skeptic_rowset16.py N 4 0`, N = 12..17 (|S| = 1410, 2386, 3767, 5648, 8134, 11340; words 20, 35, 56, 84, 120, 120; all-in 5, 10, 13, 16, 19, 19; all-out 9, 16, 31, 53, 83, 83; split 6, 9, 12, 15, 18, 18; about 12 s each); n = 16 split words 2224 2245-2249 2334 3344 3444 3445-3449 4556 4667 4778 4889.** · *workshop round 018, skeptic, refereed by theorist*

**Limits.** The quantifier is restricted to letters <= 9 and words with at least 4 placements; the 120 plateau is the number of such words once n >= 16, not a bound on a sum. The n = 16 and n = 17 split lists are the same family, so they are not independent evidence. The one placement in the orbit is the last or second-to-last for 17 of the 18 (`3344` the exception): descriptive, no rule derived. The walked n = 16 run covers 55 words (about 46% of the list) and agrees with the fast run; the other placements of split words are not classified. The 446-orbit of `3334` is identified with E-088's `235`-type orbit by size only. "Merged" words have all placements in one orbit, so "all in S or none" is true for them by definition (E-086's content is the in-S count).

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/018/skeptic_rowset16.py 16 4 0` (12 s); `.venv/bin/python workshop/rounds/018/skeptic_partial_where.py 16`.

---

## E-090 — With the full-reduction `reduceAgainstPivots` the E-084 walks give the same counts as the E-084 files in every class that could be re-run (n = 7 classes 0-2, n = 8 classes 0-1, n = 9 classes 0-3); the 10 n = 8 class 2 "key moved" steps lie past the part of depth 8 that fits in 10 minutes, so the loose end is neither reproduced nor refuted by the walk
*2026-10-01* · **`workshop/rounds/017/experimentalist_walk.py` (the E-084 walk `workshop/rounds/014/scholar_walk.py` run with the E-085 full reduction monkeypatched in, same arguments, `--stop-on-reject`, 540-585 s budgets): rejections (gate True, `tiltingPlus` False) 8/4/6 at n = 7, 4/14 at n = 8 classes 0/1, 8/4/12/12 at n = 9; the "guard-admitted and failing" row absent everywhere (about 2.2e5 guard-admitted steps in the re-run walks); summary lines identical to the E-084 files (referee diffed the saved outputs and re-ran n = 7 class 0). n = 8 class 2: 1 rejection, 0 key-moved steps, but stopped at 17 058 of 20 899 depth-8 expansions, before the 10 steps of E-084 (past expansion 16 500).** · *workshop round 017, experimentalist, refereed by scholar*

**Limits.** Only the comparison against the saved E-084 files counts: the "unpatched control" of n = 8 class 2 (cap 16 500) and the "1.3 x time" figure are void, since the library was patched by the toolsmith in the same round (E-089) and the control ran the patched library (referee). The "3.2e5 guard-admitted" of the submission is about 2.2e5 by its own table. Every walk is stopped on rejection or budget, as in E-084; nothing at n >= 10. The loose end of E-084 still rests on E-085's replay of the 10 saved parents. A full depth-8 run (about 15 min, one process) is not made; the experimentalist proposes it for `OVERNIGHT.md` (not added: it fits one 10-minute-limit command only with a checkpoint, toolsmith first).

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/017/experimentalist_walk.py 7 --class 0 --stop-on-reject --budget-sec 540` (seconds); `... 8 --class 0 --stop-on-reject --depth 8 --budget-sec 570` (about 1 min); n = 8 class 2: `... 8 --class 2 --stop-on-reject --depth 8 --budget-sec 585` (partial); outputs `experimentalist_walk_n{7,8,9}_c*.txt`.

---

## E-089 — `arrowPaths.reduceAgainstPivots` is now a normal form (unit test that fails on the old code), and `MONO=1` at n = 8, L = 6 finds no monomial cord member for the six walked LNAs that carry the non-monomial ones
*2026-10-01* · **Library patch: the loop of `reduceAgainstPivots` moves a non-pivot head to the residue and continues with the tail, so the residue has no pivot column (E-085's defect). `tests/test_procedure.py::test_reduce_against_pivots_is_a_normal_form` fails on the old code (referee stashed the patch); touched test files 42 passed, 1 xfailed, six importing test files 552 passed (referee); the E-085 pair has equal residues (`workshop/rounds/015/theorist_step7.py`); n = 8, L = 5, LNA 4 non-`MONO` gives the same 300 cord members and histogram as round 015. `MONO=1 workshop/rounds/015/toolsmith_cords.py 8 6 1 1 I I --plan` for LNA index I = 4, 9, 10, 11, 12, 13: 0 members each, walks 200-315 s with six in parallel (the referee re-ran LNA 11, 61 s alone).** · *workshop round 017, toolsmith, refereed by skeptic*

**Limits.** Canonicity (congruent implies equal residue) assumes the pivots span the whole ideal in the block; the test covers an abstract case and a commutative square, and the real-ideal half passes on the old code too. The `MONO` negative is over six LNAs among those walked (LNAs 16..428 not walked at L = 5 or 6), a plan over visited quivers and not a search from a candidate, and no run shows `MONO=1` returning a nonzero count anywhere (no positive control); "why none exist" (a cord needs a sum relation) is a heuristic. The E-084 counts under the patch are E-090. E-087's wording "MONO=1 not run" at n = 8 is superseded for these six LNAs.

**Reproduction.** `timeout 10m .venv/bin/python -m pytest -q tests/test_procedure.py -m "not slow"`; `MONO=1 timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 8 6 1 1 11 11 --plan`; outputs `workshop/rounds/017/toolsmith_mono_n8_L6_lna*.txt`, `toolsmith_nomono_n8_L5_lna4.txt`.

---

## E-088 — Every interior placement of `3334` and `2455` is one double mutation from `35` (lemma R, `(a,b,b,d) -> (a-1,b,d+1)`), so they lie in the class of `333@1`, never of `333@0` (the `444` orbit); at n = 12..17 the `444` orbit holds `333@0` and the small `235/255/455` orbit holds `333@1`
*2026-10-01* · **`workshop/rounds/017/theorist_rrule.py N` against `doubleMutation.rewritesOf`: four-run 250 of 250 and three-run 119 of 119 at n = 14 (a in 3..b, b <= 7 or 8, d <= 9), 77 of 77 at n = 12; `a = 2` not tested (not an LNA row). `theorist_label.py N` (n = 12..16) labels each closed orbit by the offsets `J` of `333@o'` it holds: `J = {0, n-6}` for `34, 444, 2234, 2444, 4445`, `J = {1, n-7}` for `35, 55, 235, 255, 455, 2455, 3334` (sizes 447/320/763/516 at n = 13/14/15/16; n = 13: 2386 vs 447); referee reproduced n = 13, 15 and added n = 17 (`3334` size 1191, `J = {1,10}`, all 10 placements; `444` 11340, `J = {0,11}`) and the fold of class `3x` (x - 4) at x = 8, 9.** · *workshop round 017, theorist, refereed by skeptic*

**Limits.** Lemma R is verified by exhaustion over the stated bounded range against the rewrite engine, not proved (the submission's "proof" kind overstates it). "Class 0 differs from class 1" is a statement about this move set only, not an invariant or a derived-inequivalence proof. `k(33x) = 2x` is E-065's argument; the new part is the end link `3x@0 -> 33(x-1)@0`, seen on one labelled path (n = 13, x = 5). Letters >= 6, x >= 10 and the parity shadows are not explained. Explains E-083/E-086's `3334`, `2455` outside the `444` orbit.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/017/theorist_rrule.py 14`; `for n in 12 13 14 15 16; do timeout 10m .venv/bin/python workshop/rounds/017/theorist_label.py $n; done`; `timeout 10m .venv/bin/python workshop/rounds/017/theorist_class.py 13 34 35 36 37 55 235 455 2455 3334 444 2244`.

---

## E-087 — An n = 8 control with cords (8 and 9 arrows) and commutativity relations finds its source LNA at depth 6 in 2 of 2 and at depth 5 in 0 of 2, at 5.7e4-6.2e4 nodes, the size of the n = 9 negatives; no monomial cord member was found at n = 6, 7
*2026-10-01* · **`workshop/rounds/015/toolsmith_cords.py` (a raw visitor over the walk from an LNA, instead of `reachedQuipuAlgebras`, which keeps only quipu trees with monomial relations and so had no cords at all, E-082). Members from LNA 000300 (index 10; member A: 9 arrows, 2 relations) and LNA 000030 (index 4; member B: 8 arrows, 1 relation) at path length 6; search both directions, canonical-key node count. Depth 6: found in both (A 56 886 nodes / 6 157 distinct, 246 s; B 61 772 / 5 225, 247 s). Depth 5: not found in either (12 112 / 2 199; 12 851 / 1 990). Cord members at n = 8, L = 5 exist only for LNAs 4, 9, 10, 11, 12, 13 (those with a "3" early in the sequence), 9-arrow members only from 10, 11, 13; at L = 6 LNA 10 has 736 members, LNA 4 has 770. n = 6 smoke test (`BOTH=1 ... 6 4 1 2 4 4`): found at depth 4 in 2 of 2; at depth 3 one not found, one found through a relabelled copy that is nearer. `MONO=1` at n = 6 (13 LNAs, L = 5) finds 0 cord members; in the raw visitor at n = 7 (LNAs 0-6, depth 5) no node with arrows >= n has a monomial relation.** *Workshop round 015, toolsmith, refereed by maverick*

**Limits.** Every cord member the walk produces has a relation that is a sum of two paths (commutativity squares), while the n = 9 candidates are monomial quipus with 2-3 cords, so this controls the search on sum-relation cords, not on the class of the E-076 candidates; "no monomial cord member" is checked at n = 6, 7 only (not n = 8, `MONO=1` not run). Two members of 736 and 770, two LNAs of 429; one 9-arrow member only (a second from LNA 10 hit the 10-minute cap). Members are keyed by labelled algebra, so "found at L and not at L - 1" is not the shortest-path statement: the n = 6 smoke test has a relabelled copy found one step nearer. The referee reproduced the n = 6 test and the `MONO=1` plan, and checked member A's n = 8 output file against the table; the member B and n = 8 depth-5 searches were not re-run. Replaces the "no cords" limit of E-082: that was a property of the enumerator, not of the walk. Does not show that the E-076 negatives are not explained by blindness to monomial cords.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 8 5 1 1 0 15 --plan` (about 5 min); `BOTH=1 timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 8 6 1 2 10 10` (member A, about 8 min); `BOTH=1 timeout 10m ... 8 6 1 1 4 4` (member B); `BOTH=1 timeout 8m ... 6 4 1 2 4 4` (5 s); `MONO=1 timeout 4m ... 6 5 1 1 0 12 --plan`; `timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_survey.py 8 6 0 0` (3 min). Outputs `toolsmith_cords_*.txt`, `toolsmith_survey_lna0.txt`.

---

## E-086 — At n = 12..15 every 3- and 4-letter word the earlier scans called merged has all its placements inside the row set of the `444` orbit or none (0 partial), and E-075's "20 of 25" at n = 14 is not reproducible (11 of 25)
*2026-10-01* · **`workshop/rounds/015/skeptic_rowset.py N` for N = 12..15 (about 1 min together): S = row set of the closed `444` orbit (1 410 / 2 386 / 3 767 / 5 648 rows); for each word listed merged by `rounds/013/skeptic_orbscan_n*.txt` and `rounds/014/experimentalist_orbscan4_n*.txt`, the start row of each placement is in S (then its orbit is S) or not (then the orbit is walked and intersected with S). 3-letter words in S / outside: 7/6, 9/10, 11/14, 13/18; 4-letter: 5/0, 10/2, 13/2, 16/4; partial: 0 at every n. At n = 14 the scan has 25 merged words, 11 in S (`234 244 346 444 445 446 447 448 456 467 478`), 886:3, 491:3, 11820:3 and five singleton orbits. The 4-letter words outside S (`2455 3334` at n = 13, 15; `2457 2466` at 14; `2457 2477` at 15) have orbits of the same size as small 3-letter merged orbits (447/763, 1636/2290, 886, 881): same size only, not compared as sets.** *Workshop round 015, skeptic, refereed by experimentalist*

**Limits.** "Every merged word" means every word the round-013/014 scans called merged; the lists were not re-derived. "Start row in S implies orbit = S" rests on the `444` orbit being closed under the limit 300 000 (asserted, held). Closes E-083's "row sets not compared" for n = 12, 14, 15 (and confirms n = 13); n = 17 is still sizes only. One orbit per n, so no rate of a letter effect (as E-079). E-075's "20 of 25" cannot be recovered from saved data: the round-010 output records the orbit size (column 5), not an id, and gives 11 for 3767 too; the only natural 20 is 25 minus the 5 single-word orbits, which is a guess. The referee reproduced all four output files byte-identically and recounted n = 14.

**Reproduction.** `for n in 12 13 14 15; do timeout 10m .venv/bin/python workshop/rounds/015/skeptic_rowset.py $n; done` (about 1 min wall in parallel); `awk '$3=="merged"{c[$6]++} END{for(k in c)print k,c[k]}' workshop/rounds/013/skeptic_orbscan_n14.txt`; outputs `skeptic_rowset_n{12..15}.txt`.

---

## E-085 — The n = 8 class 2 loose end of E-084 is a defect of the mutation rewrite, not of `tiltingPlus`: `arrowPaths.reduceAgainstPivots` is not a normal form, and the Cartan congruence fails on all 11 replayed rejecting parents (n = 7, 8, 9), agreeing with `tiltingPlus`
*2026-10-01* · **`workshop/rounds/015/theorist_cartan.py N C` (Cartan matrix R C R^T predicted for the mutated vertex against the child's, from the replayed parents of E-084): the predicted matrix has an entry -1 and the child's key differs in all 11 rejecting parents (n = 7 class 0: 5, n = 8 class 2: 1, n = 9 class 0: 5). The 10 n = 8 class 2 steps with gate True, `tiltingPlus` True and key moved (7 distinct parents; 6 of 10 lines have one pair of parallel arrows, 4 none): `tiltingPlus` is right (R C R^T is the Cartan matrix of End(T)); the child is one dimension too big at one or two Cartan entries, because step 7 of `procedure.mutateAtVertex` solves `_kernelOverIdeal` on the residues of `reduceAgainstPivots`, which reduces only while the head is a pivot. A congruent pair in n = 8 class 2, parent 1: a = 1>2>7 and b = -1>5>6>7 have `isInIdeal(a - b)` True (a - b = 1>5>7') but residues {-1567, +157'} and {-1567}; the kernel at target 7 is empty. With a full reduction (monkeypatch, `theorist_fix.py`, library untouched) all 10 steps give a congruent Cartan matrix and the parent's key, and the rejections of n = 7, 8, 9 class 0 stay rejected.** *Workshop round 015, theorist, refereed by skeptic*

**Update (round 017).** The fix is in the library (E-089) and the E-084 walks re-run under it give the same counts where comparable (E-090); the n = 8 class 2 walk was too short to reach the 10 steps.

**Limits.** Shown on the 10 steps of one saved file and the 3 rejection sets; how often the defect fires elsewhere is unknown (the other 13 sampled classes of E-084 recorded no such step). E-084's counts (rejections, about 1.3e6 guard-admitted steps) were produced with the unpatched rewrite, so children that depended on a missing relation may have been wrong; not re-run. "Does not change anything else" is checked only on the fast `arrowPaths`/`procedure` tests (referee: 23 passed, 1 xfailed, identical patched and unpatched; the author wrote 24 passed) and these parents; the slow mutation tests and BFS were not run under the patch. Whether the rank step of `tiltingPlus` also uses non-canonical residues is untested. That R C R^T is the Cartan matrix of End(T) in the vacuous case rests on Ladkani 2.3(c), not independently verified. No unit test added and no library change made: the fix is a request to the toolsmith. The referee reproduced the congruent pair, the n = 8 class 2 output and the three rejection sets; `theorist_cartan.py`, `_detail.py` and `_child.py` were not re-run.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/015/theorist_cartan.py 8 2` (also `7 0`, `9 0`); `.venv/bin/python workshop/rounds/015/theorist_step7.py` (5 s); `.venv/bin/python workshop/rounds/015/theorist_fix.py 8 2` (also `7 0`, `8 0`; 4 s); `theorist_detail.py`, `theorist_child.py`; outputs saved beside them.

---

## E-084 — Guarded walks from LNAs reach gate-admitted parents where `tiltingPlus` is False, at n = 6 (distance 8) and in each of 10 sampled classes at n = 7, 8, 9 (distance 5-7); the Coxeter guard refuses every such step, and none of about 1.3e6 guard-admitted steps fails `tiltingPlus`
*2026-10-01* · **`workshop/rounds/014/scholar_walk.py N --class C --stop-on-reject` (BFS over distinct algebras by `canonicalKey`, from the LNAs and relation duals of length n, one Coxeter-key class at a time, steps taken only when the gate admits and the key is unchanged) at n = 6, class 0 (5 616 algebras, 23 s): 8 rejecting parents at distance 8 (gate True, `tiltingPlus` False, child key != base), 9 993 guard-admitted steps all `tiltingPlus` True, 0 admitted-but-failing. n = 7 classes 0-2 (first rejection at distance 5, 6, 6), n = 8 classes 0-2 (6, 6, 7), n = 9 classes 0-3 (5, 6, 6, 6): rejections in all 10, every one also refused by the guard. n = 5 closes (11 700 algebras, 2 classes): no rejection. n = 6 classes 1-3: none to depth 11, 12, 15 (open, 119k-125k algebras). `scholar_replay.py 7 0` rebuilds the n = 7 path (4, 7, 5, 7, 5): every prefix step admitted with the key unchanged, then gate True, `tiltingPlus` False. The reached rejecting parents carry the shape of E-078 (commutative square into a vertex with one outgoing arrow); E-078's own algebra is not reached. New test `tests/test_gate_without_tilting.py` (2 tests).** *Workshop round 014, scholar, refereed by skeptic*

**Limits.** Referee reproduced the n = 6 rejection, the n = 7 replay and the test; the 540 s runs and the n = 7..9 sweep were read from saved outputs. Supersedes E-057's "there is none" and E-055's depth-limited "gate = guard": those stopped at depths short of 8 (n = 6) and 5 (n = 7). Stopped at the first rejecting level, so no count of rejections per n; n = 9 is 4 of 19 classes. "All have the A5 shape" was read from the rejecting parents inspected, not checked by script for all classes, and non-tilting rests on `tiltingPlus` alone (no independent Cartan congruence for the replayed parents). "Guard refuses every rejection" is nearly forced, since a non-tilting step usually moves the Coxeter polynomial; the evidence that matters is the 0 of about 1.3e6 guard-admitted failures (that total sums runs of different length; the table row at n = 6 uses the first-level run). *[Round 015, E-085: this loose end is a defect of `reduceAgainstPivots` in the rewrite, not a gap in `tiltingPlus`; the counts of this entry were made with the unpatched rewrite and are not re-run.]* Loose end: n = 8 class 2 has 10 steps with gate True, `tiltingPlus` True and the key moved (parallel arrows or duplicated relations), so `tiltingPlus` is not shown complete there and "the guard is redundant" is not established. The gate is unsound at n <= 9, not the guarded walk. `isTilting` stays unpromoted (STEERING round 011, q1: a gate-admitted rejection from a *guarded walk* is what was asked; this is gate-admitted, guard-refused).

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/014/scholar_walk.py 6 --class 0 --stop-on-reject` (23 s); `.venv/bin/python workshop/rounds/014/scholar_replay.py 7 0` (1.4 s); `.venv/bin/python -m pytest -q tests/test_gate_without_tilting.py`; sweep `workshop/rounds/014/scholar_sweep.sh`; outputs `scholar_walk_n*_c*.txt`.

---

## E-083 — Saved n = 17 run: `5046` and `5056` each have two closed orbits (122 673 rows at even offsets, 54 266 at odd); the 4-letter words with a 4 at n = 12..15 add members to the one big merged orbit, whose size equals that of the `444` orbit
*2026-10-01* · **`batch.py orbits 17 --cores 5046` and `--cores 5056` (5 m 36 s and 5 m 32 s, run in parallel, shared ledger): both words have orbit {0,2,4,6} of 122 673 rows and {1,3,5,7} of 54 266, no pairs, each closed under mirror, limit 1 500 000; saved as `experimentalist_n17_5046.txt`, `_5056.txt`. 4-letter scan (letters 1..9, nondecreasing, containing a 4, >= 4 placements, middle-offset orbit, reduced walk, limit 300 000): words 20 / 35 / 56 / 84 at n = 12..15, merged 5 / 12 / 15 / 20, one orbit holding 5 / 10 / 13 / 16 of them, of size 1 410 / 2 386 / 3 767 / 5 648 (the `444` orbit sizes of E-079); the rest are small orbits (`3334`, `2455` at n = 13, 15; `2457`, `2466`, `2477`).** *Workshop round 014, experimentalist, refereed by theorist*

**Limits.** Confirms E-080 at n = 17 (and lifts its "`5056` not re-run, no saved output"); a slice of E-075/E-079 for 4-letter words. The two n = 17 runs share one ledger, so the agreement of the words is not two independent computations; row sets were not compared, only sizes. "Same orbit as `444`" is by row membership at n = 13 only (referee: `234 2234 2244 2444 2445 4445 2346` inside, `3334` and `2455` outside); at n = 12, 14, 15 it rests on equal size. "Rigid" covers both words split over several orbits and words whose middle orbit lacks some placements; not separated. Referee re-ran the scan at n = 12, 13 byte-identical and recounted n = 12..15; did not recompute n = 17. As in E-079, the count is one event per n, not a rate of a letter effect. Not run: n = 16, words without a 4, `--max-word 5`.

**Reproduction.** `timeout 10m .venv/bin/python batch.py orbits 17 --cores 5046,5056 --plan`; `timeout 10m .venv/bin/python batch.py orbits 17 --cores 5046` (5.5 min; also `5056`); `for n in 12 13 14 15; do timeout 10m .venv/bin/python workshop/rounds/014/experimentalist_orbscan4.py $n; done` (9 s, 16 s, 58 s, 3 m 10 s); `.venv/bin/python workshop/rounds/014/experimentalist_orbstats4.py`.

---

## E-082 — An n = 8 control with relation-bearing sources finds its source LNA in every run at its own depth and in none one step short, at 7e3-4e4 nodes at depth 6; no member has cords
*2026-10-01* · **`workshop/rounds/014/maverick_control8.py` (extends the round-013 control with a minimum-relations argument; members rebuilt from the certificate): depth 5, 12 members with 1 relation (two per LNA 0-5, the deterministic head of the sort by fewest relations then fewest arrows): source found 12 of 12 (4 052-8 880 nodes), 0 of 12 at depth 4. Depth 6, LNA 0, 1, 2 (1 relation): found, 19 417 / 26 237 / 38 516 nodes (46-97 s); LNA 0 with 4 relations (one member, `HIGH=1`): found, 6 857 nodes; the 1-relation member of LNA 0 at depth 5: not found (4 573 nodes). Member counts at L = 6 (relations >= 1): LNA 0: 835 (1 rel 214, 2 rels 383, 3 rels 224, 4 rels 14), LNA 2: 423.** *Workshop round 014, maverick, refereed by toolsmith*

**Limits.** The members are the deterministic head of a sort, not a sample: every one has 7 arrows on 8 vertices, so none has cords, which the n = 9 candidates (2-3 cords) do. So it controls inverse-move handling for relation-bearing starts (E-069), not whether depth 6 reaches a class with cords. Node counts use the same semantics as `toolsmith_verify.py` (referee checked). Size comparison: n = 9 negatives (50 476-62 888 nodes, E-081) are 1.3-3.2 times the 1-relation n = 8 depth-6 controls and 7-9 times the one 4-relation control; the relation-richest control is the smallest (one member, so "relations restrict mutation" is a guess). Extrapolations (depth 7 about 1e5 nodes at n = 8; about 30 h for all 429 LNAs) are not results. The referee re-ran LNA 0 at L = 6 identically (3 m 19 s). Only LNAs 0-5 (L = 5) and 0-2 (L = 6) run. Narrows E-081's "non-hereditary source" limit; the cords limit stays open.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/014/maverick_control8.py 8 6 6 1 1 0 0` (3.5 min); `SHORT=1` same for depth 5; `HIGH=1` for the 4-relation member; L = 5: `... 8 5 5 1 2 0 5`; outputs `maverick_control8_L*.txt`.

---

## E-081 — The n = 9 depth-6 negatives of E-076 are full walks of 50 476 to 62 888 nodes (4 of the 16 candidates measured); an n = 7 depth-6 control finds its source LNA 16 of 16 (0 of 4 at depth 5)
*2026-10-01* · **`workshop/rounds/013/toolsmith_verify.py` (the round-010 script with a node count and a count of distinct E-042 keys on each `reached` line) gives, at n = 9, depth 6, K = 4, candidates 0, 5, 9, 13: `reached []`, nodes 50 476 / 62 888 / 62 165 / 55 247, distinct 4 437 / 7 074 / 6 395 / 5 683 (263-398 s with other jobs running). Candidate 0 at depth 3, 4: 345 / 166, 1 853 / 518 (reproduced by the referee). Per-step node growth about 5.2 (1 853 -> 50 476 over two steps is about 27x), so depth 7 is about 2.6e5 to 3.4e5 nodes, 25-30 min per candidate, consistent with E-076. Control at n = 7, depth 6: members of classes at recorded path length exactly 6 from an LNA (0 to 2 relations, mostly hereditary) return their source LNA in 16 of 16 searches (3 540 to 16 894 nodes), and in 0 of 4 at depth 5 (1 049 to 4 559 nodes).** *Workshop round 013, toolsmith, refereed by theorist*

**Limits.** A node count shows the search ran, not that it was sufficient: the n = 9 candidates (1-2 relations, 2-3 cords) are not matched in structure to the cheap control members, branching grows with n, and no n = 9 class member is known within 6 steps, so the control tests inverse-move handling (E-069) and does not show a depth-6 ball can meet an n = 9 class with no hereditary member. The comparison "3 to 4 times the nodes of the larger n = 7 controls" does not carry over to coverage. 4 of 16 candidates measured. The first control run (LNAs 0-3, 4 of the 16) has no saved output; those four rows are covered by the depth-5 negative on the same members and the depth-4 regeneration only. `distinct` counts keys and is a lower bound on classes if keys collide (not audited). LNAs 10+ at n = 7 not run.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/013/toolsmith_verify.py 9 6 4 -1 --cand 0` (263 s; also `--cand 5 9 13`); `.venv/bin/python workshop/rounds/013/toolsmith_verify.py 9 4 1 -1 --cand 0` (8 s); `timeout 10m .venv/bin/python workshop/rounds/013/toolsmith_control.py 7 6 6 2 4 9` (about 9.5 min, 12 of 12); `SHORT=1 timeout 10m .venv/bin/python workshop/rounds/013/toolsmith_control.py 7 6 6 1 0 3` (depth 5, 0 of 4).

---

## E-080 — `3a`@o is joined to `3a`@(o+2) by a shortest path of a - 1 moves for a = 5..12, and `4@0 -> 3@1` is a single anchored table rule; `5046`, `5056` have two closed orbits (even and odd offsets) of different sizes at n = 17; no invariant separating P from Q is found
*2026-10-01* · **In the reduced guarded walk (BFS over the table rules, edge and double moves, with a 2-arrow spectator added then stripped) `3a`@0 -> `3a`@2 takes a - 1 moves for a = 5..9 (`35` at n = 12, 14, 16; `36` n = 13; `37` n = 15; `38` n = 16; `39` n = 17; `35`@1 -> `35`@3 at n = 12), and a referee found 9, 10, 11 moves for a = 10, 11, 12 (n = 18, 19, 20). The path: add the spectator, one move `3a@0 + 2 -> a@1 + 3@(a-1)` (table rule for a = 5, 6; double `[2,2]` for a >= 7), then doubles `[t,t]` for t = a-1, ..., 3; for `35` at n = 12 this is `3500000000, 0500300000, 0403330000, 0333400000, 0035000000`. So P (orbit of `3@2`/`3@3`) and Q (of `5@0`/`6@0`) of E-077 are each closed under offset shifts by 2. `4@0 -> 3@1` at n = 12 is one move, the anchored rule `(4,((0,4),),((1,3),),(1,-3,-4))`; the only table rules with a lone relation of more than two arrows on the left are four of width 4, so at n = 12 `k@0` (k = 5, 6, 7) has 2 neighbours and `4@0` has 5. `5@0` has a closed orbit of 148 rows (n = 12) without `3@2`. For `5046`, `5056` at n = 13, 15, 17: two closed reduced orbits (even / odd offsets, same row sets for both words, one Coxeter key), sizes 4217 / 2116, 18157 / 8794, 122 673 / 54 266; the sizes differ, so the mirror cannot swap them and the two words are key-coarser at n = 17 (for these words only). `5046`@0 -> `5046`@2 at n = 13 takes 7 doubles.** *Workshop round 013, theorist, refereed by skeptic*

**Limits.** The paths are shortest paths in the author's own move set; no independent check that each double `[t,t]` is a valid equivalence beyond the move table. "a - 1" is read off a = 5..12, not proved. No proof, and no invariant, that a shift by 1 is impossible: no relation count, letter range, letter sum, GF(2) functional or Smith form separates P from Q, so "parity class" names two orbits. The neighbour counts of `k@0` are for n = 12 only. The n = 17 orbit sizes of `5046`/`5056` come from one run each (about 3.5 min); no output file was saved. The referee later re-ran `5046` at n = 17 (3 m 55 s) and matched it (122 673 rows at offsets 0,2,4,6; 54 266 at 1,3,5,7; one key class); `5056` at n = 17 was not re-run. The n = 13 and 15 sizes and the 7-move path were reproduced. The other 137 cores at n = 17 were not run.

**Reproduction.** `.venv/bin/python workshop/rounds/013/theorist_path.py 12 35 0 35 2` (also `13 36 0 36 2`, `15 37 0 37 2`, `17 39 0 39 2`, `13 5046 0 5046 2`, `12 4 0 3 1`, `12 5 0 3 2`; seconds each); `.venv/bin/python workshop/rounds/013/theorist_k0moves.py 12 3 7`; `.venv/bin/python workshop/rounds/013/theorist_orbitstats.py 12`; `timeout 10m .venv/bin/python workshop/rounds/011/theorist_word.py 17 5046` (about 3.5 min); `.venv/bin/python workshop/rounds/013/theorist_size.py 15 5046 0 100000`.

---

## E-079 — Counted by orbit, the 3-letter merged words at n = 12..15 lie in one big orbit per n (the `444` orbit, 7, 9, 11, 13 merged words), which holds both `34`-words (`234`, `346`) and 4-no-`34` words, so the data cannot separate "letter 4" from "collapse to `34`"
*2026-10-01* · **Rerunning E-075's scan (nondecreasing 3-letter words, letters 1..9, at least 4 interior offsets, n = 12..15; 291 word-cells after dropping the 4 size-1 cells `222`; 84 merged) with each word's orbit recorded: at each n exactly one big merged orbit, the one of `444`, holding 7/8, 9/12, 11/16, 13/19 of the merged words containing a 4 (the `34`-words `234`, `346` 2/2 at each n; the rest `244 444 445 446 447 448 449 456 467 478 489`). The other merged orbits are small (1 to 3 words) and mostly hold a 2 (`{236,266,466}`, `{235,255,455}`, `{237,277,477}`, `{238,288,488}`, `{239,299,499}`) or carry the three no-4-no-2 merges of E-075 with 4-words as orbit-mates (`{458,468,568}`, `{459,479,679}`); `233 255 266 277 288 457` are single-word merged orbits. `344 345 347 348 349` are rigid at n = 14, 15. Pooled over the four n (an orbit counted once per n), merged words / total: contains `34` 8/26, 4 without `34` 47/74, no 4 with a 2 26/70, no 4 no 2 3/121.** *Workshop round 013, skeptic, refereed by experimentalist*

**Limits.** This sharpens E-075; it does not refute its observation that `444` is the only merged `aaa` and that merged words with no 2 are almost all in the `444` orbit. The 55/100 against 3/121 contrast is largely one orbit's size, not a rate of a letter property. The orbit of a word is the orbit at its middle offset only, so "orbits hosting a word" is a convention for rigid words; the big orbit is counted once per n, so pooled orbit rows are not independent. The count here is 11 of the 25 merged words at n = 14 in orbit 3767; E-075's Limits say 20 of 25 and the two are not reconciled (the referee re-ran the scan with identical output; E-075 may have counted something else). `34` is neither sufficient nor necessary for merging. The `{4aa, 2aa}` pattern is listed, not tested. 3-letter words only, n = 12..15, no zeros, no n >= 16.

**Reproduction.** `for n in 12 13 14 15; do timeout 10m .venv/bin/python workshop/rounds/013/skeptic_orbscan.py $n; done` (about 5 min); `.venv/bin/python workshop/rounds/013/skeptic_orbstats.py` (seconds); outputs `skeptic_orbscan_n12..15.txt`, `skeptic_orbstats_out.txt`.

---

## E-078 — A hand-built 5-vertex algebra with a length-4 commutativity relation `abde = acde` is gate-admitted at `d` while `tiltingPlus` is False and the Cartan congruence fails; the same holds for all 6 padded versions at n = 5..7, and the length-3 square and the monomial control behave otherwise
*2026-10-01* · **For n = 5..7 (chain prefix + square a>b, a>c, b>d, c>d + d>e + chain suffix; 18 algebras), at every vertex with an outgoing arrow: relation `abde = acde` (6 algebras) has exactly one gate-admitted vertex with `tiltingPlus` False, always `d`, and the repo's rewrite fails Cartan(child) = R C R^T there; relation `abd = acd` (6) is (True, True) everywhere; relation `abde = 0, acde = 0` (6) has the gate refuse `d`. `tiltingPlus` and the congruence agree on all admitted steps (0 disagreements; gate-refused with `tiltingPlus` True: 0). Mechanism: `e_a A e_d = <abd, acd>` has dimension 2 and right multiplication by `d>e` sends both to the one element `abde = acde`, so the one map of Aihara-Iyama 2.32(b) / Ladkani 2.3(c) has kernel `abd - acd`; with the relation at length 3 the space has dimension 1 and the map is injective.** *Workshop round 011, scholar, refereed by theorist*

**Limits.** This answers E-066's request for a smaller instance of the step-7 shape: the shape occurs at n = 5, so n = 10 is not special. The algebras are hand-built, not shown reachable from an LNA by guarded steps, so it does not touch H-015's LNA scope, and it is not an extra gate-admitted rejection on a walk (STEERING round 002, question 1, still unmet). "Same shape as E-032 step 7" is not checked as an isomorphism (E-066's relation reads `+ = 0`, the script uses equality). The Cartan failure is only as independent as the repo's rewrite: it shows the rewrite is not End(T) of the one-map approximation, not that the child is not derived equivalent; no Cartan matrices are printed in the note and End(T) was not computed independently. The one-map identity itself (AI = Ladkani = `tiltingPlus`) is not derived. CHZ arXiv:2509.12983 (Cor 3.6 "monomial" wording) is **unread**: arxiv.org is refused by the proxy.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/011/scholar_square.py` (2 s).

---

## E-077 — The key-coarser lists A (even n) and B (odd n) of E-074 are the words whose placements alternate between the two single-relation orbits of `3@2` and `5@0` (even n) or `3@3` and `6@0` (odd n); the census rule gives list A at n = 12, 14, 16, 18 and list B without `5046 5056` at n = 13, 15, 17 (a fit at 12..16 and a consistency check at 17, 18, not a prediction)
*2026-10-01* · **For n = 12..20 the reduced orbits P (of `3@2` / `3@3`) and Q (of `5@0` / `6@0`) are closed, disjoint, each holds its own mirror, with `35` (even n) / `36` (odd n) at even offsets in P and odd offsets in Q (sizes n = 12..20: P 178, 449, 320, 743, 516, 1123, 774, 1597, 1102; Q 148, 224, 272, 390, 446, 614, 678, 904, 976); for n = 14..20 these four rows are the only single-relation rows with that key (n = 12, 13 have a third pair `9@..`). The rule "every placement lies in single-relation orbits, alternating between exactly two orbits that share a key and each hold their own mirror" reproduces list A at n = 12, 14, 16, 18 and list B minus `5046 5056` at n = 13, 15, 17; `406` alternates between two single-relation orbits with one key but is not in B because the mirror swaps them. Observation: the orbit of `4@0` holds `3@1` at every n = 12..16, while `k@0` for k >= 5 joins a 3-orbit only when k and n have the same parity (n <= 16).** *Workshop round 011, theorist, refereed by skeptic*

**Limits.** The rule was tuned on lists A, B at n = 12..16 (the mirror clause was added to remove `406`), so 12..16 is a fit. At n = 17, 18 no key-coarser list was ever computed: the comparison is against the hard-coded n <= 16 lists, so those two n show only that the rule's output is stable, not that the words are key-coarser there. The rule misses 2 of the 10 words of list B (`5046 5056`, placements in orbits holding no single relation) at every odd n. No selectivity null. No move sequence `35@o -> 35@(o+2)` was exhibited; P != Q is a reachability fact of the guarded walk. The letter-4 observation is for n = 12..16 and does not explain why `444` merges (E-075). Null results: no GF(2) functional of the rows, none of nine integer statistics mod 2 or 4, and no Smith normal form of five Coxeter-matrix polynomials separates P from Q. A referee re-ran `predict.py 18`, `census2.py 18`, `triples.py 8 20` with matching output.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/011/theorist_predict.py 18` (6 s); `timeout 10m .venv/bin/python workshop/rounds/011/theorist_census2.py 18` (about 1 min); `.venv/bin/python workshop/rounds/011/theorist_triples.py 8 20`; `.venv/bin/python workshop/rounds/011/theorist_k0.py 12 16`.

---

## E-076 — All 16 K = 4 candidates at n = 9 reach no quipu-class member at depth 6 (the 13 not yet searched, 308-567 s each under a 600 s cap)
*2026-10-01* · **`workshop/rounds/010/toolsmith_verify.py 9 6 4 -1 --cand I` for I in 1, 2, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15 finishes with `reached []` and rc 0 each; with indices 0, 4, 8 (E-072, E-073) all 16 candidates are negative at depth 6. Wall seconds with four shards at once: 308, 359, 524, 453, 488, 314, 524, 518, 567, 474, 322, 365, 309 (cand 11 only 33 s under the cap). A referee re-ran index 15 alone: `reached []`, 308 s, the same as in the shard, and checked all 16 indices against `--list`.** *Workshop round 011, experimentalist, refereed by toolsmith*

**Limits.** A bounded negative: by E-069 the search finds a class only if a member lies within its depth, and no depth-6 control exists at n = 9 (E-072's L = 5 control at n = 6, 7 is the nearest). The completion test is rc 0 plus a `reached` line (a timeout leaves neither). The output has no node count, so an early exit cannot be told from a full search. Only the 16 candidates of K = 4 were searched, not the 160 of K = 100. Depth 7 is about 5.5 times depth 6 (two data pairs), so it does not fit a 10-minute shard (Menu 4); `--budget-hours` is checked only between candidates.

**Reproduction.** `workshop/rounds/011/experimentalist_shards.sh I` (300-570 s each; it runs the `toolsmith_verify.py` command above under `timeout 10m`); outputs `workshop/rounds/011/experimentalist_cand{I}.txt`, times `experimentalist_shard_times.txt`.

---

## E-075 — A scan of every nondecreasing 3-letter word at n = 12..15 finds `444` the only merged `aaa` (a = 3..9), but words containing a 4 are over-represented among merged words (55/100 against 3/121 with neither 4 nor 2), so "a = 4 special" does not single out the collapse to `34`
*2026-10-01* · **Over 291 word-cells (letters 1..9, at least 4 interior offsets), "merged" = the reduced-walk orbit of the word at one interior offset holds it at every offset. Among `aaa`, a = 3..9, only `444` merges at n = 12..15 (333, 555, 666, 777, 888, 999 rigid in all 22 cells; a referee's probe at n = 16 agrees for 333, 444, 555; `222` is a size-1 orbit, so a = 2 is untestable). Pooled merged rates: all 84/291; contains a 4 55/100; no 4 29/191; no 4 and no 2 3/121 (`568` at 14, 15, `679` at 15); no 4 with a 2 26/70. Words with a 4 and no `34` merge (`445 446 447 448 456 457 467 468 478`) while `344 345 347 348 349` do not (`346` does).** *Workshop round 010, skeptic, refereed by experimentalist*

**Limits.** *[Round 015, E-086: "20 of 25" is not reproducible; the scan gives 11 of 25.]* The cells are heavily dependent: at n = 14 the orbit 3767 holds 20 of the 25 merged words, so the 0.55 rate measures largely one big orbit and no p-value is given; the referee asked for the count of merged words in the `444` orbit against other orbits, which was not made. The all-offsets test is easier for short (large-letter) words. The result supports "words with a 4 are over-represented in the large orbit", not "merging is a property of the letter 4". `34x` lies outside E-071's stated criterion, so its rigidity does not by itself refute the `34` route; the data cannot tell the `444 -> 34` mechanism from any other. One offset per word; n = 12..15 only; words with zeros and 4-letter words not scanned.

**Reproduction.** From the repository root: `timeout 10m .venv/bin/python workshop/rounds/010/skeptic_scan.py 14 4` (about 1 min; n = 12, 13, 15 likewise); `.venv/bin/python workshop/rounds/010/skeptic_stats.py` (output `skeptic_stats_out.txt`); `timeout 10m .venv/bin/python workshop/rounds/010/skeptic_probe.py 14 44 333 444 555`.

---

## E-074 — The key-coarser cores of the `--max-word 4` catalogue are the same 9 words at n = 12, 14, 16 and the same 10 at n = 13, 15, and the size-20300 orbits of `348`, `349` at n = 16 are the `4056` orbit and its mirror (n <= 16, word length <= 4)
*2026-10-01* · **At n = 14, 15, 16 all 139 cores close (limit 1500000); orbit-plus-mirror equals the key in 130, 129, 130, the key is coarser in 9, 10, 9 and finer or incomparable in none. List A (n = 12, 14, 16): `35 455 3334 3336 5003 5055 5504 5505 5506`; list B (n = 13, 15): `36 405 466 3335 5004 5006 5046 5056 5066 5605`; in each the orbits are the two parity classes and the key merges them. At n = 16 the walks of `4056`@1, `348`@2, `349`@3 have identical row sets (20300 rows), likewise `4056`@2, `348`@3, `349`@1; the two sets are disjoint and each holds the mirror of the other's start row.** *Workshop round 010, experimentalist, refereed by skeptic*

**Limits.** The equality of the lists is of word lists over three even and two odd n (n = 10 not rerun here; E-064 has its counts); nothing for n >= 17 or `--max-word 5`. `46 3355 3445` at n = 16 not rechecked against the 20300 orbit (E-064 identified them with `4056`). Why these words are the parity classes is unexplained. Each n ran in resumed windows of at most 9 minutes (n = 16 about 60 minutes in all).

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py N --jobs 4` for N = 14, 15, 16, repeating until it prints `n = N: 139 cores` (the ledger resumes; never two windows of one n at once); `timeout 10m .venv/bin/python workshop/rounds/010/experimentalist_same20300.py` (2 to 6 min); outputs `experimentalist_keycoarser_out.txt`, `experimentalist_same20300_out.txt`.

---

## E-073 — `toolsmith_verify.py` shards the n = 9 H-017 search by candidate and budget; the E-072 depth-5 negative reproduces and a second K = 1 candidate at depth 6 reaches nothing in 434 s (one candidate fits a 10-minute shard)
*2026-10-01* · **`workshop/rounds/010/toolsmith_verify.py` is `rounds/004/maverick_verify.py` plus `--list`, `--cand I[,J..]` and `--budget-hours H` (exit 2 before the next candidate once H hours are spent; checked between candidates only, so it does not interrupt a running one); search code unchanged (a referee diffed it). Candidate numbering depends on K: K = 1 gives 4 candidates, K = 4 gives 16, K = 100 gives 160, and the K = 1 candidates 0..3 are K = 4 indices 0, 4, 8, 12. At n = 9, K = 1: depth 5 reaches [] for all four (261 s; 42/66/79/68 s); K = 1 candidate 2 at depth 6 reaches [] in 434 s (referee 440 s), 5.5 times its depth-5 time. With E-072's candidate 1 (280 s) that is two depth 5 to 6 pairs, 5.4 to 5.7.** *Workshop round 010, toolsmith, refereed by theorist*

**Limits.** A reproduction of E-072 with tooling; the new data are one depth-6 negative with saved output. K = 4 index 8 is already done at depth 6 (so is index 4 by E-072's candidate 1): a 16-shard plan should skip them. 12 of the 16 candidates remain untimed at depth 6, so a candidate over the 10-minute cap is possible. The `--budget-hours 0.002` test (candidate 0 runs, rc 2) depends on machine speed; `--budget-hours 0` expecting all skipped does not. Under `--cand` the budget message lists only the selected candidates not searched.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/010/toolsmith_verify.py 9 5 4 -1 --list`; `... 9 5 1 -1 --budget-hours 1` (261 s); `... 9 6 1 -1 --cand 2` (440 s); outputs `toolsmith_verify_n9_d5.txt`, `toolsmith_verify_n9_d6_cand2.txt`.

---

## E-072 — The H-017 search passes an L = 5 control (42 of 42 at n = 6, 8 of 8 near-trivial LNAs at n = 7, 0 of 50 at depth L - 1); at n = 9 the four K = 1 below-diagonal candidates reach nothing at depth 5, and depth 6 costs about 5.7 times depth 5
*2026-09-30* · **`maverick_control5.py` rebuilds the longest-relation-set member of each LNA reached by a depth-5 walk (path length L = 5) as a path algebra and searches from it with the Coxeter guard on: n = 6, every LNA (42): source returned 42/42 at depth 5, 0/42 at depth 4; n = 7, the first 8 LNAs by sort order: 8/8 and 0/8. That fills E-069's missing L >= 5 case for n = 6. At n = 9 (`maverick_verify.py 9 D 1 -1`, the 4 measured K = 1 candidates): nothing reached at depth 4 (9-17 s each), nothing at depth 5 (49-98 s each, 321 s in all), nothing at depth 6 for candidate 1 (280 s). Depth 6 for the four is about 30 min, for all 16 candidates about 2 h (one candidate per 10-minute shard); depth 7 about 27 min per candidate** → H-017, E-069, E-063 · *workshop round 009, maverick, refereed by toolsmith*

**Limits.** The depth-5 negative excludes only members within 5 steps; the forward walks of E-063 needed depth 6 for `3033030` and `4444400`. The "recorded paths are shortest" check (50/50) runs the same exhaustive DFS and is a consistency check only. The n = 7 sample is 8 of 132 LNAs, all near-trivial; the n = 9 depth 5 and 6 outputs were not saved; the other 12 of the 16 candidates were not sized; the growth factor 5.7 rests on one depth 5-to-6 pair. A class with no hereditary member is still uncontrolled.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/009/maverick_control5.py 6 5 5 1` (42/42, about 3 min); `... 7 5 5 1 8` (8/8, 111 s); `timeout 10m .venv/bin/python workshop/rounds/004/maverick_verify.py 9 5 1 -1` (321 s); outputs `maverick_control5_n6.txt`, `maverick_control5_n7.txt`.

---

## E-071 — The drift families `aax` close into chain pairs with `k = 2x + 3 - a` for a = 3, 5, 6 (and 7 at n = 15) but `44x` joins the big orbit, matching a computed criterion: the seed `aaa` collapses to `34` only for a = 4
*2026-09-30* · **For `aax` the double mutation drifts `aax@o -> aa(x-1)@(o+1)` and the chain ends at the seed `aaa`; the family closes into pairs `{c, n-c}` unless the seed's span-2 collapse reaches the sliding word `44`. The collapse of `aaa` is `(a-1)a`: `23`, `45`, `56` do not reach a slider and the families are rigid, `34` does (7-step path `34@5 <- 403@6 ... 44@5` at n = 14) and `444@o` is in the orbit of `333@0` at every offset. n = 14 (`55x` x = 5..8, `66x` x = 6..8, predictions printed before the run): all rigid, `s = 6, 4, 2, 0` for `55x`, `5, 3, 1` for `66x`, all orbits closed; `33x` rechecked; `44x` one orbit of size 3767. Referee's extra runs: `55x` n = 12, 13; `66x` n = 13; `44x` n = 12 (merged, 1410); `77x` n = 15 (rigid, pairs `[0,5],[1,4],[2,3]`). The broader criterion "the orbit holds a span <= 3 word at every offset" fails (10 of 24 cells agree)** → H-021, E-065 · *workshop round 009, theorist, refereed by experimentalist*

**Limits.** A computed criterion with one positive datum (a = 4): not proved. Why `34` alone reaches the slider (`3333 <-> 3403` has no analogue for `23`, `45`, `56`) is unexplained. The drift for x > a is evidenced by the pair structure only, not displayed. a = 2, 8, 9 and n >= 16 not run; `34x`, `45x` are outside the criterion (no drift for `45x`).

**Reproduction.** `.venv/bin/python workshop/rounds/009/theorist_nbrs.py 44 4-6`; `.venv/bin/python workshop/rounds/009/theorist_path34.py`; `timeout 10m .venv/bin/python workshop/rounds/009/theorist_closure.py 14 55 5-8` (about 4 min; outputs `theorist_closure_*_n14.txt`); `timeout 10m .venv/bin/python workshop/rounds/009/theorist_translator.py 14 33 4-8` (the failed broader criterion).

---

## E-070 — The equal-size singleton pairs of `344`, `348`, `349` at n = 15..17 and `4046` at n = 14..16 are each one orbit and its mirror, and orbit-plus-mirror equals the key in all 12 cells; the key-coarser cores of E-064 at n = 12, 13 are disjoint from the 7 cores of E-059
*2026-09-30* · **All 12 orbit walks closed (limit 400000). Every unmerged pair of equal size is joined by the mirror (e.g. `348@16` `{2,3}` of size 20300, `349@17` `{1,4}`, `4046@16` `{1,4}{2,3}`), so the E-068 size-paired singletons are one orbit plus its mirror. The key-coarser cores are 9 at n = 12 (`35 455 3334 3336 5003 5055 5504 5505 5506`) and 10 at n = 13 (`36 405 466 3335 5004 5006 5046 5056 5066 5605`); neither meets the 7 of E-059 (`344 366 4044 4403 4404 4405 4605`). The key-coarser cores are parity-class orbits that each hold their own mirror; the 7 pair at n = 12 and mirror-join at n = 13, but the joins also occur at even n (`348@16`, `349@16`, `4046@14`, `4046@16`), so "pair at even n, mirror-join at odd n" is not supported** → H-021, E-064, E-068 · *workshop round 009, experimentalist, refereed by skeptic*

**Limits.** Only `n = 17` lies outside E-064's range; the n = 14..16 joins duplicate E-064 rows. The key-coarser lists were printed for n = 12, 13 only. Whether size 20300 of `348`, `349` at n = 16 is the orbit of `4056` (E-064) was not tested. Nothing for n >= 18 or x >= 10.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/009/experimentalist_mirrorjoin.py 17 344 348 349 --limit 400000` (about 7 min; n = 16 about 6 min; `4046` at 16 about 3 min; output `experimentalist_mirrorjoin_out.txt`); `.venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 12 --jobs 4`, same at 13.

---

## E-069 — The H-017 mutation search returns the source LNA in 273 of 273 round trips at n = 7 (91 of 132 LNAs) and 84 of 84 at n = 6 from depth L, and in 0 of 84 (n = 6) and 0 of 132 (n = 7) at depth L - 1; so a depth-4 negative means "no member within 4 steps" and no more
*2026-09-30* · **Positive control for the round-004 search (`search.linesReachedFrom` via `families.verify`, Coxeter guard on): quipu-with-relations members found by a depth-4 walk (the 3 longest-path members per LNA, path length L = 4 at n = 7, 3 at n = 6) are rebuilt as path algebras and searched at depth L: the source LNA is returned every time. At depth L - 1 it is returned never. Consequence: round 004's depth-4 "reached nothing" for 16 below-diagonal candidates at n = 9 excludes only class members within 4 steps; the forward walks needed depth 6 for `3033030` and `4444400`** → H-017, E-063 · *workshop round 007, maverick, refereed by experimentalist*

**Limits.** The round trip tests inverse-move handling, not independent discovery. The negative control relies on L being a shortest path (not checked). n = 7 covers 91 of 132 LNAs (10-minute cap; the tail is untested); no L >= 5 case; at n <= 8 every LNA lies in a quipu class, so a class with no hereditary member (n = 9) is not controlled. Referee's extra run: `SHORT=1 maverick_control.py 7 4 4 1`, 0/132.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/007/maverick_control.py 6 3 3 2` (84/84, 28 s); `... 7 4 4 3` (cap; output `maverick_control_n7.txt`); `SHORT=1 timeout 10m .venv/bin/python workshop/rounds/007/maverick_control.py 6 4 4 2` (0/84).

---

## E-068 — For `34x` at n = 14..17 (x = 4, 5, 7, 8, 9) the offset orbits pair as `o <-> hi - o`, i.e. `k = x + 3` and not `2x`; `346` is one orbit; `45x` has no reflection; `4046` is a size-paired reflection (`k = 11`), the parity translation belongs to `5046`/`5056` at odd n
*2026-09-30* · **20 of 20 cells (x = 4, 5, 7, 8, 9; n = 14..17) close their orbits with pair sum `s = hi = n - x - 3`; `k = x + 3` is this restated (n-independence adds nothing beyond 0 and hi sharing an orbit). The middle singletons are paired by equal orbit size only, not shown to be mirror images. `346` is one orbit at all n (the `333@0` orbit of E-065 not rechecked). `456`, `457` one orbit; `455`, `458`, `459` split by offset parity at some n only. `4046` at n = 12..16: `s = n - 11`, one lone large orbit at `hi`, size-paired; `5046` splits by parity at n = 13, 15. The E-060 statement that `4046` gives `{0,2},{1,3}` at n = 13 did not reproduce: `{0,2},{1},{3}`** → H-021, H-020, E-060, E-061, E-065 · *workshop round 007, experimentalist, refereed by skeptic*

**Limits.** x >= 10 and n >= 18 not run; `44x` not run; low-power cells `349@14`, `348@14`, `349@15` (one or two pairs); the mirror-join check of E-064 was not run on the singleton pairs.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/007/experimentalist_kd.py 14 34 3-9` (55 s; n = 16 about 8 min, n = 17 per x), `.venv/bin/python workshop/rounds/007/experimentalist_table.py`, `timeout 10m .venv/bin/python workshop/rounds/007/experimentalist_core4046.py 13 4046`.

---

## E-067 — Null test for the `|R| <= 4` reflection fits at n = 13: 39 of 109 fits are vacuous (one orbit) and of the 70 informative ones a shuffled partition fits with mean probability 0.73; the interior-core centre formula `s = first + last outside offset` (13/13) survives, the n = 15, 16 fits of the 12 E-060 cores survive
*2026-09-30* · **Null: keep each core's offsets and orbit block sizes, reassign offsets to blocks uniformly, run `fit()` (`d <= 6`). n = 13: 7 of 70 informative fits have P(null reaches `d <= d_obs`) < 0.05 (32 of 70 < 0.20). n = 15: 10 of 12, n = 16: 9 of 9 E-060 cores < 0.05. Interior cores: chance expects 4.1 of 13 centre hits, observed 13; naive product of the P values 1e-7, effective about 1e-3 to 1e-4 after collapsing the correlated families. E-061's "allI 45/62, allO 10/13" are padded by the 39 one-orbit cores (informative: 9/26 and 7/10; the allO 7/10 is a statement against null A, about 2.3 expected)** → H-021, E-060, E-061 · *workshop round 007, skeptic, refereed by theorist*

**Limits.** The uniform-label null ignores that neighbouring offsets tend to share orbits; no neighbour-aware or contiguous null is committed. The interior class and the 12 n = 15/16 cores were chosen after fitting at 13. Nothing about `k(c)` itself was tested.

**Reproduction.** `.venv/bin/python workshop/rounds/007/skeptic_null.py 300` (output `skeptic_null_out.txt`, seconds), `.venv/bin/python workshop/rounds/007/skeptic_null2.py 1000` (centre formula, n = 13; output `skeptic_null2_n13.txt`).

---

## E-064 — Over the 139 placed cores of the `--max-word 4` catalogue at n = 10, 12..16, orbit-plus-mirror classes refine the Coxeter-key classes in every core; the 20300 pairs of `4056`, `46`, `3355`, `3445` at n = 16 are one orbit and its mirror
*2026-09-30* · **For each n in 10, 12, 13, 14, 15, 16, all orbits of the 139 placed cores (484 words, 345 without a placement) close (limit 1500000). The partition of offsets by orbit-plus-mirror (orbits joined when one holds the mirror of an offset of the other) refines the key partition in every core: key finer than orbit+mirror in 0 cores, incomparable in 0, equal in 132 (n = 10), 130 (n = 12, 14, 16), 129 (n = 13, 15). The exceptions are key-coarser cores (7 at n = 10, 9 at even n >= 12, 10 at odd n), where the orbits are the two parity classes `{0,2,..}{1,3,..}`, each holding its own mirror, and the key is one class. At n = 16 the eight orbits of size 20300 (`4056` offsets 1, 2; `46` 3, 4; `3355` 4, 5; `3445` 2, 3) are two orbits X and its mirror X*: X holds `4056`@1, `46`@3, `3355`@5, `3445`@2; 28 pairs have equal or disjoint row sets (12 share all rows, 16 none). So E-062's unmerged middle pair and E-058's `{1,2}` of `4056` are one phenomenon, and "key pairs them, orbits do not" is resolved by the mirror join, not by an error of the key** → H-021, E-052, E-058, E-062, F-053 · *workshop round 006, toolsmith, refereed by skeptic*

**Limits.** Cores of at most 4 letters, n <= 16, closed orbits only. Key equality is not proved to imply orbit+mirror equality; different keys implied different classes in all 139 x 6 cases, so it is a prefilter candidate, not one to skip walks with. The key-coarser sets are not compared with the 7 cores of E-059. The ledger's `mirrors` is the loose reading (an orbit holding the mirror of any offset); that is exactly the orbit of the mirror, so the join is sound (referee). Referee reproduced the counts from the ledgers and the `--same-orbit` run (12 shared, 16 disjoint pairs).

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/006/toolsmith_orbitclass.py 14` (598 s); n = 15 (2 slices), 16 (4 slices): rerun until exit 0 (ledger resumes); `... --same-orbit` (about 6 min). n = 12 takes 48 s.

---

## E-065 — `k(33x) = 2x` and `d = x - 3` follow from a drift of the double mutation (`33x@o -> 33(x-1)@(o+1)`, label `c = x + o` conserved) and the self-dual seed `333`; the lower bound is derived, the upper bound only computed, and the argument fails for `44x`
*2026-09-30* · **(1) For x >= 5 the floating rule table does nothing to `33x` in the interior (orbit of size 1); for x = 4 only `334 <-> 333`. The interior moves come from the double mutation (F-032). (2) The double mutation on the middle relation sends `33x@o` to `33(x-1)@(o+1)` (and its dual the reverse), so `c = x + o` is constant along a chain `S_c = {33y@(c-y)}`; checked against `doubleMutation.rewritesOf` for all 33 placements at n = 14. `333` is self-dual and joins the chain `S_c` to the mirror of `S_{n-c}`, so `33x@o` and `33x@o'` share an orbit iff `{o + x, o' + x} = {c, n-c}`, i.e. `o' = o` or `o' = n - 2x - o`: `s = n - 2x`, `k = 2x`; offsets `o > s` have no partner and there are `x - 3` of them. (3) Test: for every `33x@o`, x = 3..8, at n = 14, 15, 16 (135 placements, all orbits closed), the `33y@p` rows in the orbit are exactly those with `y + p in {x + o, n - x - o}`: 135/135. Referee added n = 17, x = 7, o = 0..2: closed, matches (sizes 1363, 1449, 1449). (4) Extensions: `44x` has the same drift and self-dual seed `444`, so the argument would predict `k = 2x - 1`, `d = x - 4`, but every `44x`, x = 4..8 at n = 14..16 lies in one orbit with `333@0`: **the argument does not extend**. `345` is explained (`334@o -> 4444@o <- 345@o`), `34x` for other x is observed only (`k = x + 3`, x = 5 at n = 14..16; x = 7, 8 have no power), `45x` has no drift move** → H-021, H-020, F-032, F-051, E-061 · *workshop round 006, theorist, refereed by experimentalist*

**Limits.** The end link (the mirror chain joined back to the 33 rows by anchored/edge moves) and the upper bound "nothing else is in the orbit" are not derived. The test compares only the `33y` rows held in the orbit, not the whole orbit; "conserved label along a drift chain" is what is shown, not a conservation law of the orbit, which leaves the family. n = 14..16, x <= 8 (x = 9 excluded by a script cap), plus the referee's n = 17, x = 7. The 34x/44x/45x link runs at n = 16 were partly unfinished. Forward orbits of the reduced walk, not derived equivalence.

**Reproduction.** `.venv/bin/python workshop/rounds/006/theorist_chain.py 14` (53 s; also 15, 16); `theorist_rule.py` (1 s); outputs `workshop/rounds/006/theorist_chain_n{14,15,16}.txt`; extensions `theorist_link.py`, outputs `theorist_link_{34,44,45}x_n{15,16}.txt`.

---

## E-066 — The rejection at E-032 step 7 is located: the parent has a commutativity element `c = [8,6,4] + [8,10,4]` of `e_8 A e_4` with `c * (4->9) = 0`, where the map `p |-> (p beta)_beta` of Aihara-Iyama 2.32(b) / Ladkani 2.3(c) loses injectivity
*2026-09-30* · **Replaying `[4,6,4,6,9,4,4,6]` from the relation dual of `03033030`: the step-7 parent has 10 vertices, arrows 1>2>3>5>7>8, 8>6, 8>10, 6>4, 10>4, 4>9, the commutativity relation `8,6,4,9 + 8,10,4,9 = 0`, and the only arrow out of 4 is 4>9. The map on `e_i A e_4` has full rank for every i != 4 except i = 8 (dim 2, rank 1). Every step has gate = True; the side test is right=True / left=False at steps 1, 3, 6 and right=False at step 7. The rejection is not new (E-055, E-032, F-038); new is the explicit kernel element and its location** → H-015, H-010, E-032, E-055, E-057 · *workshop round 006, scholar, refereed by skeptic*

**Limits.** AI 2.32(b), Ladkani 2.3(c) and `tiltingPlus` being one map is stated by the author without a derivation (appears to be, not shown), so agreement is not independent corroboration of the code. "n = 10 is the first size with a commutative square into a vertex with one outgoing arrow" is untested (only n <= 7 was covered by E-057; a smaller instance at n = 5..7 was not built). The remark that CHZ Cor 3.6 in its path-wise wording would pass step 7 on this non-monomial parent rests on the repo's summary only (arXiv blocked), so **UNVERIFIED**; a flag is added to `literature/2509.12983`.

**Reproduction.** `timeout 10m .venv/bin/python workshop/rounds/006/scholar_step7.py` and `.../scholar_sides.py` (seconds each).

---

## E-063 — H-017 survives to depth 6 at n = 9, but the Coxeter polynomial cannot see (cords, relations); the Euler-form signature separates "outside every quipu class" for n = 8..11
*2026-09-30* · **n = 9: the quipus proved in the classes of the nine LNAs outside a quipu class have (cords, relations) exactly the ten recorded pairs at depth 4 and 5 (all nine LNAs) and depth 6 (`3033030`, `4444400`); none has relations <= cords. Per-LNA minima depend on depth (`4444400`: 2 at depth 5, 1 at depth 6). Coxeter trace is -1 for every tree algebra with monomial relations, and `c_{n-2}` is not a function of (cords, relations). Signature of `C + C^T`: for every LNA of length 8..11, `pos <= n-2` iff the LNA lies in no quipu class (n = 8: none outside; 9: 9; 10: 262; 11: 2647). Quipu => `pos >= n-1` holds by computation for n <= 12; the converse at n = 12 is unrun; at n = 13 the proof method breaks (`P^(1,0,3,0,1)_(1,1,1,1)`, one bad quipu of 241), but no LNA in that class is exhibited, so the criterion is not shown false there** → H-017, F-034, F-045, F-048 · *workshop round 004, maverick, refereed by skeptic*

**What it says and does not say.** Speculation level: tested on small cases. H-017 is not refuted and not proved. One relation (gldim <= 2) is a rank-2 perturbation of the form and lowers `pos` by at most 1, so a quipu in such a class needs at least 1 or 2 relations, independent of cords; it cannot give "relations > cords". The trace statement for hereditary trees is Happel's `a_1 = 1` (`literature/2509.02375-coxeter-coefficients-trees.md`, Theorem 1.1); its extension to monomial relations (arrows contribute `n-1`, a zero relation contributes 0 to `sum C_ij (C^-1)_ij`) is the round's own; the referee checked it on 3000 random monomial trees, n = 3..11, 0 failures. The signature criterion at n = 9 restates F-045's corank-2 set; new at 10 and 11 and small. **Not evidence:** the search from 16 below-diagonal polynomial-match candidates to depth 4 reached no LNA, with no positive control, so it neither supports nor refutes the sharp prediction. Not reproduced by the referee: n = 12 (over 10 min), n = 11, n = 14, 15 counts, and the depth 5/6 walks.

**Reproduction.** `.venv/bin/python workshop/rounds/004/maverick_corank.py 9` (also 8, 10, 11 about 4 min); `maverick_quipu_pos.py 13` (47 s); `maverick_coeff.py 6`; `maverick_reached.py 9 4` (107 s), `maverick_reached.py 9 6 3033030` (under 10 min); `maverick_census.py 9`, `maverick_euler.py 9`, `maverick_profile.py 9`, `maverick_verify.py 9 4 4 -1`. Referee reproduced the signature tables n <= 10, the n = 13 bad quipu, and `c_{n-1} = 1`.

---

## E-062 — For the 12 cores of E-060, the centre `s = n - k` and `d` are unchanged at n = 15 and 16; at 16 three cores (`46 3355 3445`) lose the fit through an unmerged middle pair of equal size
*2026-09-30* · **n = 15: 12 of 12 fit, every core has the `k` and `d` of n = 13 and 14, all orbits closed. n = 16: 9 of 12 fit with the same `k`, `d`; `46 3355 3445` have no fit under the round-002 rule (`|d| <= 6`, merged middle pair required): every orbit is closed under the centre `16 - k` except the two middle offsets, which are two singleton orbits of the same size 20300 (`46 {3}{4}`, `3355 {4}{5}`, `3445 {2}{3}`). So the 12 pair at 13, 14, 15 and mostly at 16: the parity split of E-059 does not carry over to them. Slide rule `s` = first + last outside offset at 15: 8 of 12 (fails `4045 4506 4556`, and `334` all inside); interior blocks 5 of 5 (`504`, m = 6, one case with m >= 6)** → H-021, E-058, E-059, E-060, F-053 · *workshop round 004, experimentalist, refereed by scholar*

**Limits.** "No fit" means the middle pair is unmerged under the round-002 rule; the equal sizes are the only sign of the reflection, so this is not a loss of the symmetry as such. The same size 20300 and the same kind of unmerged middle pair are recorded at n = 16 for `4056` in E-058, so the phenomenon is not new; new is that it hits 3 of the 12 pairing cores at 16 with the centre still `n - k`. Whether 20300 is one shared orbit across these cores was not checked. n = 17, 18 unrun; the other 127 cores unrun at 15 and 16; no mechanism. Referee re-ran `46` at 16 (178 s) and the shift table: identical.

**Reproduction.** `.venv/bin/python workshop/rounds/004/experimentalist_shift.py workshop/rounds/004/experimentalist_census_12cores_n15.jsonl workshop/rounds/004/experimentalist_census_12cores_n16.jsonl` (seconds); the censuses: `workshop/rounds/002/experimentalist_census.py 15 --cores <core>` (about 3 min per core), 16 (160 s alone); run at most 4 at once.

---

## E-061 — The interior / end-touch split is a committed column and explains none of the 7 failures; for the family `33x`, `k = 2x` and `d = x - 3` at n = 13..17, x = 3..7
*2026-09-30* · **`theorist_shortfall.py` now prints the class (allI, allO, int, end0, endhi, endboth), the prediction `s` = first + last outside offset and every slide-consistent centre. n = 13, 109 cores with a fit: int 13/13, end0 7/10, endhi 10/11, allO 10/13, allI 45/62 (E-060's 13/13 and 17/21). The fitted centre is slide-consistent in 109/109 cores, so the slide rules centres out and never selects one. End-touching failures: the orbits take the smaller consistent centre for `4045 3556 4556` and the larger for `4506`, so there is no min/max rule. All-outside failures `4046 5046 5056`: the slide is silent; orbits `{0,2},{1,3}` look like a translation of period 2, not a reflection (unverified, no null). `33x`: `k = 2x`, `hi = n - x - 3`, `d = x - 3`, the `d` top offsets singleton orbits; 20 fits at n = 13..17 for x = 3..6, no deviation; chair added `337` at n = 15, 16, 17 after the referee's request: `k = 14`, `d = 4` at all three (at 15 the fit is one pair, low power)** → H-021, H-020, E-060, F-053 · *workshop round 004, theorist, refereed by skeptic*

**Limits.** `k = 2x` is a line through five families (x = 3..7), each checked at several n, not a law; it is a description, and the mechanism is not found. At n = 13, `|R| <= 4` fits (one adjacent pair) have almost no power, which weakens both the pass counts and the reading of the 3 all-outside cores. The fitted centre is unique for 59 of the 109 cores (the rest have several closing centres, so `d` is not a fit artefact only for those 59); slides are n = 13 only. Referee re-ran all three scripts: byte-identical output.

**Reproduction.** `.venv/bin/python workshop/rounds/004/theorist_shortfall.py`, `theorist_k.py`, `theorist_33x.py` (seconds, from the committed data); the `33x` censuses: `workshop/rounds/002/experimentalist_census.py $n --cores 33,333,334,335,336 --out ...` for n = 14..17 (about 2 min in parallel); `337`: same script with `--cores 337` for n = 15, 16, 17 (`workshop/rounds/004/chair_337_n{15,16,17}.jsonl`, about 5 min).

---

## E-060 — For 13 of 13 cores whose outside block is interior, the reflection centre is the first plus the last outside offset (`d = t - h`), at n = 13; 3 of 3 at n = 14
*2026-09-30* · **Slides computed for the 395 orbit walks of E-056 (n = 13, 139 cores; no undecided, no orbit with mixed verdicts). Cores with a reflection fit and a mixed slide whose outside block touches neither end: 13, all 13 have fitted `s` = first outside offset + last outside offset, i.e. overhang `d = hi - s = t - h` (tail minus head, H-020). Mixed slide with the block touching an end: 17 of 21 fit (fails `4045 3556 4506 4556`). All outside: 10 of 13 (fails `4046 5046 5056`). All inside (62 cores): the slide is silent about `d`, not tested. n = 14, 12 chosen cores: interior 3/3 (`45 46 504`), end-touching 4/7, as at 13** → H-021, H-020, F-053, E-052, E-056 · *workshop round 003, theorist, refereed by skeptic*

**The result and its limits.** If a centre `s` pairs offsets and the reflection maps the outside block onto itself, then `s` = first outside + last outside and `d = t - h`; that step is a short argument from "the verdict is constant on an orbit" (E-052) and is close to a tautology. The content is that the fold happens in 13 of 13 interior cases. Failures occur only among end-touching blocks (4 of 21 fail), where several centres are slide-consistent and the orbits choose one differing by 1; this is one direction only. F-053 already states the case of `45` (`d = 1`, n = 12..17). Interior blocks at n = 13 have `m` = 2 to 4 outside offsets (5 only for `45` at 14); no interior-block core other than `45` has been tested with `m >= 5`. The n-independence of `s = n - k(c)` and `d` (the restated H-021') rests on 12 cores over one step, 13 -> 14 (all 12 fit again, `s` up by 1, `d` unchanged); the 14 sample is chosen, not a census. The restatement covers only cores that pair; 30 of 139 at n = 13 have no fit and are outside it. Not explained: why the outside block folds onto itself; the formula for `k(c)`.

**Referee's required change not done:** the interior / end-touching split is a described snippet, not a committed column; `theorist_shortfall.py` prints the slide string but not that flag. The counts above were recomputed independently by the referee and agree.

**Reproduction.** `.venv/bin/python workshop/rounds/003/theorist_shortfall.py` (seconds, from `theorist_slides_n13.jsonl`); `.venv/bin/python workshop/rounds/003/theorist_rule.py workshop/rounds/003/theorist_census_n14_sample.jsonl workshop/rounds/003/theorist_slides_n14_sample.jsonl`. The slides: `theorist_slides.py` in 4 shards (about 5 min each).

---

## E-059 — The 7 cores `344 366 4044 4403 4404 4405 4605` pair by a reflection at every even n from 12 to 18 and at none of n = 13, 15, 17 (strict mirror instead)
*2026-09-30* · **For these 7 cores only: n = 12, 14, 16, 18: 7 of 7 pair (fit `d` = 0 for `344 366 4403 4605`, 1 for `4044 4404`, 2 for `4405`, the same at every even n) and hold no strict mirror; n = 13, 15, 17: 7 of 7 have no fit and hold a strict mirror. n = 8..11: mixed, few offsets, weak. All 289 orbit walks closed** → H-021, E-056, F-053 · *workshop round 003, experimentalist, refereed by theorist*

**The run.** E-056 left the 7 as "mirror without reflection" at n = 13. Census (`experimentalist_census.py`, 3 s per core per n) and the round-002 fit, on the 7 cores for n = 8..18 (the n = 13 rows are E-056's). At n = 14, `344` has orbits `{0,7} {1,6} {2,5} {3,4}` (four clean pairs); at n = 13 and 15 the offsets that would fold split into singleton orbits of equal size (n = 13: 50 and 50; n = 15: 64 and 64). That equality is an observation; that the singletons "would pair" is interpretation, and no mechanism is offered. So the 7 are not a family that never pairs; the n = 13 defect is an odd-n effect **for these 7, chosen because they show it at 13**. Nothing is known for the other 132 cores at any n other than 13 (the n = 14 census is parked for the overnight); "no fit" is a statement about the fit rule (`d <= 6`), not the walks. The 8 strict-mirror-true cores of E-056 are these 7 plus `406`, which pairs; the "eight cores pair at n = 13 and 14" of E-052 is a different eight.

**Reproduction.** `.venv/bin/python workshop/rounds/003/experimentalist_table.py workshop/rounds/003/experimentalist_census_7cores_n8_18.jsonl` (seconds; per-core `d` and `s`: `.venv/bin/python workshop/rounds/002/experimentalist_fit.py` on the same file); one core: `.venv/bin/python workshop/rounds/002/experimentalist_census.py 14 --cores 344`. Reproduced by the referee for `344` and `4405` at n = 14, 15, 16.

---

## E-056 — H-021's mirror clause is no "exactly when" on any of three readings, at n = 13 over all 139 cores
*2026-09-30* · **Loose: P => mirror (109/109), mirror => P fails for 20. Strict (mirror of `c@p` in an orbit not holding `p`): true for 8 cores, false for 108 of the 109 pairing cores (by construction: the mirror acts inside the orbit); strict => P fails for 7 (`344 366 4044 4403 4404 4405 4605`). strict2: P => strict2 holds, strict2 => P fails for the same 7** → H-021, F-053, E-053 · *workshop round 002, experimentalist, refereed by skeptic*

**The run.** Round 001's census re-walked (395 orbit walks, all closed, none capped) and the fit
rerun on it: P = "orbits pair by a reflection `o <-> s - o`, overhang `d <= 6`" (109 cores,
`d` histogram 62/28/17/2, 30 no fit). Rows (P, loose, strict, strict2): 108 (y,y,n,y), 1 `406` (y,y,y,y),
10 `30xy 330x` (n,n,n,n), 13 (n,y,n,n) `3033 3034 3044 3303 3346 3456 3466 3566 3606 4056 4406 506 6005`,
7 (n,y,y,y). `3346` and `4056` hold only their own mirrors, so under the strict reading they agree with
"no pairing" and H-021's text about `3346` is right there; round 001's claim to the contrary held only
for the loose reading.

**What it does not show.** Anything at n >= 14 (the 7 survivors are not checked there); that they are not
"pairing with a defect" (`344` pairs `{0,6}`, `{1,5}` and swaps `2`, `4` by mirror only); two-cluster
words. The fit's `d <= 6` is generous, so only "no fit" is strong. Referee: "108" is a tautology of
the definition, not an empirical result.

**Reproduce.**
```
.venv/bin/python workshop/rounds/002/experimentalist_fit.py workshop/rounds/002/experimentalist_census_n13.jsonl   # seconds
# census: 4 shards of workshop/rounds/002/experimentalist_census.py 13 --shard k/4 --out ..., ~5 min wall
```

---

## E-057 — On non-monomial parents Ladkani 2.3(c) rejects exactly the gate-refused vertices (n = 5, 6, 7)  *(annotated round 014: the "none" in this entry is depth-limited; E-084 finds gate-admitted rejections at n = 6 (distance 8) and n = 7..9 (distance 5-7), all guard-refused)*
*2026-09-30* · **n = 5 depth 2: 8 refused / 0 gate-admitted rejections; n = 6 depth 3: 214 / 0 (714 accepted, all Cartan-congruent); n = 7 depth 2: 432 / 0 (1,008 accepted); 0 skipped parents; negative control at starts 42/42, 168/168, 660/660 refused-and-False, 70, 252, 924 admitted-and-True** → H-015, E-055 · *workshop round 002, scholar, refereed by skeptic*

**The run.** Unguarded BFS over distinct algebras from the LNAs, all vertices with an arrow out, on the
non-monomial parents (commutative squares `b1 b2 = b3 b4`; 12, 194, 240 of 98, 910, 1,188 algebras).
"Gate" = `mutationIsPossibleAtVertex`; "guard" = Coxeter key of the reduced child equals the parent's.
Skipped parents (directed cycle or no base key) are now counted: 0.

**What it does not show.** A second gate-admitted rejection: there is none, so E-032's ALARM step 7
remains the only one, and `tiltingPlus` adds nothing beyond the gate in this range. Nothing for the
guard (H-015), or for n >= 10 where the gate and 2.3(c) can differ (F-038). Round 001's own runs did not
count skipped parents (expected 0, unverified).

**Reproduce.**
```
N=6 DEPTH=3 timeout 10m .venv/bin/python workshop/rounds/002/scholar_nonmono.py   # about 20 s
N=7 DEPTH=2 timeout 10m .venv/bin/python workshop/rounds/002/scholar_nonmono.py   # about 46 s
```
(needs `workshop/rounds/001/scholar_h015.py`.)

---

## E-058 — At n = 16 the offsets `{1,2}` of `4056` are two mirror-image orbits, though the Coxeter key pairs them
*2026-09-30* · **`4056`, n = 16: walked orbits `{0,3}` 19798 rows, `{1}` 20300, `{2}` 20300, `{4}` 77735, `{5}` 8134, `{6}` 1416, all closed; key classes `{0,3}{1,2}{4}{5}{6}`. Orbit(1) holds the mirror of the start of 2 and vice versa, neither holds its own mirror** → F-053, F-026, H-021, E-054 · *workshop round 002, skeptic, refereed by theorist*

**The run.** The round-001 walk repeated at n = 16 (252 s; theorist reran in 237 s with identical output)
plus a mirror test (`skeptic_mirror16.py`). Since `mirrorRow` is the relation dual (F-026) and the
mirror keeps the derived class, offsets 1 and 2 are derived equivalent by inference from one mirror
membership, not by a walk joining them. So "orbit = key class" (five cores at n = 13..15) is false at
n = 16, and F-053's "each pair is one self-dual orbit" fails for `{1,2}` there. This is the first
counterexample, from one core at one length. Also corrected: the `6600066` key sums are n - 13, and the
n = 24 scan of 585 words splits 309 one centre / 155 key class of size >= 3 / 121 no equal keys
(saved: `workshop/rounds/002/skeptic_scan_n24.txt`); E-054's "7 / 123" was wrong. "Onset" remains a
description with no mechanism.

**What it does not show.** Pairing for any core at n > 16 (key only); whether other cores split like this.

**Reproduce.**
```
timeout 10m .venv/bin/python workshop/rounds/002/skeptic_cmp.py 16 4056     # 252 s
timeout 10m .venv/bin/python workshop/rounds/002/skeptic_mirror16.py        # about 2 min
.venv/bin/python workshop/rounds/002/skeptic_scan.py 24                     # about 3 min
```

---

## E-055 — Every step the gate admits at n = 6, 7 is a tilting mutation by Ladkani's exact criterion
*2026-09-29* · **61,718 admitted steps (n = 6 depth 6, n = 7 depth 4, unguarded walk): all pass Prop. 2.3(c) of arXiv:1001.4765 and Cartan(child) = r C rᵀ; the E-032 ALARM step 7 fails both** → H-015, F-038, E-032 · *workshop round 001, scholar, refereed by skeptic*

**What prompted it.** H-015 asks for a second invariant along guarded paths.
Ladkani's Prop. 2.3(c) is an iff for a step to be a tilting complex; the
literature note recommended it in place of `isMutable`, and no entry had run it.

**The run.** For each algebra reached, every gate-admitted step, test 2.3(c) (a
rank computation, `tiltingPlus` in `workshop/rounds/001/scholar_h015.py`) and the Cartan
congruence. n = 6 depth 6: 9,476 algebras, 29,822 steps; n = 7 depth 4: 8,988
algebras, 31,896 steps; 0 failures. "Guard" in the script means Coxeter key
equal after the rewrite; it never fired, so guard-passing = gate-passing here.
Parents with a directed cycle or no base key are skipped, uncounted. The ALARM
path (`RD=1`, step 7, vertex 4): gate admits, 2.3(c) False, congruence False,
key moves. Referee negative control: at every vertex of every n = 4, 5, 6 LNA and
its relation dual, `tiltingPlus` is False at all gate-refused vertices with an
arrow (10/42/168) and True at all admitted ones (20/70/252).

**What it does not show.** Nothing about H-015: the guard is inert at these
sizes (F-038), so this says the gate is already sufficient at n <= 7, depth <= 6.
Only one non-monomial negative case (the ALARM step) tests `tiltingPlus`; the
n = 7 counts were not re-run by the referee; n = 8 depth 2 was cut unfinished.
`tiltingPlus` is not in the library and untested. The region where the guard
matters (n >= 9, depth >= 6) is untouched.

**Reproduce.**
```
timeout 10m .venv/bin/python workshop/rounds/001/scholar_h015.py 6 --depth 6 --unguarded   # 3 min 40 s
timeout 10m .venv/bin/python workshop/rounds/001/scholar_h015.py 7 --depth 4 --unguarded   # 6 min 10 s
RD=1 timeout 10m .venv/bin/python workshop/rounds/001/scholar_h015_f038.py                  # seconds
```

---

## E-054 — The Coxeter key separates the offsets of `3346` at every length, so its five orbits are real
*2026-09-29* · **`3346` never pairs (key: n = 12..20, 30, 40; orbits: 13, 14, 15); `4056` pairs offsets `o ↔ n-13-o` with the last three alone (orbits to n = 15, key beyond)** → F-053, H-021, H-020, E-052 · *workshop round 001, skeptic, refereed by theorist*

**The run.** Reduced-walk orbit partition against the partition of offsets by
`coxeterTables.lnaCoxeterKey`, for `3346` (n = 13, 14, 15), `4056` (13, 14, 15),
`45` (13, 14) and `350066`, `6600066` (13): all orbits closed (largest 17905 rows),
orbit partition = key partition in every row. Distinct keys prove distinct
orbits (the moves preserve the key), so the exceptions are not caps, walk or
gauge; equal keys only bound an orbit above. Key-only partitions: `4056` pairs
`n-13`-sum with a tail of three alone for n = 14..20, 30; `350066` has nothing
to 16 and `{4,5}` at 17; `6600066` is `{0}` at 13, `{0,1}` at 14, then pairs
of sum n - 13 (the submission said n - 14; corrected by the referee). Scan at
n = 24 of 585 words over digits 0, 3..9 (<= 4 letters): 309 with one pairing
centre, 7 with a key class of three offsets, 123 with no equal keys.
*Corrected in round 002 (E-058): the true split at n = 24 is 309 one centre,
155 with a key class of three or more offsets, 121 with no equal keys; the "7"
was the number examined, not a count.*

**What it does not show.** That the pairing of `4056` beyond n = 15 is an
orbit fact (key-only, a necessary condition); a mechanism for the onset length
(13 for `4056`, 8 for `45`, 6 for `3344`), which is description, not cause; that
the other four H-020 failures are the same (only two were reachable). The
scan's summary output was not saved.

**Reproduce.**
```
timeout 10m .venv/bin/python workshop/rounds/001/skeptic_cmp.py 14 3346 4056 45   # about 95 s
.venv/bin/python workshop/rounds/001/skeptic_cox.py 3346 12 14 20 40               # seconds
timeout 9m .venv/bin/python workshop/rounds/001/skeptic_scan.py 24                 # about 3 min
```

---

## E-053 — Census of all 139 single-cluster cores of `--max-word 4` at n = 13 under the reduced walk
*2026-09-29* · **all orbits close; "the orbit holds a mirror" is true of 129 of 139 cores, so it does not discriminate; 20 cores hold a mirror and fit no reflection** → H-021, F-053, E-052 · *workshop round 001, experimentalist, refereed by skeptic (minor revision outstanding)*

**The run.** Per core, `freeMoves.orbitReport(free = REDUCED, limit = 1500000)`
at each offset, recording offsets and mirrors held (536 orbit walks, 0 caps,
largest 4217 rows). Best-fit reflection `o ↔ s - o` allowing an overhang `d`:
d = 0: 62 cores, 1: 28, 2: 17, 3: 2, no fit (<= 6): 30. Of the 30, ten (`30xy`,
`330x` words) hold no mirror and pair nothing; twenty hold a mirror and no
reflection (`344 366 506 3033 3034 3044 3303 3346 3456 3466 3566 3606 4044 4056
4403 4404 4405 4406 4605 6005`). Every core that pairs holds a mirror. At n = 14
only seven cores were run (`45 506 3033 3035 3344 3346 4056`, 111 s): `45`, `3344`
match E-052; `3346` is five singleton orbits each holding its own mirror.

**Caveat found by the referee.** "Holds a mirror" was scored including a
mirror at the orbit's own offset. Under H-021's literal wording that refutes the
"exactly when" and contradicts its "no orbit holding a mirror" for `3346`. Under
the stricter reading (the mirror of `c@p` lies in the orbit of `c@q`, `q != p`)
`3346`, `4056` and `3033` would agree with "no pairing" and only cross-offset
cases such as `344 4044 4403` survive as counterexamples; that count has not been
run. The overhang fit is loose (d up to 6 on 4 to 7 offsets) and its code is not
in the repository; only the raw orbits reproduce from the cited script.

**Reproduce.** Script `t1.py` was left in a scratchpad (not saved); the recipe is
`batch._singleCores(4, 6, False)`, `batch._rowFor(13, c, o)`,
`freeMoves.orbitReport(13, row, free = freeMoves.REDUCED, limit = 1500000)`,
`freeMoves.mirrorRow`. About 25 min on 4 processes.

---

## E-052 — The outside band of a core is its reflection pairs
*2026-09-23* · **under the reduced walk the outside offsets of `45` fall into closed orbits `{o, n - 8 - o}`, one per pair, at every length from 12 to 17; the same pairing holds for seven more cores and fails for `3346`** → F-053, H-021, H-020, F-051

**What prompted it.** Last night's shared census at `n = 16` called `45` at
offsets 2, 4 and 6 undecided and 3, 5 and 7 outside, and the outside ones after
the first closed in a few rows each -- as though some offsets shared an orbit
and others did not. F-051 had `45` at 1 and 5 in one plain orbit at `n = 14`
and read it as "the interior is one orbit".

**The run.** For each offset of a core, the reduced walk alone
(`freeMoves.orbitReport(free = REDUCED)`) to closure or a certificate, and which
other offsets' rows the orbit holds. With the faster rule lookup of E-051 this
is seconds at 13 and minutes at 17:

| n | `45` orbits of the outside band (rows) | inside |
|---|---|---|
| 12 | {1,3} 740 · {2} 1766 | 0, 4, 5 |
| 13 | {1,4} 1127 · {2,3} 4217 | 0, 5, 6 |
| 14 | {1,5} 1636 · {2,4} 11820 · {3} 2179 | 0, 6, 7 |
| 15 | {1,6} 2290 · {2,5} 18157 · {3,4} 18416 | 0, 7, 8 |
| 16 | {1,7} 3114 · {2,6} 77735 · {3,5} 16174 · {4} 26919 | 0, 8, 9 |
| 17 | {1,8} 4135 · {2,7} 122673 · {3,6} 57828 · {4,5} 62128 | 0, 9, 10 |

Every orbit closed, so `45` is `i o^(n-9) i i` at every length to 17 -- H-020's
head 1, tail 2 -- under the reduced walk, where the census had it undecided at
16 and 17. (At `n = 18` the census trial of E-051 closed `45@2`, `@3` and `@4`
at 355328, 44100 and 157281 rows.) Each orbit holds exactly two offsets, `o`
and `n - 8 - o`, and **its own mirror**: the row of `45@o` and the row of
`504@(n - 7 - o)`. The inside offsets pair the same way (`0` with `n - 8`), and
the one offset the reflection cannot reach, `n - 7`, is the extra inside at the
sink. That is why `45`'s tail is its head plus one.

**Other cores, at `n = 13` and 14:**

| core | pairs sum to | orbits at 14 |
|---|---|---|
| `46`, `56` | `n - 9` | {1,4} 11820 · {2,3} 7758 |
| `505` | `n - 8` | {1,5} 11820 · {2,4} 790 · {3} 2471 |
| `555` | `n - 8` | {1,5} 1636 · {2,4} 11820 · {3} 2179 |
| `556` | `n - 10` | {0,4} 1636 · {1,3} 11820 · {2} 2179 |
| `3344` | `n - 6` | {2,6} 1636 · {3,5} 11820 · {4} 2179 |
| `3345` | `n - 9` | {0,5} 213 · {1,4} 3988 · {2,3} 3586 |
| `666` | `n - 9` | {1,3} 2116 · {2} 962 *(13 only)* |
| `4056` | -- | {0,1} 7758 · {2} 11820 |
| `3346` | **none** | five orbits, one per offset |

So the pairing is not special to `45`. The same few orbits recur across cores
of one length -- 1636, 11820 and 2179 rows at `n = 14` hold `45`, `555`, `556`,
`3344` and part of `505` and `46` -- which is F-051's "many cores share it" seen
one reflected pair at a time. `3346` is the exception: every offset its own
orbit, none holding another offset. `4056` pairs 0 with 1 and leaves 2 alone.

Reproduce: `freeMoves.orbitReport(n, batch._rowFor(n, core, o), free =
freeMoves.REDUCED, limit = 1500000)` for each offset, and test which other
offsets' rows lie in the returned set. `tests/test_reduced_walk.py` pins `45`
at 13.

---

## E-051 — Two nights, one of them mostly wasted, and a walk thirty to eighty times faster
*2026-09-23* · **the shared census stalled from `n = 15` on walks capped at 20000 rows; the rule lookup was 90 percent of a walk and is now 59x faster; the census reuses other ledgers' verdicts and never walks a capped class twice** → H-020, H-019, H-018, F-053

Two `overnight.py` lines from `OVERNIGHT.md`, both run to their nine hours.

**The night of 2026-09-22** ran the line E-050 called superseded (plain walk):
`cores 13` and `cores 15 --max-word 4` finished (1512 and 2343 placements; 0
and 29 undecided), `cores 17 --max-word 2 --pair-word 2 --gaps 5,6` finished
(762 placements, 18.8 core-hours on three workers), and the two samples went on.

**The night of 2026-09-23** ran E-050's "tonight" line: ten shared censuses and
two samples, one worker each.

| census (shared) | done / placements | inside | outside | undecided | where the 9 h went |
|---|---|---|---|---|---|
| `13 --max-word 4` | 1457 / 1457 | 955 | 502 | 0 | 0.55 h, done |
| `14 --max-word 4` | 1850 / 1850 | 1135 | 715 | 0 | 2.0 h, done |
| `14 --max-word 5 --gaps ,` | 3330 / 3330 | 2048 | 1282 | 0 | 1.8 h, done |
| `15 --max-word 4` | 1682 / 2244 | 1071 | 583 | 28 | 4.1 h on the 28 undecided |
| `16 --max-word 4` | 410 / 2637 | 338 | 36 | 36 | 7.4 h on the 36 undecided |
| `16 --max-word 5 --gaps ,` | 410 / 4634 | *the same 410* | | | *all of it duplicated* |
| `17 --max-word 4` | 371 / 3031 | 318 | 23 | 30 | 6.8 h undecided |
| `18 --max-word 4` | 172 / 3424 | 131 | 3 | 38 | 8.8 h undecided |
| `17 --max-word 2 --pair-word 2 --gaps 5,6` | 286 / 762 | 242 | 19 | 25 | 6.0 h undecided |
| `18 --max-word 2 --pair-word 2 --gaps 5,6` | 186 / 963 | 148 | 4 | 34 | 7.8 h undecided |

**What went wrong, three ways.**

* **The cap, again.** From `n = 15` up, 75 to 97 percent of each census went on
  placements that hit the 20000-row cap -- 520 to 860 s each -- and got no
  answer. The shared walk walks plain *then* reduced, each to the cap, so an
  undecided placement costs two capped walks, the second at the reduced walk's
  higher price per row. It is the reduced phase that caps: `45@2` at `n = 16`
  closes its plain orbit in 133 rows and its reduced one at 77735 (E-052). The
  plain census had these as `outside` at 20000; the shared one could not.
* **The same capped class, twice.** `45@2` and `45@6` at `n = 16` are one class
  (E-052), and each was walked to the cap in turn: a capped walk left nothing
  behind for a later start to find.
* **One census inside another.** A `--max-word 5` catalogue begins with the
  whole `--max-word 4` one, sorted by length, so the two `n = 16` jobs did the
  same 410 placements side by side, verdict for verdict, and two of twelve
  cores did nothing new all night. At `n = 14` the `--max-word 5` run redid
  most of what the `--max-word 4` run did beside it.

**The fix that mattered: the rule lookup.** A profile of a reduced walk at `n =
17` put 90 percent of the time in `lnaMoves.rewritesOf`, which tried each of the
1844 rules at every window: about 5400 `matchesAt` calls per row, each building
two `set`s. A rule's window holds exactly its left-hand side, so its first
relation sits on one of the row's own relations. `lnaMoves.matchingWindows`
indexes the rules by that first relation's length, tries only the windows the
row's relations anchor, and checks a match by comparing the run of relations
from there (starts and ends both increase, so the relations meeting a window
are a run of the list). Same pairs, same order, so every walk visits its rows
in the order it always did.

| | before | after |
|---|---|---|
| `rewritesOf`, 4500 random rows at `n = 4..18` | 19.4 s | 0.33 s (**59x**), identical lists |
| plain walk, `45@3` at `n = 17` | 67 rows/s | **5600 rows/s** |
| reduced walk, same start | 88 rows/s | **2500 rows/s** |
| shared census `13 --max-word 4`, one core | 1987 s | **32 s**, 0 of 1457 verdicts differ |
| shared census `14 --max-word 4`, one core | 7221 s | **117 s**, 0 of 1850 verdicts differ |

**Three changes to the census, with it.**

* **A capped walk's class is remembered** (`SharedWalk.isCapped`): a later start
  an earlier capped walk passed through is `undecided` by `shared cap` at no
  cost. Its rows are the earlier start's class, so the answer can only be that
  class's.
* **Verdicts are reused across the ledgers of one length** (`--no-reuse` turns
  it off). A census looks in every other `cores-n<n>-*` ledger for the row or its
  mirror: `inside` from any of them (a certificate is a certificate), `outside`
  only from the same walk family (shared and reduced agree, plain is weaker),
  `undecided` only from the same family at limits at least as large. Under the
  shared walk every reused `inside` also seeds the walk, so a later walk stops
  at its class -- most of what a resumed census used to lose with its memory.
  Rows it takes say `"by": "reused"` and name the ledger; the other ledgers are
  re-read every ten minutes.
* **`--min-word`**, a filter: `--max-word 5 --min-word 5` asks the five-letter
  words only, into the same ledger as the whole `--max-word 5` census.

Each row now carries `maxRssMB`, since the nights run on WSL's share of the
laptop's memory (about 6.8 GB of 13.7).

**Sized by a trial** of twelve minutes per census, one core each, at
`--orbit-limit 500000`, reusing last night's ledgers: `16 --max-word 4` reached
1476 of 2637 placements (458 walked, 1018 reused) with **no** undecided; the
`n = 17` pairs at gaps 5,6 walked 150 with none undecided; `n = 18` reached
`45@2`, whose reduced orbit **closes at 355328 rows** in 510 s. Peak memory
235, 227 and 280 MB. Every inside answer on record for a single-cluster
placement came within 2300 rows (10700 for a pair), so the cap buys outside
verdicts only, and 500000 is enough for every `45` orbit to `n = 18`.

**What the nights say, read with the rest.**

* **H-020 under the reduced walk.** At 13 and 14 both shared censuses are
  complete, and no single-cluster slide holds an inside in its interior. Over
  every pair of lengths from 13 to 18, **1192 comparisons** of complete slides,
  1186 hold; the 6 that fail are all at `n = 13`, for words of six or seven
  letters whose slide there is one or two offsets long (`350066` `ii` → `iio`,
  `6600066` `i` → `oo`). `35`, `36` and `404` are inside at every offset at every
  length from 13 to 18 -- the reduced walk removes their plain outside bands, as
  E-049 found at 11 and 12. `4056` is `ooii`, `oooii`, `ooooii` at 13 to 15:
  F-052's `oii` at 12 was below the length where the law starts.
* **Five-letter cores.** `14 --max-word 5 --gaps ,` (last night) and `13
  --max-word 5 --gaps , --orbit-limit 200000` (a check of the reuse today: 2680
  placements in under two minutes, 1087 of them reused) are both complete, and
  neither has a slide with an inside in its interior.
* **Two clusters with a free arrow** (H-018, E-047). Plain walk, `n = 17`, gaps
  5 and 6, complete: when both halves are inside alone the pair is inside **330
  of 330** times. The plain walk also rescued 62 pairs with an outside half --
  and **every rescued half is a `35` or a `36`** (36 and 26), which the reduced
  walk places at every offset. Under the shared walk at 17 and 18 so far, 217
  pairs with both halves inside are all inside and **no half is rescued**: each
  pair not inside has a half that is outside or undecided alone. E-047's
  rescues look like an artifact of the one-way free move.
* **The leftover rate** (H-019). `sample 15 --orbit-limit 100000`, 1444 draws:
  **57.4% +- 1.3**, every leftover orbit **closed** (the largest at 60852 rows,
  128 of the 942 closed walks past 20000). E-045's 48.8% closed plus 9.2%
  capped is 58.0%, so at `n = 15` the capped draws were all outside. `sample
  17` at 20000, now 1131 draws: 49.6% closed, 22.9% capped, 72.5% in all.

Reproduce:

```bash
python batch.py cores 16 --max-word 4 --summary
python batch.py cores 17 --max-word 2 --pair-word 2 --gaps 5,6 --walk plain --summary
python batch.py sample 15 --count 20000 --orbit-limit 100000 --summary
```

The half-by-half table and the cross-length comparison join the ledgers on row
names; neither is printed by `--summary`.

---

## E-050 — A census that shares what each walk settles
*2026-09-22* · **the reduced walk's verdict on every placement at `n = 11` and 12, at a thirtieth of its cost and a ninth of the plain census's; sharing does not survive being split over workers** → F-052, H-020, F-051

**The idea.** Asked in response to E-049: once one placement is settled, all
its aliases and its mirror are settled too, so the census need not keep to the
reduced cores; it should use every move to the fullest. Every move, the free
move and the relation dual is an equivalence, so what a walk learns is about a
*class*. "Inside" spreads both ways: every row a walk passes through is in its
start's class, so a certificate for any of them is a certificate for all.
"Outside" is directional, but a closed orbit's rows are each known to reach
nothing, so a later walk need not expand them (F-051's cache).

**First, the question as asked: are some aliases faster?** From E-049's
per-placement data, grouping each alias by where its two-arrow relations sit
relative to the first relation of its core:

| where the `2`s are | aliases | inside with the bare core outside | both inside: alias faster |
|---|---|---|---|
| one vertex before the core (`-1`) | 285 | 15 | 8 of 218 |
| two before (`-2`) | 70 | 2 | 0 of 63 |
| two and one before (`-2,-1`) | 70 | 2 | 4 of 63 |
| three after (`+3`) | 13 | 0 | 0 of 13 |

So no decoration makes an answer the bare core already gets come faster -- 12
of 357, and on average the alias walks slightly *more* rows. What a `2` right
in front of the core does is reach answers the bare core does not reach at all.
The reduced walk's slowness is not the aliases' doing either: on the classes
that come out inside it spent 430 s at `n = 11` and 1450 s at 12, where the
fastest plain member of each class took 94 s and 166 s in all. It was its order
of search and the price of each step.

**The shared walk.** `freeMoves.SharedWalk`: a union-find over classes keyed by
`min(stripped row, its mirror)`; a set of rows known to lie in closed orbits
without a certificate, one per phase; each unit walks plain, then reduced only
if the plain walk closed. A walk stops at the first row whose class is inside
and skips every row already in a closed orbit. A unit whose class a later unit
shows to be inside is listed in that unit's `promotes`, and `--summary` applies
it. One process, catalogue order, every placement including aliases, mirror
pairs asked once:

| n | plain census | reduced census | shared, plain phase | shared, reduced phase | agrees with reduced |
|---|---|---|---|---|---|
| 11 | 464 s | 950 s | 18 s | 14 s, on 146 left | **709 of 709** |
| 12 | 1516 s | 5628 s | 75 s | 92 s, on 313 left | **1063 of 1063** |

The plain phase alone already gives the reduced walk's verdict on all but 5 and
14 placements; the reduced phase closes those, finding 3 and 7 certificates and
closing the rest. Putting the aliases first changes nothing but a few cached
answers. No promotion happened at either length: every placement's verdict was
settled by the time it was asked.

**Sharing does not survive being split.** Through `jobs.runTask` at `n = 12`:
one worker 170 s; four workers, units dealt one at a time, 147 s of wall clock
and 585 s of work; four workers dealt contiguous runs of the catalogue, 158 s
and 601 s. The orbits are shared *across* words -- F-051's 484-row orbit holds
`45`, `504` and `555` -- so every worker walks the big ones again whichever way
the units are dealt. A shared census is one worker per census, and the cores go
to running more censuses.

**At `n = 14`**, one worker, the whole `--max-word 4` catalogue (1850
placements after the mirror): **2745 s** of wall clock and 187 MB at most;
1135 inside, 715 outside, none undecided. The plain census of the same length
(2148 placements, both halves of each mirror pair) took 9.9 core-hours on the
laptop, which is about 5.8 on this machine, so the shared census is about
**7.6x** cheaper at 14 against 8.9x at 12 -- it grew sixteenfold from 12 to 14
where the plain one grew fourteenfold. The `45` slide is `ioooooii`, head 1 and
tail 2, as under the plain walk. Whether every verdict equals the reduced
walk's was not checked at this length: the reduced census alone would take the
better part of a day here.

Reproduce:

```bash
python batch.py cores 12 --max-word 4 --jobs 1            # the shared census
python batch.py cores 12 --max-word 4 --walk reduced --no-mirror --jobs 4
```

The comparison joins the two on each placement's class. `tests/test_reduced_walk.py`
pins it at `n = 11` for the `--max-word 3` catalogue.

---

## E-049 — A core and the same core with a two-arrow relation, walked both ways
*2026-09-22* · **the plain walk gave one derived class two verdicts 19 times at `n = 11` and 12; walked as one state, 60 placements move from outside to inside and none the other way; the census also asks every mirror pair twice** → F-052, H-020, F-051

**The question.** F-051 lists `45` at offset 1 and `245` at offset 0 as two of
"four different cores" in one orbit. They are one LNA up to a relation of two
arrows, so one derived class by the free move (F-028). Does the census treat
them as one, and what does asking both cost?

**How often it happens.** Under `--max-word 4 --gaps 1,2,3` a word holding a
`2` is its stripped word at another offset, and that stripped row is always
another placement of the same catalogue:

| n | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|
| placements | 859 | 1262 | 1705 | 2148 | 2591 | 3034 | 3477 |
| holding a two-arrow relation | 189 | 249 | 309 | 369 | 429 | 489 | 549 |
| stripped row also in the census | 189 | 249 | 309 | 369 | 429 | 489 | 549 |

Five of E-045's six undecideds -- `245`, `2045`, `2245`, `2555`, `2556` -- are
such words: `45` at one or two offsets further on, `555` and `556` at one.

**The two did not get the same verdict.** Every alias and its stripped
placement were run through `_verdictFor` at `n = 11` and 12 with the default
limits. **6 of 189 and 13 of 249 disagree**, always the same way round: the
alias `inside`, the stripped row `outside`. The smallest: `404` at offset 1 of
`n = 11` has **no move at all** (a closed orbit of one row), while `2404` at
offset 0 reaches an almost separate row in nine. The cause is that
`freeMoves.movesFrom(free = True)` only deletes two-arrow relations. The free
move is symmetric and the walk made it one-way, and the added relation is
exactly the spectator a rule needs (F-023).

**Walked as one state.** `freeMoves.REDUCED`: every state is a reduced row, and a
step is a move out of it or out of it with one two-arrow relation added,
stripped again. Every reduced placement at `n = 11` and 12 was run under both
walks and compared with the best verdict any of its aliases had under the plain
walk:

| n | placements, plain | reduced | reduced worse than the best alias | outside → inside | core-seconds, plain | reduced |
|---|---|---|---|---|---|---|
| 11 | 859 | 670 | **0** | 20 | 464 | 950 |
| 12 | 1262 | 1013 | **0** | 40 | 1516 | 5628 |

So the reduced walk loses nothing any alias had and places 60 placements the
plain walk called outside. The `45` and `504` slides are unchanged at both
lengths (`iooii`, `ioooii`; `iiooi`, `iioooi`), so F-042 stands. The ones that
move include `404` (`ioooi` → `iiiii` at 11, `iooooi` → `iiiiii` at 12),
`405`, `5004`, `36`, `6006` (`ooo` → `iii` at 12), and `4056` (`oio` →
`oii` at 12), the one shape H-020 had to amend its statement for.

**It costs more, and that is not a saving.** A reduced step tries every vertex
a two-arrow relation can be added at, so a row costs about twice as much, and
the outside orbits it has to close are the dearest part: 307 outside placements
took 4145 core-seconds at `n = 12`, against 1311 for 389 under the plain walk.
Walking plain first and reduced only where plain did not find a certificate
barely helps (5454 against 5628), because it is the outside orbits and not the
inside answers that cost. The factor went from 2.0 at 11 to 3.7 at 12;
nothing longer has been measured.

**The mirror is a second duplication.** The relation dual sends a row to a row
of the same class and the move set is closed under it (F-026), so the verdicts
of a row and its mirror should agree. Over the same runs: **300 of 300** mirror
pairs agree at `n = 11` under the plain walk, 288 of 288 under the reduced one,
398 of 398 and 384 of 384 at `n = 12`. Asking one of each pair removes 12 to 17
percent of a census (3034 → 2637 placements at `n = 16` under the plain walk,
2545 → 2159 under the reduced one).

**What changed in the code.** `batch.py cores` and `batch.py sample` take
`--walk plain|reduced`; the default is `plain`, so every command in
`OVERNIGHT.md` means what it meant, and `reduced` writes a ledger whose name
ends `-reduced`. Under `reduced` the catalogue drops every word holding a `2`.
Under both, the census asks only the first of each mirror pair in catalogue
order and `--summary` reads the other half of each slide in the mirror;
`--no-mirror` turns that off. The mirror is a filter, not an instrument, so
the ledger is the same one.

Reproduce: the comparison is `_verdictFor(n, word, offset, 20000, 6000, free =
...)` over `batch._placements` with `--walk plain` and `--walk reduced`, which
at `n = 12` is about 2 core-hours. `tests/test_reduced_walk.py` pins the `404`
case, the catalogue and the mirror.

---

## E-048 — The leftover rate at `n = 13` and `n = 17`, by one instrument
*2026-09-22* · **`n = 13` is 39.9% +- 0.7 with no cap in it; `n = 17` is somewhere between 47% and 70%** → H-019

Two depth-0 samples from `OVERNIGHT.md`'s first "tonight" line, run on the night
of 2026-09-20 alongside the censuses of E-046 and E-047. Both stopped on the
9-hour budget, which for a sample is a prefix and still uniform.

| n | drawn | core-h | s/draw | theorem | moves | **leftover, orbit closed** | leftover, orbit capped |
|---|---|---|---|---|---|---|---|
| 13 | 4575 | 18.0 | 14.1 | 13.5% | 46.6% | **39.9% +- 0.7** | **0** |
| 15 *(E-045, depth 4, old rows)* | 879 | 62.7 | 257 | 7.7% | 34.2% | **48.8% +- 1.7** | 9.2% |
| 17 | 214 | 9.0 | 152 | 3.7% | 26.6% | **47.2% +- 3.4** | 22.4% |

**`n = 13` is the first clean point above `n = 11`.** Every one of its 1825
leftovers had a forward orbit that *closed*, so the rate carries no cap slack at
all -- it is exactly "the moves do not reach an almost separate row", the same
statement as F-032's exhaustive 0.6, 5.4 and 16. The series by one definition is
now 0.7, 5.4, 16, 28 (60 draws), **39.9**.

**Above that the cap decides the answer.** The closed-orbit fraction, which is a
lower bound on the leftover rate, reads 48.8% at `n = 15` and 47.2% +- 3.4 at
`n = 17`: flat. The capped fraction, which could go either way, goes from 9% to
22%. So "still rising past 15" and "levelling off near a half" are both
consistent with what is on disk, and only a larger `--orbit-limit` separates
them. The `n = 15` row is E-045's, whose rows predate the closed/capped split;
its 81 capped leftovers were counted there from the orbit sizes.

**What it costs, which the menu had wrong.** `OVERNIGHT.md` said a depth-0 draw
is about a millisecond below `n = 13`. At `n = 13` it is **14 s**, and 17.5 of
the 18 core-hours went on leftovers -- closing an orbit is the whole cost, and
at `n = 13` that is most draws. At `n = 17` a draw is 152 s, and the 48 capped
ones took 4.7 of the 9 core-hours between them. One worker at `n = 17` is 214
draws a night; the `--count 20000` on the menu was two orders of magnitude past
what a night does, which is harmless (the budget cuts it) but should not be read
as a plan.

**By overlap, at `n = 13`**, the leftovers are 254 at overlap 2, 696 at 3, 509
at 4, 251 at 5 and 115 above. Overlap 2 -- F-042's cell `(2, 1)` and its
neighbours -- is one leftover in seven, as it should be if no bound on the overlap
cuts the leftovers off.

**A correction to E-044 found on the way.** Its table gives `n = 12` 208012
LNAs. That is Catalan(12), the count at `n = 13` (this run's summary prints it
for `n = 13`); `n = 12` has Catalan(11) = 58786. The rate there is unaffected,
since the draw is by index and the count is only printed.

Reproduce:

```bash
python batch.py sample 13 --count 50000 --summary
python batch.py sample 17 --count 20000 --summary
```

---

## E-047 — Two heavy clusters in one line, at lengths with room for them
*2026-09-22* · **outsiders with a free arrow between the clusters exist from `n = 14`; every one has a half that is outside on its own** → H-018, F-040

F-040 found no LNA outside a quipu class with two heavy clusters and a free arrow
between them, at every length up to 11, and noted the shape barely fits there.
Three runs from `OVERNIGHT.md` asked the same of the move closure at lengths
with room:

| run | placements | core-h | inside | outside | undecided | free-gap rows |
|---|---|---|---|---|---|---|
| `cores 15 --max-word 2 --pair-word 2 --gaps 1,2,3,4` | 1548 | 18.1 | 776 | 751 | 21 | 130 |
| `cores 17 --max-word 2 --pair-word 2 --gaps 1,2,3,4` | 2256 | 41.9 | 1017 | 1029 | 210 | 210 |
| `cores 16 --max-word 2 --pair-word 3 --gaps 1,2,3 --core-limit 900` | 859 | 17.9 | 391 | 401 | 67 | **0** |

The single-core censuses of E-046 hold pair words too (their default `--gaps
1,2,3`), and contribute 4, 10, 30 and 50 free-gap rows at `n = 11, 12, 14, 16`.

**The count F-040 made, at lengths it could not reach.** Over all 434 rows with
two heavy clusters and a free arrow between them: **388 inside, 32 outside, 14
undecided.** None at `n = 11` or 12 is outside. The first is at `n = 14`:
`330004500000`, a `33` against the source and a `45` five arrows in, whose
forward orbit closes at 1895 rows without an almost separate member.

**Every outsider is explained by one of its halves.** Looking each half up alone,
at the same offset in the same length:

| left half alone | right half alone | the pair | rows |
|---|---|---|---|
| inside | inside | inside | **364** |
| inside | inside | outside | **0** |
| inside | outside | outside | 32 |
| inside | outside | **inside** | **24** |
| inside | outside | undecided | 14 |

So with a free arrow between them, two placeable clusters make a placeable pair
every time (364 of 364), and each of the 32 outside pairs has a half that is
outside where it sits. The `45` in `330004500000` sits in its own outside band
(`45` at `n = 14` is `ioooooii`; it is at offset 5). What looked like the case
that would break the single-cluster reading is a single-cluster outsider with a
harmless neighbour.

**But the halves are not independent, and the 24 say so.** A `33`, `34` or `44`
at or near the source, three or four arrows before a `35` or `36`, carries the
pair inside where the `35` or `36` alone is outside: `3300035` at offsets 0 to 2
at `n = 17`, `33000035` at 0 and 1, and so on through `34` and `44`, at `n = 15,
16, 17`. The rescue reaches less far as the gap widens -- three offsets at gap 3,
two at gap 4 -- and every right half involved is one whose own slide has an
inside head (`35`: `ii…`, `36`: `i…`). One reading is that the left cluster,
pushed off the source, lets the right one reach its own head; nothing here
checks that.

**Without a free arrow the interaction goes the other way too.** Among the rows
where the two clusters touch through a one-arrow overlap, eight have both halves
inside alone and the pair outside -- all of them `35 0^g xx` at offset 1
(`3500033` at `n = 11`, `3500034` and `3500044` at 12, `3500036`, `3500046`,
`3500056`, `3500066` at 14, `3500035` at 15), where the `5` reaches over the gap
into the second cluster.

**The run that answered nothing, and why.** `--pair-word 3 --core-limit 900` at
`n = 16` has no free-gap row at all. The catalogue is sorted by word length, so
its first 900 words are the gap-1 pairs -- and **a gap of one or two never has a
free arrow** when every relation is two arrows or more: the last relation of the
first half reaches over it. Measured on the catalogue at `n = 16`:

| pair-word 2, gap | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| placements | 406 | 496 | 500 | 400 | 300 | 200 |
| with a free arrow | 0 | 0 | 50 | 120 | 180 | **200** |

and for `--pair-word 3`, 0 of 6135 at gap 1, 170 of 7520 at gap 2, and all 1447
at gap 6. The night's 859 placements at `--pair-word 3` are real rows about
touching clusters and are kept, but the free-gap question needs `--gaps 5,6`.

Reproduce:

```bash
python batch.py cores 15 --max-word 2 --pair-word 2 --gaps 1,2,3,4 --summary
python batch.py cores 17 --max-word 2 --pair-word 2 --gaps 1,2,3,4 --summary
python batch.py cores 16 --max-word 2 --pair-word 3 --gaps 1,2,3 --summary
```

The half-by-half table is not something `--summary` prints; it came from joining
each pair row to the single-core rows of the same length on `(core, offset)`.

---

## E-046 — The core census at lengths 11 to 17, and the six undecideds resolved
*2026-09-22* · **H-020 holds for every single cluster from `n = 13` to 17; the undecideds were all outside; the census walks the same orbits over and over** → H-020, F-051

Four `overnight.py` invocations between 2026-09-20 19:15 and 2026-09-21 23:17,
all from `OVERNIGHT.md`: the first "tonight" line, the "run of lengths", the
"two-cluster shapes" (E-047) and the "six undecideds, sharpened". The
`--max-word 4` ledgers now stand at:

| n | placements | cores | inside | outside | undecided | complete | core-h |
|---|---|---|---|---|---|---|---|
| 11 | 859 | 337 | 675 | 184 | 0 | yes | 0.3 |
| 12 | 1262 | 403 | 873 | 389 | 0 | yes | 0.7 |
| 13 | 367 | 62 | 287 | 80 | 0 | no (E-045) | -- |
| 14 | 2148 | 443 | 1252 | 896 | 0 | **yes** | 9.9 |
| 15 | 932 | 125 | 666 | 260 | 6 | no (E-045) | -- |
| 16 | 3034 | 443 | 1599 | 1278 | 157 | **yes** | 34.6 |
| 17 | 856 | 89 | 636 | 196 | 24 | `--core-limit 250` | 11.3 |

**H-020, asked properly.** Every core with one heavy cluster that is decided at
two lengths from 13 to 17 was compared between them: **975 comparisons over 283
cores, no failure** of "the longer slide is the shorter one with the interior
verdict repeated". The `45` core, F-042's example, reads `iooii`, `ioooii`, …,
`iooooooooii` from `n = 11` to 17 -- head 1, tail 2 at every length -- and its
opposite `504` is head 2, tail 1 at every length. At `n = 14` the 343
single-cluster cores split 118 inside everywhere, 97 outside everywhere, 127
inside at the ends and outside between, and **one** with a shape the head/tail
reading does not cover: `4056`, `oio` / `oooio` / `oooooio` at 12, 14, 16, which
is inside one offset from the sink and outside at the sink itself. Its suffix
`io` is as fixed as any tail; it is a word, not a count. No single-cluster slide
is inside in its interior and outside at an end.

**Below `n = 13` the law is not visible yet, as F-042 predicted.** Twelve
single-cluster cores change shape between 11 and 12 or 12 and 14 -- `44066`
`io` → `oooo`, `6600044` `oi` → `oooo`, `550066` `i` → `ooo` -- and every one of
them is a long word whose slide at the shorter length is three offsets or
fewer, where a placement is near both ends at once.

**Two-cluster cores do not obey it, and should not.** 31 incompatibilities
between lengths, all from pair words: `3500035` is `iio` at 14 and `ooooo` at 16,
because its right `35` moves away from the sink as the line grows while its left
one stays at the source. E-047 is what they obey instead.

**The six undecideds of E-045, at three times the orbit limit.** `--orbit-limit
60000 --join-limit 40000 --cores 245,2045,2245,2555,2556,3344` at `n = 15`, 43
placements, 2.3 core-hours: **all six are outside**, each a closed orbit of
**21709 rows** -- just past the 20000 the census had. The heads and tails E-045
read off are unchanged. And across the whole census, the 92 slides that hold an
undecided are every one consistent with the other lengths when `?` is read as
`o`, and two of them also when read as `i`; none only as `i`. The undecideds
are the largest outside orbits, which sit next to the inside end (F-051), and
they are almost certainly outside.

**Where the time went -- and it is the opposite of E-045.** The inside answers
are now nearly free: 1.8 of `n = 16`'s 34.6 core-hours. Outside took 12.5, and
**the 157 undecided took 20.3**, 7.8 minutes each, for placements every other
line of evidence says are outside. The cost table in `OVERNIGHT.md` estimated
8.7 core-hours for `n = 16`; the real figure is four times that, and it is the
cap that the estimate missed.

**The census walks the same orbits again and again.** Placements with the same
closed-orbit size are frequent -- at `n = 14`, the 896 outside placements have
**107** distinct orbit sizes between them, and one size, 7393, was walked 72
times for 3.7 core-hours. Walking a handful of them again and comparing the sets
(F-051): same size is same orbit, or a mirror pair of orbits. If each distinct
orbit were walked once, the outside rows would cost 0.5 core-hours instead of
9.2 at `n = 14`, 1.0 instead of 12.5 at `n = 16`, and 2.1 instead of 13.9 for
the `n = 17` pairs.

**A catalogue word is not a core.** `--core-limit 250` at `n = 17` gave 89
cores, not 250: the catalogue lists every word the digit ranges allow, and most
four-letter words are not LNAs at all (their relations' ends do not increase)
and have no placement at any length. The cut is still the same at every length,
which is what it is for, but it is a cut of the catalogue, not a count of cores.

Reproduce:

```bash
python batch.py cores 14 --max-word 4 --summary
python batch.py cores 16 --max-word 4 --summary
python batch.py cores 15 --max-word 4 --orbit-limit 60000 --join-limit 40000 --summary
```

The cross-length comparison is not in `--summary`: it reads every
`cores-n*-w4p2a6g123-o20000j6000.jsonl` ledger, builds each core's slide per
length, and checks that each pair of decided slides differs by a run of the
interior verdict.

---

## E-045 — The length-15 night: two censuses and a sample, all three cut off
*2026-09-20* · **no census finished; the instrument was the bottleneck, not the machine** → H-018, H-019, H-020

The night `overnight.py --hours 9 --only sample15 cores15 cores13` was left to
run. All three jobs used their whole budget and all three stopped on it, with
the work they had done on disk:

| job | units done | of | core-hours | what it was for |
|---|---|---|---|---|
| `cores15` | 932 | 2591 | 53.8 | H-018, the census at a length with room |
| `cores13` | 367 | 1705 | 18.0 | the control length F-042 already covered |
| `sample15` | 879 | 6000 | 62.7 | H-019, the leftover rate at `n = 15` |

The machine was not the problem: 134 core-hours came back from 15 workers in 9
hours, which is a saturated machine. **The instrument was.** Three things came
out of reading the ledgers, and the first is much the largest.

**1. The walk was paying for the whole orbit to answer a membership question.**
Of `cores15`'s 932 placements, the 666 that came back `inside` cost 49 of the
53.8 core-hours; the 260 `outside` ones cost 3.8. That is backwards -- `inside`
is the *easy* answer -- and the reason was that `_verdictFor` called
`freeMoves.orbitOf`, which enumerates to its 20000-row cap, and only then looked
through the result for an almost separate row. Re-running a random twelve of
those `inside` placements with the walk stopping at the first almost separate
row it meets: every one of them found its certificate within 351 rows, nine of
the twelve within 100, and the twelve together took **5.8 seconds against the
2818 seconds they cost on the night** -- 486x. `sampling.probe` had the same
shape and the same 15.5 core-hours of `moves` rows to show for it.

`freeMoves.orbitReport` is the fix, with a `stopWhen` predicate and, because a
walk that stops early must not be mistaken for one that closed, an explicit
`stoppedBy` of `closed`, `found` or `cap`. The verdicts are unchanged; what
changed is the price. Measured afterwards, single core, over a random sample of
each catalogue:

| census | placements | median | mean | one core |
|---|---|---|---|---|
| `cores 13 --max-word 4` | 1705 | 0.18s | 2.35s | **1.1 h** |
| `cores 15 --max-word 4` | 2591 | 0.28s | 5.27s | **3.8 h** |
| `cores 16 --max-word 4` | 3034 | 0.61s | 10.35s | **8.7 h** |
| `cores 17 --max-word 4` | 3477 | 2.66s | 38.82s | **37.5 h** |

`cores15` was four nights of work and is now half an hour on eight workers. The
catalogue grows slowly and the price of a placement does not: `n = 17` is still
a night, `n = 18` at this width is several.

**2. The two censuses walked their catalogues at different speeds, so only
their overlap could be read.** `cores13` reached 62 cores and `cores15` reached
125, both from the front of the same catalogue. A census read against another
length is the entire point of running one, and the budget, not the question,
decided where each stopped. `--core-limit` and `--cores` now cut the catalogue
deliberately instead, and neither is in the ledger's name, so two nights can
split one census between them.

**3. The shapes the length was chosen for got no units at all.** `--gaps` puts
two-cluster cores in the catalogue after every single core, and at the rate the
night ran neither length reached them: "rows with two heavy clusters and a free
arrow between them: 0" in both summaries. F-040's question -- the one that needs
`n >= 13` to be askable -- was not asked. It needs its own run, which
`--max-word 2 --pair-word 2` gives, and that is now in `OVERNIGHT.md`.

**What the partial data looks like, recorded and not concluded from.** Of the 62
cores done at both lengths, 25 have an `outside` somewhere in their slide. Every
one of those 25 has the same shape at both lengths: some number of inside
offsets at the source end, some number at the sink end, and outside everywhere
between, with **both counts identical at `n = 13` and `n = 15`** and the outside
middle absorbing the two extra offsets. Not one of the 25 has an inside in its
interior. H-020 is that observation written down as something to test; it rests
on two lengths, one of them a third finished, and is not a finding.

Two smaller things worth having on the record:

* **Six placements came back `undecided`, and five of the six sit exactly at the
  boundary** between the outside middle and the inside tail (`ooooo?ii`,
  `oooo?ii`, `ooooo?i`), the sixth at the head boundary (`ii?ooooo`). The hard
  cases are where the answer changes, which is where a sharper limit is worth
  spending on: `--cores 245,2045,2245,2555,2556,3344 --join-limit 40000`.
* **81 of the 510 sampled leftovers, one in six, had walks that hit the cap**
  rather than closing. `settledBy == 'leftover'` was covering both, so a
  leftover *rate* at `n = 15` is part rate and part budget. `probe` now records
  `orbitClosed` and the summary splits them, and `--orbit-limit` has joined the
  sampler's ledger name, which it should always have been in -- the same
  mistake `--pair-word` made in `cores` before the first night.

Reproduce, or continue:

```bash
python batch.py cores 15 --max-word 4 --summary
python batch.py cores 13 --max-word 4 --summary
python batch.py sample 15 --count 6000 --depth 4 --summary
```

The ledgers from the night are kept and the new runs resume from them. Their
rows carry no `stopped` or `orbitClosed` field, which is how a row written
before the walk learned to stop early can be told from one written after. The
sample ledger was renamed from `sample-n15-s0-d4.jsonl` to
`sample-n15-s0-d4-o20000.jsonl`, because `--orbit-limit` was missing from the
name and two runs under different limits disagree about what a leftover is.

`OVERNIGHT.md` is the menu this run produced: what each parameter changes, what
a census of each length costs, and the runs worth making next.

---

## E-044 — The sampler, calibrated against the lengths whose answer is known
*2026-09-19* · **agrees at `n = 9`, 10 and 11; `n = 12` is 28% leftover** → H-019

`batch.py sample` is only worth running at `n = 16` if it gives the right answer
at `n = 10`, where the answer is known exhaustively. This is that check. The
draw is uniform over all Catalan(n - 1) LNAs of the length; each is put through
the cheap pipeline one row at a time (`sampling.probe`), where F-032 measured the
same thing a length at a time.

| n | LNAs | drawn | theorem | moves | **leftover** | exhaustive leftover |
|---|---|---|---|---|---|---|
| 9 | 1430 | 60 | 40% | 58% | **1.7% +- 1.7** | **0.70%** (F-032: 0.6%) |
| 10 | 4862 | 400 | 31.2% | 64.0% | **4.8% +- 1.1** | **5.4%** (F-032) |
| 11 | 16796 | 6 | 33% | 50% | 17% +- 15 | 16% (F-032) |
| 12 | 208012 *(sic: 58786 -- E-048)* | 60 | 27% | 45% | **28% +- 5.8** | not known |

**`n = 10` is the load-bearing row**: 400 draws put the leftover rate at
4.8% +- 1.1, and the exhaustive answer is 5.4%. `n = 11` drew only 6 before the
run's time limit — the orbit walk is much slower there — so its agreement is
suggestive and no more.

**The exhaustive `n = 9` pass was run here too**, over all 1430 rows, as a check
that `probe` asked one row at a time agrees with `overlaps.py` asked a length at
a time: 610 by the theorem, 810 by the moves, **10 leftovers, 0.70%**. It does.

**`n = 12` is the new number.** 28% +- 5.8, against 0.7% at `n = 9` — the rate
is rising steeply and the sampler is the only way to see it above `n = 11`.
Preliminary: 60 draws, one seed. H-019 is the hypothesis this feeds and says
what would settle it.

**A leftover rate from a sample is an upper bound, not an estimate of the
truth.** `freeMoves.orbitOf` walks forwards only and stops at `--orbit-limit`,
so a row it fails to place may still be placeable — by a longer walk, or by a
walk from the other direction, which is the asymmetry E-037 recorded. The rates
above are therefore "not placed by this much walking", and the exhaustive
figures they are checked against carry the same caveat.

**The run that was cut off resumed correctly**, which was not the point of the
experiment but is worth recording: `n = 11` stopped at 6 of 400 draws on a time
limit, and the ledger holds those 6, so the same command continues from draw 7.

Reproduce:

```bash
python batch.py sample 10 --count 400 --jobs 4
python batch.py sample 10 --summary
```

---

## E-043 — The deduplicated walk, checked against the plain walk
*2026-09-19* · **same answers everywhere it was checked; 2.2x to 6.1x faster** → F-049, F-050

The dedup of E-042 is only worth having if it changes no answer. Checked by
running both walks over **every LNA of a length** and comparing three things at
once: the set of lines collected, the set of hereditary forms collected, and the
set of canonical keys the visitor was shown.

| length | depth | LNAs | mismatches | nodes plain | nodes deduped | ratio |
|---|---|---|---|---|---|---|
| 5 | 4 | 14 | **0** | 847 | 674 | 1.3x |
| 6 | 5 | 42 | **0** | 14701 | 7427 | 2.0x |
| 7 | 5 | 132 | **0** | 94498 | 39280 | 2.4x |
| 8 | 4 | 429 | **0** | 143068 | 78295 | 1.8x |

**It failed first, and the failure was real.** Before the sign gauge was
quotiented out, `n = 7`, `30300`, depth 5 lost one node: the walk reached one
algebra under two presentations differing only in the sign of a zero relation,
which are the same ideal, so the dedup was right to identify them — and the
*procedure* then gave two different answers from them. That is F-050, and it was
found by this check rather than by reading.

**Wall clock, on the searches the merge hunt runs** (`n = 9`, one container, one
process):

| start | depth | plain | deduped | speedup |
|---|---|---|---|---|
| `3345000` | 4 | 8.0 s | 3.6 s | 2.2x |
| `3345000` | 5 | 38.4 s | 10.7 s | 3.6x |
| `3345000` | 6 | 184.9 s | 30.5 s | **6.1x** |
| `3033030` | 4 | 4.9 s | 2.4 s | 2.1x |
| `3033030` | 5 | 18.4 s | 5.9 s | 3.1x |
| `3033030` | 6 | 77.1 s | 15.2 s | **5.1x** |

The speedup is below the node ratio because the key costs something at every
node; it grows with depth for the same reason the node ratio does.

**What the dedup does not preserve is the number of *visits*.** At `n = 9`,
`3345000`, depth 6 the plain walk collects 54 line records and the deduped walk
13 — and the **set** of LNAs reached is identical, 4 either way, at depths 4 and
5 as well. Every consumer in the repo builds a set, so this is invisible to all
of them; a caller that counted records would be counting routes, which was never
a meaningful number.

**And end to end, through `merges.py` itself.** The same run at `n = 9`, depth 4,
`--all-groups`, once with the dedup and once with `--no-dedupe`: all 9 members
searched, **every `lnasReached`, `reachedOrbits`, `coveredOrbits` and `alarms`
identical**, 121.5 s against 215.1 s of search time — 1.8x. That is the shallow
end; the depths the overnight job runs at are 5 to 8, where the single-search
measurement above is 3.6x to 6.1x.

`merges.py` now passes a `fingerprint.Visited` by default and records what it
skipped in the checkpoint's `walk` field, so a run's cost can be read back
rather than re-measured. `--no-dedupe` restores the old walk, for measuring what
the dedup changes rather than for producing answers.

Reproduce: `tests/test_fingerprint.py::test_dedup_reaches_exactly_what_the_plain_walk_reaches`
pins the n = 5 to 7 cases at the depths above.

```bash
python merges.py 9 --depths 4 --all-groups --checkpoint logs/a.jsonl
python merges.py 9 --depths 4 --all-groups --no-dedupe --checkpoint logs/b.jsonl
```

---

## E-042 — How much of a mutation search is repeated work
*2026-09-19* · **3.6x at depth 4, 6.3x at 5, 11.4x at 6, and compounding** → F-049

`mutationSearchDepthFirst` walks mutation *sequences* and has never deduplicated
the algebras they reach. This counts what that costs, over the two leftover
orbits at `n = 9` that `merges.py` searches from, by keying every node the
visitor is shown with `search.quiverKey`.

| start | depth | nodes | distinct | ratio | parallel-arrow nodes |
|---|---|---|---|---|---|
| `3345000` | 4 | 779 | 216 | 3.6x | 0 |
| `3345000` | 5 | 3887 | 615 | 6.3x | 0 |
| `3345000` | 6 | 19483 | 1708 | **11.4x** | 40 (0.2%) |
| `3033030` | 4 | 350 | 108 | 3.2x | 4 (1.1%) |
| `3033030` | 5 | 1435 | 262 | 5.5x | 44 (3.1%) |
| `3033030` | 6 | 6088 | 639 | **9.5x** | 284 (4.7%) |

**The ratio roughly doubles per level**, which is the number that matters: the
waste is not a constant overhead to be shrugged at but the dominant term at the
depths the merge hunt wants. Extrapolating the two columns, depth 8 is 30x to
40x, and depth 8 at `n = 10` is exactly what `overnight.py` is running.

**Why there are so many routes.** A mutation is invertible and mutations at
distant vertices commute, so the number of sequences reaching a given algebra
grows with the depth faster than the number of algebras does. Nothing about
this is specific to these two starts.

**The parallel-arrow census, which decides whether a canonical form is hard.**
Over the two depth-5 walks, 44 of 5322 nodes have parallel arrows, and **every
one of them has exactly one bundle of exactly two arrows**. So the arrow-key
ambiguity that NOTES.md has recorded as an open want since F-039 costs a
minimisation over 2 relabelings, not a search. There was no node with two
bundles, and none with a bundle of three.

**Vertex labels do not move under mutation**, so none of this is graph
isomorphism testing: two algebras reached from one start are equal on the nose
or not at all. That is why the key is exact and cheap, and why the probabilistic
fingerprint this measurement was meant to justify turned out not to be needed
for identification at all — only, in `fingerprint.digest`, for storing a very
large visited set in less memory.

Reproduce:

```python
from quivermutation import nakayama as nk, search, fingerprint
alg = nk.LinearNakayamaAlgebra(9, [3, 5, 0, 5, 0, 0, 0])
visited = fingerprint.Visited()
search.mutationSearchDepthFirst(alg, 6, [], 'x', printOutput = False)   # plain
print(visited.summarise())
```

---

## E-041 — The two separators of the sweep, checked
*2026-09-19* · **both hold; one reproduces the quipu boundary exactly, the other certifies 619 rows at `n = 11`** → F-047, F-048

E-040 ended by naming two criteria as its highest-value unverified output, and
warning that a criterion quoted out of a paper should be assumed to be missing a
hypothesis until it has been run against a known answer. This is that run.

### (a) The `Z`-congruence invariant (math/0610685 Cor. 3.15) → F-047

Profile: the Smith normal form of `g(Φ)` for each irreducible factor `g` of the
Coxeter polynomial.

* **Soundness, on every member of every orbit** rather than a sample: 1430 LNAs in
  48 orbits at `n = 9`, 4862 in 113 at `n = 10`. **Zero orbits split.**
* **Sharpness:** splits 1 of the 11 cospectral orbit-groups at `n = 9`, 3 of 25 at
  `n = 10`.
* **The `n = 9` split was then named by the quipu theorem**, which is the check
  that turns a refinement into a separation: all four orbits on one profile are
  `P^(1,4)_(1,0,1)`, both on the other are `P^(1,2)_(1,1,2)`, no orbit on the
  wrong side. The invariant reproduces the true class boundary on the exact pair
  the Coxeter polynomial cannot see.

Cost: about three minutes for `n = 9`, fifteen for `n = 10`, in sympy. The Smith
normal form is the expensive part and it is per irreducible factor, so this scales
with the factorisation rather than with `n`.

### (b) The periodicity criterion (math/0611201 Thm. 3.4) → F-048

`Φ` periodic and the Euler form indefinite certifies not piecewise hereditary.

* Fires on **0** LNAs at `n = 5` to `9`, **3** at `n = 10`, **638** at `n = 11`.
* At `n = 10` the three are `34504030`, `50505000`, `45050400` and
  `piecewiseHereditary` certifies **none** of them — the rows backlog 20 wants
  `lemma:taupathimpliesnotpwh` for. At `n = 11`, 619 of the 638 are new.
* **Falsification test:** an almost separate LNA is piecewise hereditary by the
  quipu theorem, so the criterion must never fire on one. Over 30648 rows at
  `n = 5` to `11` it fires on 641, **not one almost separate**.

This is the criterion E-040's warning was about, in its usable form. The version
that misfired there was de la Peña's, which needs "not of Dynkin module type"
supplied separately; Ladkani's asks for an indefinite Euler form instead, and
indefiniteness rules the Dynkin and Euclidean cases out by itself. **The same
mathematics, and one statement of it is safe to implement while the other is
not** — which is the concrete lesson, rather than the general caution.

### What is still not done

Both are verified as *criteria*; neither is in the code. F-047's profile has been
checked at two lengths and F-048's at seven, so porting them is now a matter of
writing them into `invariants` and `piecewiseHereditary` rather than of deciding
whether they are true. Backlog 32 is updated.

---

## E-040 — A literature sweep aimed at merging classes of LNAs
*2026-09-19* · **21 summaries; one criterion that classifies the tame half outright, and three merges at `n = 11` no move of ours makes** → F-043, F-044, F-045, F-046

The question put to the literature was narrow on purpose: **what merges two
LNAs?** Not what classifies Nakayama algebras, not what invariants exist — what
would let two of our orbits be joined, or be proved distinct. Anything that could
not be tied to that in a sentence was rejected.

### How it was searched

Three passes, and the third is the one that paid.

1. **Backwards**, through the reference lists of the three papers the project
   already uses (arXiv:2112.08129, 2305.06642, 2310.08346). Seventeen distinct
   references between them — a small enough set to read in full. This is what the
   `literature/README.md` candidate list was built from.
2. **Outwards**, from those into the authors' own corpora: Ladkani's 27 arXiv
   papers, the Happel school, the silting-mutation line.
3. **Forwards**, by citation. Semantic Scholar's graph API on the three papers'
   arXiv ids, plus `export.arxiv.org` title and abstract search on "Nakayama" ×
   "derived equivalence". **This found the two best papers in the sweep, and no
   backward reference list could have**: arXiv:2302.02880 (Ueda) and
   arXiv:2203.15735 (Dong–Lin–Ruan) are both later than everything we cite, and
   Brüstle came in as a reference of arXiv:1910.01494, which itself was only found
   forwards. *Do the forward pass first next time.*

Roughly 60 papers screened on abstracts, 21 read closely enough for a file, 22 of
Ladkani's rejected with a recorded one-line reason each so they are not re-screened.

### What came back, sorted by what it does

**Merges.** `research/literature/` now holds four sources that produce derived
equivalences between LNAs: Ueda's Cor. 1.3 (F-043), the Happel–Seidel symmetry and
its extension in Lenzing–Meltzer–Ruan Prop. 4.1, Dong–Lin–Ruan Prop. 4.5, and
Ladkani's `A(mn, m+1) ≃ kA_m ⊗ kA_n` (0911.5137 Cor. 1.2). Every one of them is
about **radical powers** `kA_n/rad^r` or a tensor of two lines — the thinnest
family of rows we have. Nothing found speaks about an LNA with relations of
several different lengths, which is the open case.

**A classification.** Brüstle's Theorem 1.2, which decides the derived class of
any LNA with non-negative Euler form from three numbers: F-045.

**Separators.** Ladkani's Cartan-matrix-up-to-`Z`-congruence (math/0610685
Cor. 3.13), reported to split cospectral quipu groups the Coxeter polynomial
cannot; and two periodicity criteria (math/0611201 Thm. 3.4; de la Peña,
arXiv:1310.1557) certifying non-piecewise-heredity, one of which is reported to
certify `34504030`, `50505000` and `45050400` at `n = 10`, which our own criteria
miss. **Neither has been re-verified here** — they are the obvious next thing to
check.

**Two lines closed.** `HH*(A) = k` for every LNA (arXiv:2312.14699), so idea 22
is dead; and no extension of the Avella-Alaminos–Geiß invariant to string algebras
exists, so R-008's line stays closed — arXiv:1910.01494's skewed-gentle conclusion
does not reach us, because its hypothesis is "no simple projective module" and
every LNA has one (`P_n = S_n`, from the sink).

### What was checked against the code, and what it cost

Everything below was run here rather than taken on trust, which is the only reason
the findings above are findings.

| check | result |
|---|---|
| Ueda Cor. 1.3, 15 instances at `n ≤ 16`, against `movesJoin` | all 15 joined — but by the **double mutation**, not the rule table (F-043) |
| Happel–Seidel Table 1, 11 star types, against `treeCoxeterKey` | 11/11, and an invented 12th row correctly fails (F-044) |
| Happel–Seidel Table 1, 12 sheaf types, against `canonicalWeightType` | 12/12 (F-044) |
| Brüstle Thm. 1.2, reimplemented, at `n = 9, 10, 11` | reproduces F-011's tame partition; 1 new merge per length; 0 contradictions (F-045) |
| LMR Prop. 4.1, 21 instances, against `movesJoin` | 15 joined, 6 not — all six with both orbits **exhausted**, not capped, at `n = 11, 13, 15` (F-046) |
| the three `n = 11` merges, against `derivedOrbits(11)` | four of our orbits merge into two (F-046) |
| de la Peña's periodicity criterion, naive reading, at `n = 9` | **certifies 273 LNAs including `A_9` itself** — the Dynkin exclusion is the whole criterion, caveat recorded in the file |

That last row is the one to remember: a criterion quoted out of a proof, applied
without its exclusions, certified the hereditary line as non-piecewise-hereditary.
Every summary in this sweep that states an implementable criterion should be
assumed to be missing a hypothesis until it has been run against something whose
answer we already know.

### What was not done

* The two separators above are unverified here (Ladkani's congruence invariant and
  the two periodicity criteria). They are the highest-value follow-up, because a
  separator is what F-010 has wanted since R-008.
* `n = 12` and up for Brüstle: the Smith normal form is cheap but `derivedOrbits`
  is not, so there is nothing to compare against past `n = 11`.
* Three papers are summarised **from secondary sources** — Happel–Seidel, Rickard,
  Assem–Happel — because they are journal-only and pre-arXiv. Happel–Seidel's has
  been checked (F-044); the other two have not.
* The cellular-automaton sweep of `literature/README.md` (H-009) was not touched;
  this sweep was about merging, not about the move rules as a rewriting system.

---

## E-039 — Ueda's radical-power equivalence, checked against the move table
*2026-09-19* · **15 instances at `n <= 16`, all 15 joined by the moves alone** → F-043

Prompted by the literature sweep: arXiv:2302.02880 Corollary 1.3 gives a triangle
equivalence `per N(n, l+1) -> per N(n, l)` for `n = p(p+1)q + p(p-1)r`,
`l = (p+1)q + pr`. Enumerating `p` in 2..5, `q` in 1..4, `r` in 0..4 plus the
half-integer case `p = 2`, and keeping `4 <= n <= 16`, gives 15 distinct `(n, l)`.

For each: build both LNAs, compare `coxeterKey`, and ask `freeMoves.movesJoin`
whether the moves join them, with the meeting row then checked reachable from
**both** ends by `orbitOf(..., target = meet)` -- a `movesJoin` result alone is
one walk from each side and worth confirming when the conclusion is that a paper
adds nothing.

* Coxeter keys agree in all 15, as they must.
* All 15 joined with `free = True` (derived equivalence, the same relation Ueda
  proves).
* All 15 joined again with **`free = False`**, so they are joined by *mutation*
  moves alone -- a stronger statement than the paper's, for these instances.
* Both-ends check passed in all 15.

Cost: under a minute for the whole sweep, `limit = 20000` rows per walk, never
approached.

**What was not done.** `n > 16` was not tried: the parameter grid thins out fast
(the next instances are at `n = 18` and `n = 20`) and the orbits grow, and the
point was made. The *other* corollary of the paper -- an equivalence from every
`N(n,l)` to an algebra of global dimension at most 2 -- was not checked, because
the target is not an LNA and the pipeline has nothing to compare it against.

Reproduce: `ueda2.py` in the merge session's scratchpad.

---

## E-038 — The branch's free-move table, re-measured under the guarded search
*2026-09-19* · **both rows reproduce exactly, and the merge changes no number** → E-037, F-041

E-037 was measured on `claude/quiver-nakayama-investigation` before F-038 and
F-039 landed on `main`, so it was measured with a search that took steps which
are not derived equivalences and with a Cartan matrix that counted a parallel
pair of arrows as one path. Merging the branch puts its measurements on top of
three changes at once:

* `search.mutationSearchDepthFirst` now refuses a step whose Coxeter key moves
  (`coxeterGuard`, on by default), so the walk is strictly shorter-reaching;
* `invariants.integerCartanMatrix` counts **arrow** paths and takes the rank of
  the ideal where the cheap count is not provably right, so the guard is reading
  a different number than it would have;
* `search.quiverKey` returns `None` at a quiver with parallel arrows, so those
  nodes are no longer meeting points (below).

Any of the three could have moved the table. Re-run of E-037 section 2, same
depth, `alsoDual = True`:

| | `n = 8` | `n = 9` |
|---|---|---|
| single deletions | 572 | 2002 |
| joined by the known moves | 562 | 1937 |
| joined by meeting at depth 3 | 10 | 52 |
| left | **0** | **13** |

**Identical to E-037 in every cell**, and the 13 left at `n = 9` are the same 13
pairs: `2302230/2300230`, `2302300/2300300`, `2302330/2300330`,
`2302302/2300302`, `2302030/2300030`, `2303020/2303000`, `3302302/3300302`,
`3022302/3020302`, `3002302/3000302`, `0230300/0030300`, `0230302/0030302`,
`0302302/0300302`, `0303020/0303000`. So F-041 stands as written: the guard
removed nothing the branch had counted, which is the same answer F-039 got for
the sweeps it re-ran — at these lengths the unguarded search was not reaching
anything the guarded one cannot.

`n = 10` was **not** re-run: 36s at `n = 8` and 296s at `n = 9` extrapolates past
what this was worth, and the two rows that were re-run are the two the finding's
claim rests on. The `n = 10` row of E-037 is therefore still a pre-guard
measurement and is marked as such there.

**One soundness fix the merge did need.** `quiverKey` keyed a node on
`pathAlg.rels`, and F-039 established that `rels` is a *lossy projection* once a
quiver has parallel arrows -- a path named by its vertices no longer says which
of two arrows it runs along, so two different algebras write down the same
`rels`. Before the merge this could not bite, because `procedure.isMutable`
refused to mutate a quiver with a parallel pair at all and the search never
descended past one. R-013 removed that refusal. Two searches could then have
"met" at a key that is not an algebra, which would be a fabricated proof of
derived equivalence -- the precise failure F-038 is about, reintroduced by a
different door. `arrowRels` cannot serve as the key either: the procedure hands
out arrow keys in the order it builds them, so the same algebra reached by two
routes carries different ones and there is no canonical form to compare across
routes (F-039 says this in as many words). So a parallel-arrow node is now not a
meeting point in either direction; the search still walks through it. No meeting
recorded at `n = 8` or `n = 9` was at such a node, which is why the table above
did not move.

Reproduce: the script is `reverify.py` in the merge session's scratchpad, and
the two rows it prints are `tests/test_other_families.py`'s
`test_every_single_two_arrow_deletion_at_length_eight_is_a_mutation` generalised
to a second length. Full suite after the merge: **3015 passed, 1 xfailed**.

---

## E-037 — What the short quivers can show, and two questions asked properly
*2026-09-17* · **the evidence base is narrower than it looks; the free move settles at `n = 8`; the barricade does not hold and bounded overlap is not the pattern** → F-040, F-041, F-042, H-017

*Renumbered at merge from `E-032`, which was taken on `main` first by an unrelated entry while this branch was open. Session logs and commit messages from the branch use the old identifier.*

Prompted from outside the code, and the prompt was the useful part: everything
known about the classes the quipu theorem misses comes from lengths 9 to 11,
where an LNA has very little room -- at `n = 11` no vertex is more than five from
an end -- so patterns found there may be patterns of the small cases rather than
of the problem.

### 1. How narrow the evidence is (F-040)

Counting heavy clusters -- runs of relations linked by overlaps of two arrows or
more, which is what puts an LNA outside the theorem:

| | `n = 9` | `n = 10` | `n = 11` |
|---|---|---|---|
| LNAs outside a quipu class | 9 | 262 | 2647 |
| of those, with two heavy clusters | 0 | 2 | 42 |
| with two clusters and a **free arrow between them** | **0** | **0** | **0** |
| with the cluster touching an end | 8 | 221 | 2119 |

So every unclassifiable LNA at every length worked on here is **one overlapping
cluster, usually against an end**. Two clusters with a free arrow between them
first fit at `n = 10`, and every LNA that has them is in a quipu class. The
*barricade* -- two clusters walling in a two-arrow relation -- needs 12 arrows and
first fits at `n = 13`.

### 2. The free move, asked one relation at a time (F-041)

*Measured before the Coxeter guard of F-038 existed. Re-run after the merge at
`n = 8` and `n = 9`, both rows identical; the `n = 10` row below has not been
re-measured. E-038.*

H-012 asks whether deleting a two-arrow relation is a mutation equivalence. Two
changes to how it is asked:

* **one deletion at a time.** The whole strip is a composition of single
  deletions, and each of those is a much shorter journey.
* **meeting in the middle.** `search.meetingPoints`: two searches that reach the
  same quiver have joined their algebras, at twice the depth for the same cost.
  Nothing in the pipeline did this -- `resolveMergeCandidates` keeps only the
  lines a search lands on and throws the rest of the tree away.

| | `n = 8` | `n = 9` | `n = 10` |
|---|---|---|---|
| single deletions | 572 | 2002 | 7072 |
| joined by the known moves | 562 | 1937 | 6768 |
| joined by meeting at depth 3 | 10 | 52 | 206 |
| left | **0** | 13 | 98 |

`n = 8` is settled outright. Of the 13 left at `n = 9`, eleven have both sides
almost separate with the same quipu, so the theorem already calls them one
derived class; only `2302330 / 2300330` and `3302302 / 3300302` are outside
everything, and they do not meet at depth 4 either. A depth-4 pass over the 98
left at `n = 10` joins none of them, at 1060s: depth 4 buys nothing over depth 3
here, on either length, which says the remaining pairs are either much further
apart than 4 + 4 mutations or not joined at all.

### 3. The barricade, built on purpose

The shape H-012's doubts are about, built at the lengths where it first fits: a
heavy cluster, a free arrow, a two-arrow relation, a free arrow, another heavy
cluster. `freeMoves.orbitOf` walks the moves out of one row, which is what makes
a length-13-to-16 question affordable at all -- `derivedOrbits` would partition
208012 rows to answer it about one.

**The moves get it out, in all 95 shapes tried.** Every pairing of six heavy
clusters on the left and right, at gaps of one to three arrows on each side and
lengths 13 to 16: `33000200330` at `n = 13` reaches its strip in an orbit of a
thousand rows, and so does every wider version. **The barricade does not trap a
two-arrow relation.**

**Two false starts, and both were the measurement rather than the mathematics.**

* Run first with the rule table left out -- F-032 having found that the double
  mutation subsumes it at `n <= 10` -- the same barricades came out **not**
  joined, with orbits of 53 to 89 rows against 1073 with the table. The table is
  not subsumed at these lengths.
* With the table in, 49 of the 95 still came out not joined -- and every one of
  them had an orbit of *exactly the 20000-row cap*, so what was measured was the
  budget. Walking from both ends instead (`freeMoves.movesJoin`, the move-level
  twin of `search.meetingPoints`) joins **all 49 in 45 seconds**, most of them
  instantly. A one-way walk that runs out of budget says nothing whatever, and
  the orbit size is the tell: if it equals the cap, there is no result.

So the shape H-012's doubts are about does not hold the relation in, at the
lengths where it first exists. What is untested is the user's fuller version --
four clusters at `n = 30` to `50` -- and the shapes here are the minimal ones.

### 4. Bounded overlap is not the pattern (F-042)

The generalisation asked for was a weaker version of "almost separate": overlaps
of at most two, or at most so many overlaps above one. Crossing those coordinates
against membership of a quipu class kills both. The cell `(max overlap 2, exactly
one of them)` -- the smallest possible step past the theorem -- already holds
outsiders at `n = 9` (`3033030`), `n = 10` (12 of them) and `n = 11` (84).

Chasing what those outsiders have in common instead: the ones with the most free
arrows are the same two little cores at every length, `45` and `504`, and sliding
`45` along the quiver gives a clean law -- inside when it touches the source or
comes within one arrow of the sink, outside everywhere in between, with the band
growing by one place per vertex. `0450000` at `n = 9` is inside, `04500000` at
`n = 10` is not. The condition that decides it is **where the cluster sits**, and
no condition on the relations by themselves can see the difference. F-042.

### 5. What the quipu members look like (H-017)

Walking every LNA outside a quipu class at `n = 9` to depth 5 and grouping the
quipu algebras reached by (cords, relations): the pairs that occur are (1,2)
(1,3) (1,4) (1,5) (2,3) (2,4) (2,5) (2,6) (3,4) (3,5), and **never one with
relations at most cords**. Reading `relations - cords` as a defect, a quipu class
is defect `<= 0` and these are all `>= 1`; the minimum over a class is 2 for
`3033030` and 1 for the eight-member class. That is the beginning of a normal
form, and it makes a sharp prediction about the polynomial-only candidates of
F-034 that sit below the diagonal. H-017.

---

## E-036 — What conditional deeper probing costs, and what it reaches at n = 9
*2026-09-18* · **the parallel-arrow region is reached by exactly one of the nine leftover members at n = 9; giving its branches two more mutations costs 2-4% of the run, and nine mutations into the region it still reaches nothing but its own dual** → H-016

A depth-bounded search gives every branch the same budget. `search.DeeperWhen`
gives extra mutations to the branches that reach a quiver meeting a condition,
with a per-branch budget so the walk still terminates; `merges.py --deeper-on`
asks for it from a command line. This measures the two things that decide
whether it is worth using: what it costs, and whether the extra depth it buys
reaches anything.

### 1. Cost and gain at n = 9

Every member of every leftover orbit at `n = 9` (9 members: the eight of
`3345000` and `3033030`), searched from itself and from its relation dual, plain
against `parallel-arrows:2` — two extra mutations for a branch that reaches a
quiver with a parallel pair, once per branch. Four cores; the seconds are the
sum over the nine members, not wall clock.

| depth | firings | grants | LNAs gained | LNAs lost | plain | probed |
|---|---|---|---|---|---|---|
| 4 | 58 | 4 | 0 | 0 | 89s | 91s |
| 5 | 600 | 32 | 0 | 0 | 418s | 435s |

**Every firing is `3033030`.** The other eight members never reach a quiver with
parallel arrows at all at these depths, so the probe is free for them and the
whole cost is one member's: 19s to 36s at depth 5, where raising the depth for
everyone from 5 to 7 would cost about twenty-five times the run. That ratio is
the case for the mechanism, and it is the one thing here that is not about
parallel arrows in particular.

**Nothing is gained, to depth 5.** The same LNAs are reached either way. This is
a weaker negative than E-035's: `3033030` is alone in its orbit *and* alone in
its Coxeter polynomial group, so the only outcome that would show here is it
reaching a **seeded** LNA — a leftover turning out to be in a quipu class after
all — and at depth 5 it reaches one LNA, its own dual.

### 2. One pass against two

Whether granting depth inside one search differs from recording the interesting
quivers and searching from those afterwards. The second contains the first, and
the containment can only be strict for a branch that leaves the region and comes
back, since a firing *below* the one that bought the depth is already inside the
subtree the grant paid for, at exactly the remaining depth a re-search would
give it.

| start | depth | extra | plain | one pass | two passes |
|---|---|---|---|---|---|
| `3030` | 4 | 3 | 50 | 105 | 105 |
| `3030` | 5 | 2 | 118 | 195 | 195 |
| `30330` | 4 | 2 | 67 | 86 | 86 |
| `30330` | 5 | 2 | 158 | 237 | 237 |

Nodes reached, counted by the quiver with its arrow names and its relations.
**Equal in every case**: no branch at these sizes leaves the region and returns
within the depth searched. So the one-go run is not the weaker of the two in
practice, which is the practical answer — the two-pass route's advantage is that
the count of firings is visible before the second round is paid for, not that it
reaches more.

### 3. The one member that reaches the region, pushed to depth 9 inside it

`3033030` is the only member of either leftover orbit whose walk ever reaches a
quiver with parallel arrows, so it is the whole of the `n = 9` test and it is
cheap. From it and its relation dual, plain against `parallel-arrows:2`:

| depth | firings | grants | deepest firing | reaches | seconds |
|---|---|---|---|---|---|
| 6 | — | — | — | itself | 98s |
| 6 + 2 | 4322 | 154 | 8 | itself | 221s |
| 7 | — | — | — | itself | 350s |
| 7 + 2 | 29122 | 826 | 9 | itself | 1200s |

The four ran together on four cores and the last two shared the machine with a
test run, so the seconds are an upper bound and the ratio between them is the
part worth reading.

"Deepest firing" is the length of the longest mutation path at which the
condition still held, so the last row walked **nine** mutations into the region.
It reaches nothing but its own relation dual, which is what depth 5 already
reached.

This is the `n = 9` half of what H-016 asks for, and past the depth it asks for.
It does not settle H-016: `3033030` is alone in its Coxeter polynomial group, so
the only thing it *could* show is a leftover turning out to be in a quipu class,
and one member at one length is not the hypothesis. `n = 10` and `n = 11`, where
H-013's leftover orbits sit, have not been looked at this way. But it is a
negative at a depth and a length E-035 could not reach, and it cost 20 minutes
on one core rather than the run over every member that a uniform depth 9 would
have been.

### 4. Over a whole classification

`classify.py --deeper-on` gives one probe to all three searching steps. At
`n = 8`, the shortest length whose classification needs a search at all:

| | classes | rows | firings | grants | seconds |
|---|---|---|---|---|---|
| plain | 11 | 429 | — | — | 38s |
| `parallel-arrows:2` | 11 | 429 | 44 | 4 | 39s |

**The answer does not move**, which is the check that matters: extra depth may
place a row the depth could not reach, and may never place one differently. The
published table of arXiv:2305.06642 is what both are checked against, as
`tests/test_classify_end_to_end.py` does. The condition does fire here, four
times buying depth, so this is the probe running over a real classification
rather than a no-op.

### Reproducing

Section 4 is `classify.py 8 --quiet` with and without
`--deeper-on parallel-arrows:2`.

Section 1 is `merges.py 9 --depths 4 5 --all-groups` run twice, once with
`--deeper-on parallel-arrows:2` and once without, comparing `lnasReached` and
`deeperFirings` per checkpoint record. The two runs can share one checkpoint
file: it records the condition, and a plain record does not count as covering a
probed search or the other way round.

Section 2 is the computation of
`tests/test_deeper_probing.py::test_recording_and_searching_again_contains_deepening_in_one_pass`
at the four settings in the table.

Section 3 is `merges.searchFrom((9, (3, 0, 3, 3, 0, 3, 0), depth, spec))` for
each row.

---

## E-035 — Lifting the parallel-arrow restriction, and re-measuring E-033
*2026-09-18* · **every wrong-key node at n = 6 and n = 7 to depth 5 was a mis-count; the new region reaches nothing new at these sizes** → F-039, R-013, H-016

E-033 walked the search tree with the Coxeter key in hand and split the nodes
where it had moved into "parallel arrows, harmless" and "clean, the real fault".
This is the same sweep after the relations were moved onto paths that name their
arrows, so a parallel pair can be stated, counted and mutated at.

### What was changed

* `arrowPaths` — an arrow is `(tail, head, key)`, a path is a tuple of arrows.
* `procedure` — steps 1 to 7 per arrow; step 5 divides by the arrow, step 7 reads
  each candidate's first arrow back as its relation and its tail back into the
  old quiver, step 6 and the carried-past relations name the composite they use.
* `procedure.isMutable` — no longer refuses a quiver with a parallel pair.
* `invariants` — the Cartan matrix counts arrow paths, exactly and cheaply.
* `search` — the illegal-relation check is over arrow relations.
* `arrowPaths.homDimensionByClosure` — the cheap count closes the commutativity
  relations to a fixed point, where `paths.numberOfPathsUpToRels` applied each of
  a subset once.
* `invariants.integerCartanMatrix` — and the key is **exact** wherever the cheap
  count is not provably right, which is wherever the ideal is not monomial. No
  closure makes the cheap count see a relation of three or more paths, and step 4
  produces one at every vertex with three arrows out.

### 1. The sweep, before and after

`mutationSearchDepthFirst` from every LNA of the length, `coxeterGuard = False`
so the corrupt region is visible, `coxeterKey` compared at every node against the
start's. One core.

| | | nodes | wrong key | parallel | clean | lines | seconds |
|---|---|---|---|---|---|---|---|
| `n = 6`, depth 5 | before | 14,693 | 4 | 4 | 0 | 1,789 | 28 |
| | after | 14,701 | **0** | 0 | 0 | 1,789 | 27 |
| `n = 7`, depth 5 | before | 94,446 | 79 | 75 | 4 | 7,175 | 254 |
| | after | 94,498 | **0** | 0 | 0 | 7,175 | 268 |
| `n = 7`, depth 6 | before | 336,760 | 339 | 277 | 62 | 20,683 | 866 |
| | after | 337,360 | **10** | 0 | 10 | 20,683 | 949 |

The node count rises by 8, 52 and 600: a parallel-arrow node is no longer
terminal. Lines are unchanged, exactly, at all three. The cost is 10% at depth 6,
which is the exact Cartan matrix against the cheap one, and it is not optional.

Three mechanisms, all three mis-counts:

* the *parallel* nodes -- 4, 75 and 277 of them -- had the key computed over
  vertex sequences, so a parallel pair contributed 1 to the Cartan matrix where
  it contributes 2;
* some *clean* nodes are the incomplete closure. The four at `n = 7` depth 5 are
  `40030` by `[1, 4, 2, 5, 2]`, by `[1, 2, 4, 5, 2]`, by `[1, 1, 5, 4]` and one
  more;
* the rest are relations of three or more paths, which the cheap count has no
  reading of and ignores, e.g. `34400` by `[1, 3, 4, 2, 2, 1]`.

The ten that survive at depth 6 are **one bad step and its descendants**: all ten
are reached from `30330` along the `[4, 1, 3, 1, ...]` family. Replaying `[1, 4, 2, 5, 2]` in both engines gives the **same quiver and
  the same `rels`**, and the key holds in one and moves in the other, which is
  what says it is the measurement.

### 2. What still moves the key

From the relation dual of `33030` at `n = 7` by `[4, 1, 3, 1, 3, 3]`, F-038's own
smallest clean case, the key moves at step 6 under the arrow model too, and there
the cheap and the exact Cartan matrices agree -- so it is the algebra. R-012
stands and the guard stays.

### 3. The classification is unchanged

`classifyLength` on both engines, same machine:

| | classes | rows | before | after |
|---|---|---|---|---|
| `n = 6` | 4 | 42 | 0.3 s | 0.3 s |
| `n = 7` | 6 | 132 | 1.2 s | 1.3 s |
| `n = 8` | 11 | 429 | 29.9 s | 33.4 s |

Class for class, size for size, identical at all three, and about 10% dearer --
the exact Cartan matrix where the ideal is not monomial, against the cheap count
everywhere.

### 4. What the region beyond a parallel pair looks like

From `3030` at `n = 6`, mutated at `[1, 3, 4, 1, 4]`, the quiver has two arrows
`1 -> 6` and one relation between the two parallel paths `5 -> 1 -> 6`. The old
gate refused all six vertices. The paper's criterion admits three, every one
keeps the Coxeter key, and mutating at 3 comes back out to a quiver with **no**
parallel arrows and a genuine commutativity relation
`5 -> 1 -> 3 -> 6 = 5 -> 1 -> 6` -- a quiver the search could not reach at any
depth before.

### 5. Does the new region reach anything? Not at these sizes

The point of walking through a parallel-arrow node is what lies beyond it, so:
for every LNA of the length **and its relation dual**, the set of LNAs the
guarded search reaches, with the gate allowing parallel arrows and with it
refusing them, at the same depth.

| | starts | lines reached, allowing | refusing | starts gaining | losing |
|---|---|---|---|---|---|
| `n = 6`, depth 6 | 74 | 801 | 801 | 0 | 0 |
| `n = 7`, depth 6 | 244 | 3,390 | 3,390 | 0 | 0 |

(The counts are the sum over starts of how many LNAs that start reaches, so a
line reached from two starts counts twice; what matters is that no start gained
or lost one.)

**So the region is reachable and walkable and yields nothing new here.** Two
reasons not to read that as "it never will". The comparison holds the *depth*
fixed, and entering the region and returning from it costs steps, so at depth 6
the part of it that can come back to a line at all is thin -- the one walk
looked at by hand, `3030` by `[1, 3, 4, 1, 4]` then 3, takes six mutations to
get back out to a quiver with no parallel pair, and that quiver is not a line.
And `n <= 7` is fully covered by the move rules with no search at all (F-021),
so there is nothing left for a search to find at these lengths whatever it walks
through. H-016.

Reproduce: `tests/test_parallel_arrows.py`; the sweep is a `visitor` on
`search.mutationSearchDepthFirst` comparing `search._coxeterKeyOrNone` at each
node, and the before column is the same script against `git show 78328e7`. The
reachability comparison patches `procedure.isMutable` to pass
`allowParallelArrows = False` and compares `search.linesReachedFrom` either way;
it is 84 minutes at `n = 7` depth 6 on one core.

---

## E-034 — Can the guarded search be fooled where the polynomial is known to fail?
*2026-09-18* · **not at depth 6, at the smallest collision** → H-015

F-038's guard refuses any step that moves the Coxeter polynomial. That is
necessary for a derived equivalence; H-015 asks whether it is sufficient. The
place to look is where the polynomial is known to be blind, and F-010 says
exactly where that is: cospectral quipus, the smallest collision at order 9,

    P^(1,4)_(1,0,1) = A_{9,(1,3)}^{(3,6)}   class 3060000
    P^(1,2)_(1,1,2) = A_{9,(1,4)}^{(3,4)}   class 3004000

Different trees, so **not derived equivalent**, yet one Coxeter polynomial —
confirmed here, both `(1, 1, -1, -3, -4, -4, -3, -1, 1, 1)`. If a guarded search
out of one reaches anything the other reaches, then a step the guard admits is
not a derived equivalence and H-015 falls.

Searching from each, from the member and from the relation dual:

| depth | guard | LNAs from `3060000` | from `3004000` | **shared** | time |
|---|---|---|---|---|---|
| 5 | on | 2 | 7 | **0** | 85 s |
| 5 | off | 2 | 7 | **0** | 45 s |
| 6 | on | 2 | 8 | **0** | 395 s |
| 6 | off | 2 | 8 | **0** | 190 s |

**Nothing shared, either way.** `3060000` is remarkably rigid — it reaches only
`6000030`, its own dual, at either depth — while `3004000` moves to seven or
eight. The guarded and unguarded searches return *identical* sets here, which is
consistent with F-038: at this depth the corrupt region is entered but does not
come back round to a line.

**What this is worth, and what it is not.** It is the sharpest single test
available and H-015 survives it. But a negative is a depth bound, not a proof,
and this probes **one** collision: order 10 has two collision groups and order 11
has four (F-010), none of them tried. It also cannot detect a guard-passing step
between two algebras that are cospectral for some *other* reason than being
quipus.

Reproduce: the pair is `python classify.py 9 --collisions`; the search is
`search.mutationSearchDepthFirst` from each with `coxeterGuard` both ways.

---

## E-033 — Is the mutation search sound?
*2026-09-18* · **no, and the fix costs 1.85× and changes no answer** → F-038, R-012, F-037

Prompted by the two ALARMs of E-032. The question is narrow and had never been
asked directly: **does `mutationSearchDepthFirst` stay inside one derived
equivalence class?** Every step is meant to be a tilting mutation, so the Coxeter
polynomial must be constant over the whole search tree — at every quiver reached,
not only at the lines it records.

### 1. Walk the tree with the invariant in hand

A `visitor` that computes `coxeterKey` at each node and compares it to the start,
aborting at the first mismatch. **It fails at `n = 6` in 7 seconds**: from
`3030`, the path `[1, 3, 4, 1, 4]` reaches a quiver whose key has moved from
`(1,1,-1,-2,-1,1,1)` to `(1,1,0,-1,0,1,1)`. Replaying it one mutation at a time
shows the last step producing a **parallel arrow**, `1 -> 6` twice.

### 2. Classify every wrong node

| | nodes | wrong key | parallel | cyclic | **clean** | lines | lines wrong |
|---|---|---|---|---|---|---|---|
| `n = 6`, depth 5 | 25,398 | 4 | 4 | 0 | **0** | 3,263 | 0 |
| `n = 7`, depth 6 | 609,474 | 604 | 507 | 0 | **97** | 37,911 | 0 |
| `n = 8`, depth 5 | 1,093,976 | 1,204 | 1,030 | 0 | **174** | 55,175 | 0 |

The parallel-arrow ones are a limitation of the model and are **harmless**:
`procedure.isMutable` refuses every vertex of a quiver that has parallel arrows
anywhere, so the node is terminal, and such a quiver can never be mistaken for a
line. Not one oriented cycle appeared at all.

**The clean ones are the fault** — acyclic, no parallel arrows, key moved,
nothing to stop the search descending. Smallest: `n = 7`, from the relation dual
of `33030`, path `[4, 1, 3, 1, 3, 3]`, where step 6 loses a commutativity
relation outright. Written up as F-038, retracted as R-012.

**No answer was wrong at these sizes.** All 96,349 lines collected carried the
starting key. The corrupt region exists but had not reached an answer.

### 3. The guard

`search.mutationSearchDepthFirst(..., coxeterGuard = True)`, now the default,
refuses a step whose key differs from the start's — R-005's third requirement,
applied per step rather than per rule.

| | without the guard | with it |
|---|---|---|
| wrong-key nodes, `n = 6` depth 5 | 4 | **0** |
| lines reached, `n = 6` depth 5 | 3,263 | 3,263 |
| core-seconds, `n = 6` depth 5 | 75 | 138 (**1.85×**) |
| core-seconds, `n = 7` depth 5 | 759 | 1,403 (**1.85×**) |
| LNAs losing a reached line, `n = 6`, `n = 7` | — | **0** |
| LNAs gaining one | — | **0** |

Comparing the *sets* of lines reached, from every LNA and from its dual, at
`n = 6` and `n = 7` to depth 5: not one line lost, not one gained. **Existing
results at these sizes stand unchanged.**

A cyclic quiver has no unimodular Cartan matrix, so `coxeterKey` raises there;
`_coxeterKeyOrNone` returns None and such a step is let through, because the
search does not descend from a cycle anyway. Without that, starting a search at a
cyclic algebra crashed — `tests/test_cycles.py` caught it.

### 4. The links of E-032 re-derived and replayed

A link is a positive claim, so each was found again and checked move by move:
every step admissible, no illegal relation, **no parallel arrow or oriented
cycle**, and the key held. All five pass.

| link | paths into the target | checked | steps |
|---|---|---|---|
| `n = 10` `34504030 -> 50505000` | 19 of 20 lines collected | 5 | 7 |
| `n = 10` `05040330 -> 33460000` | 9 | 5 | 7 |
| `n = 11` `030233030 -> 300330400` | 3 | 3 | 3–4 |
| `n = 11` `346004030 -> 060040400` | 4 | 4 | 3–4 |
| `n = 11` `302340030 -> 300403030` | 3 | 3 | 3–4 |

The first is F-037.

### 5. The ALARM itself, reproduced

The depth-8 search from `03033030` that raised it takes 8257 s on one core, so it
was split by first mutation and the branches run in parallel. It reproduces
exactly, and it is the clean mechanism of part 2:

* From the **member**, ten branches, **not one line reached at all** — everything
  the overnight run recorded from this side is the start itself.
* From the **relation dual**, branch `9` reaches one line correctly, and branch
  `4` reaches **four lines, all four wrong**, all of them `30233330`, which is a
  member of orbit `00330400`. No parallel arrows, no oriented cycle.

Bisecting `[4, 6, 4, 6, 9, 4, 4, 6]`, every step admissible and acyclic
throughout:

    steps 1-6   key (1, 1, -2, -3, 1, 4, 1, -3, -2, 1, 1)   held
    step 7 at vertex 4   -> (1, 1, -1, -2, -1, 0, -1, -2, -1, 1, 1)   moved
    step 8 at vertex 6   lands on the line 30233330, outside the class

So the ALARM was neither a broken orbit (H-013's guess) nor the alarm test's own
`polyOf` defect (E-032, part 4): it is one admissible-but-not-derived-equivalent
mutation at step 7 of eight, and the guard refuses it.

Re-running that one branch both ways settles it:

    coxeterGuard = False   4 lines, 4 WRONG, 534s   ['30233330']
    coxeterGuard = True    0 lines, 0 WRONG, 941s   []

1.76× here, and the false answer is gone.

**This is why the small-`n` sweeps found nothing.** The corruption needs a deep
enough tree to come back round to a line: seven clean steps, one bad one, then
one more. At `n ≤ 8` and depth ≤ 6 the corrupt region is reached but never
returns to a line, which is exactly what part 2's "lines wrong: 0" column says.

### 6. Which existing results the fix disturbs: none found

The guard only ever *refuses* a step, so it can only shrink what a search
reaches. Every **negative** in the record — "reaches nothing seeded", "stayed
apart to depth 8" — is therefore untouched. Every **positive** needed checking,
and the shallow ones are where most of them live: F-034 and H-014 rest on
depth-3 walks, F-036 and E-031 on depth-4 searches.

Comparing the sets of lines reached with the guard and without, from each LNA and
from its dual:

| | sampled | LNAs where the guard changes what is reached |
|---|---|---|
| `n = 10`, depth 3 | 300 of 4862 | **0** |
| `n = 10`, depth 4 | 150 of 4862 | **0** |
| `n = 9`, depth 4 | 300 of 1430 | **0** |
| `n = 6`, depth 5 | all 42 | **0** |
| `n = 7`, depth 5 | all 132 | **0** |

Together with parts 2 and 4 this says the damage was confined to depth 7 and
beyond, and the only wrong answer the record ever contained is the ALARM itself,
which was already excluded from E-032's unions by the alarm test. **E-032's
conclusions stand unchanged: 43–46 classes at `n = 10`, 84–115 at `n = 11`.**

### 7. The `break` that meant `continue`

`search.py` abandoned every remaining vertex at a node as soon as one mutation
there produced an illegal relation, rather than just that vertex — and the loop
runs over `reversed(vertices)`, so it lost every lower-numbered one. Fixed.
**It never fired in E-032**: `isIllegalRelation` prints when it triggers and all
three overnight logs contain zero such lines over 121 core-hours.

---

## E-032 — The overnight run: H-013 at `n = 10` and `n = 11`
*2026-09-18* · **`n = 10` settled to 43–46 classes; one prediction wrong; and an ALARM that indicts the search rather than the orbits** → E-033

`python overnight.py`, started 2026-09-17 15:46, budget 9 h, on 16 cores under
WSL. Three jobs, and all three reached a definite state:

| job | command | outcome |
|---|---|---|
| `merges10` | `merges.py 10 --depths 5 6 7 8 --jobs 7` | **finished**, exit 0, 00:38 |
| `merges11` | `merges.py 11 --depths 4 5 6 --jobs 7` | stopped on its budget, exit 2, 00:51 |
| `classify10` | `classify.py 10 --resume` | terminated at the deadline, 01:07 |

279 searches at `n = 10` (57.7 core-hours) and 3133 at `n = 11` (63.3), by depth
`{5: 122, 6: 122, 7: 26, 8: 9}` and `{4: 2415, 5: 718}`. Depth 7 and 8 are thin
because `settled(poly)` stops searching a group once its orbits have merged. The
slowest single search was 9691 s, at depth 8.

### 1. `n = 10` is finished, and H-013 was right twice and wrong once

| polynomial group | orbits | H-013 predicted | what happened |
|---|---|---|---|
| `(λ-1)²(λ+1)²(λ²+1)(λ⁴+λ³+λ²+λ+1)`, `C(2,4,5)`'s | 2 (69, 42) | merge, depth ≤ 7 | **merged at depth 6**, 9 searches found it |
| `(λ-1)²(λ+1)²(λ²+λ+1)(λ⁴-λ²+1)` | 4 (4, 2, 2, 1) | at least two stay apart | **all four stayed apart to depth 8** |
| `(λ+1)²(λ²-λ+1)(λ⁶-λ³+1)` = `T¹⁰+T⁹+T+1` | 2 (1, 1) | **no link to depth 8** | **merged at depth 7**, from both sides |

So 12 orbits fall to **at most 10** non-quipu classes, and the derived classes at
`n = 10` number between `36 + 7 = 43` and `36 + 10 = 46`, where H-013 said 43–48.

The wrong prediction is the one worth keeping. `34504030` and `50505000` are the
two the Coxeter polynomial can never separate — the pair `remark:Coxeter` of
arXiv:2310.08346 makes its point with, and the ones H-013 called singletons under
every move known. A depth-7 search links them in both directions, in 553 s from
`34504030` and 1951 s from `50505000`. **A link is a positive claim and the
search is now known to be unsound in some places, so this one is checked
separately in E-033 rather than believed here.**

### 2. `n = 11` got through depth 4 and most of depth 5

54 orbits in 20 polynomial groups, down to **at most 51** non-quipu classes, so
between `64 + 20 = 84` and `64 + 51 = 115`. Three merges, all first seen at depth
4 and all confirmed from both sides:

    030033030 <-> 300330400     28 searches
    040344030 <-> 060040400     12
    300340030 <-> 300403030     12

No group beyond those merged at depth 5, and the two largest groups — 9 orbits
each — stayed entirely apart. Rerunning the same command resumes from the
checkpoint.

### 3. `classify10` got nowhere worth keeping

1010 of 4862 rows attempted in 9 hours, 31 resolved and 5 still unresolved at
depth 6. It is an independent route to the `n = 10` count and it is far too slow
to be one; `merges.py` answered the same question in a fraction of the time
because it only searches what the moves leave over. **Do not schedule
`classify.py 10` again as a whole-range run.**

### 4. The ALARM

Two depth-8 searches, from `03033030` and `30330300`, reported reaching orbit
`00330400`, which carries a different Coxeter polynomial. H-013 said an alarm
"would refute F-032's orbits". It does not: both orbits are internally
consistent, each carrying exactly one polynomial over all its members (4 and 142
of them). What the alarm indicts is the **search**. That is E-033.

Two things about the alarm test itself, found while reading it:

* `polyOf` is built from the leftover orbits only, so `polyOf.get(other)` is
  `None` for any orbit a seed reaches. A leftover linking to a **quipu** orbit —
  which would be the most interesting result the run could produce, an LNA the
  theorem misses turning out to be in a theorem class after all — is therefore
  reported as an ALARM and **excluded from the union**, not recorded as a merge.
  It did not happen here, but the test would have hidden it.
* An alarm does not stop the run or mark the group, so the two alarms sat in the
  log for five hours.

---

## E-031 — Is reorientation a mutation, and does it help the search?
*2026-09-17* · **yes, and it is the merge step rather than the search that needed it** → F-036

Suggested from outside the code: all orientations of a relation-free tree should
be mutation equivalent, by mutating at a source -- which only flips the outgoing
arrows -- and then at each newly created source, with left mutation at sinks for
the other direction. The classification has been merging classes on a shared
hereditary form all along, which is a statement about mutation classes resting on
a fact about derived ones, so this is the step in between.

### 1. It is true, and the sequence can be written down (F-036)

Right mutation at a source of a relation-free tree quiver reverses exactly the
arrows there, creates no relations and does not renumber -- 2339 cases over every
tree and orientation at orders 3 to 7, no failures, and the same for left
mutation at a sink. `reflections.reflectionSequence` turns one orientation into
another by flipping every vertex on one side of a differing edge, in topological
order; verified against the engine at orders 4 to 10, every step admissible, 5568
reorientations, no failures.

Right mutations **alone** already connect every orientation of every tree up to
order 8 -- which matters because the search walks only those -- but the worst
distance grows: 4, 6, 9, 12, 16 at orders 4 to 8, against 3, 3, 6, 6, 10 when
sinks are allowed too.

### 2. It does *not* widen a bounded search

Allowing a search to jump to any orientation whenever it reaches a relation-free
quiver (`reflections.linesReachedThroughReflections`), against the plain search
from the same start:

| start | depth | plain | with reorientation |
|---|---|---|---|
| `kA_7` | 2 | 5 LNAs | 5 |
| `kA_7` | 3 | 8 | 8 |
| `kA_8` | 2 | 5 | 5 |
| `kA_8` | 3 | 9 | 9 |

Not one new LNA, and the reason is in part 1: a right mutation at a source *is* a
reflection, so the plain search already performs them where they are cheap. What
the lemma adds is the reflections that are **not** cheap -- and those do not lead
anywhere new within the depth either. Seeding backwards, from every one of the
128 orientations of each quipu of order 8 at depth 2, reaches 3 LNAs for the line
and 0 or 1 for every other quipu: the hereditary side is simply a long way from
any LNA.

### 3. It is the merge step that needed it

`mergeReport` merges two classes when both reach the same tree. Searching every
LNA to depth 4 and pairing the ones that reach the same tree:

| | `n = 7` | `n = 8` |
|---|---|---|
| LNAs reaching a relation-free quiver | 23 of 132 | 26 of 429 |
| pairs reaching the same tree | 79 | 80 |
| of those, pairs reaching **isomorphic quivers** | 61 | 65 |
| pairs whose orientations must be joined | **18** | **15** |

Those merges are derived equivalences until the orientations are joined, and
joining them takes 4 to 11 mutations, which is beyond any depth the pipeline runs
at. `reflections.mutationBridge` produces the whole path and the engine confirms
where it lands.

**A false start worth recording.** The first version of
`relationFreeQuiversReached` recorded what the search of the *opposite* algebra
found without carrying it back through the opposite, so half the orientations in
each list were reversed. It made the answer to part 3 come out as 4 of 79 at
`n = 7` and 0 of 80 at `n = 8`, where the true counts are 18 and 15 -- the
measurement looked like a nearly-empty result when it is a fifth of the pairs.

### 4. What it does not do: H-012

The natural hope was that the bridge would close H-012's gaps -- the 8 LNAs at
`n = 8` and 44 at `n = 9` that the known moves do not join to their stripped
form. It cannot: **none of those 52 rows reaches a relation-free quiver at all at
depth 4**, on either side of the pair, so there is no tree to bridge through.
Recorded against H-012.

---

## E-030 — Two other families, measured by Coxeter polynomial
*2026-09-17* · **no tree outside the quipu shape; quipus *with* relations carry every class the theorem misses** → F-033, F-034, F-035, H-014

The question: the quipu theorem names a class by a tree with no relations, and
every LNA it does not cover has to be classified some other way. Is there a
second family that does for those what quipus do for the almost separate ones?
Two candidates were measured, both through the Coxeter polynomial, which is a
derived invariant and forces the number of simples — so a candidate can only
match an LNA of its own length, and a difference settles it.

### The instrument

`invariants.coxeterCoefficients`. For a quiver with no oriented cycles the Cartan
matrix is unimodular, so

    det(lambda I - Phi) = det(lambda C^T + C),

which is a determinant of integers and `lambda` with no inversion in it. Taking
it at `n + 1` points and interpolating gives the polynomial as an integer
coefficient tuple: exact, hashable, and a hundred times faster than the symbolic
route — 16796 LNAs at `n = 11` in 19 s.

`coxeterTables` turns that into the three tables a search matches against: every
LNA's polynomial, every quipu's, and every LNA's **status** — in a quipu class by
F-032's moves (QUIPU), in none because no quipu of the order carries its
polynomial (NOT_QUIPU), or neither (UNPLACED). It reproduces F-032 exactly:
1421/9/0 at `n = 9`, 4600/260/2 at `n = 10`, 14149/2631/16 at `n = 11`.

### 1. Every tree, not only the quipus (F-033)

`python families.py trees 9 10 11` — and 12.

| order | trees | not quipus | cospectral with a quipu | **leads** |
|---|---|---|---|---|
| 9 | 47 | 29 | 3 | **0** |
| 10 | 106 | 70 | 0 | **0** |
| 11 | 235 | 171 | 7 | **0** |
| 12 | 551 | 424 | 15 | **0** |

F-031 asked this of the trees of maximum degree three; these are all of them,
694 non-quipu trees over the four orders. 25 of them do share a polynomial with
some LNA — but in every case every LNA under that polynomial is one the moves
place in a quipu class, and the tree is cospectral with that quipu rather than
isomorphic to it. Two tree algebras are derived equivalent exactly when the trees
are isomorphic, so those are refuted outright, with no search.

### 2. Quipus that carry relations (F-034)

`python families.py quipus 9 --min-arrows 2`, and the same at 10 and 11.

Enumerated: each quipu of the order, each orientation of its edges up to the
tree's automorphisms, each admissible monomial ideal — an antichain of directed
paths under "is a contiguous subpath of". The relation-free ideal and the
linearly oriented line are left out, being the quipu theorem's own case and the
LNAs themselves.

| order | shortest relation | ideals walked | matching | polynomials covered |
|---|---|---|---|---|
| 9 | 2 arrows | 370 483 | 3677 (2820 up to isomorphism) | **2 of 2** |
| 10 | 2 arrows | 3 411 263 | 306 624 | **7 of 7** |
| 11 | **3 arrows** | 5 465 194 | 1 246 011 | **20 of 20** |

Every Coxeter polynomial of an LNA outside a quipu class is carried by quipu
algebras with relations, in quantity — 1746 of them on `3033030` alone at
`n = 9`, over 16 different quipu shapes.

**Confirmed from the other side.** A polynomial match is necessary and not
sufficient, so the classes were also walked: every quiver a mutation search out
of an LNA reaches is in its class by construction, and `reachedQuipuAlgebras`
keeps the ones that are quipus with monomial relations. Walked from **every** LNA
outside a quipu class, and its dual, to depth 3:

| | `n = 9` | `n = 10` | `n = 11` (sample of 200) |
|---|---|---|---|
| LNAs outside a quipu class | 9 | 262 | 2647 |
| reaching a quipu with relations | **9** | **262** | **200 of 200** |
| reaching none | 0 | 0 | 0 |
| confirmed per LNA: min / median / max | 8 / 16 / 18 | 3 / 20 / 84 | 5 / 27 / 117 |
| distinct algebras confirmed | 178 (depth 4) | 3510 | -- |

Every algebra reached this way is in the enumeration — the only things reached
and not enumerated were the linearly oriented lines, which are the LNAs. The
`n = 10` walk takes 7 minutes and `n = 11` would take an hour and a half, so
`n = 11` is a random sample of 200 of its 2647 rows, seeded, and every one of
them reaches something too.

**And what they reach is not arbitrary.** `python families.py members 9` prints
the confirmed members per class, simplest first, and seven of the nine LNAs at
`n = 9` reach the *same quipu* — `P^(6)_(1,1)`, the line on eight vertices with a
pendant at the second — carrying their own relations shifted by one vertex with
one absorbed into the branch. That is the raw material for H-014's second part,
which is the part that would make this a theorem rather than a census.

**The smallest instance of the phenomenon is at order 4.** `D_4` with a single
two-arrow relation has the Coxeter polynomial of `kA_4`, and **one** mutation
takes it there. So a quipu with relations being derived equivalent to a line is
not exotic; what is new is that it happens for the lines the theorem cannot name.

### 3. A relation of two arrows is not free on a quipu (F-035)

`python families.py free 4 5 6 7 8`. Deleting every two-arrow relation and
comparing the polynomial:

| order | ideals with one | polynomial kept | **changed** |
|---|---|---|---|
| 4 | 11 | 7 | **4** |
| 5 | 72 | 48 | **24** |
| 6 | 543 | 300 | **243** |
| 7 | 4160 | 2138 | **2022** |
| 8 | 34938 | 15337 | **19601** |

The control passes: on the **linearly oriented line** the polynomial is kept
every time, which is `corollary:lengthtworelations` of arXiv:2310.08346 and
F-028. Off the line it fails immediately — the order-4 counterexample above is
the smallest. So `--min-arrows 3`, which is what makes order 11 affordable, is a
real restriction of the family and not a normalisation, and the order-11 run is a
statement about ideals whose relations all have three arrows or more.

### 4. Relation-free sightings, as instrumentation

`search.relationFreeSightings` records every quiver a search reaches with no
relations left, with whether its underlying graph is a tree, whether it is a
quipu, whether the quiver has an oriented cycle, and the path that got there;
`classify.py --sightings FILE` writes them as JSON lines. Nothing is recorded
unless a sink is open, so it costs nothing when it is not asked for. It answers a
question nobody had asked of the searches: the hereditary algebras they pass
through are collected and all but the first thrown away, and whether any of them
is *not* a tree has never been looked at.

**First measurement: `python classify.py 9 --sightings` records none at all.**
That is the expected answer and worth having written down. The quipu classes are
named by the theorem without a search, and the only rows the classification does
search at `n = 9` are the nine in no quipu class -- which contain no hereditary
algebra, so there is nothing for a search to find. The instrument will only have
something to say where a search runs into a class that *does* have one, which
means `--form-depth` at a length with an unnamed class, or the deeper resolving
runs of `merges.py`.

---

## E-029 — Reading `proposition:doubleMutation`, and what it does to coverage
*2026-09-17* · **a mechanism, found by reading; n = 9 needs no search, n = 10 is 16 orbits** → F-032

The first session run locally, with the `.tex` sources of all three papers to
hand. NOTES item 5 said to look for mechanisms before rules; this is the third,
and the largest.

1. **Read** the proposition and its proof (`main.tex` lines 331–497 of
   arXiv:2310.08346v1). Replaced the second-hand lead in the literature summary.
2. **Implemented** it as an interval rewrite, dual through F-026, `s = 1`
   extension behind `allowSource`. Reproduced `example:mutationToA11_5` steps
   `L_8`, `L_9` and `example:A13tworelations` steps `R_1`, `R_1^2` exactly.
3. **Verified** against the engine at `n = 5..10`, 14 processes: 17556
   confirmations, 0 failures, 1 min 22 s for `n = 9, 10`.
4. **Measured** coverage with and without the table, the free move, and the
   extension. The extension changes nothing (`doubles, paper only` gives the same
   numbers at `n = 7, 8, 9`). Interior-only (`s > 1` and `t < n`) gives 273 / 429
   and 770 / 1430 at `n = 8, 9`, with or without the free move.
5. **Named what is left** by Coxeter polynomial against the quipu polynomials of
   the order (34 at `n = 10`, 60 at `n = 11`, matching F-010's collision counts),
   and by the two implemented non-piecewise-hereditary certificates.
6. **H-012, by known moves only**: with the double mutation, edge moves and the
   rule table (all mutations), 8 of 429 LNAs at `n = 8` and 44 of 1430 at `n = 9`
   are not joined to their stripped form. That is a statement about the known
   moves, not a counterexample; nothing settled.
7. **H-011's sharp question**: of the 262 rows left at `n = 10`, 190 have a
   relation at both source and sink (73 %); at `n = 11`, 1743 of 2647. Not a
   characterisation — but the question has changed, since what is left is no
   longer a gap in the rules (F-032).
8. **What an interior application does** (`s > 1`, `t < n`), `n = 6..10`: at
   `n = 10`, 7164 applications; maximum overlap down in 1416, up in 1416, same in
   4320; relation count −1 in 1430, +1 in 1430. So interior moves *do* lower
   overlap when bystanders cross `r`, while on an isolated pair they only slide
   it. A claim to the contrary was written into H-010 and corrected before
   commit.
9. **Relation dual.** The first `merges.py` smoke test at `n = 10`, depth 3,
   found 18 links and every one was a relation dual, reached at depth 0 because
   the search starts from each member's dual. Closing the orbits under it is free
   and takes `n = 10` from 16 leftover orbits to **12 in 7 polynomial groups**,
   `n = 11` from 86 to **54 in 20**. So the derived classes number 43 to 48 at
   `n = 10` and 84 to 118 at `n = 11` (H-013).
10. **Reading, not yet used**: `lemma:taupathimpliesnotpwh` certifies both members
   of `example:A10double`, which A9, A13 and vertex deletion all miss.

11. **Sweep for further mechanisms, negative.** arXiv:2310.08346 §2 has exactly
    two derived equivalences (the free move and this one). arXiv:2305.06642's
    `algorithm:CRswap` needs almost separate relations and is what the seeding
    already encodes; `corollary:2Rels` is the free move restricted to that case.
    arXiv:2112.08129 is the procedure itself. No further *merge* mechanism in the
    three papers; the one unused tool is a *certificate*,
    `lemma:taupathimpliesnotpwh`.

Reproduce: `python overlaps.py 9 10 11 --free --doubles --no-rules`, and
`pytest tests/test_double_mutation.py` (the `n = 9, 10` engine checks are
`slow`).

**What this makes obsolete in the handover plan.** The raised-bound rule
discovery (`--max-arrows 7 --max-width 8` interior, `--extend` again): its purpose
was coverage, and the table adds nothing at `n = 10` on top of this. Not run.

---

## E-028 — Checkpointing the classification, and what a resumed n = 8 gives
*2026-09-17* · **the published table, out of two interrupted halves**

E-008 records that a long classification does not survive the night, and F-014
records the cost: the n = 10 run named 61 classes and kept none of them, because
only the search step wrote anything as it went. The naming and resolving steps
now write the table after every class and record what they finished in a sidecar
JSON file beside the CSV.

Checked by interruption rather than by argument. `classifyLength(8, ...)` with a
budget of zero seconds stops inside the search with 10 of the 429 rows still
unplaced and exits 2; resumed, it finishes and gives

    133, 65, 64, 64, 40, 26, 13, 10, 9, 4, 1

which is arXiv:2305.06642's n = 8 table exactly, with nothing left as a
candidate. Ten tests in `tests/test_checkpointing.py`, ~54 s.

**One thing the first version got wrong.** `stoppedEarly` was computed as "the
deadline has passed by the time the run returns", which is not the same as "a
step broke out". A zero budget at n = 7 expires before the first check and still
leaves a complete classification, because seeding places that whole length
outright -- and the run reported itself as stopped and unfinished, which would
send someone to resume a run with nothing in it. It now means a step actually
broke out of its loop, and is pinned as a test.

**What is still not checkpointed, and deliberately.** `probe.py` holds its
search in memory and prints at the end, so the deep run H-010 asks for -- seven
mutations at clearance 9, where six already took 43 minutes -- either completes
or is lost. `overnight.sh` runs it under a hard timeout for that reason rather
than pretending otherwise. Making it resumable is the obvious next piece of work
if that probe is going to be run repeatedly.

---

## E-027 — The free move, the square walked instead of searched, and the trees
*2026-09-16* · **the free move beats the whole table; two new families; no non-quipu tree reached** → F-028, F-029, F-030, F-031

Four threads, all suggested from outside the search, and the cheapest of them is
the largest result this project has had.

### 1. Relations of two arrows are free (F-028)

`corollary:lengthtworelations` of arXiv:2310.08346 has been sitting in
`research/literature/` unused since that paper was read. Deleting every
two-arrow relation keeps the Coxeter polynomial in all 4861 cases at `n = 3..10`
and keeps the quipu name in all 44320 cases at `n = 3..13`. Added to the orbit
computation it merges 20052 pairs at `n = 12` that the 1794 verified rules of the
table at the time do not, and cuts what a search must still place by 76% at
`n = 9` and 90% at `n = 8`.

The reduced space is exactly the LNAs of one fewer vertex, by shortening every
relation by one arrow — checked as a bijection, not just a count, for
`n = 4..13`.

### 2. Walking the square instead of searching for it

NOTES backlog 26. Open a relation with one mutation — F-027 says that always
gives a 2-by-k square — then mutate at consecutive vertices along the side it
opens. The cost is linear in the relation's length, where a search is exponential
in the depth, so a family's later members cost no more than its first. The
instrument was validated against F-020's lone slide, which it reproduced out to
`d = 9` in 13 seconds; discovery had reached `d = 7` at far greater cost.

**The march must be allowed to repeat its first vertex.** A first version
advanced one vertex per mutation and found nothing at all, because the two known
end families are `[1, 1]` and `[2, 2]`. That is worth recording: a walk that
cannot stand still cannot see either of them.

| planted | where | outcome |
|---|---|---|
| a lone relation, `l = 3..6` | interior | nothing — the square closes only by undoing itself, at every length, not just F-027's `l = 5` |
| a long relation and a two-arrow one, all gaps | interior | nothing but the two-arrow relation sliding past; the long one is a pure spectator |
| two relations of ≥ 3 arrows, `L, m = 3..7`, all gaps | interior | nothing, in 305 configurations |
| the same, all gaps ≥ 2 | either end | nothing |
| an overlapping pair, gap 1 or 2 | interior | **F-030**, the pair-to-triple family |
| a relation at an end | either end | **F-029**, the doubling |

So the square is a real mechanism and a cheap one, and what it finds is
concentrated exactly where F-022 and F-024 said the action was: on relations
that overlap, or against an end. A relation with room around it does nothing,
whatever its length and whatever it is next to.

### 3. The two families, and what they cost the story

F-029 is not expressible as a table rule at all — the honest conclusion is that
the encoding needs widening, which is NOTES backlog 27. F-030 is expressible, and
the table already held its first two members and none of the rest, which is
H-008's shape for the third time.

Between them and the free move, **`A_8` is fully covered with no search**: 21
orbits, nothing left. `A_9` falls from 380 orbits and 222 rows needing a search
to 77 and 37.

**A false start worth recording.** The collapse directions of F-029 were first
written by inverting the doubling's condition and sequence by hand. Both were
wrong: the sequence, because the procedure relabels and a left mutation at a
vertex is not undone by a right mutation there — the collapse is two right
mutations at the *source*; and the condition, which fired on 97 cases where 282
were available. Defining each collapse as "the LNA whose doubling is this one"
fixed both at once, and all four moves then fire 907 times apiece over `n = 7..10`
with no failures.

### 4. What is left at n = 9, now that it is small enough to read

37 rows in 77 orbits, and they have a property in common: **every one has a
relation at the source and a relation at the sink**, and every one is already
reduced. Only 1 of the 37 is certified non-piecewise-hereditary. At `n = 10` the
same two counts are 670 and 660 of 887, so the characterisation is strong there
but not complete. Recorded against H-011, whose mechanism it is the natural limit
of: an LNA with both ends occupied has no free end to walk a run to.

### 5. Trees that are not quipus (F-031)

Every tree of maximum degree three up to order 12, tested for quipu-ness and then
compared by Coxeter polynomial against every LNA of the same length. The first
non-quipu appears at order 10 and is unique — the centre with three neighbours,
each carrying two leaves, which is the tree the question was asked about. Eleven
non-quipu trees over orders 10, 11 and 12; 80444 LNAs compared; not one shared
Coxeter polynomial.

---

## E-026 — The dual as a mirror, and what a rule walks through
*2026-09-16* · **410 duals, all holding; every intermediate a square with a side of two; no shortcut survives** → F-026, F-027, R-011

Four things suggested from outside the search, all cheap, and the first of them
corrects a finding made the same day.

### 1. The proper mirror of a rule

F-025 read an asymmetry between the two ends of the quiver off a comparison
between a rule at the sink and *the same pattern* at the source. That is not the
mirror. The mirror is the relation dual: reverse every arrow **and** exchange
right mutation for left.

| | |
|---|---|
| at the sink | `(0:3) (1:6) -> (0:2) (1:6)` via `[2, 2]` -- 4 confirmations, 0 failures |
| the same pattern at the source | 1 confirmation, **3 failures** |
| its **dual** at the source | `(0:6) (4:3) -> (0:6) (5:2)` via `[-7, -7]` -- 4 confirmations, **0 failures** |

Checked at `(l, m)` = (3,6), (4,7), (3,4), (5,6). The pattern's dual is a pair
sharing an *end*, not a pair sharing a *start*, which is why comparing a pattern
with itself at the other end says nothing. R-011.

### 2. Closing the table under it

Whether the sequence's **order** reverses under the dual was the one thing not
obvious. Over 50 rules sampled from both halves of the table, order **kept**
works for all 50 (8 of them exclusively; the other 42 have sequences symmetric
enough that both work) and order reversed works for none exclusively. So the
order is kept.

Then, over the whole table: 1384 rules, **none self-dual**, **410 duals
missing**, verified at up to four lengths each --

> **410 hold, 0 fail, 0 never apply**, in 106 s on four processes.

The table is generated closed now (364 floating, 1430 anchored). Coverage barely
moves -- 7 rows at n = 10, none at n = 9 -- because those orbits were already
joined another way. The value is that it is free, that it halves what a search
has to look for, and that it is what corrected F-025. F-026.

### 3. What a rule walks through

Every multi-mutation rule in the table run one step at a time, each intermediate
quiver classified:

| | |
|---|---|
| commutative squares | **1364**, sides `2 x k` for k = 2..8 |
| squares with a short side other than 2 | **0** |
| a line again, mid-sequence | 336 |
| anything else | 5 |

So a mutation leaves the line only into a square with a side of exactly two, and
a rule is: open a relation into such a square, do something along its long side,
close it back. Exactly the structure the search has been finding by brute force.

**Tried and it does not immediately give a construction.** Opening the lone
relation of `00500000` in A_10 into its 2-by-4 square and searching four
mutations over the square's vertices finds one way back to a line: `[-3]`, the
undo. A square with nothing to interact with closes only onto itself, which is
F-023's lesson again -- the companion relation is the whole point. F-027.

### 4. Whether mixing directions shortens the rules already known

**Mixed sequences are not an unexplored region.** `localMutationSequences` tries
both signs at every vertex at every step, so every run so far has been free to
find them, and 274 of the table's 1794 rules do mix a left mutation with a right
one -- 82 floating and 192 anchored, almost all of them two or three mutations
long.

**And no rule in the table got shorter.** 150 rules of three mutations or more
were sampled; for each, one LNA it matches, and a search at one mutation fewer
over the vertices of its window. Three came back with a shorter sequence at that
LNA:

| rule | listed | shorter, at one LNA |
|---|---|---|
| `(4:2) -> (0:2)` | `[5, 4, 3, 2]` | `[-2, -7, -8]` |
| `(5:2) -> (0:2)` | `[6, 5, 4, 3, 2]` | `[-2, -8, -9]` |
| `(0:2) -> (6:2)` | `[-3, …, -8]` | `[1, 9, 8]` |

All three are the lone short-relation slide of F-020, which would have been the
interesting place to find a shortcut -- and **none of the three is a rule**.
Stated as floating rewrites they give 1 confirmation and 21 failures apiece;
anchored to the left end, 1 confirmation and 8 failures. They work at the single
LNA the search tried, where the window is flush against both ends of the
shortest quiver its width fits in, and nowhere else. Checked for every d from 1
to 7 with the same answer.

**So F-020's one mutation per arrow travelled stands**, and this is the evidence
against the obvious objection to it. A run that finds a shorter sequence for one
LNA has found nothing until `verifyMove` says otherwise.

**The limits of this, stated so it can be redone properly.** One matching LNA per
rule, and mutations only within one vertex of the window. A shortcut that needs
to reach further out, or that only applies at some positions, would not show up
here.

---

## E-025 — Discovery with the bounds raised, and probes deep and wide
*2026-09-16* · **3045 rules, n = 8 to 98%, and H-010 tested two deeper** → F-024, F-025

E-024's residue said what to do: of the LNAs at n = 9 for which no rule had the
pattern at all, 32 of 33 contained a relation of five arrows or more, and every
run so far had stopped at five arrows in a six-arrow window, where such a
relation cannot sit beside another. So: raise the bounds, and probe the
configurations the findings actually turn on.

### The discovery run

```bash
python discover.py --anchor both --max-arrows 7 --max-width 8 \
    --anchor-lengths 13,14 --jobs 3 --verify-cap 12
```

| stage | |
|---|---|
| patterns x ends x lengths | 412 x 2 x 2 = 1648 searches, 3 mutations each |
| rewrites described | 5807, in 6146 s |
| recurring at both lengths | 3826 |
| verified with no failures | **3045**, in 579 s |
| neither listed nor a floating rule restricted to an end | 2929 |
| changing the orbit partition at lengths 6 to **11** | **728** |

Two hours of search, ten minutes of verification. 1469 of the fresh rules are at
the source and 1467 at the sink, which is the consistency check the relation
dual demands.

**The criterion had to change with the bounds, and this is the trap.** A window
of nine arrows does not fit in A_9 at all, so judging by what changes the
partition at n <= 9 -- which is what E-023 and E-024 did -- *cannot* select a
wide rule however useful it is. Judged that way this run yields 209 rules, all
of window 7 or 8. Judged at lengths 6 to 11 it yields **728**, of which 208 have
a window of nine arrows or more. The earlier curations should be read with that
in mind: they were not wrong at the lengths they measured, but they could not
see past them.

**Coverage, with no mutation search at all:**

| n | LNAs | before this run | after |
|---|---|---|---|
| 7 | 132 | 100% | 100% |
| 8 | 429 | 95% | **98%** -- 10 rows left |
| 9 | 1430 | 73% | **84%** |
| 10 | 4862 | 55% | **63%** |
| 11 | 16796 | 43% | **47%** |

The measurement now reaches n = 10 and n = 11, which it had not before; at 16 s
for n = 10 and about four minutes for n = 11 there was never a reason not to.

### The probes

`probe.py`, written for this run, plants one named pattern and enumerates what
the mutations near it reach, reporting the arrows a run can actually rewrite.

**The frozen pair, deeper and honestly interior (F-024).** The earlier probes
allowed mutations at every vertex of A_13; re-run in A_21 where the ends are out
of reach, `(1:3) (2:3)` reaches 2 LNAs at three mutations, 4 at four, 4 at five
and 6 at six -- against 8, 14, 22 and 36 with an end in reach. None of them
lowers the overlap. Six mutations *with* an end in reach does lower it, to zero,
by walking the relation down to arrow 1; that sequence is not translation
invariant and fails at every shift tried.

**Long pairs, in the interior.** `(1:3) (2:6)`, `(1:5) (2:6)`, `(1:5) (2:7)`,
`(1:3) (2:7)` at four mutations, and the equal pairs `(1:6) (2:6)` and
`(1:7) (2:7)`: every one frozen, and the overlap goes *up* in a third to a half
of what they reach.

**Long pairs at the ends, which is where the new family came from (F-025).** The
same four against each end at three mutations: nothing at all moves at the
source, and at the sink every one loses an arrow off the shorter relation --
`(1:3) (2:6)` and `(1:3) (2:7)` going to overlap 0 outright.
`sinkShortRelationShrinkRules` generates that family and it is verified for
every 3 <= l < m <= 9, 21 members, 4 confirmations apiece, no failures. The
mirror at the source gives 1 confirmation and 3 failures.

**A run of three no longer dissolves when the relations are long.** F-022 had it
that a run of three heavily overlapping relations comes apart where a run of two
does not, on the evidence of `(1:3) (2:3) (3:3)` and two others. At three
mutations `(1:3) (2:6) (3:7)` and `(1:5) (2:6) (3:7)` reach two LNAs each and
neither lowers the overlap. So the dissolution of a run of three is not a
property of the run; it is a property of the *short* runs that were tested, and
what it probably costs is mutations -- F-020's one-per-arrow-travelled again.
Do not quote F-022's run-of-three line without that qualification.

### What to do next

The bounds can go up again -- `--max-arrows 9 --max-width 10`, and the interior
run at the same bounds, which this session did not get to. And the whole
judgement should now be made at lengths 6 to 11 as a matter of course, since it
is affordable and the alternative silently discards every wide rule.

---

## E-024 — Widening the rules to tolerate a bystander
*2026-09-16* · **625 verified, and A_7 needs no search at all** → F-023, H-011

F-023 said what stops a known rule from firing is a relation in its window that
it does not touch. This is that read as a construction rather than a diagnosis:
take each rule in the table, put one untouched relation -- a *spectator* --
somewhere in its window, growing the window by up to three arrows to make room,
and let `verifyMove` decide.

```bash
python discover.py --extend --jobs 4
```

| stage | |
|---|---|
| widenings generated from the 368 rules then in the table | 10609, in 13 s |
| of those, firing on an LNA no search had placed at n <= 9 | **1008** |
| verified with no failures | **625**, in 474 s |
| changing the orbit partition at n <= 9 | **270** -- 125 floating, 145 anchored |

Eight minutes, no search, no mutation run speculatively. The filter is what makes
it affordable and is worth keeping: generate freely, throw away everything that
would not fire on a row still needing a search, then verify.

**Coverage with no search at all.**

| n | theorem | before this batch | after |
|---|---|---|---|
| 6 | 81% | 100% | 100% |
| 7 | 67% | 96% | **100%** |
| 8 | 54% | 81% | **95%** -- 23 rows left of 429 |
| 9 | 43% | 60% | **73%** -- 392 of 1430 |

A classification of A_7 is now a table lookup. At n = 8 twenty-three rows need a
mutation search and at n = 9, 392.

**The number H-011 said to watch stayed small.** Re-running the blocked-rule
diagnostic -- for each unplaced LNA, the rule whose left-hand pattern is present
with the fewest extra relations in its window:

| | n = 8 before | n = 8 after | n = 9 after |
|---|---|---|---|
| blocked by a bystander | 126 | 18 | 359 |
| **no rule has this pattern at all** | 29 | **5** | **33** |

So the mechanical half is still the whole story: what is left is overwhelmingly
more of the same, and another widening pass is the obvious next run.

**And the 33 turn out to say the same thing.** Looked at by hand, they are
almost all a pair of relations one of which is *long*: `(1:3) (2:6)`,
`(1:5) (2:6)`, `(1:5) (2:7)` and the like. 32 of the 33 contain a relation of
five arrows or more and 23 contain one of six or more -- and discovery has never
been given a pattern like that. Both runs used `--max-arrows 5 --max-width 6`,
so a six-arrow relation could not appear beside another at all. This is F-013's
lesson again: *absence of a pattern from the table is evidence about the search,
not about the mathematics*. Raise the bounds before concluding anything about
these.

**What this does not say.** Every one of these rules was verified, but 355 of
the 625 are not listed, and the 270 that are were chosen for changing the
partition at n <= 9. That is a curation against the lengths measured, not a
claim about n >= 10. Re-run the command.

---

## E-023 — Discovery against the ends of the quiver
*2026-09-16* · **630 anchored rules, and coverage at n = 6 becomes complete** → F-023

The first run of `discoverAnchoredMoves`, on the framework F-022 added.

```bash
python discover.py --anchor both --max-arrows 5 --max-width 6 --jobs 4 --verify-cap 12
```

74 patterns x 2 ends x 2 lengths (A_11 and A_12) = 296 searches at three
mutations, margin 3.

| stage | |
|---|---|
| rewrites described | 1344, in 286 s |
| recurring at both lengths | 892 |
| verified with no failures | 724, in 70 s |
| a floating rule restricted to an end | 94 |
| genuinely anchored | **630** -- 315 at each end |
| changing the orbit partition at n <= 9 | **229**, and those are what is listed |

**A first attempt at the same run had to be abandoned**, and why is worth
recording. `verifyMove` enumerated the LNAs by building a path algebra for each,
which at length 12 is 58786 of them and six seconds -- per rule, and there were
886 to check. The enumeration is now cached as relation-length rows
(`nakayama.allRelationLengths`), the algebras built only where a rule actually
matches: 0.3 s instead of 6.5, and the verification of all 886 fell from hours
to 70 seconds. Anything that verifies many rules over the same lengths should go
through that function.

**What the run cost and bought.** Ten minutes end to end. Coverage with no
search at all: 100% at n = 6 (from 83%), 96% at n = 7 (72%), 81% at n = 8 (57%),
60% at n = 9 (45%).

**What is still not reached, and the next question.** 567 rows at n = 9, all of
them heavily overlapping, 306 at overlap 2. The diagnostic that pointed at the
spectators -- for each unplaced LNA, the rule whose left-hand pattern is present
with the fewest extra relations in the window -- says at n = 8 that 126 of 155
are blocked by a bystander and 29 by having no rule at all. Run it again after
the next batch: the count of "no rule has this pattern" is the one to watch,
because it is the part more discovery cannot fix.

---

## E-022 — Whether the relation dual widens the move orbits
*2026-09-16* · **it halves the orbit count and adds no coverage**

The relation dual -- reverse every arrow, renumber -- is one of the three
class-preserving operations of arXiv:2305.06642 and holds for *any* LNA, not
only an almost separate one (`nakayama.relationDual`). It is free, it is not in
the move orbit, and the obvious thought is that adding it would carry rows
across the overlap line for nothing. It does not.

| n | floating | + anchored | + anchored + dual |
|---|---|---|---|
| 7 | 95 covered, 84 orbits | 107, 71 | 107, **45** |
| 8 | 246, 310 | 274, 277 | 274, **156** |
| 9 | 644, 1106 | 726, 1019 | 726, **542** |

The orbit count roughly halves at every length and the covered count does not
move by one row. The reason is structural rather than accidental: the almost
separate condition is itself dual-symmetric, so the dual maps seeded to seeded,
and the rule table already contains the mirror of every rule it contains, so it
maps orbit to orbit. The dual therefore identifies orbits pairwise and never
joins a covered one to an uncovered one.

**Worth knowing, and worth not repeating.** Halving the orbit count is real and
would be worth having if orbits were the expensive object; they are not, the
uncovered rows are. Do not reach for the dual again expecting coverage.

---

## E-021 — What can be done to a heavily overlapping run, in the interior
*2026-09-16* · **the pair is frozen, a run of three is not** → F-022, R-010

The experiment H-003 asked for, aimed where F-021 says to aim it. Each pattern
planted in the middle of A_13 at offset 4 -- four arrows of empty quiver on the
left, six on the right -- with `lnaMoves.localMutationSequences` enumerating
every admissible sequence at vertices within the margin, and the reached LNAs
reported by maximum overlap.

**The isolated pair, at four settings.**

| pattern | start overlap | mutations | margin | reached | any lower |
|---|---|---|---|---|---|
| `(1:3) (2:3)` | 2 | 3 | 3 | 8 | no |
| `(1:3) (2:3)` | 2 | 4 | 3 | 14 | no |
| `(1:3) (2:3)` | 2 | 5 | 3 | 22 | no |
| `(1:3) (2:3)` | 2 | 4 | 6 | 34 | no |
| `(1:4) (2:4)` | 3 | 3 | 3 | 17 | no |
| `(1:5) (2:5)` | 4 | 3 | 3 | 16 | no |
| `(1:3) (2:4)` | 2 | 3 | 3 | 16 | no (3 of them go **up** to 3) |
| `(1:4) (3:3)` | 2 | 3 | 3 | 16 | no (3 go up to 3) |

About a quarter of an hour in total, the depth-5 probe a third of it. Neither depth nor margin
is the dial: doubling the margin at four mutations reaches 34 LNAs instead of
14 and not one of them has a smaller overlap.

**The margin-6 row is stronger than it was written as, and the description was
wrong.** A margin of 6 around a pattern at arrows 5 to 8 of A_13 admits the
vertices 1 to 13 -- *every vertex of the quiver*, both ends included. So that
row is not a probe of the interior at all: it says that from `00003300000`,
**four mutations anywhere in A_13** reach 34 LNAs and none of them has a smaller
overlap. That is a claim about the LNA rather than about locality, and it is the
stronger one. It was recorded here as an interior probe with a wide margin,
which it was not. `probe.py --allow-ends` is how to ask that question on
purpose; without the flag the quiver is lengthened to keep the ends out of
reach, so an interior probe stays one.

**A third relation, and which third relations count.**

| pattern | overlapping run | three mutations |
|---|---|---|
| `(1:3) (2:3) (3:3)` | 3 | down to **0**, via `[6, 5, 6]` |
| `(1:3) (2:4) (3:4)` | 3 | down to **0** |
| `(1:4) (2:4) (4:3)` | 3 | down to **0** |
| `(1:4) (2:4) (3:4)` | 3 | down to 2, from 3 |
| `(1:2) (2:3) (3:3)` | 2 | 31 reached, none lower |
| `(1:3) (2:3) (4:2)` | 2 | 31 reached, none lower |
| `(1:3) (2:3) (5:2)` | 2 | 35 reached, none lower |

The parameter is the length of the run of relations linked by an overlap of two
or more, not the number of relations present: `(1:2) (2:3) (3:3)` has three
relations and is as frozen as the bare pair, because its first shares one arrow
and not two. The rewrites that dissolve a run of three were already in the
table, found at length 8 and in E-011 -- so nothing here is a new rule, and that
is the result. **Do not run a deeper interior search for a rule that pulls an
isolated pair apart**; four probes at three settings of depth and two of margin
say there is none to find, and F-022 says where the pair does come apart.

Reproduce with `lnaMoves.localMutationSequences(13, relLengths, lo, hi, steps,
margin)` on `lnaMoves.embedPattern(13, pattern, 4)`.

---

## E-020 — What the theorem and the move orbits reach, by relation overlap
*2026-09-16* · **the gap is exactly overlap two and above** → F-021

The measurement H-003 has been asking for since 2026-09-13, now that there is a
coordinate to make it in. Every LNA of a length partitioned into orbits under
the verified rules -- applied as rewrites on the relation lengths, with no
mutation computed, which F-017 licenses -- and an orbit called covered when it
contains one the quipu theorem names.

With the 123 floating rules:

| n | LNAs | overlap 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|---|
| 6 | 42 | 16/16 | 18/18 | 1/7 | 0/1 | | | |
| 7 | 132 | 32/32 | 57/57 | 4/33 | 2/9 | 0/1 | | |
| 8 | 429 | 64/64 | 169/169 | 9/132 | 4/52 | 0/11 | 0/1 | |
| 9 | 1430 | 128/128 | 482/482 | 24/484 | 10/247 | 0/75 | 0/13 | 0/1 |

covered over total at each maximum overlap. Two readings, and both matter.
Everything at overlap 0 or 1 is covered, at every length -- which is the almost
separate set exactly, so the theorem's reach is not merely *mostly* the low
overlap rows, it is precisely them. And above the line the table reaches 34 rows
out of 820 at n = 9, none at all past overlap 3.

The leftovers' heavily overlapping runs, commonest first at n = 9: `(1:3) (2:3)`
391 times, `(1:4) (2:4)` 198, `(1:3) (2:4)` and `(1:4) (3:3)` 144 each. By the
longest run in the LNA, 434 of the 786 have nothing longer than a pair.

Seconds per length. `python overlaps.py 6 7 8 9 --cores`, and the same numbers
are pinned in `tests/test_overlap.py`.

---

## E-019 — Two more families, and whether a rule's inverse is free
*2026-09-15* · **two families confirmed, the inverse shortcut refuted** → F-020

Acting on F-020's own lesson rather than raising the search bound.

**The families.** Two candidates read straight off consecutive window widths in
the enlarged table, then verified at three lengths each with `verifyMove`:

| family | d | result |
|---|---|---|
| `(0:2) (2:2)` → `(1:2) (d+2:2)` via `[-3, -5, …]` | 1–6 | 8 confirmations each, no failures |
| `(1:2) (4:2)` → `(0:2) (d+4:2)` via `[2, -7, …]` | 1–5 | 8 confirmations each, no failures |

Discovery had found d = 1 and 2 of each; d = 3 needs four mutations and d = 6
needs seven, so the rest were out of reach of any search run so far. Minutes to
check, against the hours a four-mutation run costs. Both are generated now.

**The inverse shortcut, and it does not hold.** Family A's listed left slide is
exactly its right slide with the sequence reversed, each vertex negated, and each
then moved one step toward zero — which looked like it might be a property of the
window's numbering and so give every rule's inverse for nothing. Applied to all
96 listed rules and verified:

| sequence | inverts | fails | inverse leaves the window |
|---|---|---|---|
| all one direction | 12 | 20 | 28 |
| mixed directions | 0 | 5 | 31 |

So **12 of 96**. The transform works for the slide families because a slide's
sequence is a single uniform run; it is not a general fact about the table, and
the spreading pair -- whose sequence mixes a right mutation with left ones -- is
the counterexample closest to hand. Do not try it again.

---

## E-018 — How far the gentle condition reaches
*2026-09-15* · **one class, the hereditary one** → F-019, R-008

Idea 22's premise, checked before implementing anything.

Every LNA of lengths 5 to 9 tested for `all(arrows in (0, 2))`, then grouped by
the class the classification puts it in:

| n | gentle LNAs | 2^(n-2) | classes containing one |
|---|---|---|---|
| 5 | 8 | 8 | `P^(0)_(0,4)` only |
| 6 | 16 | 16 | `P^(0)_(0,5)` only |
| 7 | 32 | 32 | `P^(0)_(0,6)` only |
| 8 | 64 | 64 | `P^(0)_(0,7)` only |
| 9 | 128 | 128 | `P^(0)_(0,8)` only |

So every gentle LNA is in the class of the path algebra of A_n, and the other 21
classes at n = 9 contain none -- including both members of the cospectral pair,
which have 18 members each and not one gentle among them. Seconds to run, against
the days an AAG implementation would have taken.

The reason is a theorem, not a coincidence: all relations of length 2 implies
almost separate relations, and operation 2 of `cor:EquivNakayamaAlgebras` drops
such a relation without changing the class, so dropping them all leaves the
hereditary algebra. Pinned for n = 4 to 10 in
`test_every_gentle_lna_is_the_hereditary_one`.

Also established while looking: `WebSearch` reaches the literature from the
session sandbox even though `curl` and `WebFetch` to arxiv.org are refused by the
egress proxy. Enough to find and identify a paper, not to read one.

---

## E-017 — n = 9 on the corrected engine
*2026-09-15* · **22 classes, then 20** → F-018

`classify.py 9` on the engine of F-015 with the gate of F-016, default depths
(`--depth 6 --resolve-depth 6`). About 90 minutes.

The run placed all 1430 rows and left nothing a candidate, but at **22** classes:
18 quipus plus `C(2,3,5)` (46 LNAs), `C(2,2,6)` (13), `C(2,4,4)` (8) and one not
piecewise hereditary. It reported three groups "proved distinct despite sharing a
Coxeter polynomial", two of which were `C(2,3,5)` against `P^(5)_(1,2)` and
`C(2,2,6)` against `P^(1,1)_(1,3,1)`.

**Diagnosis, in this order.**

1. The `C(...)` names come from `canonicalWeightType`, which reads them off the
   class' Coxeter polynomial — so they cannot separate two classes that share
   one. Circular.
2. A direct probe: iterative deepening from every member of each of the two
   classes, and from every member's relation dual, reporting every foreign class
   reached. Both merged **at depth 2**, from the first member tried and from its
   dual as well. Seconds, against the 90 minutes of the run.
3. (2,3,5) and (2,2,6) are domestic weight types, and their extended Dynkin trees
   are `P^(5)_(1,2)` and `P^(1,1)_(1,3,1)` — the two classes they were separated
   from. Computed, not asserted.

**Re-run of the post-search half only** (the rows were sound; only the merge step
was wrong), on the same table: both merged at the first depth tried, then
`3033030` certified not piecewise hereditary and `3345000` named tubular
`C(2,4,4)`. **1430 LNAs, 20 classes, 0 candidates, 1 separated** — the separated
group being the cospectral pair of F-010, exactly F-011. Under a minute.

**Then re-run whole, from scratch, on the fixed pipeline: the same 20**, with the
same sizes class for class. The naming order is visible in the log -- the theorem
names 18 classes, `resolveMergeCandidates` then merges 11 away (`2233030` into
`P^(1,1)_(1,3,1)` and `2334400` into `P^(5)_(1,2)` among them, both at depth 6),
and only then do the fallbacks name what is left: `3033030` not piecewise
hereditary by Proposition A9, `3345000` the tubular `C(2,4,4)`. One separated
group, the cospectral pair of F-010. About 35 minutes, 10 of them the search.

---

## E-016 — Are the move rules local?
*2026-09-15* · **yes, both halves** → F-017

H-009's own caveat, checked before anything else was built on it.

1. **Applicability.** `matchesAt` against a predicate reading only the window's
   cells plus one bit (a relation covering the window's first arrow having
   started earlier), over every rule in `VERIFIED_MOVES` x every admissible LNA
   x every window position:

   | n | comparisons | matches | disagreements |
   |---|---|---|---|
   | 5, 6, 7 | 67,712 | 150 | 0 |
   | 8, 9 | 924,352 | 1084 | 0 |

2. **Legality.** The whole table re-verified where each rule fits: **1218
   confirmations, zero failures** -- every match is an admissible sequence
   landing on the predicted LNA with the Coxeter polynomial kept.

Cheap: seconds for lengths 5 to 7, a couple of minutes for 8 and 9, and about
four minutes for the legality half. Rows 1 and 2 are tests now
(`test_whether_a_move_applies_is_a_local_condition`,
`test_each_rule_holds_wherever_it_applies`), so **do not repeat them by hand**.

The result that was not the question: the state has to be the **arrow** row, not
the vertex row -- see F-017. Anyone starting the CA literature sweep should start
there rather than from `relLengths`.

---

## E-015 — Every mutation the loosened gate newly allows
*2026-09-15* · **280 of them, all Coxeter-preserving** → F-016

Before switching the search's gate from the strict reading to the paper's
criterion, every mutation the switch would newly allow was enumerated and
checked. Walking out of every LNA of the length, at every vertex of every quiver
reached, comparing the old gate against the new one and computing the Coxeter
polynomial wherever they disagreed:

| n | depth | allowed by both | newly allowed | Coxeter moved | new gate narrower |
|---|---|---|---|---|---|
| 5 | 3 | 304 | 22 | 0 | 0 |
| 6 | 3 | 1450 | 138 | 0 | 0 |
| 7 | 2 | 1938 | 120 | 0 | 0 |

The last column matters as much as the others: a criterion that was *narrower*
anywhere would have meant the switch loses a mutation the published runs used,
and it never is.

Then the classifications, which are the acceptance test: n = 6, 7 and 8 all give
the same classes with the same sizes as before, and n = 7 dropped from 38
seconds to 20.

The old criterion is kept in `tests/test_procedure.py` as `strictlyMutable`,
which is what makes the comparison re-runnable; the first two rows are a test
now. **Do not repeat the n = 7 row** — about four minutes, and it says the same
thing as the other two.

---

## E-014 — The procedure on coefficients, against the one it replaced
*2026-09-14* · **agreement everywhere but two cases, which are R-007** → F-015

Five runs, all gated on `mutationIsPossibleAtVertex` so both implementations walk
the same mutations:

1. **One mutation, every admissible vertex, every LNA of n = 4..8.** 45 + 126 +
   462 + 1716 = 2349 comparisons, **zero** differences. Minutes.
2. **Depth-3 walks, n = 5 and 6.** 1446 and 7496 step comparisons, **zero**
   differences.
3. **Depth-3 walks, n = 7.** 37470 step comparisons, **2** differences, both
   after three mutations, both a relation the old implementation did not
   produce. These are the whole of R-007.
4. **The exact cleanup on the old steps' output**, n = 5 and 6 at depth 3 and
   n = 7 at depth 2: 1446 + 7496 + 4710 = 13652 comparisons, **zero**
   differences. Worth having separately, because it says the disagreement is in
   step 7 and not in the cleanup.
5. **Coefficients against the guess**, over every quiver within depth 3 of every
   LNA of n = 5 and 6: 1239 Cartan matrices, **zero** differences.

**A mixed engine is not an option, and this is how that was learned.** Running
the old steps 1-7 with the exact cleanup passed run 4 above and then reached
*two* different hereditary forms from `A_6` `3030` at depth 7 — a degree-4 tree
alongside `P^(1,1)_(1,0,1)`, which cannot both be one class. The exact cleanup
expects the relations step 7 produces; with step 7's output missing a relation it
cuts the wrong generators. Use one engine or the other, whole.

**Timings**, n = 7 over the 462 admissible single mutations: procedure 0.26 s
against 0.94 s, admissibility 0.14 s against 1.74 s. The exact versions are
3.6x and 12x *faster*.

**Do not repeat runs 1, 2, 4 and 5** — they are `tests/test_procedure.py` now.
Run 3 at n = 7 depth 3 takes about eight minutes and is worth re-running only if
step 7 changes.

---

## E-013 — Audit of the quipu symmetry, after R-006 was challenged
*2026-09-14* · **no defect found** → F-014

Four runs, in increasing cost:

1. **Canonicalisation against `networkx.is_isomorphic`**, over every quipu
   parameter pair of orders 3–11 (12 names at order 3 up to 28656 at order 11).
   Same canonical parameters iff isomorphic graphs, both directions, zero
   exceptions. Seconds.
2. **The paper's class-preserving operations against the quipu fibres**, over
   every LNA of lengths 4–10 with almost separate relations and no length-2
   relation. Orbits equal fibres exactly at every length; largest orbit 8, the
   paper's bound. Seconds. Now `tests/test_quipu_symmetry.py`.
3. **Tree enumeration against parameter enumeration.** The old
   `generateAllQuipus` (enumerate non-isomorphic trees, test the degrees) and
   `quipuForms.allQuipusOfOrder` (enumerate the P^(m)_(k) parameters,
   canonicalise) agree on the counts for orders 4–12, and no tree the first
   accepts is rejected by `quipuForms.isQuipu`. The old function only tested the
   "degree-3 vertices lie on one path" condition when there were more than three
   of them, which looked like a hole — but three branch vertices in a tree of
   maximum degree 3 always do lie on one path, since a path through two of them
   passes through the third, so four is the smallest number that can fail.

   Kept, since it is a genuinely independent route: it is now
   `quipuForms.quipusByTreeEnumeration`, with the degree test written out as
   `isQuipuByDegrees`, and the agreement is a test rather than a note here.
4. **Hereditary form by mutation search**, from all four long-relation members of
   `P^(1,4)_(1,0,1)` (`0003030`, `3030000`, `3060000`, `6000030`) and both of
   `P^(1,2)_(1,1,2)` (`0400030`, `3004000`), at depth 6, plus iterative
   deepening 2–6 from `3060000` and `3004000`. **Nothing reached** — no
   relation-free quiver from any of them. Tens of minutes.

Run 4 is the one that would have been independent of `thm:QuipuToAn`, and it is
simply out of range here, the same way `A_{7,(2,4)}^{(3,3)}` is (see the test
`test_the_theorem_answers_where_the_search_gives_up`). **Do not repeat it at
depth 6 or less.** Depth 7+ at n = 9 was not attempted and is expected to be
hours; the cheaper route to an independent check is a derived invariant computed
from the algebra, not a deeper search.

---

## E-012 — Pair slide at relation lengths 2 to 7
*2026-09-14* · **confirmed a family**

`lnaMoves.verifyMove` on the pair-slide rewrite for each `l`, both directions,
over lengths `l+3 .. l+6`. 22 confirmations per direction per length, zero
failures throughout. → F-013, confirming H-001.

Seconds to run. Should have been the first thing tried after finding the rule at
`l = 3`.

---

## E-011 — Interior discovery, three mutations
*2026-09-14, concluded 2026-09-15* · **44 rules, and 30 false ones caught** → F-020, R-009

`lnaMoves.discoverLocalMoves`, 26 patterns of up to 3 relations spanning ≤ 5
arrows, planted at offset 4 in A_13 and offset 5 in A_14, `maxSteps=3`,
`margin=3`. Re-run as

    python discover.py --jobs 2

after the original `interior.py` turned out never to have been committed.

**Discovery.** 52 searches, 315 s on two cores. 336 rewrites described, **166
recurring across both embeddings**, 134 of them not already in the table.

**Verification, first attempt — wrong, and instructively so.** All 134 checked at
the fixed lengths 7 to 10, which E-010 had used: 74 passed. But the lengths have
to follow the window, and a window of 9 arrows fits in A_10 at exactly one
position, flush against both ends. Each window-9 rule therefore got one
confirmation from one length.

**Verification, redone per rule at `width + 1 .. width + 4`.** 44 survive; **all
30 window-9 rules fail at length 11**, where the window can sit clear of the
ends — wrong rules, not thin ones. R-009.

| window | rules | lengths checked | confirmations |
|---|---|---|---|
| 5 | 2 | 6, 7, 8, 9 | 22 |
| 6 | 16 | 7, 8, 9, 10 | 22 |
| 7 | 18 | 8, 9 (+11, 12) | 3 (+14, 42) |
| 8 | 8 | 9, 10 (+11, 12) | 3 (+2, 8) |
| 9 | 0 | 10, 11 | **all 30 failed at 11** |

The window-7 and window-8 survivors were then checked at lengths 11 and 12 as
well, since two lengths is the minimum that rules out an end effect and those had
only two: **52 checks, no failures** (14 and 42 confirmations at 11 and 12 for a
window of 7; 2 and 8 for a window of 8). All 44 are in `VERIFIED_MOVES` now,
taking the listed table from 52 rules to 96 and the table with families from 64
to 116.

**And the rule that mattered was not one of the 44.** Among them,
`(0:2) -> (3:2)` via three left mutations, next to E-010's one- and two-mutation
versions, is the third member of a family whose `d`-th member needs `d`
mutations — so discovery at any bounded depth sees only an initial segment of it.
Generating the family instead gives every member: F-020, and H-008 confirmed.
That is the return on this run, more than the 44.

Cost: roughly 30 s per (pattern, embedding) at `maxSteps=3`; the re-verification
is the expensive half, since a window of 8 wants length 12 and its 58786 LNAs.
Tests H-007, and H-007 bit back.

---

## E-010 — Whole-quiver discovery at length 8
*2026-09-13* · **34 rules**

`discoverMoves([8], maxSteps=2)`: 101 candidates, 85 not already known, 34
verified. All of window 6 — the width that first has room to sit clear of both
ends at that length.

Each was confirmed only 3 times at lengths 7–8, which is thin, so all 34 were
**re-verified over lengths 7 to 10**: all survived, 22 confirmations each. Keep
doing this for wide rules found at short lengths.

---

## E-009 — Whole-quiver discovery at length 7, three mutations
*2026-09-13* · **2 rules** — poor yield

`discoverMoves([7], maxSteps=3)`: 61 candidates, 45 new, **2** verified.

The yield is low because `describeLink` only admits a *local* rewrite — one whose
window contains every relation it touches — and at length 7 a three-mutation
sequence usually disturbs the whole quiver, so nothing recurs across positions.
This is the experiment that motivated interior embedding (H-007). **Do not repeat
at this length.**

---

## E-008 — Classification of n = 10
*2026-09-13, updated 2026-09-14* · **unfinished — resume it**

`classifyLength(10)`, several attempts, none yet complete. Furthest reached:
about 1900 of 4862 rows.

**Long runs do not survive.** Three separate causes, all worth knowing:

1. two attempts were killed by over-broad `pkill -f` patterns issued by the
   session itself — a pattern that also matches the shell issuing it kills the
   shell, and anything sharing its process group;
2. one died under `setsid` when the machine went away between sittings;
3. n = 10 takes hours, so any of the above is likely to happen at least once.

**So: resume rather than restart.** The table is written after every class
searched, and `classify.py --resume` continues from the existing CSV:

    python classify.py 10 --resume

Before assuming a long job is still running, check `ps` — a stalled row count
looks the same as a dead process.

---

## E-007 — Certificate propagation by vertex deletion
*2026-09-13* · **0 / 0 / 1 / 24 / 308**

`notPiecewiseHereditaryByDeletion` over every LNA of lengths 4 to 11 → F-012.
The zero below length 9 is the correctness check, not an absence of result.

---

## E-006 — Cospectral quipu enumeration to order 13
*2026-09-13* · **the collision map**

`quipuForms.cospectralQuipuGroups(n)` for n = 4..13, cross-checked against equal
Coxeter polynomials computed through each algebra's Cartan matrix for n = 4..11.
The two agree exactly. → F-010.

Seconds to run, no mutation search involved. `python classify.py <n> --collisions`.

---

## E-005 — Orbit verification under the move table
*2026-09-13* · **clean**

Every LNA of lengths 5 to 9: compute its orbit under the verified moves, apply
each recorded mutation sequence to check it reaches the class it claims, and check
the Coxeter polynomial is constant on the orbit. 1764 orbit members, zero
failures.

Run **before** trusting any change to the rule table — it is what caught R-005.

---

## E-004 — Reduction preserves the Cartan matrix
*2026-09-13* · **clean, after R-003**

Every legal mutation of depth ≤ 3 out of every LNA of lengths 5 to 8: 38095
reductions, zero changes. Exact and heuristic Cartan matrices agree throughout.
→ F-008.

First run gave 9 apparent failures; all were R-003, not the reduction.

---

## E-003 — Exact against heuristic Cartan matrix
*2026-09-13* · **agree everywhere tested**

All 624 LNAs of length ≤ 8, and all 8101 quivers reached by walking every legal
mutation of depth ≤ 3 out of all 188 LNAs of lengths 5 to 7. No disagreement, so
no published Coxeter polynomial moves. The shapes where the two models differ
(F-004) have not turned up in an LNA search.

---

## E-002 — Classifications of n = 5 to 9
*2026-09-12 – 2026-09-13* · **match the published table**

| n | LNAs | classes | time |
|---|---|---|---|
| 6 | 42 | 4 | ~13 s |
| 7 | 132 | 6 | ~37 s |
| 8 | 429 | 11 | ~4 min |
| 9 | 1430 | 20 | ~56 min |

n = 6, 7, 8 match arXiv:2305.06642 exactly. n = 9 → F-011. Lengths 6–8 are pinned
as `slow` tests.

---

## E-001 — Reproducing the published n ≤ 8 classification
*2026-09-12* · **the baseline**

The first cross-check of the restored code against the papers: relation-set counts
against the Catalan numbers, the worked example of arXiv:2112.08129 step by step,
Coxeter polynomials of A_n and D_n, and the class membership of the n ≤ 8 table.
Everything agreed once F-001 was fixed.
