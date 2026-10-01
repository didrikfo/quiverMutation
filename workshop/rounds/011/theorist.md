# The key-coarser lists A, B are the words that alternate between the two single-relation orbits of the key of `35` (even n) or `36` (odd n); this predicts the n = 17, 18 lists except `5046`, `5056`

author: theorist · round: 011 · kind: conjecture (tested prediction, no proof)
thread: T1/T3 · bears on: H-021, E-064, E-070, E-074, E-075

## Claim

(1) For n = 8..20 the key of `35` (even n) / `36` (odd n) is the same at every offset, and among the single-relation rows `k@o`
(one relation of k arrows at offset o) exactly two mirror pairs carry it for n >= 14: `3@2, 3@(n-6)` and `5@0, 5@(n-6)` (even n);
`3@3, 3@(n-7)` and `6@0, 6@(n-7)` (odd n) (n = 12, 13 have a third pair, `9@..`). Their reduced orbits are two small, disjoint,
closed orbits P, Q, each holding its own mirror, with `35`@even and `36`@even offsets in P and the odd offsets in Q: checked for every n = 12..20 (table below).
So the "parity classes" of lists A, B are P and Q: **P and Q are the orbits of `3@2` and `5@0` (`3@3`, `6@0`)**, and a word is key-coarser iff its placements
alternate P, Q, P, Q (or the same for another two-orbit key class with the same property).
(2) Census rule: walk the reduced orbit of every single-relation row; call a catalogue word *predicted* if all its placements lie in
those orbits, alternate between exactly two of them, the two orbits have one Coxeter key, and each holds its own mirror. This reproduces list A
exactly at n = 12, 14, 16, 18 and list B minus `5046 5056` at n = 13, 15, 17 (those two have placements in orbits holding no single relation: sizes 18157 and 8794 at n = 15).
Prediction for n = 17, 18 was the lists B (minus `5046 5056`, which lie outside the method) and A: both came out. It is a necessary-side check only: a key-coarser word at n = 17, 18 outside the
single-relation orbits would be missed, as `5046 5056` are at 13, 15, 17.
(3) The mirror clause matters: `406` alternates between two single-relation orbits with one key (n = 13, 15, 17) but is not in B, because its two orbits are swapped by the mirror
(mirror of `406`@even lies in the orbit of the odd offsets; n = 15 checked), so orbit+mirror merges them and key = orbit+mirror.
(4) Letter 4 (E-075), observation only: the orbit of `4@0` holds `3@1` (and `3@(n-5)`) at every n = 12..16 (sizes 1410, 2386, 3767, 5648, 8134), whereas
`k@0` for k >= 5 holds a `3@j` only when k and n have the same parity (`5@0`, `7@0`, `9@0` hold none at even n; `6@0`, `8@0` hold none at odd n). So 4 is the one letter whose
`k@0` is joined to the 3-orbit at both parities, and the smallest cut-off `k` is 5 (even n) / 6 (odd n), exactly the Q seeds. This explains why 4-words are over-represented in the big orbit, not why `444` merges.
Not claimed: a proof that `35`@o lies in P iff o is even; any n beyond 20 for the orbits; `--max-word 5`; an independent test of `5046 5056`.

## Evidence

Orbits of the seeds (reduced walk, all closed, `theorist_predict_out.txt`; n = 12, 14 printed separately, same form):

| n | P seed: size | Q seed: size | word | P offsets | Q offsets | disjoint, each holds own mirror |
|---|---|---|---|---|---|---|
| 12 | 3@2: 178 | 5@0: 148 | 35 | 0,2,4 | 1,3,5 | yes |
| 13 | 3@3: 449 | 6@0: 224 | 36 | 0,2,4 | 1,3,5 | yes |
| 14 | 3@2: 320 | 5@0: 272 | 35 | 0,2,4,6 | 1,3,5,7 | yes |
| 15 | 3@3: 743 | 6@0: 390 | 36 | 0..6 even | 1..7 odd | yes |
| 16 | 3@2: 516 | 5@0: 446 | 35 | even to 8 | odd to 9 | yes |
| 17 | 3@3: 1123 | 6@0: 614 | 36 | even to 8 | odd to 9 | yes |
| 18 | 3@2: 774 | 5@0: 678 | 35 | even to 10 | odd to 11 | yes |
| 19 | 3@3: 1597 | 6@0: 904 | 36 | even to 10 | odd to 11 | yes |
| 20 | 3@2: 1102 | 5@0: 976 | 35 | even to 12 | odd to 13 | yes |

Key classes (`theorist_triples.py`, keys only, n = 8..20): `35` (`36`) has one key at all its offsets for every n, and for n = 14..20 the single-relation rows with that key are exactly the four above.
Other parity-class words come from other key classes with two single-relation orbits: `3336` at n = 14 lies in orbits of `7@0` (312) and `3@4` (364); `405`, `5004` at n = 15 in orbits of `4@2` (43) and `5@1` (75).
Census (`theorist_census2_out.txt`; first block pre-mirror-condition, second block with it): even n, prediction = list A at 12, 14, 16, 18; odd n, prediction = list B without `5046 5056` at 13, 15, 17 (with the mirror
clause; without it `406` is an extra). Each n takes 3 to 60 s (single-relation orbits only, n = 18: 54 s).
Null results (not repeated): no GF(2)-affine functional of (r_i mod 2, [r_i != 0]) separates P from Q at n = 12 (`theorist_invariant.py`: only the trivial end-entries), none of nine integer
statistics mod 2 or 4 is constant on an orbit (`theorist_stats.py`), and the Smith normal forms of `Phi - 1, Phi + 1, Phi^2 + 1, Phi^2 +- Phi + 1` agree for all offsets of `35` and `3334` (`theorist_snf.py`):
the Z-conjugacy class of the Coxeter matrix does not separate P from Q (tested on those five forms only).
Weakest point: the identification is by computed orbit membership; no move sequence from `35`@0 to `35`@2 is exhibited, and P != Q is a reachability fact of the guarded walk, not shown to be a derived inequivalence.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/011/theorist_predict.py 18      # 6 s; N = 12..20
timeout 10m .venv/bin/python workshop/rounds/011/theorist_census2.py 18       # 54 s; also 12..17
.venv/bin/python workshop/rounds/011/theorist_triples.py 8 20                 # keys, seconds
.venv/bin/python workshop/rounds/011/theorist_k0.py 12 16                     # letter 4
timeout 10m .venv/bin/python workshop/rounds/011/theorist_word.py 15 406 5046 # about 1 min
.venv/bin/python workshop/rounds/011/theorist_mirror406.py 15 406
```

## Prior record

E-074 states the lists and says "why unexplained"; E-064 that the orbits are parity classes each holding its own mirror; E-051 notes `35`, `36` are placed alone at every offset in the reduced walk
(the same words). Nothing in `research/` links the lists to single-relation rows or to the orbits of `3@2`, `5@0` (grep of "single relation", "3@2", "5@0" in `research/`: no hit). E-075's merging
statistics and E-065's `34 -> 44` are consistent with item (4) but do not contain it. No retraction involved.

## Code changed

New scripts only, all `workshop/rounds/011/theorist_*.py`; no library change, no tests run. A `theorist_word.py 17 5046` run was left unfinished (over budget).

## Next

- Theorist: exhibit the move `35@o -> 35@(o+2)` in the guarded walk (orbits have only 178 rows at n = 12: BFS path and the moves used), and why `4@0 -> 3@1` exists but `5@0 -> 3@2` does not.
- Experimentalist: `5046`, `5056` at n = 17 (orbit sizes about 10^5 to 10^6: plan first) and `orbitclass` at 17 only if the census rule is to be tested for words outside the single-relation orbits.
- Skeptic: is "both single-relation orbits hold their own mirror and share a key" ever true for a pair whose word is not key-coarser (other than `406`)? Census at n = 14 found none.
