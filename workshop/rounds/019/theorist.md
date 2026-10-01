# The one in-orbit placement of a split four-letter word is the one whose last relation ends at the sink (gap 0) or one before it (gap 1); lemma R alone never reaches 333@0, it stops at a sink-touching shape

author: theorist · round: 019 · kind: result
thread: T1/T2 · bears on: E-088, E-091, H-021

## Claim

Let S be the closed `444` orbit (reduced move set, 2-arrow relations deleted). For the split words of E-091 (letters <= 9, a 4, >= 4 placements,
some but not all placements in S; 6, 9, 12, 15, 18 words at n = 12, 13, 14, 15, 16) define the right gap g = n - (end vertex of the last relation).
**(1) Rule.** The in-S placement of a word has a g that depends only on the word, not on n (n = 12..16, 60 word-n pairs, 0 exceptions):
g = 0 for `2224`, `2334`, `4556 4667 4778 4889` (the `4 b b b+1` words); g = 1 for `224x` (x = 5..9) and `344x` (x = 4..9); `3344` is the one
exception: its in-S placement is offset 1 (left gap 1; mirror of `44033` at right gap 1, so it is the g = 1 type seen from the left).
So "last or second-to-last" is exactly g = 0 (6 words at n = 16) or g = 1 (11 words), plus `3344`: 17 of 18 as E-091 has it, and now a rule that
holds at every n (offset = o_max - g, o_max = n - 4 - last letter), not a count at n = 16.
**(2) Does it reduce by R to `333@0`?** No, for all 6/9/12/15/18 words at n = 12, 13, 14, 16 (0 of 45): R followed to its end (each step checked to be a
one-step double-mutation neighbour in `doubleMutation.rewritesOf`) stops at a terminal that is itself in S but is not `333@0`: `2224 -> 4`,
`224x -> 4x`, `344x -> 4(x+1)` (`3444 -> 45`), `2334 -> 334`, `4 b b d -> 3 b (d+1)` (`4556 -> 357`, `4667 -> 368`, `4778 -> 379`, `4889 -> 38(10)`),
and `3344` is R-inert. The terminals are the right-boundary shapes. The rest of the way to `333@0` is not R: a BFS in the full reduced move set
(`theorist_split.py`) goes through the `34@k <-> 403@k-1` shuttle (the `4@0 -> 3@1` rule of E-080 slid along the line) in every case, 7..13 steps at n = 12, 14.
**(3) Why one placement.** Short shapes are in S at the boundary only (table below, n = 12..16): `4` at left gap 0 or right gap 0; `4y` (y = 5..8) at left gap 0 or
right gap 1; `334`, `335`, `357` at right gap 0; `44` anywhere. R keeps the end of the last relation (the new last interval `(s+2, s+3+d)` has the old end), so a word
inherits the right gap of its terminal and has exactly one placement with that gap. The non-S placements are interior and sit in the other drift classes
(E-065) or in shadow orbits (j = -1 in `theorist_jlabel_n12.txt`).

Not claimed: a proof. The reason that `4y`, `334`, `357` are in S only at the stated slots is read from tables (n = 12..16, y <= 8), not derived; the step
"R keeps g" is checked by the end formula and by the computed terminals, not for all words. The 5- and 6-letter versions, `a = 2` R-steps and letters >= 10 are not covered.
R is only a one-step neighbour when the shortened interval does not swallow a relation on its left: my first run applied R to `3344` and returned `4@3`, a
row not in S, which is not a neighbour; `rchain` now filters every R-step against `rewritesOf` (E-088 tested isolated runs only).

## Evidence

Right gap g of the in-S placement (from `theorist_gaps_n{12..16}.txt`; all n give the same g per word, present words only):

| word | g | R-terminal (n = 16) |
|---|---|---|
| 2224 | 0 | `4@11` |
| 2334 | 0 | `334@9` |
| 4556 4667 4778 4889 | 0 | `357 368 379 38(10)` |
| 2245..2249 | 1 | `45..49` |
| 3444..3449 | 1 | `45..4(10)` |
| 3344 | left gap 1 | R-inert (`3344@1`) |

Shapes in S (offset, right gap g), n = 12 / 16: `4`: (0,7),(7,0) / (0,11),(11,0); `45`: (0,5),(4,1) / (0,9),(8,1); `334`: (5,0) / (9,0); `357`: (2,0) / (6,0); `44`: all 7 / 11 placements.
`3335 3345 3346 3347` have no placement in S at n = 12..16. Files: `theorist_shapes.txt`, `theorist_rchain_n{12,13,14,16}.txt`, `theorist_split_n{12,14}.txt`.
Class labels (333 offsets held, orbit size) for all placements at n = 12: in-S ones read `0/1410`, others `-1` (shadow, 37..1766) or `2,3,1` (the `235` orbits):
`theorist_jlabel_n12.txt`.

## Reproduction

```
for n in 12 13 14 15 16; do timeout 10m .venv/bin/python workshop/rounds/019/theorist_gaps.py $n > workshop/rounds/019/theorist_gaps_n$n.txt; done   # about 15 s each
for n in 12 13 14 16; do timeout 10m .venv/bin/python workshop/rounds/019/theorist_rchain.py $n > workshop/rounds/019/theorist_rchain_n$n.txt; done  # about 15 s each
for n in 12 13 14 16; do timeout 10m .venv/bin/python workshop/rounds/019/theorist_shapes.py $n; done > workshop/rounds/019/theorist_shapes.txt
timeout 10m .venv/bin/python workshop/rounds/019/theorist_split.py 12 > workshop/rounds/019/theorist_split_n12.txt   # also 14; paths via full moves
timeout 10m .venv/bin/python workshop/rounds/019/theorist_jlabel.py 12 > workshop/rounds/019/theorist_jlabel_n12.txt
```

## Prior record

E-091 (the 6..18 split words, one in-S placement each, last or second-to-last "descriptive, no rule"), E-088 (lemma R; `3x@0 -> 33(x-1)@0`), E-080
(`4@0 -> 3@1`). New: the rule in terms of the right gap with n-independent values, the terminal shapes, and the negative answer to "reduces by R to `333@0`".
Not in `research/RETRACTIONS.md` as far as grep shows.

## Code changed

None in the library. New scripts `theorist_{gaps,rchain,shapes,split,jlabel}.py` and outputs in `workshop/rounds/019/`. No tests touched.

## Next

- skeptic: reproduce the g table at n = 17, and test 4-letter words with a 4 and 4+ placements that are outside S entirely (3334, 2455, ...) for their terminal and gap.
- theorist (next): derive "`4y` in S only at left gap 0 / right gap 1" from the `34 <-> 403` shuttle; then the g of an arbitrary word is read off its R-terminal.
- E-091's "last or second-to-last" could be restated as g in {0, 1}; chair's call.
