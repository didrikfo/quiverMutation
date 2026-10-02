# Vertex deletion does not transport classes by any position rule; stripping a free end of length >= 3 does, at n = 8, 9, 10

author: maverick · round: 027 · kind: result (mostly negative)
thread: new (S-1) · bears on: F-012, H-020, F-011, STEERING S-1 questions 1 and 2

## Claim

Speculation level: tested on small cases (n = 9 complete; n = 8, 10 as a check).
(a) Nothing positional works. For same-class LNA pairs (L, L') at n = 9 the chance that delta_i(L) and delta_j(L') lie in one class at n = 8 is 0.37-0.49 for every (i, j) (0.43 at i = j = middle, 0.49 at i = j = 1), against 0.17-0.19 for all pairs of LNAs. So deletion keeps classes together about twice as often as chance and far from always; no choice of (i, j) gets near 1. Choosing, per LNA, the vertex nearest the middle among those covered by fewest relations lands 78% of LNAs (size-weighted) in the modal image class of their class; a fixed vertex gives 69%.
(b) A coverage bound is trivial: in every n = 9 class (checked with the 16 unresolved cospectral LNAs dropped) some one n = 8 class is reached by some deletion from 99.8% of its members, so "exists a deletion that agrees" says nothing; only a rule says something.
(c) The one rule that works: delete a free vertex at an end (a vertex left of every relation, or right of every relation). Let K be the number of free vertices at that end. If K >= 3, the image class depends only on the source class (all classes, n = 8: 3 of 3, n = 9: 5 of 5, n = 10: 10 of 10, head and tail pooled). K >= 2 fails for 1 class at n = 9 and 2 at n = 10; K >= 1 fails for 8 of 15 classes at n = 9. A free end vertex is any vertex of the end run, and all of them give the same LNA (the word with one zero removed). Deleting a free vertex in a gap between relations also gave one image class (42 of 42 at n = 9, 184 of 184 at n = 10), but those are few cases and I did not test whether the count is by chance.
Does not claim: that the class of the image is *the same core's class* by an independent route (I only show image class is a function of source class); a proof; anything for n >= 11 (not run); that the n = 9 cospectral pair is classified (see Evidence).

## Evidence

Classes: Coxeter key, except a cospectral key is split into its quipu classes by the orbit of free + edge + double moves joined with the mirror (`maverick_classes.py`). n = 8: 11 classes; n = 9: 20 classes plus an unresolved set of 16 LNAs in the cospectral key `P^(1,2)_(1,1,2)` / `P^(1,4)_(1,0,1)` (these are in neither seeded orbit; I dropped them from rate tables, so n = 9 counts above are over 1414 of 1430 LNAs); n = 10: 176 unresolved LNAs in two cospectral keys, also dropped (their image classes at n - 1 carry a '?' label, so a few n = 10 image labels are coarser than reality). Image classes at n - 1 use the same labels.

Pair statistic, ordered same-class pairs (295 298 at n = 9), both deletions always admissible:

| n | (i, j) | same-class pairs | all pairs |
|---|---|---|---|
| 9 | mid, mid | 0.429 | 0.189 |
| 9 | 1, 1 | 0.489 | 0.172 |
| 9 | 1, n | 0.492 | 0.172 |
| 10 | mid, mid | 0.375 | 0.137 |
| 10 | 1, 1 | 0.465 | 0.141 |
| 10 | 1, n | 0.467 | 0.141 |

Rates for a class depend on the class: 1.0 for the 128-member class of the full-length `1^9` word and the 8-member one, 0.36 for the 300-member class, so the average hides a split between rigid and loose classes. In the n = 9 matrix rows 1 and 2 (and 8, 9) are identical, as they must be (deleting vertex 1 or 2 of a free head gives one LNA).

Free end strip, per source class, set of image classes over all members that have a free end with K >= K0 (head and tail pooled; `maverick_endstrip.py`):

| n | K0 = 1 | K0 = 2 | K0 = 3 |
|---|---|---|---|
| 8 | 8 of 10 classes one image | 6 of 6 | 3 of 3 |
| 9 | 7 of 15 | 9 of 10 | 5 of 5 |
| 10 | 21 of 34 | 17 of 19 | 10 of 10 |

The K0 = 2 failures are the 300-member class of n = 9 (key (1,1,0,-1,-2,-2,-1,0,1,1)) and two classes at n = 10, each with two images. At K0 = 3 the map sources -> images is injective at n = 8, 9 and not at n = 10 (10 classes into 9). Per-end breakdown of the same fact (image class given source class and K): K = 1 226/295, K = 2 160/178, K >= 3 all pure (n = 9; n = 10 the same: 739/965, 504/558, then 264/264, 108/108, 45/45, 18/18, 7/7).

Predictable change (question 2), where K <= 2: not examined beyond counting how many image classes there are (2 or 3). I did not find an invariant that moves by a rule. The Coxeter key of the image is of course not a function of the key of the source for a non-free deletion, but I did not test a key-level rule (e.g. image key given source key and the deleted vertex type) beyond the purity above. Question 3 (core with room to move, then deleted) untouched.

## Reproduction

All from the repository root; each under 1 minute (n = 10 about 30 s):
```
for n in 8 9 10; do .venv/bin/python workshop/rounds/027/maverick_classes.py $n; done   # class labels, unresolved count
.venv/bin/python workshop/rounds/027/maverick_delete.py 9      # (i,j) rate matrix per pair, per class
.venv/bin/python workshop/rounds/027/maverick_null.py 9        # same vs all pairs (also 10)
.venv/bin/python workshop/rounds/027/maverick_rules.py 9       # fixed vertex / cover-min rules, coverage bound
.venv/bin/python workshop/rounds/027/maverick_dist.py 9        # image class distributions
.venv/bin/python workshop/rounds/027/maverick_free.py 9        # free / uncovered vertices
.venv/bin/python workshop/rounds/027/maverick_freetype.py 9    # head / tail / gap split (8, 10 too)
.venv/bin/python workshop/rounds/027/maverick_endstrip.py 10   # the K0 table (8, 9 too)
```

## Prior record

F-012 and `piecewiseHereditary.py` state that deletion keeps only piecewise heredity; STEERING S-1 asks the compatibility question. I grepped FINDINGS/HYPOTHESES for head/free-end statements: H-020 (head and tail of a core) is the nearest and concerns where a core can be placed, not deletion; I found no record of the K >= 3 free-end statement. It may be a consequence of H-020 / the free move (F-032) that I did not derive; if the theorist knows a one-line reason, the finding is only a check.

## Code changed

None in `quivermutation/`. New scripts in `workshop/rounds/027/`: `maverick_classes.py` (labels), `maverick_delete.py`, `maverick_null.py`, `maverick_rules.py`, `maverick_dist.py`, `maverick_free.py`, `maverick_freetype.py`, `maverick_endstrip.py`. No tests touched or run (no library change).

## Next

- Resolve the 16 unresolved n = 9 LNAs (3030000 and 15 others in the cospectral key; free + edge + double + mirror orbits do not reach a seed) by the Brusle route (F-045) or mutation search; then the n = 10 176. toolsmith/experimentalist.
- n = 11 test of the K >= 3 statement and of K = 2 failures (classes need the same labelling; `maverick_classes.py` at 11 not run). experimentalist.
- theorist: why a free end of >= 3 vertices is transported (the K = 2 failing class at n = 9 is a good case to look at); is it the free move plus H-020's head?
