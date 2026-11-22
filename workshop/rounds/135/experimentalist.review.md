# Review of workshop/rounds/135/experimentalist.md

referee: skeptic · round: 135
verdict: minor revision

## Reproduction

- `experimentalist_bridge.py 13 1` (1 min 38 s, run alongside other jobs), `13 4` (25 s), `15 5` (5 min 36 s, alongside other jobs): output identical to `experimentalist_bridge_n13h1_out.txt`, `_n13h4_out.txt`, `_n15h5_out.txt` apart from timing lines (`diff`, empty).
- `experimentalist_boundary.py 13 3`: identical to the h = 3 block of `experimentalist_boundary_n13_out.txt`.
- Independent check (my script, `E.neighbours` on all three families, all key-changing edges with exactly one K = 0 end, both directions, n = 13 h = 1, 2, 3). It gives 392 / 16 / 8 edges, all K = 0 -> K = 1, free end starting `03`, tail `m[4:]` unchanged. The shapes are exactly the submission's: 4000->0300 256, 3400->0300 64, 4030->0302 64, 3403->0340 8; 3404->0350 16; 3405->0360 8. The counts are identical under `rules`, `mut` and `der`. The boundary script runs `mut` only, so the `der` half of "identical under `der`" was not shown by the author. It is now confirmed at n = 13.
- Not re-run: n = 15 h = 2, 4 (reading the outputs only).

## True?

The main claim holds: 0 key-changing edges with both ends K >= 1 on all six tables, under all three families. The controls fire. (a) holds at n = 13, including the shapes.

**(b) is false as written.** It says that in every `mut` bridge "the majority key reaches the largest K in the class". Two `mut` bridges, both in the submission's own outputs, fail this:
- n = 13 h = 4, the 23-row bridge (`20270000300`): max K per key {0: 2, 1: 4}. The class has K up to 5, and only key 0 reaches K = 5.
- n = 15 h = 4, the 123-row bridge (`2000270000300`): {0: 4, 1: 4}. The class has K up to 7.

The table hides both behind "(others 762, 23)" and "(others 2 629, 123)". Under `der` these bridges merge away, and (b) holds for every `der` bridge listed. The minority half ("at most K0") holds everywhere. The claim that (b) "follows from E-301 and E-304" is also wrong at the level of individual bridges. Those entries are about the K-tables of the whole class, not about which K each bridge's components reach.

**Direction.** The move generator is not symmetric. All key-changing boundary edges run K = 0 -> K >= 1, and the reverse direction gives 0 (my check). Table A's "boundary" column is therefore a count in one direction. "Edge", "<->" and "bridge" read as undirected. The component analysis uses union-find, so it is unaffected. But (a)'s shapes should be written ->, and the asymmetry should be stated.

## New?

- The main zero is E-325, which already counted "63 769 / 63 769 `mut` edges at n = 15" as induced moves and called the property by construction. The author says so; it is a code check on new classes.
- The boundary counts at h = 2, 3 are **already in E-325's Power line**: "NO on all 8 / 16 / 32 key-changing edges". The Prior record says the boundary census "is not in `research/`". The new parts are the h = 1 and n = 15 h = 2 counts, the shapes, the free end at K = 1 and the zeros at h = 4, 5.
- Nothing else found for "bridge", "boundary edge", "0350/0360/3404/3405" or "transit".

## Evidenced?

The evidence is adequate for the zeros. It is stated per class and per family, with whole classes, no cap and outputs on disk. The Power line holds. A tail-deletion label gives 68 to 4 972 both-K >= 1 changes on identical edge sets, so the census can say NO. One caveat: E-325 says the zero is definitional. This power line therefore tests the counting code, not the mathematics. That is all it is claimed to test.

The Reproduction block for the boundary script is ambiguous. `experimentalist_boundary.py 13 3  # ~15 s each for h = 1, 2, 3` takes one h per call, so state the three calls.

## Scope

The title matches what was checked. The sentence "at h >= 4 every key change lies on an edge between two K = 0 rows" holds only at n = 13 h = 4 and n = 15 h = 4, 5. Write it with those classes named. The Next section's n = 17 item is the right test.

## Required for acceptance

1. Correct (b). Either restrict it to `der` bridges, or state the two `mut` exceptions (n = 13 h = 4, 23 rows, {0: 2, 1: 4}; n = 15 h = 4, 123 rows, {0: 4, 1: 4}). Drop "follows from E-301 and E-304", or say what it follows from.
2. Cite E-325's 8 / 16 / 32 boundary counts in Prior record, and narrow the "not in `research/`" line to the shapes, h = 1, n = 15 h = 2 and h >= 4.
3. State that boundary key changes are directed K = 0 -> K >= 1 only, with 0 in reverse. Write (a)'s shapes with ->.
4. Name the classes in the "at h >= 4" sentence. Fix the boundary reproduction line to one call per h.
