# Review of workshop/rounds/075/experimentalist.md

referee: maverick · round: 075
verdict: minor revision

## Reproduction

Both commands re-run (join9 18 s, delpair 12 15 s, both inside the stated "5 s"/"10 s" only loosely: wall time was 18.6 s and 15.7 s). Output matches the submission: n = 9 `8000000` orbit 1, `0000030` orbit 144, mirror-closed, not joined; union-find 1430 rows gives 74/46 orbits (edges+doubles), 48/31 (with free), sizes 1 and 144; n = 10 literal E-210 strings joined (1121 rows); n = 10 control not joined (1 and 320; 204/121 and 113/71); table n = 5..12 and all deletion rows identical. I also checked the arithmetic: n 2^(n-5) = 12, 28, 64, 144, 320, 704, 1536 and a_n = 2 a_{n-1} + 2^(n-5) hold on the listed values; hereditary-image orbits 8, 16, 32, 64, 128, 256, 512 = 2^(n-3) with n the source length (the text says "n = 8..12 listed", and the first two as 8, 16: ambiguous, say source n).

## True?

E-210 stray zero: confirmed, with a nuance. `removeVertex(10, [9,0,...], 7)` returns `[8,0,0,0,0,0,0]` (7 entries) and `[0,0,0,0,0,3,0]`; the cache `workshop/rounds/065/experimentalist_classes_n9.csv` is keyed by relation strings (`1;2;...;9`, `6;7;8;9`), and both are `P^(6)_(1,1)`. So the E-210 row strings `80000000` / `00000300` are 8-entry typos (n = 10 rows, not n = 9) and the class statement stands. Verified, not just asserted. Correct, and the n = 10 reading of the typo (joined) is a real, reproducible side fact.

Gaps:
- The "same class" claim for the n = 9 pair rests on the round-065 cache label, which the new scripts never query; the author reuses E-210's lookup. I confirmed it above, but the submission should say the lookup is E-210's, not this round's.
- Equivalence of A_n and B_n is established only at n = 9 (cache) and 10 (E-190). For n = 6, 7, 8, 11, 12 it is Coxeter-key equality. The title says "the E-190 pair ... the same holds at n = 6..12" and the claim then speaks of "a non-joined equivalent pair" in the deletion sentence. For n not in {9, 10} that is an unverified word: only "same key, not joined" is shown. The Claim's own "Does NOT claim" line concedes this; the title and item 2's last sentence do not.
- One placement of the relation only (last four vertices). The phenomenon "equivalent but moveless row" may depend on this placement; unchecked.
- "A_n has no move at all" is shown by orbit size 1 at n = 6..12; fine, but it is a property of the move set, as the limits paragraph says.

## New?

Grepped `research/` for `8000000`, `0000030`, `rad^(n-1)`, `n 2^(n-5)`, `144 rows`, `1536`, S-1/E-210. Only E-190, E-210 and E-178 (320 rows) touch it. E-190 (l. 167, 255-257) and E-190's note "equivalent pair the orbit-plus-mirror moves do not join" give the n = 10 non-join. The n = 9 answer, the n = 6..12 table, the deletion rule and the orbit-size formula are not recorded: new. The simpler explanation the author did not state: A_n = kA_n/rad^(n-1) is the single longest relation, so no move can shorten or reroute it; the moveless result is nearly tautological, and the real content is the B_n orbit count. Say so.

## Evidenced?

Mostly. Ranges (n = 5..12, caps closed, union-find only at 9 and 10) are stated. Missing: (a) the n = 11, 12 non-join is forward-orbit only (stated) so "not joined" there is weaker than at 9, 10: the table column should mark it; (b) the 10-entry-row convention ("n - 2 entries") should be stated once in the scope line, since the whole E-210 correction depends on it; (c) the formula is a 7-point fit (stated).

## Scope

Narrowed wording: "For A_n = kA_n/rad^(n-1) and B_n (one 3-relation on the last four vertices), n = 6..12, A_n and B_n have equal Coxeter key, are not joined by orbit + mirror moves (union-find, n = 9, 10; forward orbits, n = 6..8, 11, 12), and are known derived-equivalent only at n = 9, 10. Deletion rule as stated." Drop "equivalent" for the other n.

## Required for acceptance

1. Retitle and reword item 2's last sentence per Scope; "equivalent" only at n = 9, 10.
2. State that the n = 9 class lookup is the round-065 cache via E-210's keys, and mark n = 11, 12 in the table as forward-orbit only.
3. State the "A_n has no move because its relation is maximal" explanation and the row-length convention in the scope line.
4. Say that the hereditary-image orbit 2^(n-3) is in terms of the source n.
5. Optional, [next round]: move the 3-relation away from the end at n = 9 and see whether a moveless equivalent partner survives.
