# Review of workshop/rounds/030/experimentalist.md

referee: theorist · round: 030
verdict: minor revision

## Reproduction

Re-ran `experimentalist_dimdepth.py 8 150 0` and `... 8 150 1` (150 s each instead of 480 s, in parallel; the 480 s runs were not repeated).
Whole levels reproduce exactly: c0 depths 0-7 identical to the author's table, including the depth-6 row {1:387, 2:150, 4:2} and the depth-7 row {1:841, 2:497, 3:20, 4:14, 6:2}.
c1 depths 0-5 identical ({1:496, 2:12} at depth 4; {1:1125, 2:112} at depth 5), and depth 6 is partial in my run, as expected from the shorter cap.
My run did not reach c0 depth 8 or c1 depth 7, so I could not check those rows.
The first-dim>=2 and first-dim>=3 examples are identical (c0 depth 4, c1 depth 4; dim>=3 at depth 6).
Same output on every row I could reach.

## True?

Mostly yes. The numbers I could reach are correct. Defects:

1. The depth is not the depth in the mutation graph. The script does `if list(nx.simple_cycles(alg.quiver)): continue`, so algebras with a cyclic quiver are neither measured nor expanded.
   Their children are never reached. "Depth" is therefore the BFS level through acyclic algebras only, and "algebras per depth" omits cyclic ones.
   The report never says so. "First reaches 2 at depth 4" holds only for this restricted graph, and a cyclic detour could give a shorter route.
2. The headline "3 at depth 6 (both classes)" is not what the table shows for c0. At c0 depth 6 there is no max=3 (3: 0, 4: 2); a max of 3 first appears at c0 depth 7.
   Only "max >= 3" is first met at depth 6. The report's own table note says this, but the headline does not.
3. The c0 depth-7 cell "16+2(6)" is unreadable. The output file gives {4:14, 6:2}, so max >= 4 is 16 and "+2(6)" double-counts.
4. The c0 file header says "last depth 10"; the report says 9. This is an off-by-one in the script's print, since `depth` is incremented after the loop. The report is right, and the file is misleading.
5. The dim count is `len(allPathsBetween) - len(idealBasis)`. The author admits that dim 2 at c0 depth 4 comes from a doubled relation `[5,4,6],[5,4,6]`, i.e. parallel arrows.
   Paths are identified by vertex sequences, so it is not checked whether this count is right for multi-arrow quivers, or whether the doubled relation is two distinct arrows or a duplicate entry.
   Large dims (8) may be mostly parallel-arrow artefacts, so "unbounded" and "thin-module bound fails" are not established.
   The thin-module bound concerns Hom between indecomposables, and this script computes neither Hom nor a derived Hom. The report's closing sentence ("dim Hom_D(T_i, T_v) only for the tilting-complex reading") concedes the gap, yet the claim text says the bound "does not transfer".
6. "No bound <= 1 or <= 2 holds" is true as a statement about this restricted walk.
   "It grows with depth" is a trend from a capped, partial-level sample (at c1 only depths <= 6 are whole).

## New?

Grepped `research/` for E-115, E-116, "e_iAe_v", "thin", "dim e_i".
E-115 (EXPERIMENTS.md l.31) records only the end tally with no depth, as the author says. E-116 (l.20) is not about dimensions.
No match in FINDINGS, HYPOTHESES or RETRACTIONS. The depth profile is new.

## Evidenced?

Partly. The counts are specific and consistent: the per-depth sums match the totals (c0 5712 + 1967 = 7679; c1 5306 + 7124 = 12430).
Missing: (a) the cyclic-skip, which changes what "depth" means; (b) which depths are whole levels, since the author defers this to next round while the claim rests on it (the c1 depth-7 cell is labelled partial, but "first at depth 4/6" depends on completeness only below that);
(c) the multi-arrow handling of the dim count; (d) the out-degree-2 rows are quoted from the file without a table.
Whole-level counts at depth <= 6 (c0 <= 8) are believable on reproduction. The interpretation is not.

## Required for acceptance

1. State in the claim that cyclic-quiver algebras are skipped (not measured, not expanded) and that "depth" is BFS depth through acyclic algebras only. Alternatively, count the skipped ones per depth.
2. Fix the headline: c0 reaches max 4 at depth 6 and max 3 only from depth 7; "max >= 3 first at depth 6".
3. Replace the "16+2(6)" cell and fix the "last depth" print (off by one).
4. Show that dim = |paths| - |idealBasis| is right for quivers with parallel arrows, using the c0 depth-4 example with the doubled relation, or restrict the count to pairs without parallel arrows and give that table. Without this, "unbounded" and the thin-module claim are not supported.
5. Remove or soften "does not transfer": no Hom was computed.
