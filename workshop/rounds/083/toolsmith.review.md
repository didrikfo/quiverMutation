# Review of workshop/rounds/083/toolsmith.md

referee: skeptic · round: 083
verdict: accept

## Reproduction

Re-ran in a fresh dir (/tmp/sk83), with the author's scripts and unchanged commands.
- `build 3060000`: 178 s. Ball sizes 1, 13, 83, 386, 1584, 6354; 4770 depth-5 parents. Matches.
- `scan` in the three shards (38 s, 35 s, 37 s): the same 10 KEPT indices (1275, 1433, 1746, 1818, 1915, 1932, 1970, 2003, 4551, 4722) with the same paths. Per-shard moved/kept 272/2, 261/6, 97/2, so 640 candidates in total. Matches.
- `toolsmith_d5pkl.py` + `skeptic_join.py side R / side C / meet` for the five joins the author read only from the terminal (1275, 1433, 1746, 1818, 1915): totals 8, 8, 7, 7, 7; all "Hom-verified True meet iso True". Matches the table. The root-side chains differ slightly from the prose (1275 meets via F1 F5 R9 R8 R7, not F4; 1746 via F1 F4 F7 F6 F5), which does not matter: the totals agree.
- Hom output for all 10 is in `toolsmith_d5hom_out.txt` (read, not re-run): tiltingPlus False and totals 2,2,2,2,3,3,4,4,3,3, as stated.
- I did not re-run the Hom tests or the control (`ctl`).

## True?

I found no error. Checked: the counts add up (272+2+261+6+97+2 = 640 flagged; 10/640 = 1.6%; 10/4770 = 0.21%). Hom(T,T[-1]) != 0 for all 10 matches the file. The claim is hedged correctly (upper-bound totals, J = 0 premise and generation by T assumed, "none leaves" would be refuted by a non-joining child).

Not checked by the author, and not required:
- Distinctness of the 10 children. The author says so: 1915 and 4722 share the E-233 child's ball sizes (1, 14, 103, 579), and I confirm 1915's child ball is 1, 14, 103, 579. So "10 steps" does not mean "10 distinct classes", and the text does not claim it.
- Only the S = -1 and S = 0 placements of the one-vertex R step are tested, as in E-233.

## New?

Grep of `research/` for E-233, E-229, E-204, key-kept, perI, depth 5: E-229 (depth <= 4, one kept) and E-233 (that one Hom-tested and joined; "depth-5 parents" left open) are the only prior records. The depth-5 count, 10 kept, is new. Nothing in RETRACTIONS covers it.

## Evidenced?

Yes. The scope line gives root, depth, move set, placements and what is not covered. Each of the 10 has its path, the Hom entries, and the join total and meet type. Small gap: the join outputs for the first five indices are not in `toolsmith_d5join_out.txt` (the author says they were read from the terminal). I re-ran them and they agree.

## Scope

Title and claim match what was checked (n = 9, root 3060000, parents at depth 5, perI-flagged steps only). The remark "E-204's no-J != 0 sentence is now false at depth 5" is fair for perI-flagged steps. Multi-vertex J != 0 steps not flagged by perI are excluded, which the text says.

## Required for acceptance

None. Optional, for the record:
1. Save the first five join outputs into `toolsmith_d5join_out.txt`.
2. Distinctness of the 10 children (up to relabelling) stays open, as the author lists.
