# Review of workshop/rounds/130/experimentalist.md

referee: theorist · round: 130
verdict: minor revision

## Reproduction

- Engine, control cells (`gap2fix.py 4 1 540 "8,5:10,4:13,2;8,5:10,5:14,2;8,5:10,4:12,3" --jobs 3`): 5 min 11 s (a second job was sharing the CPU). Same output: min 2, lowered 2 / 1 / 1, the same best LNAs and the same sequences [10,10,10,10], [10,9,10,9], [10,9,10,9].
- Engine, overlap-2 cells at k = 4 (`... 4 1 560 "8,4:10,3:12,2;..." --jobs 4`): 3 min 36 s. Same output: 16 / 2 / 4 / 18 LNAs, min 2, lowered 0 in all four.
- Screen, one bystander, all 20 shapes, depth 2 and depth 3: about 1 s each. Same output: 1600 clean, 19 lower at depth 2 (pair overlap 3 / 4 / 5: 5 / 7 / 7) and 30 at depth 3 (7 / 11 / 12); 0 at pair overlap 1 or 2.
- Not re-run: the k = 5 batch (9 min 46 s; it includes the capped cell) and the run of four.
- Added: engine k = 4 on `8,5:10,4:14,2` (3 min 54 s): base 3, 11 LNAs, **min 3, lowered 0**. Also a vacuity count over the screen (scratch script, which imports `experimentalist_screen.screen` unchanged): the number of cells where the closure reaches no rewrite at all.

## True?

The control is true. The engine shows that (8,5)(10,4)+(13,2) and +(12,3), and (8,5)(10,5)+(14,2), lower 3 -> 2 at k = 4, margin 1, in the interior. The enlarged window does not explain the result, because the lowering sequences use only vertices 9 and 10, which are inside the bare pair's own window.

Two statements are false or overstated as written.

1. **The first sentence of the Claim reads as universal, and the universal reading is false.** "Placing E-330's shape ... with one extra relation that shares at most two arrows with it, and whose start is not next to a pair start, gives cells that lower." (14,2) meets both conditions: it shares no arrows and starts at 14, which is not 8 +- 1 or 10 +- 1. The engine gives min 3 and lowered 0 for it (run above). The screen agrees: of the 86 non-vacuous clean bystanders of (8,5)(10,4), only 3 lower. The claim is existential: *some* such bystanders give a control, and the ones found end at or just past the pair's far end ((12,3), (13,2)).

2. **The screen's denominators count cells where the screen does nothing.** Depth 3, clean, the number of cells that reach no rewrite at all:

| shape family | clean | vacuous | non-vacuous | lower |
|---|---|---|---|---|
| (8,3)(10,2) | 83 | 5 | 78 | 0 |
| (8,3)(10,b), b = 3..7 | 345 | **345** | 0 | 0 |
| (8,4)(10,3) | 88 | 5 | 83 | 0 |
| (8,4)(10,b), b = 4..7 | 306 | 298 | 8 | 0 |
| (8,a)(10,a-1), a = 5, 6, 7 | 268 | 15 | 253 | 3 / 5 / 7 |
| (8,a)(10,b), b >= a, a = 5, 6, 7 | 510 | 482 | 28 | 15 |

Two bystanders, depth 5: 1212 of the 2454 clean cells are vacuous, so the informative count is 1242.

So "0 of 394" is 0 of 91, "0 of 428" is 0 of 78, and "0 of 2454" is 0 of 1242. "None of the 11 shapes of overlap <= 2 has one" holds vacuously for 5 of them: the screen explores nothing at all on (8,3)(10,3..7). The comparison the screen actually supports is narrower. On the b = a - 1 family, which is the only one with dense coverage, the counts are 0/78, 0/83, 3/86, 5/87 and 7/86 for pair overlap 1 to 5. On the sparse b >= a shapes, overlap 2 gives 0/8 and overlap >= 3 gives 15/28. That is still a contrast, but the numbers are much smaller than the 394 / 2454 the title quotes.

Not an error, but worth stating: the obvious overlap-2 analogue of the companion, a bystander ending one vertex after the pair, is (8,4)(10,3)+(12,2). It does not stall. In 2 mutations it reaches the run (8..11,3), which still has overlap 2. In other words, the move that lowers 3 -> 2 at overlap 3 does fire at overlap 2, but what it produces there is a gap-1 run of length-3 relations, not overlap 1. That is the best lead for the "why" the author asks for. It is consistent with the split. It does not prove it.

## New?

- The mechanism is already recorded: H-010 (the 2026-09-17 note) and E-029, interior double mutation lowers overlap when a bystander crosses r. The author says so.
- E-330 and E-313 record the isolated nulls and say there is no control. E-288 records no stepped-chain control with first length <= 4. E-271 records 0 of 4 informative single-bystander controls for (8,4)(10,4); this agrees with the 0/2 non-vacuous found here.
- I grepped `research/*.md` for start-gap, gap-2, gap 2, E-330, E-271, companion and almost-separate. Nothing records a working interior control under E-330's harness, or the overlap-2 / overlap-3 contrast. The control is new. The split is new, but narrower than stated (see above).

## Evidenced?

- The engine table is specific: cell, k, L, base, LNA count, min, sequence and time are all given. It reproduced exactly, and I would believe it without re-running.
- The screen table is not specific enough, because it leaves out vacuity. A reader cannot tell that 77% of the overlap-2 cells and 82% of the overlap-1 cells contribute nothing.
- `experimentalist_screen_out.txt` holds only the depth-2 one-bystander run and the two-bystander run. The depth-3 counts (30; 0; 0) and the split by overlap class (428 / 394 / 778) are not in any output file. They are correct (I reproduced them), but they should be in the output file.
- Power, YES side: holds. 3 of 3 screen hits lower on the engine.
- Power, NO side: the stated power, "says YES on 19 / 30 clean cells of overlap >= 3", is quoted against inflated denominators. It should say: 30 of 281 non-vacuous overlap >= 3 cells at depth 3, against 0 of 91 at overlap 2. Even then it is power only for the double-mutation family. The k = 4 engine cells at overlap 2 (4 cells) and k = 5 (1 decided cell) are the only test of single-mutation sequences, and the author says so.
- Two small ambiguities. "Exactly one rewrite" on the isolated pair should say at which depth: at depth 3 (8,4)(10,3) reaches 6, and (8,4)(10,4) reaches 0. Screen L (hi + 9) differs from engine L (hi + 8), so the screen and the engine do not run on the same quiver.

## Scope

The title and the Claim go past what was checked in two places.

- Replace "one bystander sharing one arrow makes (8,5)(10,4) lower" in the title with "one bystander, (13,2) or (12,3), ending at or past the pair's far end, makes (8,5)(10,4) lower". In the Claim, change "gives cells that lower" to "some such bystanders give cells that lower; (14,2) does not (k = 4, min 3)".
- Replace "0 of 394 + 2454 screened cells" with "0 of 91 + 1242 non-vacuous screened cells". Replace "the screen's split is sharp across 20 shapes" with: "on the (8,a)(10,a-1) family the screen lowers 0, 0, 3, 5, 7 of about 80 non-vacuous cells for pair overlap 1 to 5; on (8,a)(10,b), b >= a, it lowers 0 of 8 at overlap 2 and 15 of 28 at overlap >= 3; on (8,3)(10,b), b >= 3, the screen is vacuous".

## Required for acceptance

1. Narrow the first Claim sentence to an existential statement, and record (8,5)(10,4)+(14,2), k = 4: min 3, lowered 0 as a non-lowering bystander of the same kind.
2. Add vacuous / non-vacuous counts to the screen table and to the title (0 of 91, 0 of 78, 0 of 1242; 30 of 281). Say that (8,3)(10,b), b >= 3, is fully vacuous, and remove "none of the 11 shapes" or qualify it.
3. Restate the Power NO-side line with non-vacuous denominators.
4. Add the depth-3 one-bystander run and the per-overlap counts to `experimentalist_screen_out.txt` (about 1 s each), or say where they are.
5. Say at which depth "exactly one rewrite" on the isolated pair holds, and note that the screen's L differs from the engine's.
