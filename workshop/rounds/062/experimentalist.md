# Floating H-020 rules of width 6, 7, 8 also hold at the lengths E-165 skipped: 0 failures in 8248 applications

author: experimentalist · round: 062 · kind: result
thread: T4 · bears on: H-020, E-165, E-159
scope: `lnaMoves.VERIFIED_MOVES` floating (unanchored) rules of width 6, 7, 8 (86, 130, 114 rules) at lengths 11 (w=6), 10, 11 (w=7), 11, 12 (w=8); every rule at every one of these lengths, no cap, no sample. Not covered: anchored rules, widths 9..11, lengths above 13.

## Response to referee
Referee (theorist): minor revision. Items 1 (claim wording: lengths w+1..12, 13 complete only for w = 8) and 2 (three `verifyMove` checks named; commit 4f5a46d) done in the Claim by the chair; optional item 3 (zero-application count) not done.

## Claim

Every floating rule of width 6 (86 rules) at length 11, width 7 (130) at lengths 10 and 11, and width 8 (114) at lengths 11 and 12 passes `verifyMove` (predicted LNA, every step admissible, Coxeter polynomial unchanged): 0 failures in 8248 filed applications. Together with E-159 and E-165 each width 6..8 floating rule is now checked at every length from w+1 up to 12 for w = 6, 7 and up to 13 for w = 8 (E-165 for 12, 13; `lengthsToCheck` for the low end); length 13 is complete only for w = 8. `verifyMove` runs its three docstring checks (predicted LNA, every step admissible, Coxeter polynomial unchanged); run at commit 4f5a46d (package code unchanged this round). A failure at any (width, length) below would refute it.

## Evidence

| width | length | rules | applications | failures |
|---|---|---|---|---|
| 6 | 11 | 86 | 3612 | 0 |
| 7 | 10 | 130 | 650 | 0 |
| 7 | 11 | 130 | 1820 | 0 |
| 8 | 11 | 114 | 570 | 0 |
| 8 | 12 | 114 | 1596 | 0 |

Each (width, length) was cut into 4 shards (`START::4` of the rule list), total wall time about 45 s, 8 s, 32 s, 15 s, 75 s per cell on 4 cores. No FAIL lines in any shard log (`grep -c FAIL`).

Note on the tool: `theorist_rulelen.py` filters rules by `L >= w+5`, which silently excludes w=7 at L=10 and w=8 at L=11, 12 (w+3, w+3, w+4); a pass of that script at those lengths would run zero rules. The new script drops that filter. It also cannot be asked for `--plan` in the old script (it reads argv as integers); the new one has `--plan`.

## Reproduction

```
for c in "11 6" "10 7" "11 7" "11 8" "12 8"; do set -- $c; for s in 0 1 2 3; do timeout 10m .venv/bin/python workshop/rounds/062/experimentalist_gaps.py $1 $2 $s 4 > workshop/rounds/062/gaps/g_$1_$2_$s.txt & done; wait; done
```
Add `--plan` to see rule counts only. Total about 3 minutes. Shard logs are in `workshop/rounds/062/gaps/` (small).

## Prior record

E-165 (lengths 12, 13 only; names the gap), E-159 (w <= 5 at w+5..11). This closes the E-165 length gap for floating rules; the only remaining unchecked floating cases are widths 9..11 beyond w+2 and lengths above 13 (E-165 scope line).

## Code changed

New `workshop/rounds/062/experimentalist_gaps.py` only. No tests (no package file touched).

## Next

Chair: E-165's scope line can drop "lengths 11 (w=6), 10-11 (w=7), 11-12 (w=8) unchecked". Experimentalist: the 54 rules of width 9..11 at w+3.. lengths, and the same length gaps for anchored rules, if the anchored table matters.
