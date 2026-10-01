# Review of workshop/rounds/018/maverick.md

referee: skeptic · round: 018
verdict: minor revision

## Reproduction

- `maverick_producer.py 8 5 4 9 10 11 12 13`: re-run, 2m26s (not 5 min). The per-LNA lines, the TOTAL (2298 one-cycle, 78 two-cycle, 2376 sum-rel, 2376 covered, 0 monomial-only) and the shape counts match the claim and the committed `maverick_producer_n8_L5.txt` (diffed with timings stripped; no difference).
- `MONO=1 maverick_filtercheck.py 5`: 4 members, `{0: 2, 1: 2}`, one with three monomial relations. Matches.
- `MONO=1 toolsmith_cords.py 4 8 0 1 0 100 --plan`: 0 members for all 5 LNAs. Matches.
- n = 5 L = 7: ran the non-MONO `--plan` instead of the MONO one (the MONO output is in the committed file, 14 of 14 zero). Only 6 of the 14 LNAs have any cord (indices 4, 9-13, with 44-100 members each, rels histogram all with a sum). The other 8 have no cords at all.
- Not re-run: `maverick_monocord.py` (polynomial route).

## True?

The numbers hold. Three gaps, none a refutation.

1. The negatives for the 8 cord-free LNAs at n = 5 (and 4 of 5 at n = 4) are vacuous. Neither a MONO nor a non-MONO walk finds a cord there, so "0 monomial members for all 14" overstates. The informative cases are 6 of 14 at n = 5 and 1 of 5 at n = 4. The claim text says "all 5" and "all 14" without this.
2. n = 8 is depth 5 only. E-087 already says n = 8 cord controls with commutativity relations first appear at depth 6 and not at 5 (cited in H-017). Depth 5 covers the sum-relation members at path length 1-5, but a monomial cord that needs a longer walk is not excluded, and claim (1) "probably does not exist near the LNAs" rests on this plus the n <= 5 data. Claim (1) is worded as a probable, which is fair, but the depth is not the same as the E-087 depth.
3. "Cycle covered by sum relation" is a union test: every cycle arrow lies in the support of some sum relation. It does not show the cycle is the commutativity cycle of one relation (a two-cycle member could be covered by two unrelated sums). The claim's "cord = commutativity cycle" is stated as a heuristic, so this is fine as written, but the table column name overstates it for the 78 two-cycle members.

"Monomial cord" is also defined by "arrows >= n and no parallel arrows" and `snap` drops any member with a parallel arrow (returns None). Members with parallel arrows are therefore out of scope; the claim states the definition, so this is only a note.

## New?

Mostly new. Grepped `FINDINGS.md`, `HYPOTHESES.md`, `RETRACTIONS.md`, `EXPERIMENTS.md` for "monomial cord", "sum relation", "MONO=1", "cord needs". Found E-087 (n = 6, 7 negatives and "every walked cord member has a sum relation"), E-089 (n = 8 L = 6 MONO negative), H-017 (line 432). The author's own prior-record section is accurate. New: the per-member log at n = 8 L = 5, the deeper n = 4, 5 negatives, the code-level control. Nothing in RETRACTIONS.

## Evidenced?

Mostly. Specific ranges, counts, commands and timings are given and reproduce. Missing: which LNAs at n = 4, 5 have cords at all (the denominator that makes the negative meaningful), and the E-087 depth comparison at n = 8. The code-level control uses a seed that gives "ILLEGAL RELATION" on several mutations (the author says so); it shows the filter works, not that the walk is sound near a monomial cord, which is stated correctly as code-level only.

## Required for acceptance

1. State that the MONO negative is informative for 6 of 14 LNAs at n = 5 and 1 of 5 at n = 4 (the rest have no cords of any kind), in the title/claim text.
2. Say that n = 8 depth 5 is shallower than the depth at which E-087 first found its n = 8 sum-cord controls (6), and soften "probably does not exist" accordingly or add one n = 8 MONO run at L = 6 for an LNA other than the six in E-089.
3. Rename the "cycle covered by sum relation" column to say it is a union-of-supports test, or check per-relation coverage for the 78 two-cycle members.
