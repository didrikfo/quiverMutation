# Experimentalist notebook (rewritten round 107)

## What I now believe
- T7/H-010 (107): no lowering control for a start-gap-2 pair at k <= 4 (margin 1-2). Consecutive-start runs (8,9,10) lower: len 3: 2 -> 0, len 4: 3 -> 2, len 5: 4 -> 3 (k=3). Gap-2 chains (8,10,12), lengths 3-5, reach 0 LNAs (vacuous) or stay at base (1-6 LNAs). So (8:4)(10:4) two-bystander cells stay uninformative; k=6 not run (hours).
- harness `localMutationSequences` has power on runs, none on gap-2 chains; "reached 0" is a cap, not a verdict.
- T6/H-017 (105): defect not constant per key class at n=9 depth 5 (E-283).
- 103 literature: gate part published for gentle LNAs only; Fosse Cor. 4.6 A_n ~ B_n n>=5, orbit(B_n) proper subset.
- Practical: pipe through `tail` buffers all output; write to a file instead. Cost per added step ~2-3x.

## What I tried
- 107: `rounds/107/experimentalist_chain.py` (arbitrary relation list), output `experimentalist_chain_out.txt`.
- 105: `rounds/105/experimentalist_defect_census.py`. 103: `_arxivq/_arxivabs/_gentle_gate/_fosse_check`. 099: `_pairs/_single/_single4`.

## What I would do next
1. T7 overnight: chain (8,4)(10,4)(12,3) and (12,4) at k=5 margin 1; only if a lowering appears, k=6 two-bystander cells.
2. T6 depth-6 defect minima per LNA (overnight shards).
3. Gate check on rows with a length-3 relation (Ladkani 1001.4765); S-1: A_n from B_6 outside the move set.
