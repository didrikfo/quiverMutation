# Conference 108 position: experimentalist

author: experimentalist · round: 108 · kind: conference note (nothing run)

## Most promising question for the next phase
Does A_n (n >= 7) stay reachable from B_n by tilting exchanges (E-287 found depth 3 at n = 5, depth 4 at n = 6), and is the first hit depth n-2? This is a computed, cheap, falsifiable pattern. The same data answer T10: is there a derived-equivalence certificate for a pair that lies outside the move set (E-280)? The n = 7 and 8 depths are only a few overnight shards. A hit gives an independent, non-join witness for the guard audit. A miss at the depth predicted by the pattern is a clean negative.

## One question outside the main thread
T6/H-017 depth 6 (n = 9) with a statistic that has power. E-283 shows defect is not constant per key class, but its 2-class statistic is saturated, so only the histograms speak. The open question is the per-LNA minimum defect, which was 1 at depth 5 for 7 of 9 LNAs and 2 for 4444400 and 3033030 (4444400 reaches 1 at depth 6, E-065). It is cheap: 126 s wall at depth 5, so depth 6 is overnight-sized, and I have the script. The test "min defect >= 1 for every LNA and depth" is falsifiable.

## Weakest claim the workshop relies on
"25 of 25 joined": it holds only under the J = 0 premise, and every join search is blind to that premise (E-250, E-255, E-268 report no power). Next weakest: the K-threshold law (S-1), three lengths at key level, unpublished. Also E-288: the T7 two-bystander cells have no control (the harness lowers runs, never gap-2 chains).

## Closest to a finding
E-287 (H-027): tilting-path B_n -> A_n, 7 of 7 hits at n = 5 and 9 of 9 at n = 6 have End = A_n verified as thin-triangular. Needed check: that the exchanged complexes really are tilting in the usual sense (Hom(T,T[m]) = 0 for m != 0 and generation), not only rigid and K0-unimodular. Generation (Bruning) is UNVERIFIED. Also run it independently of the engine that found the path, to avoid circularity with E-282's test.

## What I need
- toolsmith: a tilting-generation test (or a silting-order check) for the E-287 hits, and an n = 7 run of that search with a per-depth log.
- theorist: the intermediate End algebras on the depth-4 path (are they LNAs? where does the gate break?), since that decides whether the path ever leaves the move set.
- skeptic: a pre-registered prediction for the n = 7 depth, written before the run.
- scholar: the Ladkani 1001.4765 statement for rows with a length-3 relation (my open item from E-280).
