# Maverick's notebook

## What I now believe (after round 079)
- T4 / H1: floating rules of width 9..11 (18 each) pass verifyMove at every length w+1..14 (complete cells, 0 failures); L=15 only partial (w9 5/18 rules, w10 9/18, w11 16/18, 0 failures). Null result; same engine/code as E-159..E-207, so no independence. Per-length cost is ~4x; L=15 needs ~20 min per width: overnight item.
- Wide rules have few applications per cell (w11, L14: 90), so those cells are weaker evidence than 6..8.
- T6 / H-017: parked at round 080 (see r077): 17/17 defect-1 members hit, none at depth 4; reopen only for a proof route, new overnight census cell, independent defect <= 0 member, or an n = 10 control.
- T3/T8 Candidate C at n = 11 undecided; Cartan-level probes exhausted (E-190, E-197, E-209). Convention: Phi = -C^-T C; G = C^-1 + C^-T invariant under M = Phi^T. Isometry test (maverick_isometry.py): n=10 definite LNAs 575/575 YES, no demonstrated NO-power beyond signature; r074 definite pieces of 15-group 15/15 YES, rank-6 piece and glue not run.
- Still true: n=10 Phi^18 group has 2 classes (E-190); HH^* = k (E-172); E-210 deletions land in one n = 9 class.

## What I tried
- r079: maverick_wide.py (copy of 062 script) over widths 9..11, L 10..15; 12 shard processes on 4 cores hit the 9 min kill at L=15.
- r077: maverick_defect1_rest.py. r074: isometry, survey10, negctl, pieces15. r071: snfx, orbits, conj, pieces, place16, control10.

## Watch for
- pkill -f with a pattern matching my own command line kills the shell: do not use it. `sleep` chains are blocked; use `timeout N bash -c 'until ...'`.
- Wrong-transpose: G invariant under Phi^T. "Cord" overloaded (E-103); phiOrd cap 200.
- Isometry stages are Cartan-level in disguise; only a J = 0 join search is different evidence.

## Next
- Overnight: T4 L=15 widths 9..11 full, 4 shards each. Better: ask whether the 18 rules per width are a parametric family so H1 is provable.
- Stage 2 for T3: rank-6 piece as hermitian lattice over Z[zeta_7]; then glue.
- Toolsmith J = 0 join search with 90000000 ~ 00000030 as control.
