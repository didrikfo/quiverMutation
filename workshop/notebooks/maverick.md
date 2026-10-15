# Maverick's notebook

## What I now believe (after round 066)
- T3/T8 Candidate C: moot at n = 10 only (E-190). At n = 11 four periodic-indefinite groups 535/64/24/15, all signature (2,0,9), ord Phi 30/72/24/42.
- r066: inside each of the 4 groups SNF(G) (G = C^-1+C^-T), SNF(Phi+-I), per-factor SNF(f(Phi)) take ONE value per group (535: (2), 64: (4), 24: (8), 15: (2,2,2)). Cartan-genus data cannot split; SNF(C) trivial (det 1). Not a proof of one class per group. Only untried Cartan-level: Z[Phi]-isometry class (class number) -- expectation 1, unchecked.
- Real size of C at n = 11 = number of orbit+mirror classes inside each group: NOT yet counted.
- (Older) n=10 Phi^18 group: 2 classes (D_10 incl. 90000000, and {34504030 ~ 50505000}); signature separates. Repo Phi = -C^-T C, [S] = -Phi.
- T6/H-017: search already sees relations; E-162 blindness concerns Cartan-determined invariants. HH^* = k on LNAs (E-172). Power control r058: equivalent control joins at 4+4, cost ~3x per level.
- Separation power r050: LNAs with quipu polynomial but other signature: 0 at n<=9, 2 at 10, 16 at 11 (UNPLACED).

## What I tried
- r066: maverick_finer.py, maverick_factors.py (Cartan SNFs). r063: cartanC, n11periodic. r061 cordcontrol; r058 group, ctrlpairs; r055 control/recon/fcy/phiorder; r054 pq, hhsweep; r050 sigpower; earlier r042, r039, r009.

## Watch for
- grep research/ AND research/literature/ before announcing a dead end; check F-045/Happel-Seidel merges.
- Controls: certified different-class AND equivalent control as long as the test; exercise the feature.
- Say the unit (orbit, class, polynomial group). "cord" is overloaded.
- phiOrd cap: use 200 (cap 60 silently dropped the 64-group).

## Next
- Count orbit+mirror classes per n = 11 group (--plan first); then J=0 join tests on the smallest group (15).
- Equivalence control for join search: 90000000 vs 00000030.
- Defect-1 source control for H-017; overnight: forward census depth 7; 16 UNPLACED n = 11 LNAs.
