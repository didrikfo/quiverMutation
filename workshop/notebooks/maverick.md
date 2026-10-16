# Maverick's notebook

## What I now believe (after round 069)
- T6/H-017: defect-1 source control PASSES (r069, 3345000, L=5: 3/3 hit, 0/3 at L-1; cells (2,3),(2,3),(1,2)). With E-185 (defect 2, L=5,6) the search's positive controls are complete. I proposed closing T6 for the ledger; only the overnight depth-7 census (needs human approval) reopens it. Untested: (5,3,4) cell (3 cords), other LNAs' defect-1 members (3455000 has 3 of (2,3); 4444400 none).
- T3/T8 Candidate C: moot at n = 10 only (E-190). At n = 11 four periodic-indefinite groups 535/64/24/15, all signature (2,0,9), ord Phi 30/72/24/42.
- r066: inside each n=11 group SNF(G), SNF(Phi+-I), per-factor SNF(f(Phi)) take ONE value per group (535:(2), 64:(4), 24:(8), 15:(2,2,2)). Cartan-genus data cannot split; not a proof of one class per group. Untried: Z[Phi]-isometry class number (expect 1).
- Real size of C at n = 11 = orbit+mirror classes per group: NOT yet counted.
- n=10 Phi^18 group: 2 classes (D_10 incl. 90000000; {34504030 ~ 50505000}); signature separates. Repo Phi = -C^-T C, [S] = -Phi.
- HH^* = k on LNAs (E-172). Power control r058: equivalent control joins at 4+4, ~3x per level.
- r050: LNAs with quipu polynomial but other signature: 0 at n<=9, 2 at 10, 16 at 11 (UNPLACED).

## What I tried
- r069: maverick_defect1.py (defect-1 control, 517 s). r066: finer, factors. r063: cartanC, n11periodic. r061 cordcontrol; r058 group, ctrlpairs; r055; r054; r050; earlier.

## Watch for
- grep research/ AND research/literature/ before announcing a dead end; check F-045/Happel-Seidel merges.
- Controls: certified different-class AND equivalent control as long as the test; exercise the feature.
- Say the unit (orbit, class, polynomial group). "cord" is overloaded (E-103).
- phiOrd cap: use 200. Forward walk n=9 depth 5 costs 50-90 s per LNA; each back-search member ~150 s.

## Next
- Count orbit+mirror classes per n = 11 group (--plan first); then J=0 join tests on the smallest group (15).
- Equivalence control for join search: 90000000 vs 00000030.
- Overnight: forward census depth 7; 16 UNPLACED n = 11 LNAs.
