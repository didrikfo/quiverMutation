# Maverick's notebook

## What I now believe (after round 063)
- T3/T8 Candidate C is moot at n = 10 ONLY (r063 response: at n = 11 the 638 periodic-indefinite LNAs form 4 key groups 535/64/24/15, all signature (2,0,9), C live; maverick_n11periodic.py, 37 s; do not recommend closing beyond n = 10): the Phi^18 group has Phi^9 = -I, [S] = -Phi, so Cartan level allows only S^9 ~ [b], b even (a = 18 also). Spectrum = kD_10 (h = 18); S^9 = [8] there.
- The group has TWO derived classes: D_10 = {class 0 (320 rows), 90000000 = kA_10/rad^9 (Happel-Seidel, F-045: G pos. def., det 4)} and W = {34504030 ~ 50505000} (G signature (8,2)). The "3 candidates / 2 certified inequivalences" of E-178 is really 2 classes / 1 inequivalence. Signature already separates; C separates nothing.
- Recommended: close C at n = 10 only; T3/T8 stays done-at-catalogue-level (E-066). Sign convention: repo Phi = -C^-T C, [S] = -Phi, so Phi^a = (-1)^(a+b).
- T6/H-017: search already sees relations; E-162 blindness is about Cartan-determined invariants. Refutation-only; proposed park (not close). Cord control r061 (n = 9, 3033030): L = 5 3/3, L = 6 1/1; no defect-1 source yet.
- HH^* = k on LNAs n <= 10 (E-172); "key" = Coxeter polynomial. Power control r058: equivalent control 34504030 ~ 50505000 joins at 4+4; cost ~3x per level.
- Separation power (r050): LNAs with a quipu polynomial but other signature: 0 at n <= 9, 2 at n = 10, 16 at n = 11 (UNPLACED). S-1 lone 3: free-end K threshold n = 11, 13, 15.

## What I tried
- r063: maverick_cartanC.py, maverick_n11periodic.py (Cartan-level C, signatures, D_10 comparison). r061: cordcontrol. r058: group, ctrlpairs. r055: control, recon, fcy, phiorder. r054: pq, hhsweep. r050: sigpower. Earlier: r042, r039, r009.

## Watch for
- grep research/ AND research/literature/ for a closing theorem before announcing a dead end; also check whether a "separate candidate" is already merged by F-045 / Happel-Seidel.
- A control must be certified different-class, an equivalent control as long as the test, and exercise the feature (J != 0, relations, right cord kind).
- Denominators: say the unit (orbit, class, polynomial group). Words overload: "cord".

## Next
- Equivalence control for the join search: 90000000 vs 00000030 (equivalent, different orbit-plus-mirror classes).
- n = 11: count derived classes inside the 4 groups (535/64/24/15); done: C is live in each at the Cartan level.
- Defect-1 source control for H-017 (cells (5,1,2), (5,3,4)); L = 6 on the other 8 LNAs. Overnight: forward census depth 7; 16 UNPLACED n = 11 LNAs.
