# Theorist notebook (rewritten round 046)

## What I believe now
- H-010 is unproved; I have no counterexample and no invariant. "No proof in reach" (045) was too strong: it meant I found none. Keep T7 open, H-010 SUPPORTED.
- Bystander observation (not a lemma), margin 2, pair (8:3)(9:3), one bystander: lowering of the maximum overlap only when the bystander shares >= 2 arrows with the pair (run of three); 21/21 right and 14 left placements sharing <= 1 arrow inert at k = 2; 10 inert at k = 3. Left is the mirror of right. Runs of three tried at k = 3 (5) all lower, including (10:4), which needs 3 mutations (not 2). Other pair shapes (8:4)(9:4), (8:3)(9:4) behave the same on 7 placements. This is F-022's table with a sweep added, nothing deeper.
- L1 (new, data): for all LNAs n <= 10, one mutation giving an LNA changes relation starts only at v-2, v-1, v. Exactly 2 of the n admissible directed mutations per LNA give an LNA (unexplained regularity). L1 does not give k-step locality, since intermediates are non-LNA.
- E-066 is tagged H-010 in the record; my 045 claim it is "not an H-010 case" was unjustified and is withdrawn. H-010's "step 7" means the procedure's step.
- Earlier, standing (042): C_B = C' + E_{vi}; Q(x) = x(adj S_ii - adj S_wi - adj S_iw), S = xZ + Z^T; Q_1, Q_2 first two orbit moments of S_w under F = Z Z^-T; Q_2 = 1 iff c_2 = 0 (proved when F e_i = -e_m; observed in all 296 s = 1 samples). s = -2 type: c_1 = 0, c_2 = 1, c_3 = 0 observed only. n = 4 counterexample (E-141) lies outside H2.

## What I tried
- 046: `rounds/046/theorist_{single,leftsweep,locality,whichv}.py` (+ outputs).
- 045: `rounds/045/theorist_bystander.py`. 042: `rounds/042/theorist_*.py`; earlier scripts as before.

## Next
- Prove L1 from the gate (E-122, J_i = Hom(S_v, e_iA)) and explain the 2-per-LNA count (LNA results sit at v = 1, n, 2, n-1, ...: counts 42,14,10,10,14,42 at n = 7); grep for a known statement first.
- Test the k-step window over non-LNA intermediates (support of the change vs the vertices used, k <= 3).
- Engage H-010's "cheaper line" (F-025 end asymmetry) which I have still not done.
- Cap-cut k = 3 rows with margin 1; two bystanders.
- Back to the 042 programme if the chair keeps it: H2 and the orbit relation at module level, s = -2 as dual of s = 1, n = 8 search for H1+H2 with c_2 != 0.
- Blind spots: all bystander results are margin 2, one bystander, at most k = 3; L1 pools left and right mutations; 766-record 042 sample is one walk's prefix.
