# Theorist notebook (rewritten round 045)

## What I believe now
- T7/H-010 (045): no proof in reach. The statement quantifies over all mutation sequences, whose intermediates are non-LNA algebras; I know no invariant on those that separates LNAs of different overlap (Coxeter data cannot: the derived class contains overlap-0 LNAs, F-022). Locality is provable-looking (gate J_i = Hom(S_v,e_iA), E-122: a right mutation at v is blocked iff a relation of >= 2 arrows ends at the arrow out of v; step 7 touches only paths through v) so depth-k interior orbits are translation covariant, but that does not give all k.
- Step 7 of 2112.08129 (E-066, commutative square) is a rejection example, not an H-010 case. Do not use it for H-010 again.
- Lemma L (data, k <= 3, one bystander, m <= 4): an interior heavy pair is lowered only when a bystander shares >= 2 arrows with it (run of three); bystanders sharing 1 or 0 arrows are inert (`rounds/045/theorist_bystander.py`, output `theorist_bystander_3_2.txt`, partial: 10-min cap).
- Recommended closing T7 with H-010 staying SUPPORTED (search to depth 6, F-024).
- Earlier, still standing (042): C_B = C' + E_{vi}; with H1, H2, Q(x) = x(adj S_ii - adj S_wi - adj S_iw), S = xZ + Z^T; Q_1, Q_2 are first two orbit moments of S_w under F = Z Z^-T; Q_2 = 1 iff c_2 = 0 (proved when F e_i = -e_m, observed in all 296 s = 1 samples); s = -2 type (540 at n = 6) moments c_1 = 0, c_2 = 1, c_3 = 0 observed only. Not a matrix identity of Z: needs realisability. n = 4 counterexample (E-141) lies outside H2.

## What I tried
- 045: `rounds/045/theorist_bystander.py` (pair + one bystander, gap scan, k = 2 full, k = 3 to m = 3, g = 1).
- 042: `rounds/042/theorist_{dump_off,reduce,lemma,moments,orbits,exceptions,cbcheck,n4,zrandom,explore,blocks}.py`; 041, 039, 037, 034, 031 scripts as before.

## Next
- Locality lemma L1: for every LNA n <= 9 and every admissible v, if the mutation yields an LNA then the relation lengths change only within distance 1 of v (and the exact change is the slide). Cheap, and is the half of an H-010 proof that is within reach.
- Finish the k = 3 bystander table with margin 1; add left bystanders and two bystanders.
- Back to the 042 programme if the chair keeps it: prove H2 and the orbit relation at module level (AR theory, eAe at n = 6), explain s = -2 as dual of s = 1, search n = 8 for H1+H2 with c_2 != 0.
- Blind spots: Lemma L is a search result at 3 mutations with margin 2; proof sketch of locality (radius ~ k) is not written out carefully; 766-record sample in 042 is one walk's prefix.
