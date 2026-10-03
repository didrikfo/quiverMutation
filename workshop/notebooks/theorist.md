# Theorist notebook (rewritten round 031)

## What I believe now
- J_i = ker(P_v -> sum_b P_tb) = H^{-1}(cone)_i = socle copies of S_v in e_iA (round 031, derived; not in research/). Circuit in Gamma_i <=> the left mutation is silting but not tilting. A k-circuit is a socle element supported on k paths.
- So the "no nn / no circuit >= 3" fact on walks cannot come from the shape of g_i: layered hand algebras (i -> a_1..a_m -> v -> t1,t2, partitions of paths, m = 2..4, n = 6..8, scalar 1) realise every circuit shape with J != 0 and no single path in J. The restriction must come from derived equivalence to an LNA.
- Hand examples T3 (circuit 3, n = 7) and D (nn, n = 6) have no LNA key even with 1-2 pendants (0 hits in ~300 extensions); W + pendant does (control works). 900 + 48 non-W layered members at m = 3, 4: 0 hits. Key obstruction is evidence, not a proof, and the layered family is a thin slice (and a weak W control for m >= 3).
- Earlier: E-110 lemma (monomial + differences, scalar 1); walks n = 8 c0/c1 have all Gamma_i components <= 2 edges (E-113); dim e_iAe_v bound by cone estimate E-116 is too weak (d reaches 8). Hom-dimension route is dead for proving the bound.

## What I tried
- 031: `rounds/031/theorist_t3.py` (T3, D, W + pendants), `theorist_layers.py m npend` (exhaustive layered family). Relations on blocks of size >= 3 must be written as chains of 2-term relations (a 3-term list means a sum, caused a bug: J = 0 false positives).
- 029: `theorist_circuit.py`, `theorist_keys.py`. 026: W exact match on 17 802 rows.

## Next
- Find an invariant reason D / T3 keys are not LNA keys (Coxeter polynomial roots / Euler form signature, cf. E-063 signature pos <= n-2); try a statement "socle element with support >= 3 or nn forces Euler form not of LNA type".
- Check Aihara-Iyama 2.31/2.32 for the socle/H^{-1} reading before anyone promotes it.
- Broader hand search: add relations on the a_k side and non-layered sources, scalar != 1 (the 4 coefficient-2 pairs), 3 pendants (needs OVERNIGHT-style time).
- Blind spot: a key miss in a thin family was read as "cannot occur"; do not say that.
