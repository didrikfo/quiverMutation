# The kernel J_i is H^{-1} of the mutation cone (the right socle of e_iA), so a circuit means "silting but not tilting"; no non-W circuit algebra with an LNA key was found among 962 layered hand algebras (n = 6, 7, 8) plus pendants, but there is no proof

author: theorist · round: 031 · kind: negative (reduction + search, no proof) · thread: T5 · bears on: E-110, E-113, E-116

## Claim

1. **Reduction (derived, short).** For a gate-admitted v (no single path in J) the minimal left add(A/P_v)-approximation of P_v is B = sum over out-arrows b of P_tb, because every radical map out of P_v factors through an arrow. So J = ker(P_v -> B) as a module, and J_i = H^{-1}(C)_i = Hom(T_i, C[-1]) for the cone C = (P_v -> B). Equivalently J_i = {x in e_iAe_v : x * rad A = 0} = Hom(S_v, e_iA): the copies of S_v in the socle of the right projective e_iA. So "Gamma_i has a circuit" means exactly "the left mutation A/P_v + C is silting but not tilting", and a circuit is a socle element x in e_iAe_v whose support in the path basis is the circuit's edge set (k paths for a k-circuit). W is the socle element p1 - p2 with p1 b1 = p2 b1 != 0, both killed by b2.
2. **Consequence for the assignment.** The two-term kernel structure alone does NOT forbid nn or circuits >= 3: for the layered quivers i -> {a_k} -> v -> {t1,t2} with monomial and commutativity relations (scalar 1), every circuit shape occurs with J != 0 and no single path in J (table). Whatever excludes them on walks is the derived-equivalence-to-LNA hypothesis, not the shape of g_i. This is the answer to the T5 request "does the mechanism come from the two-term kernel structure": no, it cannot (negative result, in the family tested).
3. **Hand examples against the key test.** Smallest circuit-3 example T3: vertices i, a1, a2, a3, v, t1, t2 (n = 7), relations a1 v t1 = 0, a3 v t2 = 0, p2 b2 = p1 b2, p2 b1 = p3 b1 (x = p1 - p2 + p3 in J, dim e_iAe_v = 3). Its Coxeter key is in no LNA/dual-LNA key set at n = 7, and no pendant extension (one or two extra vertices, no new relations, not at v) lands on an LNA key at n = 8, 9 (180 extensions tried). The nn algebra D likewise (10 and 130 extensions: 0 hits). Control: the W algebra (n = 6) has no LNA key itself, but W plus one pendant vertex does (6 of 10 one-pendant extensions at n = 7, 30 of 130 with two), so the test can say yes.
4. **Family search.** Layered family L_m (n = m + 4): for each sink t_j an arbitrary set partition of the m paths i a_k v t_j with a zero block. Over all members with J != 0 and no single path in J, with 0 or 1 pendants (m = 4) or 0, 1, 2 (m = 2, 3; the m = 3 runs used the same pendant code) the only members with an LNA key are W-shaped (m = 2, with pendants). Members with a non-W circuit: m = 2: 1 (nn), m = 3: 48 (9 nn, 39 long), m = 4: 900 (93 nn or double nn, 12 W+nn, 795 long = 212 long3 + 583 long4); all 0 hits. Pure-W members at m = 3, 4 (6, 36) also have 0 hits, so for m >= 3 the family says nothing about W either: it is a weak control there.

Not claimed: any theorem that nn or circuits >= 3 cannot occur on LNA-derived algebras; the layered family is a thin slice (one source, one v with two sinks, no relations on the a_k side, no scalar != 1); a key obstruction is not an explanation, and I have no invariant-level reason why the nn key differs from LNA keys while W's does not.

## Evidence

Step 1 weakest point: that B has no summand other than the P_tb (needs minimality and that no non-radical map is used; the library's gate and `kerdim` use exactly this g_i, so the identification J = {x : x * alpha = 0 for all out-arrows alpha} is by definition in `scholar_longsquare.py`; the H^{-1} reading is the standard cone computation and is not tested numerically here).

| family | n | J != 0, no single path in J: W / nn / long | LNA-key hits among non-W |
|---|---|---|---|
| T3, D, W + pendants (hand) | 6-9 | 1 each + 180/130/130 extensions | 0 / 0 (W: 6 of 10, 30 of 130) |
| layered m = 2 (pendants <= 2) | 6 | 2 / 1 / 0 | 0 |
| layered m = 3 (pendants <= 2) | 7 | 6 / 9 / 39 | 0 |
| layered m = 4 (pendants <= 1) | 8 | 36 / 93 (+12 W+nn) / 795 | 0 |


## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/031/theorist_t3.py 2        # 12 s
timeout 10m .venv/bin/python workshop/rounds/031/theorist_layers.py 3 2  # seconds
timeout 10m .venv/bin/python workshop/rounds/031/theorist_layers.py 4 1  # 82 s
```

## Prior record

E-110 defines J and the circuit lemma and says why nn / long circuits do not occur on walks is open; E-113 gives the walk data; E-116 uses the same mutation triangle for the dim bound. `grep -i socle` and `H^{-1}` in `research/` find nothing, so the socle / non-tilting reading of J is not recorded (it is elementary, so check Aihara-Iyama Prop 2.31/2.32 before promoting). D is E-103's example; the T3 example and the pendant test are new, and the family search is new. Nothing in RETRACTIONS bears on it.

## Code changed

None in the library. New: `workshop/rounds/031/theorist_t3.py`, `theorist_layers.py`. No tests run (no library change).

## Next

- Skeptic: refute the Claim 1 reading (J = socle, silting-not-tilting) on one explicit mutation, and extend pendants to the T3 / D examples at 3 extra vertices (needs more time, proposal not run: 130 -> ~1500 extensions per example).
- Theorist (next): an invariant reason; compare Euler forms / Coxeter polynomials of D and T3 with LNAs, aiming at a statement like "nn forces a root of the Coxeter polynomial off every LNA's"; the pending test is whether the polynomial alone already separates.
- Experimentalist: whether on walks Hom(T_i, C[-1]) != 0 forces dim Hom(S_v, e_iA) = 1 (compute dim J_i on all 155+ J != 0 rows; E-110 shows 1 in each Gamma shape).
