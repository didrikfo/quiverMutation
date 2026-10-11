# The rewrite's (k,i) Cartan entry is dim coker g_i because it counts paths into k* modulo steps 4 and 6; dually its (i,k) entry is dim ker psi_i, so the congruence fails exactly when some dim ker g_i != 0

author: scholar · round: 021 · kind: proof · thread: E-095/E-097 (tilting <=> Cartan congruence) · bears on: H-015, E-087, E-095, E-097

## Claim

Write k for the mutated vertex, k* for the new vertex, A = kQ/I acyclic, alpha: k -> t(alpha) the arrows out of k,
g_i : e_iAe_k -> (+)_alpha e_iAe_{t alpha}, p |-> (p alpha)_alpha (paths written in walking order, as in the code).
Let B be the algebra presented by steps 1-7 (unreduced; reduction preserves the algebra). Assuming step 7 is complete
(it returns generators of the whole ideal out of k*), for every i != k:

 (a) dim e_i B e_{k*} = dim coker g_i, from steps 4 and 6 alone (step 7 does not enter);
 (b) dim e_{k*} B e_i = (+)_alpha dim e_{t alpha}Ae_i - dim e_kAe_i = dim ker psi_i, psi_i : (+)_alpha e_{t alpha}Ae_i -> e_kAe_i, (q_alpha) |-> sum alpha q_alpha,
     from step 7 (this is the part E-095 called "read off the data").

With R = R^+_k (the code's `rplus`), X = R C R^T has X[k,i] = chi(P_i, T^+_k) = sum_alpha dim e_iAe_{t alpha} - dim e_iAe_k = dim coker g_i - dim ker g_i,
and X[i,k] = (b) exactly (psi_i is onto for i != k, so no ker/coker correction). So with Y the child's Cartan matrix, X - Y = -dim ker g_i at [k,i] and 0 at [i,k]:
E-095's "row k off the diagonal, equal to -dim ker g_i". The congruence holds iff all dim ker g_i = 0, i.e. g injective, i.e. `tiltingPlus`.
Not claimed: diagonal and rows/columns i,j != k (not derived; the scripts show no difference there, E-095), and nothing about whether B is the true End(T).

Second claim (E-097's condition): every rejecting parent met by the guarded walks at n = 5..7 has the shape "v has exactly one outgoing arrow v -> e and
a commutativity relation between two paths from one vertex a ending x, v, e (x distinct)"; the strict A5 square (arrows a>b, a>c, b>v, c>v, v>e) is
only 767 of 1123 at n = 6, so "A5-shaped" must be read as the long-sided square.

## Evidence

Proof of (a). A path i ~> k* in the new quiver ends with a flip alpha*, so is (path p: i ~> t alpha) alpha*, p in the carried/composite arrows.
Ideal elements ending at k* come from relations ending at k*, which are step 4 (for beta: h -> k, sum_alpha (alpha beta) alpha* = 0) and step 6
(a relation into k extended by alpha, hence p in I there); step 5 relations (rbar alpha* = r/alpha) and step 7 relations start at t(alpha) or k*
and contain k* as a source or middle vertex, and B is acyclic at the gate so no path returns to k*. Reading p back in the old quiver (composite = beta then alpha),
step 4 is the image of g_i (it relates (p alpha)_alpha for p ending with beta), step 6 is I. So e_iBe_{k*} = (+) e_iAe_{t alpha} / im g_i.   [QED, uses only steps 4, 6.]

Proof of (b). Paths k* ~> v are rbar_r s, s an old path k_r ~> v (a tail that passes k does so as one composite; a cycle is excluded). The code's map
Phi(eps) = (sum eps_P (r_P/alpha) s_P mod I)_alpha. Since sum_alpha alpha (r/alpha) = r, Phi lands in ker psi_v. Conversely, an element (q_alpha) of ker psi_v is
a path combination sum alpha q_alpha in I_{k v}; I_{kv} is spanned by u r s with r a minimal relation, and if u is nontrivial u = alpha u' contributes u' r s in the
alpha-component, which is already in I (zero in e_{t alpha}Ae_v). Only u trivial survives, i.e. rbar s. So im Phi = ker psi_v. The ideal of B out of k* is exactly
ker Phi (step 7 is by construction its generating set, "forced by nearer" included), hence dim e_{k*}Be_v = rank Phi = dim ker psi_v = sum_alpha dim - dim e_kAe_v
(psi_v onto for v != k: every path k ~> v starts with an arrow). [QED given completeness of step 7; this is the step I am least sure of. It is also where E-087's defect
would have shown: a non-normal-form reduction makes rank Phi too small.]

Check by script (n = 5 closed, n = 6 and 7 class 0 under a 500 s cap, E-080 family). Child = reduced rewrite of gate-admitted step; entry [k,i] vs dim coker g_i,
entry [i,k] vs (b), for all i != k, every step:

| set | steps tilting | steps not | [k,i] = coker | [i,k] = (b) | violations |
|---|---|---|---|---|---|
| E-080 family n = 5..7 | 75 | 6 | all | all | 0 |
| n = 5 class 0 / class 1 (closed) | 16 620 / 13 680 | 0 | all | all | 0 |
| n = 6 class 0 (46 584 algebras, cap) | 99 267 | 1 123 | all | all | 0 |
| n = 7 class 0 (30 420 algebras, cap) | 56 365 | 156 | all | all | 0 |

Shape of rejecting parents (distinct algebras by canonicalKey, at the mutated vertex v):

| set | distinct rejecting parents | strict A5 square + tail | long-square (relation of 2+ paths a ~> x, v, e) | neither |
|---|---|---|---|---|
| E-080 n = 5..7 | 6 | 6 | 6 | 0 |
| n = 6 class 0 | 1 123 | 767 | 1 123 | 0 |
| n = 7 class 0 | 156 | 156 (strict test only) | not run | 0 |

n = 6 example outside the strict shape: relation (3,4,6,1) = (3,5,2,6,1), v = 6, e = 1: two paths 3 -> 4 -> 6 and 3 -> 5 -> 2 -> 6 of unequal length, then 6 -> 1.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/021/scholar_step7_entries.py --e078                  # 2 s
timeout 10m .venv/bin/python workshop/rounds/021/scholar_step7_entries.py 5 --class 0              # 60 s (class 1: 45 s)
timeout 10m .venv/bin/python workshop/rounds/021/scholar_step7_entries.py 6 --class 0 --budget-sec 500   # 500 s
timeout 10m .venv/bin/python workshop/rounds/021/scholar_step7_entries.py 7 --class 0 --budget-sec 500   # 500 s (strict test only; saved output scholar_step7_n7c0.txt)
```
Counts at n = 6, 7 depend on the cap and load (E-097 got 907 / 143; here 1 123 / 156 with a longer walk, run alone).

## Prior record

E-095 states `chi(P_i,T^+_k) = dim coker - dim ker` as a sketch and the rewrite entry "read off the data"; E-097 gives the dim ker histogram and
says A5-shaped "no shape check made". This submission supplies the derivation of the rewrite entry from steps 4, 6 (a) and step 7 (b) and the shape check.
Not in RETRACTIONS. The new content is (b), the dual entry, and that "A5-shaped" is the long-sided square: strict A5 fails for 356 of 1 123 at n = 6.
The one-map identity AI 2.32(b) = Ladkani 2.3(c) = `tiltingPlus` is still the authors' and not derived here.

## Code changed

None in `quivermutation/`. New `workshop/rounds/021/scholar_step7_entries.py` (copy of round 018's cross-check, with the two entries and shapes).

## Next

- theorist: referee (b) -- is "step 7 returns generators of the full ideal out of k*" true, or only up to the "forced by nearer" test? A counterexample would be a (k*,v) entry above dim ker psi_v.
- toolsmith: a cheap assertion in `mutateAtVertex`/`reduce` that dim e_iBe_{k*} = dim coker g_i would check (a) with no Cartan matrix; and a general-shape test (distinct tails) replaces the strict A5 reading in `tests/test_gate_without_tilting.py` wording.
- skeptic: whether the long-square shape (not strict A5) also satisfies E-080's "length-3 square is fine" -- an admitted-and-failing parent with the relation at the short level would break H-015's mechanism.
