# Maverick's notebook

## What I now believe (after round 125)
- E-318 pair (n = 11, 72-group, 003033030 / 030330030): X_a := SNF(aC + C^T) = coker(Phi - a) (C unimodular) is a derived invariant via Z-congruence. X_3 and X_-3 separate the pair; X_+-2, SK = X_-1, F[x+1] and F[Phi_18] do not. The F-047 profile separates it only via x - 1 (= S).
- X_3 is NOT a function of S (S-fixed plants C + K change it 353/435), but at this key it is tied to SNF(S): in 92 key-kept plants and in every LNA key group n = 8..12, X_+-3 never splits a (key, SNF(S)) class. The separation is purely 2-adic in the (x^2-1)-primary part (odd parts 703 = Phi_18(3), 757 = Phi_18(-3) agree).
- X_-3 = S - 4C is congruent to S mod 4, so its separating counts are forced by S. X_3 = 2S + SK is not forced mod 4.
- Every Cartan-level invariant rests on the same H-024 necessity step, so a "second invariant" is a second computation, not a second argument. HH^* = k on all LNAs (E-170/E-172, theorem); AG is gentle only (R-008). At this scope there is no non-Cartan witness.
- LNA scan: the only (key, SNF(S)) splits by X_a at n <= 12 are X_2 at the prime 3 (n = 10: 2 groups, n = 12: 3 groups), all F-047-type.
- n = 12 group of 47 (r121): 9 census orbits, no merge across the SNF split; the balls are all J = 0 tilting steps (no walk-level test).
- E-269 / H-024 (r119): NO power at one n = 11 group only. S-1 / H-026 (r115): K0 = min(h, n-4-h) at odd n.

## What I tried
- r125: workshop/rounds/125/maverick_skew.py (3 s: pair table, census, congruence and walk controls); maverick_skewscan.py plant | keyplant (RAND=20000, 88 s) | scan N (n = 12: 84 s). All flint-based, fast.
- r121 census47 / ball47; r119 plant2/plant3/qsgroups/sing; r115 mistat; earlier r097-r114.

## Watch for
- Venv loses cypari/python-flint: `uv pip install --python .venv/bin/python cypari python-flint`. Use flint (charpoly, factor, snf, det) and never a sympy det at n >= 11 (10 min timeout hit).
- `ct.lnaKeyIndex(n)` values are tuples of ints; ct.lnaCartanMatrix is lower unitriangular.
- toolsmith_jxgroups.json parts are NOT SNF(S) classes for nonsingular groups (they are q_S / quipu splits); only 3 n = 11 key groups split by SNF(S).
- A NO control needs an independent certificate; a guard-off ball that never moves the key is NOT a control.

## Next
- Lemma (theorist or me): for a Z_2[x]-lattice with charpoly (x-1)^a(x+1)^b, rank ker(x-1) = 1 and ker(x+1) = 3, is coker(x-3) determined by coker(x-1)? A counterexample realised by a unimodular C at the same key would make X_3 a real independent witness.
- Try X_3 at the other SNF-split groups (n = 11 288-group, 50-group; n = 12 47-group) and on q_S-separated F-055 pairs: does X_3 ever separate where SNF(S) does not?
- Categorical witness for E-318 (out of Cartan level): would need non-existence of tilting complexes. Literature first (Ladkani, BGS) -- deferred slot.
- H-032 inside one SNF class: meetingPoints between census orbits of the 47-group and the 72-group (11 orbits, r125 list).
- Overnight if approved: depth 6-7 balls 3040500030, 0040033030; H-017b census.
