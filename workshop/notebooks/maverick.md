# Maverick's notebook

## What I now believe (after round 125)
- E-318 pair (n = 11, 72-group, 003033030 / 030330030): X_a := SNF(aC + C^T) = coker(Phi - a) (C unimodular) is a derived invariant via Z-congruence. X_3 and X_-3 separate the pair; X_+-2, SK = X_-1, F[x+1] and F[Phi_18] do not. The F-047 profile separates it only via x - 1 (= S).
- X_3 is NOT a function of the matrix S on unimodular forms (S-fixed plants change it 353/435), but none of the 82 S-and-key-fixed plants changes it. On LNAs it is tied to (key, SNF(S)): 92 key-kept plants and every key group n = 8..12, X_+-3 never splits a class. NOT a general fact: on random unimodular unitriangular forms n = 4..6 X_3 splits 1/22/234 (key, SNF(S), SNF(SK)) classes (r125 referee + maverick_formctl.py), many 2-adically. The separation here is 2-adic (x^2-1)-primary (Phi_18 splits off: Res = 3 odd, x^6+x^3+1 irreducible mod 2).
- X_-3 = S - 4C is congruent to S mod 4, so its separating counts are forced by S. X_3 = 2S + SK is not forced mod 4.
- Every Cartan-level invariant is a function of the Z-congruence class and rests on the same H-024 necessity step: for one pair all Cartan witnesses are logically equivalent ("second computation, not second argument"); "carries nothing beyond SNF(S)" only has content as a partition statement over groups. HH^* = k on all LNAs (E-170/E-172); AG is gentle only (R-008).
- LNA scan: the only (key, SNF(S)) splits by X_a at n <= 12 are X_2 at the prime 3 (n = 10: 2 groups, n = 12: 3 groups), all F-047-type.
- n = 12 group of 47 (r121): 9 census orbits, no merge across the SNF split; the balls are all J = 0 tilting steps (no walk-level test).
- E-269 / H-024 (r119): NO power at one n = 11 group only. S-1 / H-026 (r115): K0 = min(h, n-4-h) at odd n.

## What I tried
- r125: workshop/rounds/125/maverick_skew.py (3 s: pair table, census, congruence and walk controls); maverick_skewscan.py plant | keyplant (RAND=20000, 88 s) | scan N (n = 12: 84 s; 8..11: 22 s); maverick_formctl.py (8 s, form-level control). All flint-based, fast.
- r121 census47 / ball47; r119 plant2/plant3/qsgroups/sing; r115 mistat; earlier r097-r114.

## Watch for
- Venv loses cypari/python-flint: `uv pip install --python .venv/bin/python cypari python-flint`. Use flint (charpoly, factor, snf, det) and never a sympy det at n >= 11 (10 min timeout hit).
- `ct.lnaKeyIndex(n)` values are tuples of ints; ct.lnaCartanMatrix is lower unitriangular.
- toolsmith_jxgroups.json parts are NOT SNF(S) classes for nonsingular groups (they are q_S / quipu splits); only 3 n = 11 key groups split by SNF(S).
- A NO control needs an independent certificate; a guard-off ball that never moves the key is NOT a control.
- A "never splits" result needs a shown-to-fire control; N cross pairs of one partition are one test, not N. Give a command for every number.

## Next
- Lemma, now LNA-only (the general form-level version is false): for LNAs, is coker(Phi - 3) determined by (key, SNF(S))? Would need the Nakayama shape of C. Try n = 13 scan and an LNA-shaped (Kupisch-constrained) form-level search.
- Do not recommend X_3 to toolsmith without a case: on LNAs n <= 12 it never splits an SNF(S) class.
- Categorical witness for E-318 (out of Cartan level): would need non-existence of tilting complexes. Literature first (Ladkani, BGS) -- deferred slot.
- H-032 inside one SNF class: meetingPoints between census orbits of the 47-group and the 72-group (11 orbits, r125 list).
- Overnight if approved: depth 6-7 balls 3040500030, 0040033030; H-017b census.
