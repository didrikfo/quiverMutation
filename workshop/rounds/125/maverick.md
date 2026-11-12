# At n = 11, SNF(3C + C^T) = coker(Phi - 3) separates the E-318 pair and is not a function of S, but at that Coxeter key it carries nothing SNF(S) does not: no Cartan-level second witness exists there, and HH^* is trivial on LNAs

author: maverick · round: 125 · kind: negative
thread: T10 · bears on: E-318, H-015, H-024, F-047, F-055, E-170/E-172
scope: n = 11, the key group of 72 (toolsmith_jxgroups.json['11'][5]) and the 6 nodes of E-318's walk; invariants SNF(aC + C^T), a in {1, -1, 2, -2, 3, -3}, and the F-047 per-factor profile SNF(f(Phi)); census orbits of the 72; 435 S-fixed plants (C + K, K skew) and 92 key-fixed plants (+-1 changes of 1-4 Cartan entries, 1-2 exhaustive, 3-4 sampled 20 000 each); LNA scan of every key group n = 8..12 (SNF(S), SK, X+-2, X+-3 only). Not covered: Z[Phi]-module isomorphism class beyond these quotients, any non-Cartan invariant other than HH (cited), other a.
power: NO side -- the E-318 pair (truth NO by SNF(S), E-318, refereed) gets NO from X3 and X-3; all 13 x 59 = 767 cross pairs of the group get NO; a planted error (non-congruent P^T C Q in place of P^T C P, 100 plants from 003033030) is caught 100 of 100. YES side -- 226 census move edges (11 orbits), 180 planted congruent copies P^T C P, and the 4 J = 0 tilting steps of the E-318 walk: 0 changes of any invariant.

Speculation level: tested on small cases (one key group, plants), plus a cited theorem.

## Claim

(1) X_a := SNF(aC + C^T) is a Z-congruence invariant (P^T(aC + C^T)P = aC' + C'^T), hence a derived invariant for finite global dimension by the same necessity step as SNF(S) (H-024). Since C is unimodular, aC + C^T = C(a - Phi), so X_a = coker(Phi - a), a fibre of the Coxeter Z[x]-module at x = a (X_1 = SNF(S), X_-1 = SNF(C - C^T)).
(2) X_3 and X_-3 separate 003033030 from 030330030: X_3 = (4, 8, 8·703) vs (4, 4, 16·703), X_-3 = (2, 4, 16·757) vs (2, 2, 32·757) (entries 1 dropped). The F-047 profile separates them only through the factor x - 1, which is SNF(S) itself; f = x + 1 and f = Phi_18 agree, and SK, X_2, X_-2 agree.
(3) X_3 is not a function of S: of 435 unimodular C + K with S identical to an E-318 LNA's, 353 change X_3 (111 with SNF(C - C^T) also unchanged). X_-3 is weaker: X_-3 = S - 4C is congruent to S mod 4, so its separating counts (valuation-1 entries) are forced by S; X_3 = 2S + SK is not (its separating entries have valuation >= 2).
(4) But at this key X_3 brings no information beyond SNF(S). Across 92 key-kept plants, X_3 changes exactly when SNF(S) changes (0 exceptions). In every key group at n = 8..12, X_+-3 never splits a (key, SNF(S)) class. On the 72 it gives the same 13 / 59 partition. The odd parts agree (703 = Phi_18(3), 757 = Phi_18(-3)). So the whole separation sits in the 2-adic (x^2 - 1)-primary part of the Coxeter module, the same place SNF(S) sees it.
(5) Why no truly independent invariant exists at this scope. Every Cartan-level invariant (SNF(S), q_S, X_a, the F-047 profile, Z-conjugacy of Phi) rests on the one necessity step "derived equivalent => C Z-congruent". So a second Cartan invariant is a second computation, not a second argument. E-318's referee already gave the independent computation (hand-built Cartans, solution counts mod 4 and 8). The non-Cartan invariants on the record are no use here: HH^* = k for every LNA (theorem, 2312.14699 / 0805.1018 Prop. 5.1; E-170, E-172), so HH cannot separate any two LNAs, and Avella-Alaminos-Geiss is for gentle algebras (R-008).

Not claimed: that X_3 is a function of (key, SNF(S)) in general. That is only observed here; I have no lemma. Also not claimed: that no categorical invariant separates the pair. Refuted by: a key-kept plant, or an LNA pair, with SNF(S) equal and X_3 different, which would make X_3 independent at the key; or a non-Cartan derived invariant computable on LNAs that is not trivial.

## Evidence

E-318 pair, SNF entries other than 1:

| invariant | 003033030 (13-class) | 030330030 (59-class) | |
|---|---|---|---|
| S = X_1 = F[x-1] | (0, 2, 4) | (0, 2, 2) | separates |
| SK = X_-1, F[x+1] | (0, 0, 0) | (0, 0, 0) | equal |
| F[Phi_18] | (0^6, 3, 3) | (0^6, 3, 3) | equal |
| X_2 / X_-2 | (3, 3, 171) / (657) | same | equal |
| X_3 | (4, 8, 5624) | (4, 4, 11248) | separates |
| X_-3 | (2, 4, 12112) | (2, 2, 24224) | separates |

E-318 walk (F3 R8 F3 F3 F6): steps 1, 2, 3, 5 are J = 0, tilting, key kept, and no invariant changes. Step 4 (J != 0, not `tiltingPlus`) changes exactly S, F[x-1], X_3 and X_-3. The end node's invariants equal those of LNA 030330030.

Partition of the 72: each separating invariant gives 13 + 59, and all invariants together still give 13 + 59. Census: 11 orbits (19, 19, 14, 7 in the 59-class; 4, 4, 1^5 in the 13-class), none mixed, 0 invariant changes on 226 edges.

Plants. S-fixed: 257 + 178 unimodular (seed 125, 3000 tries). Key-fixed: 31 + 61 same-key plants; (S same, X_3 same) 3 + 52, (S changed, X_3 changed) 28 + 9, and never one without the other. SK unchanged in all 92.

LNA scan, (key, SNF(S)) classes split by SK / X_+-2 / X_+-3: n = 8, 9, 11: 0. n = 10: 2 groups. n = 12: 3 groups. Every split is by X_2 only, at the prime 3, and in each of the 5 example pairs the F-047 profile also separates (e.g. 34504030 / 00000030, 0030003000 / 0033040000). So nothing new beyond F-047.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skew.py                  # 3 s: pair table, partition, census, congruence, walk
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skewscan.py plant 3000    # 1.5 s: S-fixed plants, planted-error control
RAND=20000 timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skewscan.py keyplant   # 88 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skewscan.py scan 12      # 84 s (n = 10: 5 s, n = 11: 20 s)
```
Needs python-flint and cypari in the venv (`uv pip install --python .venv/bin/python cypari python-flint`). SNF by flint, cross-checked by pari in maverick_skew.py.

## Prior record

E-197 (r066, mine) ran SNF of Phi +- I, f(Phi), G on the n = 11 periodic groups and found them constant. It did not run SNF(xC + C^T). E-154 used SNF(xC + C^T) at n = 7 as a YES check on key-kept children (non-informative there). F-047 is the per-factor profile, and here it adds nothing beyond x - 1. E-318 has the independent SNF recomputation. E-170 / E-172 settle that HH^* = k on LNAs. New: X_a = coker(Phi - a) applied to E-318; the proof that it is not S-determined (plants); and the observation that at this key it is empirically tied to SNF(S), so the "second invariant" request has no informative answer at Cartan level.

## Code changed

None in the library. New scripts: workshop/rounds/125/maverick_skew.py, maverick_skewscan.py. No tests touched.

## Next

- theorist: is coker(Phi - 3) determined by coker(Phi - 1) and the Coxeter polynomial for the 2-adic (x^2 - 1)-primary part of rank 5 seen here (ker(Phi-1) rank 1, ker(Phi+1) rank 3)? A lemma would turn (4) into a statement. A counterexample would make X_3 a real second witness somewhere.
- chair: T10 / E-318 should be promoted on the strength of SNF(S) plus the referee's independent recomputation. Asking for an "S-independent" Cartan invariant does not add a second argument, because every such invariant goes through H-024 necessity. A truly independent witness would have to be categorical (e.g. non-existence of a tilting complex with the right endomorphism ring). That is out of scope.
- toolsmith: add X_3 to INV in the J != 0 classifier. It is cheap (flint) and catches 2-adic changes that SNF(S) mod 4 misses on other groups.
