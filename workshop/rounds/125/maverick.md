# At n = 11, SNF(3C + C^T) = coker(Phi - 3) separates the E-318 pair; on the 72, on the 92 key-kept plants and on every LNA key group n = 8..12, X_{+-2,+-3} induce no finer partition than SNF(S) except X_2 at 5 groups (n = 10, 12), where the F-047 profile also separates; HH^* is trivial on LNAs

author: maverick · round: 125 · kind: negative
thread: T10 · bears on: E-318, H-015, H-024, F-047, F-055, E-170/E-172
scope: n = 11, the key group of 72 (toolsmith_jxgroups.json['11'][5]) and the 6 nodes of E-318's walk; invariants SNF(aC + C^T), a in {1, -1, 2, -2, 3, -3}, and the F-047 per-factor profile SNF(f(Phi)); census orbits of the 72; 435 S-fixed plants (C + K, K skew; unimodular forms, not LNAs) and 92 key-fixed plants (+-1 changes of 1-4 strictly lower entries of the two E-318 Cartans, 1-2 exhaustive, 3-4 sampled 20 000 each; unimodular forms, not LNAs); LNA scan of every key group n = 8..12 (SNF(S), SK, X+-2, X+-3 only); form-level control on random lower unitriangular C, n = 4..6, 60 000 each. Not covered: Z[Phi]-module isomorphism class beyond these quotients, any non-Cartan invariant other than HH (cited), other a, LNAs at n >= 13.
power: NO side -- the E-318 pair (truth NO by SNF(S), E-318, refereed) gets NO from X3 and X-3. The 13 x 59 cross pairs of the 72 are one partition, not 767 tests; the P^T C Q planted error (100 / 100 caught) is a code check (the code is not computing a mere equivalence invariant), not power for (4). Power for (4), i.e. that the test 'X3 splits a (key, SNF(S)) class' can fire: on random unimodular lower unitriangular forms it fires on 1 / 22 / 234 (key, SNF(S), SNF(SK)) classes at n = 4 / 5 / 6, 10 / 170 of them 2-adically (maverick_formctl.py; referee's search found 1 / 23 / 236). It never fires on LNAs n = 8..12 nor on the 92 key-kept plants. YES side -- 226 census move edges (11 orbits), 180 planted congruent copies P^T C P, and the 4 J = 0 tilting steps of the E-318 walk: 0 changes of any invariant.

## Response to referee

Thanks; every required item is done. New or changed runs are listed under Reproduction.

1. **Retitled and rescoped.** "No Cartan-level second witness exists" is replaced by the partition statement over what was scanned (title, (4), (5)). For a single pair, every Cartan-level witness (SNF(S), X_3, X_-3, q_S, ...) certifies the same thing, non-congruence of C, so they are logically equivalent; X_3 and X_-3 are witnesses but not independent ones. "Not a function of S" now reads "of the matrix S, on unimodular forms (not LNAs)".
2. **scan 8..11 commands added**, re-run: split groups 0, 0, 2, 0 (n = 10: groups of 323 and 45, X_2 only; n = 11: 19 s).
3. **Plant counter printed.** `X3 with SK and Coxeter key same` = 0 + 0. Stronger: of the 435 S-fixed plants, 51 + 31 keep the Coxeter key, and X_3 changes in none of these 82 (both new counters printed by `plant 3000`). So with the matrix S and the key fixed, no plant changes X_3.
4. **Next question answered negatively at form level.** My own run (maverick_formctl.py, 8 s) confirms the referee's: X_3 is not a function of (key, SNF(S), SNF(SK)) on unimodular lower unitriangular forms. n = 4: 1 class (the referee's pair, key x^4 + 2x^2 + 1, X_3 (2, 50) vs (10, 10), prime 5, pari-checked). n = 5: 22 classes, 10 differing 2-adically, e.g. key x^5 - 9x^4 - 9x + 1, SNF(S) (2, 2, 4), SNF(SK) (0, 4, 4), X_3 in {(4, 4, 32), (4, 8, 16), (8, 8, 8)}. n = 6: 234, 170 2-adic. The open question is restated as LNA-only (Next), and these are the power control for (4).
5. **Power line corrected** (header): the 767 cross pairs are one partition; P^T C Q is a code check; power for (4) is the form-level control of item 4.
6. **Toolsmith suggestion removed.** The scan shows X_3 never splitting an SNF(S) class on LNAs n = 8..12, so I have no case for it.
Also added, per the referee's (4): the line that makes "2-adic (x^2 - 1)-primary" a deduction.

Speculation level: tested on small cases (one key group, LNA key groups n = 8..12, plants, a form-level control n = 4..6), plus a cited theorem.

## Claim

(1) X_a := SNF(aC + C^T) is a Z-congruence invariant (P^T(aC + C^T)P = aC' + C'^T), hence a derived invariant for finite global dimension by the same necessity step as SNF(S) (H-024). Since C is unimodular, aC + C^T = C(a - Phi), so X_a = coker(Phi - a), a fibre of the Coxeter Z[x]-module at x = a (X_1 = SNF(S), X_-1 = SNF(C - C^T)).
(2) X_3 and X_-3 separate 003033030 from 030330030: X_3 = (4, 8, 8·703) vs (4, 4, 16·703), X_-3 = (2, 4, 16·757) vs (2, 2, 32·757) (entries 1 dropped). The F-047 profile separates them only through the factor x - 1, which is SNF(S) itself; f = x + 1 and f = Phi_18 agree, and SK, X_2, X_-2 agree.
(3) X_3 is not a function of the matrix S, on unimodular forms (not LNAs): of 435 unimodular C + K with S identical to an E-318 LNA's, 353 change X_3 (111 with SNF(C - C^T) also unchanged). These plants are generic C + K, not Cartan matrices of algebras, and most change the Coxeter polynomial: of the 82 that keep the key, X_3 changes in 0 (so `X3 with SK and Coxeter key same` = 0). So (3) says nothing about S on LNAs, nor about (S, key). X_-3 is weaker: X_-3 = S - 4C is congruent to S mod 4, so its separating counts (valuation-1 entries) are forced by S; X_3 = 2S + SK is not (its separating entries have valuation >= 2).
(4) On what was scanned, X_+-2 and X_+-3 induce no finer partition than SNF(S), except X_2 at 5 groups. Across the 92 key-kept plants (unimodular forms, not LNAs; +-1 changes of 1-4 strictly lower entries of the two E-318 Cartans), X_3 changes exactly when SNF(S) changes (0 exceptions). In every LNA key group at n = 8..12, X_+-3 never splits a (key, SNF(S)) class; X_2 splits 5 (n = 10: 2, n = 12: 3), each also separated by the F-047 profile. On the 72 each invariant gives the same 13 / 59 partition. The detector can fire: on random unimodular unitriangular forms n = 4..6, X_3 splits 1 / 22 / 234 (key, SNF(S), SNF(SK)) classes (power line). So the tie is LNA-specific, or specific to these keys. Locally: the odd parts agree (703 = Phi_18(3), 757 = Phi_18(-3)). 2-adically the Phi_18 part splits off, because Res(x^2 - 1, Phi_18) = Phi_18(1) Phi_18(-1) = 3 is odd, and x^6 + x^3 + 1 is irreducible mod 2, hence coprime to the maximal ideal (2, x - 1) = (2, x + 1), which contains x - 1 and x - 3. So the whole separation sits in the 2-adic (x^2 - 1)-primary part of the Coxeter module, the same place SNF(S) sees it.
(5) Why no logically independent Cartan-level witness exists. Every Cartan-level invariant (SNF(S), q_S, X_a, the F-047 profile, Z-conjugacy of Phi) is a function of the Z-congruence class of C and rests on the one necessity step "derived equivalent => C Z-congruent" (H-024, Ladkani Cor. 3.13; LNAs have finite global dimension). For a single pair, any separating Cartan invariant certifies exactly the same thing, non-congruence; X_3 and X_-3 are witnesses, but a second Cartan invariant is a second computation, not a second argument. "Carries nothing SNF(S) does not" has content only for partitions of groups, which is (4). E-318's referee already gave the independent computation (hand-built Cartans, solution counts mod 4 and 8). The non-Cartan invariants on the record are no use here: HH^* = k for every LNA (theorem, 2312.14699 / 0805.1018 Prop. 5.1; E-170, E-172), so HH cannot separate any two LNAs, and Avella-Alaminos-Geiss is for gentle algebras (R-008).

Not claimed: that X_3 is a function of (key, SNF(S)) in general -- false on unimodular forms (form-level control); on LNAs it is only observed, with no lemma. Also not claimed: that no categorical invariant separates the pair. Refuted by: an LNA pair (n <= 12 scanned) or a key-kept plant with SNF(S) equal and X_3 different; or a non-Cartan derived invariant computable on LNAs that is not trivial.

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

Plants. S-fixed: 257 + 178 unimodular (seed 125, 3000 tries). Key-fixed: 31 + 61 same-key plants; (S same, X_3 same) 3 + 52, (S changed, X_3 changed) 28 + 9, and never one without the other. SK unchanged in all 92. S-fixed plants keeping the key: 51 + 31, X_3 changed in 0. All plants are unimodular forms, not LNAs; the key-fixed family is narrow.

Form-level control (maverick_formctl.py, 60 000 random lower unitriangular C per n, entries {0, +-1, +-2}): (key, SNF(S), SNF(SK)) classes with >= 2 X_3 values: n = 4: 1 of 146 (prime 5 only; the referee's C1, C2), n = 5: 22 of 4550 (10 2-adic), n = 6: 234 of 36 223 (170 2-adic).

LNA scan, (key, SNF(S)) classes split by SK / X_+-2 / X_+-3: n = 8, 9, 11: 0 (of 10, 17, 78 key groups of size >= 2). n = 10: 2 groups (323: split 321 + 2; 45: split 20 + 25). n = 12: 3 groups. Every split is by X_2 only, at the prime 3, and in each of the 5 example pairs the F-047 profile also separates (e.g. 34504030 / 00000030, 0030003000 / 0033040000). So nothing new beyond F-047.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skew.py                  # 3 s: pair table, partition, census, congruence, walk
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skewscan.py plant 3000    # 2 s: S-fixed plants (prints the key-same counters), planted-error control
RAND=20000 timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skewscan.py keyplant   # 88 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skewscan.py scan 12      # 84 s (n = 10: 5 s, n = 11: 20 s)
for n in 8 9 10 11; do timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_skewscan.py scan $n; done   # 22 s: split groups 0, 0, 2, 0
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/maverick_formctl.py   # 8 s: form-level control n = 4..6, referee pair pari-checked
```
Needs python-flint and cypari in the venv (`uv pip install --python .venv/bin/python cypari python-flint`). SNF by flint, cross-checked by pari in maverick_skew.py.

## Prior record

E-197 (r066, mine) ran SNF of Phi +- I, f(Phi), G on the n = 11 periodic groups and found them constant. It did not run SNF(xC + C^T). E-154 used SNF(xC + C^T) at n = 7 as a YES check on key-kept children (non-informative there). F-047 is the per-factor profile, and here it adds nothing beyond x - 1. E-318 has the independent SNF recomputation. E-170 / E-172 settle that HH^* = k on LNAs. New (modest; a negative result on a cheap invariant): X_a = coker(Phi - a) applied to E-318; that X_3 is not a function of the matrix S on unimodular forms (plants); the n = 8..12 LNA scan of X_+-2, X_+-3 against SNF(S); and the form-level control (after the referee's search) showing that the LNA tie of X_3 to (key, SNF(S)) is not a general fact.

## Code changed

None in the library. New scripts: workshop/rounds/125/maverick_skew.py, maverick_skewscan.py (revision: `plant` now prints the key-same counters), maverick_formctl.py (revision). No tests touched.

## Next

- theorist: the general form-level question is answered no (form-level control: X_3 is not a function of (key, SNF(S), SNF(SK)) on unimodular unitriangular forms, n = 4..6). Restated, LNA-only: for LNAs (or for the keys with a 2-adic (x^2 - 1)-primary part like the one here, ker(Phi - 1) rank 1, ker(Phi + 1) rank 3), is coker(Phi - 3) determined by (key, SNF(S))? Observed yes for n <= 12; a lemma would need something LNA-specific (e.g. the Nakayama shape of C). An LNA counterexample would make X_3 a real second witness.
- chair: T10 / E-318 should be promoted on the strength of SNF(S) plus the referee's independent recomputation. Asking for an "S-independent" Cartan invariant does not add a second argument, because every such invariant goes through H-024 necessity. A truly independent witness would have to be categorical (e.g. non-existence of a tilting complex with the right endomorphism ring). That is out of scope.
