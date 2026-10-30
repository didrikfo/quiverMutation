# q_S has no power beyond signature at n = 10: its class is consistent with the equivalent periodic pair 90000000 / 00000030 (miss; one key group)

author: skeptic · round: 099 · kind: negative
thread: T3/T8 · bears on: E-265, E-259, E-190, E-178, E-210
scope: n = 10 only, the one Phi^18 key group (323 rows); exact Cartan arithmetic; one pair tested directly, the 321/2 partition is a det/signature partition with `disc_iso` run only inside the signature-equal class. A test with almost no power at n = 10; the agenda-2 item stands for n >= 11. Says nothing about n = 11, 12 or about the 140 group.

## Response to referee

1. Done (wording). The 321/2 partition is fixed by det and signature (both 4; signatures (10,0) vs (8,2)). `disc_iso` was called only against a class representative of equal det and signature; it returned True on all 320 other rows of the 321 class (so it was exercised there, but never as a splitter). "q_S constant within the 321 class" therefore means: det and signature constant, and an isometry exists to the representative. Evidence-section sentence on this stands.
2. Done. E-178 records 34504030 ~ 50505000 equivalent (7-step move, joins at 4 + 4, control in that entry), and these are exactly the 2-row (8,2) class. Cross-check: that class is one derived class by explicit moves and q_S puts both rows in one class (consistent). It also supports E-190's separation of D_10 from {34504030, 50505000}.
3. Done by grep of `research/` (E-190, E-178, FINDINGS l. 595). What is recorded for the 321 rows: E-178 gives class 0 as 320 rows forming one orbit under the move set ("key-and-orbit classes", 320 rows, one orbit, includes 00000030), i.e. joined by explicit orbit moves. The 321st row, 90000000, is NOT joined by moves: E-190 identifies it as kA_10/rad^9 = the tame D_10 class (Cartan-level / identification argument) and states the pair is "an equivalent pair the orbit-plus-mirror moves do not join" (also E-224: no ball meeting at radius 5+4; E-178: J0 join empty through 6+6). So "all 321 explicitly joined" is false: 320 by orbit moves, 1 (90000000) by identification with D_10. The HIT criterion's "any two rows known equivalent" is then met by the 320-orbit pairs and by the 90000000/00000030 pair; I tested the latter directly and a representative-vs-320 test covers the former. The q_S-constancy over the 320-orbit is the weak-power signature fact above.
4. Done. "Closes agenda 2 for n = 10" is replaced by: a test with no power at n = 10; the open item stands for n >= 11 (and the 140 group is untouched).

Additional, not requested: independent check of the isometry without `disc_iso` (`workshop/rounds/099/skeptic_qs_check.py`, output `skeptic_qs_check.out`, <1 s). Brute-force enumeration of D = S^-1 Z^10 / Z^10 for both rows (order 4 each) and all bijections preserving q = x^T S x mod 2 and addition: 2 isomorphisms exist, matching the isometry found by `disc_iso` (q value multisets {0, 1/2, 1/2, 1} on both sides). Power check: q_S(90000000) vs q_S(34504030) has 0 q-preserving isomorphisms (q values {0, 3/2, 3/2, 1}), so q_S does separate the two classes, but this is Milgram: q on this group is determined by signature mod 8 (10 vs 6), so it is the signature again. Why the test is nearly vacuous: the discriminant group is (Z/2)^2 and there are very few forms of given det and signature mod 8.

## Claim

Stated before running. HIT = q_S(C + C^T) non-isomorphic for 90000000 and 00000030 (equivalent per E-190), or for any two rows of
that key group known to be equivalent: that would refute q_S as a derived-equivalence invariant (for finite global dimension).
MISS = isomorphic. Result: MISS. q_S(90000000) is isomorphic to q_S(00000030) (explicit isometry found; |det S| = 4, order-4
group, both signature (10,0), i.e. D_10). So E-265's invariant survives this test. It does NOT show q_S is a complete invariant, nor that
the pair is joined by moves.

The test has little power here, and this should be read as such: the key group of 323 rows splits under q_S into exactly two
classes, 321 rows (signature (10,0), the D_10 class of E-190) and 2 rows (34504030, 50505000, signature (8,2)). That is the same
split signature already gives, and it agrees with E-190's "2 derived classes". An inequivalent pair is thus separated, but by a
cheaper invariant. Where q_S earned its power (E-265: n = 11, 12) the groups are non-periodic.

## Evidence

- key(90000000) = key(00000030): True. det S = 4 and 4; signature (10,0) and (10,0).
- `disc_iso` (E-265's own brute-force isometry search, imported unchanged): True, with generator images printed in the .out file.
- All 323 rows of the key group compared to class representatives (det, signature, then `disc_iso`): 2 q_S classes, sizes 321 and 2;
  90000000 and 00000030 both in the 321-class. Within the 321 class q_S is constant, consistent with E-265's within-quipu control
  (and E-178's 320 rows + 90000000 being one class, F-045 / FINDINGS l. 595).
- Shared with the audited thing: `disc_iso` and `cartan` come from workshop/rounds/097/maverick_cartanZ.py (not re-derived). A
  false-positive in `disc_iso` (answering True too easily) is excluded by its 7-of-8 splits in E-265, and here it returns False
  on the (10,0) vs (8,2) pair only through the det/signature pre-filter, so it was not exercised as a splitter in this group.
  The sole independent check is therefore det/signature. Not re-verified: the Nikulin-type statements (class = genus), irrelevant
  here since only isomorphism of q_S is asserted.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/099/skeptic_qs_periodic.py     # 16 s; output in skeptic_qs_periodic.out
```
The script stubs the `cypari` import (not installed in this environment; only the module-isomorphism search needs it, unused here).

## Prior record

E-265 states the n = 10 pair was "not matched" (E-259: no periodic group holds two quipu-named LNAs); E-197 showed q_G constant on
the n = 10/11 periodic groups (a different form, G the symmetrised Euler form of the inverse Cartan side, no split). E-190 proved
the pair equivalent (90000000 is D_10). This does not close STATE agenda 2: it is a test with no power at n = 10 and the open item stands for n >= 11. New content is only the explicit
q_S iso and the 321/2 partition; the partition is not new information beyond E-190 + signature.

## Code changed

None (one new script).

## Next

- A refutation would need an equivalent pair in a non-periodic group at n = 11, 12 that q_S separates; E-265's 0 of 513 split orbits
  is on a small subset of orbit members: a skeptic round could enlarge that control using actual E-187 moves on rows q_S-separated in
  the 15 groups at n = 12.
- Maverick/toolsmith: the 140 group (same genus) still needs the Z[Phi] computation (cypari absent here; install before any run).
