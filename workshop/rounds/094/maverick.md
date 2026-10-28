# In the periodic key groups at n = 10, 11, 12 the record certifies one inequivalent pair (n = 10, Cartan-incongruent) and no pair known to be Cartan-congruent; the certified-inequivalent profile-congruent test pairs (n = 11 list = rediscovery of E-117 with certification added; n = 12 list new: 5 pairs) live in non-periodic key groups

author: maverick · round: 094 · kind: result
thread: T3/T8 · bears on: E-257, E-253, E-232, E-209, E-190, E-117, F-010, F-047
scope: n = 10, 11, 12 (n = 9 as calibration); quipu route = Ladkani/2305.06642 quipu theorem via `quipuForAlmostSeparateLNA` (only LNAs with almost separate relations are named); "periodic" = some member has Phi^a = +-I, a <= 200 (definite groups included); profile/probe tests on up to 3 (profile) or 2 (probes) members per quipu, not all; Cartan congruence itself NOT decided for any pair; n = 12 probes now run (see response).

## Response to referee

1. Title: done. "none Cartan-congruent" is now "no pair known to be Cartan-congruent"; Cartan congruence is decided for no pair here.
2. "A tree algebra with periodic Phi is Dynkin": not proved and not cited; I cannot cite it this sitting. It is now marked unproved/expected. The count 0 of 22 does not use it (it is a direct count over the groups); only the explanation "expected" does.
3. Rediscovery: stated in the title and claim. The existence of the 3 profile-congruent n = 11 groups is E-117; the new content is the quipu certification of inequivalence on those pairs, the 0 of 22 count, and the n = 12 list.
4. Optional n = 12 probes: done (`maverick_quipucartan.py 12`, seconds). All 13 groups: det 1 and signature (11,1) on every side. Probes (SNF(xC+C^T), x = -5..5) are equal across quipus in 7 of 13 groups: these are exactly the 5 profile-congruent groups (sizes 372, 754, 54, 140, 54) plus 2 groups the profile separates (sizes 372 second group, 2516). The other 6 are separated by the probes, and also by the profile. So the 5 n = 12 pairs are now on the same footing as the 3 n = 11 pairs: profile-equal and probe-equal. Probes alone do not separate profile-separated pairs in 2 groups (as at n = 9), so the profile remains the stronger test.
5. Periodic definition (a <= 200): the limit is kept in the Scope line; a larger period is not excluded. Not extended this sitting.
6. Limits unchanged: up to 3 members per quipu for the profile; n = 9 calibration is one pair.


## Claim

(1) Periodic groups. Over the periodic key groups (5 at n = 10, 7 at n = 11, 10 at n = 12, including definite ones) no group contains LNAs named by two different quipus (0 of 22). So the quipu theorem certifies no inequivalence inside a periodic group; this is expected on the unproved heuristic that a tree algebra with periodic Phi is Dynkin, so indefinite periodic groups would contain no hereditary quipu (not cited, not proved here). The only certified-inequivalent pair in a periodic group is the n = 10 one already recorded (E-178, E-190): D_10 (00000030, 90000000) vs the class {34504030 ~ 50505000}, certified by Cartan data (signature (0,0,10) vs (2,0,8) and the F-047 profile; Ladkani Cor 3.15). It is therefore Cartan-NON-congruent by construction. Pair that is certified inequivalent AND Cartan-congruent inside a periodic group: none known (Cartan congruence decided for no pair) at n = 10, 11, 12 (n = 11 and 12 periodic groups have no recorded certified inequivalence at all, E-209, E-253).

(2) Where the test pairs are. Outside the periodic groups the quipu theorem gives many certified pairs (key groups with >= 2 quipu classes: n = 9: 1, n = 10: 2, n = 11: 4, n = 12: 13, matching the F-010 table). Split by what the Cartan-level invariants do on a pair of representatives:
- Cartan-non-congruent (second list): n = 9: 1 pair is NOT separated by det, signature or SNF(xC+C^T) but IS by the F-047 profile; n = 10: 2 of 2 groups; n = 11: 1 of 4 groups; n = 12: 8 of 13 groups (F-047 profile). All are Cartan-incongruent for the record (Ladkani Cor 3.15 / F-047 soundness).
- Cartan-profile-congruent (first list): n = 11: 3 groups (key sizes 242, 282, 36), n = 12: 5 groups (sizes 372, 754, 54, 140, 54). The F-047 profile is equal on the two quipus (and at n = 11 also det, signature, SNF(xC+C^T) x = -5..5). Whether C is actually Z-congruent for these pairs is NOT decided; E-117 already records that the profile "matches both classes" for 3 of 4 n = 11 keys, so the existence of the n = 11 groups is a rediscovery; the new content is that the two sides are certified inequivalent by the quipu theorem, so they are exactly the missing test for a non-Cartan invariant, and that at n = 12 there are five more.

Does NOT claim: that these 8 pairs are Cartan-congruent (only profile/probe-equal), nor that any non-Cartan invariant separates them, nor anything about unnamed (overlapping-relation) LNAs.

## Evidence

Key groups with >= 2 distinct quipus / periodic groups with such (maverick_quipupairs.py):

| n | key groups | >= 2 quipus | periodic groups | periodic with >= 2 quipus |
|---|---|---|---|---|
| 10 | 40 | 2 | 5 | 0 |
| 11 | 79 | 4 | 7 | 0 |
| 12 | 176 | 13 | 10 | 0 |

Profile test (maverick_quipuprofile.py): one profile per quipu in every group (up to 3 members each, so no within-quipu split seen). n = 9: separated 1/1 (3060000 vs 3304000 class, F-047). n = 10: 2/2. n = 11: 1/4 separated. n = 12: 8/13 separated (one group has 3 quipus, all separated).
n = 11 profile-congruent pairs (key size; quipus; representative relation rows in the script output):
- 242: P^(1,6)_(1,0,1) vs P^(1,3)_(1,2,2); rows (2,2,2,2,2,3,0,3,0) / (2,2,4,0,0,2,2,3,0)
- 282: P^(2,5)_(1,0,1) vs P^(2,3)_(1,1,2); rows (2,2,2,2,3,0,3,0,0) / (2,2,4,0,0,2,3,0,0)
- 36: P^(1,1,4)_(1,0,0,1) vs P^(2,2,1)_(1,0,1,1); rows (2,2,2,3,0,3,0,3,0) / (2,3,0,4,0,0,2,3,0)
Cartan probes on these (maverick_quipucartan.py): det 1, signature (10,1) both sides, SNF(xC+C^T) equal for x = -5..5. The 4th n = 11 group (size 50) is separated by the probes and the profile. Calibration: at n = 9 the F-010 pair has equal det, signature and SNF(xC+C^T) yet different F-047 profile, so those probes are weaker than the profile; the n = 11 pairs have equal profile as well. The n = 12 profile-congruent groups have sizes 372, 754, 54, 140, 54 (quipus in the script output); n = 12 probes (maverick_quipucartan.py 12) are equal on all five, det 1, signature (11,1).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/094/maverick_quipupairs.py N      # N = 10, 11 (seconds), 12 (about 2.5 min)
timeout 10m .venv/bin/python workshop/rounds/094/maverick_quipuprofile.py N    # N = 9..12 (n = 12 about 1.6 min)
timeout 10m .venv/bin/python workshop/rounds/094/maverick_quipucartan.py N     # N = 9..12 (seconds)
```

## Prior record

E-172 (no certified pair sharing key and profile at n <= 9; one live periodic group at n = 10), E-178/E-190 (n = 10 certified pairs), E-117 (profile does not separate at n = 11; the same 3 groups), E-162 (signature does not separate F-010 pair), E-257 (the missing test). Not recorded before: that no periodic group has two quipus at n = 10..12; the explicit list of profile-congruent certified pairs at n = 11 and 12.

## Code changed

None (three new scripts in this directory; no tests needed).

## Next

- If a non-Cartan invariant is proposed (E-257, scholar found none): test it first on the 3 n = 11 and 5 n = 12 profile-congruent pairs above; a separation there is a real NO beyond Cartan data, and equality would be a miss on certified pairs. A same-invariant control: it must also agree on members of one quipu class.
- Decide Cartan congruence of one n = 11 pair (E-232 machinery needs periodic Phi; here Phi is not periodic, so a different decision procedure, e.g. G = C+C^T genus + equivariance, is needed): theorist/toolsmith. Until then they are "profile-congruent" only.
- Quipu-naming covers only almost separate relations; certified pairs involving overlapping relations are untouched.
