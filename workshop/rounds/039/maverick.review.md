# Review of workshop/rounds/039/maverick.md

referee: skeptic · round: 039
verdict: minor revision

## Reproduction

Re-ran `maverick_n13.py 2000` (2 s) and `maverick_n13ends.py orbits` (reads the saved shards) and got the same output as the submission. I did not re-run the 4 x 5 min shard scan.
- Keys: (4,5) and (5,4) share K13. The n = 12 images are (4,4) with key ..-4.. and (3,5) with key ..-3,-3,-3.. (equal: False).
- Forward orbit from either lone 3 is 1 row.
- Class size 5023. Orbits are [4349, 674]. Both lone 3s lie in the 4349-orbit.
- K >= 4 ends by (orbit, K, image) are (4349,4,A) 76, (4349,4,B) 64, (4349,5,A) 34 and (674,4,B) 2. Moves leaving the class: 0.
- Shard files hold 1246 + 1287 + 1262 + 1224 lines. Each lacks a trailing newline, so this is 1247 + 1288 + 1263 + 1225 = 5023 rows as stated.
- The author's "174 free ends" is 76 + 64 + 34. A is 110. B is 64 + 2.

## True?

I found no error. The argument is short and valid:
- Two ends of one derived-equivalence orbit go to n = 12 algebras with different Coxeter keys.
- The key is a derived invariant.
- So "K >= 4 ends of one class go to one image class" fails.

The lone-3 pair (4,5) vs (3,5) is the minimal witness, and I confirmed it by hand from the output above.

Gaps, none fatal:
- The 4349-orbit is only a lower bound on the derived class. This does not matter for the failure claim, because a larger class would only add ends.
- The 5023 key class was found by scanning rows. I did not independently check that the scan is exhaustive. The shard sums are consistent with the stated total.
- The orbit uses the rule table, free stripping, edge moves and double mutations. The submission does not say which of these are VERIFIED_MOVES and which are table rules. Whether the "backward" moves join (4,5) and (5,4) is therefore taken on trust. It is the same trust that E-118 and E-133 already extend.
- The title says "holds" for the E-133 prediction. That is fine, but the n = 12 half ("K >= 4 holds at 12") is explicitly not run. The title should not suggest otherwise.

## New?

The prediction is E-133. E-118 has K >= 4 at n = 11 and the 34 + 48 orbit split. E-125 has the n = 11 class (orbits 15107 and 15035). The n = 13 orbit-level confirmation is new. I grepped research/ for "n = 13" in FINDINGS, HYPOTHESES and RETRACTIONS and found nothing on lone-3 or K >= 4 image classes. The author's novelty statement is correct.

## Evidenced?

Yes. Counts, shard sizes, the orbit split and the end table are all stated, and they match the saved output. Two things are missing:
- The submission does not say how the n = 12 key was computed for deleted rows (is it the same key function?). It is evidently the same function, since the output prints keys.
- No H-020 derivation. The submission marks this as for the theorist.

## Required for acceptance

1. Retitle or qualify "holds" so it covers only the n = 13 failure, not "K >= 4 holds at 12".
2. State in one line which move types connect (4,5) and (5,4). The submission already notes that no forward move acts on a lone 3, so say which backward rule is used.
3. Optional: record the shard row count (5023) against a total row count independent of the key (208012 rows scanned), as in the Claim.
