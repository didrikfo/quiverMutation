# Theorist notebook (rewritten round 011)

## What I believe now
- Lists A/B (E-074) are not a parity-of-n accident: the key of `35` (even n) / `36` (odd n) contains, among single-relation rows, exactly two mirror pairs for n >= 14
  (`3@2`,`5@0` even; `3@3`,`6@0` odd). Their reduced orbits P, Q are small (n = 18: 774, 678), disjoint, self-mirror; `35`@even in P, @odd in Q (n = 12..20, all checked).
- Census rule "word alternates between two single-relation orbits with one key, each self-mirror" reproduces A at n = 12, 14, 16, 18 and B minus `5046 5056` at 13, 15, 17
  (about 1 min per n). `406` alternates but its two orbits are mirror images, so it is not key-coarser. Other parity words (`3336`, `405`, `5004`) sit in other two-orbit key classes.
- Letter 4: `4@0` holds `3@1` at every n = 12..16 (orbit of 3767 at n = 14); `k@0` (k >= 5) holds a `3@j` only if k = n mod 2. Smallest cut-off k = 5 (even n) / 6 (odd n) = the Q seeds.
  Explains over-representation of 4, not `444` itself. Observed, not proved.
- Null: no GF(2) functional of rows, no integer statistic mod 2/4, and no SNF of Coxeter-matrix polynomials separates P from Q. Do not retry these.
- Older (round 009): drift family `aax` (E-065/E-071): rigid for a = 3, 5, 6, 7; a = 4 merges via `444 -> 34 -> 44`; k(33x) = 2x explained; "orbit holds a word at all offsets" does not predict merging.

## What I tried
- Round 011: theorist_{invariant,stats,snf,single,keyclass,triples,predict,census2,word,mirror406,k0,extra}.py in rounds/011. Fast route: walk single-relation orbits only.
- Not done: `5046 5056` at n = 17 (orbit about 10^5+, run was cut), any proof of `35@o -> 35@o+2`, n > 20.

## Next
- BFS path in P from `35@0` to `35@2` at n = 12 (178 rows): list the moves; same for `4@0 -> 3@1`, and show `5@0` has no such move (parity obstruction in the move set?).
- Use the single-relation orbit method on the 7 of E-059 (`344 366 4044 ...`) and on `--max-word 5`: cheap census instead of hours-long orbitclass.
- Blind spots: the criterion was fitted to lists A/B; `5046 5056` break it; keys of single relations are probably the Happel-Seidel tree classes (not checked, would explain 2 vs 3 orbits).
