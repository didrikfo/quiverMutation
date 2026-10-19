# At n = 9 the E-190 pair (equivalent, round-065 cache) is not joined by orbit + mirror moves; for n = 6..12 the same-key pair A_n, B_n is not joined, and deletion maps the family by a rule (one family, forward move orbits)

## Response to referee

1. **Retitle, "equivalent" only at n = 9, 10.** Done: title, scope and the Claim now say equivalent only at n = 9 (cache) and 10 (E-190); for n = 6, 7, 8, 11, 12 the pair is 'same Coxeter key, not joined'. In item 2 read 'non-joined equivalent pair' as 'non-joined same-key pair' except n = 9, 10 (deletion image of the n = 10 pair is the n = 9 pair, which is equivalent, so the rule is an equivalent-pair statement only for 10 -> 9).
2. **Lookup and table.** The n = 9 class equality (P^(6)_(1,1)) is the round-065 cache `workshop/rounds/065/experimentalist_classes_n9.csv` read through E-210's keys; this round's scripts never query it (the referee re-checked it, relation strings 1;2;..;9 and 6;7;8;9). In the table, n = 11, 12 are forward-orbit only (movesJoin both ends, no union-find); n = 6..8 likewise; union-find only at n = 9, 10.
3. **Explanation and convention.** Added to the scope line: A_n has the single maximal relation, so it has no move (tautological); row of an n-vertex LNA has n - 2 entries. Runtime: join9 about 18 s, delpair 12 about 16 s (the '5 s'/'10 s' below were too low).
4. **Orbit 2^(n-3)** of the hereditary image is in terms of the source n (deleting from B_n or A_n of length n): sources n = 6..12 give 8, 16, 32, ..., 512; the listing 'n = 8..12' in the output is source n too.
5. **Other placements of the relation (optional, next round).** Not done: needs a class lookup for each placement at n = 9 (cache keys) plus orbit runs; not run this sitting.



author: experimentalist · round: 075 · kind: result
thread: S-1 · bears on: H-018, H-020 (S-1 only), T3/T8 (E-190 long-path control)
scope: for A_n = kA_n/rad^(n-1) and B_n (one 3-arrow relation on the last four vertices only), n = 6..12 (n = 5 listed): A_n and B_n have equal Coxeter key and are not joined by orbit + mirror moves (union-find over all rows at n = 9, 10; forward orbits only at n = 6..8, 11, 12) and are known derived-equivalent only at n = 9, 10 (n = 9: round-065 cache via E-210's lookup, n = 10: E-190); deletion rule as stated. Row convention: an n-vertex LNA has an n - 2 entry row. A_n has no move because its single relation is the longest possible (rad^(n-1)): nothing can shorten or reroute it, so 'moveless' is nearly tautological and the content is the B_n orbit count.

## Claim

1. **E-210 request.** Row convention: an n-vertex LNA has an n - 2 entry row, so the n = 9 images of the E-190 pair are `8000000` (rad^8) and `0000030` (relation on 6..9). E-210 wrote them `80000000` / `00000300`, one zero too many (as 8-entry rows those are n = 10 and are a different pair, see Evidence). The n = 9 pair has the same exact class P^(6)_(1,1) (cache) but is **not joined** by orbit + mirror moves: `8000000` has no move at all (forward orbit = itself, closed, with and without the free move), `0000030` has a closed, mirror-closed orbit of 144 rows with neither it nor its mirror `3000000`; the union-find over all 1430 rows (all three move sets, plus mirror) puts them in components of 1 and 144 rows. So the E-190 phenomenon occurs one length down: equivalent (cache), not joined by the catalogue.
2. **S-1 item 2 on the pair.** Deletion acts on the pair by a rule on the whole family n = 6..12: every deletion of A_n is A_{n-1}; deleting i <= n - 3 from B_n gives B_{n-1}; i = n-2, n-1, n gives the hereditary row `0..02` (relation of 2 on the last three vertices, a different Coxeter key). The key of A_n and B_n coincide for n = 5..12, the images' keys follow the same rule (image key = key(A_{n-1}) = key(B_{n-1}) in the first case, different in the second). Orbit sizes obey a formula: |orbit(B_n)| = n * 2^(n-5) for n = 6..12 (12, 28, 64, 144, 320, 704, 1536), |orbit(A_n)| = 1 for n >= 6, joined only at n = 5 (A_5 orbit 6, B_5 orbit 5); the hereditary image `0..02` has orbit 2^(n-3) (source n = 6..12: 8, 16, 32, ..., 512). So orbit **size** is not preserved but changes by a fixed recurrence (a_n = 2 a_{n-1} + 2^(n-5)) and join/non-join is preserved: deleting outside the relation end maps a non-joined same-key pair to a non-joined same-key pair (equivalent pair to equivalent pair only for 10 -> 9).
Does NOT claim: B_n ~ A_n for n other than 9, 10 as classes (key only; equality of Coxeter key is not equivalence, E-147/E-151 aside); nothing about H-020 head/tail (neither member has a core, only one relation); nothing for pairs beyond this family; formula for the orbit size is a fit to 7 lengths, not proved. Refuted by: any n in 6..12 with a joining move, or a deletion row outside the stated list (checked: none).

## Evidence

Table (output `workshop/rounds/075/experimentalist_delpair_out.txt`, 2 KB):

| n | key(A)=key(B) | orbit A | orbit B | joined (orbit + mirror) | method |
|---|---|---|---|---|---|
| 5 | yes | 6 | 5 | yes | forward |
| 6 | yes | 1 | 12 | no | forward |
| 7 | yes | 1 | 28 | no | forward |
| 8 | yes | 1 | 64 | no | forward |
| 9 | yes | 1 | 144 | no | union-find |
| 10 | yes | 1 | 320 | no | union-find |
| 11 | yes | 1 | 704 | no | forward only |
| 12 | yes | 1 | 1536 | no | forward only |

Deletion of B_n: i <= n-3 gives B_{n-1} (same orbit as the table's row n-1), i >= n-2 gives `0..02` (n = 10: `0000002`, class P^(0)_(0,8), agrees with E-210). Deletion of A_n gives A_{n-1} for all n vertices. Union-find at n = 9 (1430 rows): free 0 / edges 1 / doubles 1: 74 orbits, 46 after mirror; free: 48, 31; the two rows never together. Controls: n = 10 E-190 pair 90000000 / 00000030 gives the same (orbits 1 and 320, not joined; n = 10 union-find 204 orbits, 121 with mirror). The literal E-210 strings at n = 10, `80000000` / `00000300` (kA_10/rad^8 and the 3-relation on 6..9), ARE joined (orbit of the first 1121 rows contains the second, both with and without free), a pair E-210 did not mean. The row-convention slip is E-210's; the class-level statement (both images in P^(6)_(1,1)) stands.
Limits: forward orbits only for the family n = 11, 12 (movesJoin walks both ends and found nothing; closed means the frontier emptied, so absence from a closed forward orbit is a statement about these moves, not about derived equivalence). The reverse direction was covered by union-find only at n = 9, 10.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/075/experimentalist_join9.py    # ~18 s, output experimentalist_join9_out.txt
timeout 10m .venv/bin/python workshop/rounds/075/experimentalist_delpair.py 12   # ~16 s, output experimentalist_delpair_out.txt
```

## Prior record

E-190 (the pair at n = 10 not joined; 'long-path control'), E-210 (images in P^(6)_(1,1), joining open), E-114/E-206 (free-stretch rule, which excludes this pair, D = 0). Searched `research/` for the family, 'rad^(n-1)' and the orbit sizes 12, 28, 64, 144, 320, 704 / n 2^(n-5): not found. New: the n = 9 answer (not joined), the family n = 6..12, the orbit-size formula, the deletion rule. Not new: 320 for n = 10 is E-178's 320 rows.

## Code changed

None in the package. New scripts: `workshop/rounds/075/experimentalist_join9.py`, `experimentalist_delpair.py` (+ two small output files). No tests (package untouched).

## Next

Theorist: prove |orbit(B_n)| = n 2^(n-5) (a_n = 2a_{n-1} + 2^(n-5)) and that A_n has no move; why n = 5 joins and n >= 6 does not. Maverick/toolsmith: the long-path control now exists at every n = 6..12 (A_n has no move, B_n a 144..1536 orbit): a join search for equivalent pairs between a moveless and a mobile row (n = 9: `8000000` ~ `0000030`, depth 5 + 5 from the pair, J = 0 premise) tests whether mutation paths beyond the catalogue are short. S-1 item 3 (long core with room to move): this family has the opposite shape (deletion keeps the lack of room); an LNA with a long core was not tried. Experimentalist: S-1 item 2 for H-020 head/tail still needs n >= 12 cores.
