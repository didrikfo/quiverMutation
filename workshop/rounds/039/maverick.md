# The n = 13 prediction of E-133 holds: the key class of the lone 3 at (4,5)/(5,4) is one move orbit of 4349 LNAs whose K >= 4 free ends go to two different n = 12 key classes, so "K >= 4 transports" fails first at n = 13

author: maverick · round: 039 · kind: result
thread: S-1 · bears on: E-118, E-125, E-133, H-020

## Claim

Speculation level: tested (key level for the images; move orbits computed on the whole key class).

At n = 13 the Coxeter key K13 = (1,1,0,-1,-2,-3,-4,-4,-3,-2,-1,0,1,1) is shared by exactly 5023 LNAs (full scan of the 208012 rows). The union-find of the rule table, free (length-2 stripping), edge moves and double mutations restricted to these rows (no move leaves the class) gives two orbits, of 4349 and 674 LNAs; the lone 3 at h = 4 (`00003000000`, K = 5) and at h = 5 (K = 4) lie in the same 4349-orbit. This orbit has 174 free ends with K >= 4 (K = 4: 76 + 64; K = 5: 34, both head and tail counted); their deletions land in exactly two n = 12 key classes, A = (..,-3,-4,-3,-2,..) (110 ends) and B = (..,-3,-3,-3,-2,..) (64 K = 4 ends, plus 2 in the 674-orbit). So "all ends of K >= 4 of one class map to one image class" fails at n = 13, at the length E-133 predicted from keys, and the failing source is a single derived-equivalence orbit, not a key coincidence. Smallest witness: the lone 3 itself, `0000300000000` (h,K) = (4,5): deleting the head gives (3,5) = B, deleting the tail gives (4,4) = A.
Not claimed: that A and B differ as derived classes beyond "their Coxeter keys differ" (that suffices if the key is a derived invariant, which it is); that the move orbit is the whole derived class (it is a lower bound; the key class is the upper bound, and they differ, 4349 + 674 vs 5023, only by the mirror-free split); that n = 12 K >= 4 holds (not run; prediction only); non-lone cores in this class (the table below mixes them; not interpreted).

## Evidence

| quantity | value |
|---|---|
| LNAs at n = 13 with key K13 | 5023 (shards 1247 + 1288 + 1263 + 1225) |
| move orbits inside (rules, free, edges, doubles) | 4349, 674; each is its own mirror orbit; moves leaving the class: 0 |
| lone 3 at h = 4, h = 5 | both in the 4349-orbit; forward `orbitReport` from either is 1 row (closed): no forward move acts on a lone 3, they are joined only through backward moves from other rows |
| n = 12 images of the two lone-3 ends | (4,4) key ..-4.. and (3,5) key ..-3,-3,-3..: different |
| K >= 4 ends by (orbit, K, image) | (4349,4,A) 76, (4349,4,B) 64, (4349,5,A) 34, (674,4,B) 2 |
| K = 3 ends (info) | A 244, B 142 |

Reading: the mirror pair (4,5),(5,4) at n = 13 behaves exactly like (3,4),(4,3) at n = 11; the shift K0 3 -> 4 for n 11 -> 13 is the key-level rule of E-133 and is now confirmed at the orbit level (the failing class is one orbit, as at n = 11 where orbit 15107 split 34 + 48).

## Reproduction

```
.venv/bin/python workshop/rounds/039/maverick_n13.py 2000        # 2 s: keys, forward orbit of the two lone 3s (1 row each)
for i in 0 1 2 3; do timeout 10m .venv/bin/python workshop/rounds/039/maverick_n13class.py $i 4 & done; wait   # 4 shards, ~5 min each, writes maverick_n13class_<i>.txt
.venv/bin/python workshop/rounds/039/maverick_n13ends.py orbits   # 5 s, reads the shard files
```

## Prior record

E-133 / r037: key-level prediction only (K >= 4 fails first at n = 13 at (4,5)); E-118: K >= 4 holds at n = 11 (10/10 classes), K0 = 3 fails. Grep of research/ for n = 13 lone-3 results and of RETRACTIONS.md: none. New: the full key class at n = 13 (5023), its orbits, the confirmation. Caveat inherited: image classes at n = 12 are compared by key, not by E-115 labels (not available at n = 12).

## Code changed

None in the library. New: `maverick_n13.py`, `maverick_n13class.py`, `maverick_n13ends.py` (+ four shard txt files). No tests touched.

## Next

- experimentalist: n = 12 for completeness (predicted K >= 3 fails at (3,5), K >= 4 holds; key class scan ~5 min with the same script, n = 12); and n = 15 (K >= 5) needs a 4-shard scan of a much larger row set: size with `--plan` first, probably overnight.
- theorist: why no forward move acts on a lone 3 but backward moves from other rows join (4,5) and (5,4); derive the {h,K} dependence of the key from H-020.
- S-1 Q1 answer so far: compatible deletion = from the shorter free side, tested only on lone relations.
