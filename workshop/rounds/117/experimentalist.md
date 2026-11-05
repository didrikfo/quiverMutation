# H-026 K-law at even n (n = 10..16, key level, head end, lone-3 classes): non-central classes obey min(h, n-4-h), the central class (h = n-4-h) does NOT (pre-registered prediction K0 = 5 at n = 14 falsified); at n = 15 h = 3 the image key is a function of the move component of the source row

author: experimentalist · round: 117 · kind: backlog
thread: S-1 · bears on: H-026
scope: key level (Coxeter polynomial of the vertex-1 deletion), head end, lone-3 key classes (key = key of the lone-3 row, class = all rows with that key, not checked to be one orbit); n = 14 K >= 1 (h = 5,4,3,2,1), n = 10, 12 K >= 0, n = 16 K >= 3 (h = 6,5,4), refits n = 13, 15, 17 at K >= h-1; class-level only for n = 15 h = 3, K >= 2 rows (4597) closed under forward rewrites (+4780 rows). Not run: tail end, q_S / derived-class test, non-lone classes, n >= 18.
power: (a) planted K0-1 miss (image keys at K = K0 merged) is caught (law test NO) on 8 of 8 refits; (b) rival formulas R1 (floor((n-1)/2)-2), R4 (floor(n/2)-2), R5 (h-1), R3m (long-side reading) are refuted on 7, 7, 12, 12 of the 12 odd-n fits; (c) class-level purity test: shuffling the image labels gives 159 (K=2) / 52 (K=3) mixed components on average (min 138 / 38), the real data 0.
backlog: item 3 settled (at move-component level; derived class not tested) · item 4 settled (with a finding: R2 and R3-small-side are the same as the law wherever h < n-4-h) · item 5 settled (prediction falsified at the central class)

## Claim

(1) Even n. The prediction logged before the run (`experimentalist_prediction.md`: K0 = 5, 4, 3, 2 and none at K >= 2 for (5,5), (4,6), (3,7), (2,8), (1,9) at n = 14) held for four of five classes and FAILED for the central class (5,5): its rows have ONE image key at every K >= 1 (1903 rows). The same holds at n = 16 (6,6) (K >= 3, 1469 rows, one key) and at n = 10 (3,3), n = 12 (4,4) (splits only at K = 0). Non-central even classes follow min(h, n-4-h): n = 10 h=2,1; n = 12 h=3,2; n = 14 h=4,3,2,1; n = 16 h=5,4 (K0 = h, 10 classes; (14,1) has K0 = 1 with 5 keys at K=1). So H-026 as worded (K0 = min(h, n-4-h) for all lone-3 classes) is false at the self-mirror class h = n-4-h (even n); what the data support is: K0 = min(h, n-4-h) if h != n-4-h, and no split for K >= 1 if h = n-4-h (n = 10, 12, 14, and n = 16 for K >= 3). The refinement was made after seeing n = 14; n = 10, 12, 16 are post-hoc confirmations, not pre-registered. Refuted by: a self-mirror class with 2 image keys at some K >= 1, or a non-central class with K0 != min(h, n-4-h). Not claimed: a derivation; why the symmetric class is different (the lone row is its own mirror: head and tail deletion coincide).
(2) Power of the rivals. All 12 odd-n fits (n = 13, 15, 17 refit here with a fresh pass, n = 19, 21 recorded) fix h = the short side. On that reading R2 = min(h, n-5-h) and R3 = min(h, n-3-h) are algebraically equal to the law wherever h < n-4-h, so no existing or new non-central data can separate them; they split from the law only at the central even class, where all three fail the same way (observed 0). R1, R4 are refuted by 7 of 12 odd fits (non-central ones), R5 and the long-side reading of R3 by 12 of 12; on even non-central data (10 classes) R1, R4 fail on 6 and 10. So "the rival min(h, n-3-h) fails" was not a real rival; the law's discriminating content is the short-side reading. The planted K0-1 miss was caught on all 8 refits.
(3) Class level, n = 15 h = 3. The key class (K >= 2: 4597 rows; with the 4780 rows reached by forward rewrites) splits under the move rewrites (`lnaMoves.rewritesOf`, union of forward edges) into 1677 components. Every component has a single image key at K = 2 (373 components with >= 2 rows) and at K = 3 (117 such; 1213 rows; image keys 956 / 257). Shuffled labels would give 52 +- ~6 mixed components at K = 3, the data 0. My logged prediction ("not separated by move orbit") is therefore wrong at the component level: the split IS a function of the move component. Caveat: components of forward edges are sub-orbits (reverse edges and rows outside K >= 2 reached only backwards are missing), so "different orbits" is not shown, only "not one component under the available rewrites"; derived equivalence classes (q_S, Z[Phi]) not tested, so the key class could still be one derived class with a non-invariant deletion.

## Evidence

n = 14 table (image key counts per K; full text `experimentalist_n14_out.txt`): (5,5) 1091/486/222/71/33 for K = 1..5, always key 0; (4,6) K=1..4 two keys (1503/856, 661/412, 285/136, 122/65), K = 5, 6 one; (3,7) two keys K = 1..3 (2183/780, 956/265, 410/129), one from K = 4; (2,8) two keys K = 1, 2 (3216/522, 1412/257), one from K = 3; (1,9) five keys at K = 1, one from K = 2. The key-0 columns of n = 14 are the n = 15 h=3 key-0 columns shifted by one K (2183, 956, 410, 171, ...), as are the 780/265/129 of the n = 13 lone sequence: the table is the same lattice of numbers at every n (not explained). Other tables: `experimentalist_even_extra_out.txt` (n = 10, 12, 16), `experimentalist_power_out.txt`, `experimentalist_orbits_out.txt`.

| class | n=10 | n=12 | n=14 | n=16 |
|---|---|---|---|---|
| central (h = n-4-h) | 0 | 0 | 0 (pred. 5) | none K >= 3 |
| h = central-1 | 2 | 3 | 4 | 5 |
| next | 1 | 2 | 3 | 4 |

(entries K0 = last K with 2+ image keys; 0 = splits at K = 0 only.)

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/117/experimentalist_kpower.py table 14 1 5,4,3,2,1   # 46 s
timeout 10m .venv/bin/python workshop/rounds/117/experimentalist_kpower.py table 12 0 4,3,2        # 12 s (also n=10 K1=0 h 3,2,1; n=16 K1=3 h 6,5,4, 48 s)
timeout 10m .venv/bin/python workshop/rounds/117/experimentalist_kpower.py power                   # about 6 min
timeout 10m .venv/bin/python workshop/rounds/117/experimentalist_kpower.py orbits                  # about 70 s
```
The enumerator, Cartan builder and Coxeter polynomial are copied from `rounds/115/experimentalist_kloc.py`; the library cross-check (pwh.removeVertex, 0 mismatches) was done in E-304, not repeated here.

## Prior record

E-296, E-299, E-301, E-304 (odd n only; all ten fits odd), H-026, R-015 (digit-joined writer artefact). Even-n central failure, the odd/even parity of all fits, and the component-level result are not in `research/` (grep "self-mirror", "even n", "n = 14" in EXPERIMENTS.md: nothing for H-026).

## Code changed

None in the library. New: `workshop/rounds/117/experimentalist_kpower.py`. No tests needed.

## Next

Skeptic/theorist: the self-mirror exception should be derived (head deletion = tail deletion by the mirror symmetry of a (h,h) row, so images are mirror pairs and keys coincide?); H-026 should be reworded (referee the refinement; it was fitted after n = 14). Converted H-entry for cycle 2: "image keys of deletion split along move components of the source key class" (n = 15 h = 3; needs q_S / derived-class check and a full orbit, including reverse edges). Experimentalist next: tail end at even n; non-lone classes at even n.
