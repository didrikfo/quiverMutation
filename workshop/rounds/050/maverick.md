# H-017's Euler signature is a function of the Cartan matrix, hence a derived-class invariant (E-152, E-063) that is blind to (cords, relations); as a separator from a quipu's Coxeter polynomial it acts on 2 LNAs at n = 10 (known, F-048) and 16 at n = 11 (not F-048 rows)

author: maverick · round: 050 · kind: negative
thread: T6 · bears on: H-017, E-063, E-152, F-045, F-048
scope: n = 6..11 all LNAs (42..16796), hereditary quipus n = 6..11, one cospectral quipu pair at n = 9. Signature of C+C^T by floating eigenvalues (tolerance 1e-9). Separation is tested against every quipu sharing the LNA's Coxeter key (checked at n = 9, 10, 11). No search run, no n = 12.

## Response to referee

1. **Prior record (done).** Items 1-2 restate what E-152 ("one signature per class", invariants "cannot separate anything the key does not") and E-063 ("signature cannot give relations > cords") already conclude; only the explicit Sylvester step is mine. The counts 262 = 260 + 2 and 2647 = 2631 + 16 are E-063's; my only addition is the NOT_QUIPU / UNPLACED split, a bookkeeping refinement. Item 5's "UNPLACED on the non-quipu side" is item 3 restated for `pos <= n-2`. Prior record below updated.
2. **F-048 on the 16 (done, `maverick_sixteen.py 11`).** None of the 16 is an F-048 row (F-048 fires on 638 at n = 11; Phi periodic and Euler form indefinite, tested on the integer Coxeter matrix `-C^T C^-1`, power up to 600). So at n = 11 the signature excludes 16 LNAs from their quipu's class that F-048 does not certify, and `lnaStatus` leaves UNPLACED. At n = 10 both of the 2 are among F-048's 3 rows, so known. F-048 is the stronger separator at n = 10, 11 overall (3 and 638 rows vs 2 and 16); the signature is a different, not nested, one. Caveat: "not in the class of that quipu" is what the signature proves; I did not test other certificates (A9/A13 criteria, E-040) on the 16, but they are UNPLACED, i.e. the repo's existing certificates did not fire.
3. **Several quipus sharing a polynomial (done).** The first script kept one signature per key (last writer). Now every sharer is used: n = 9, 10, 11 have 1, 2, 4 keys carrying more than one quipu, and in every such key all quipus have the same signature (0 keys with differing signatures). So "unique signature per key" holds and the count is unchanged: signature differs from all sharers for 0, 2, 16 LNAs at n = 9, 10, 11, and from some sharer for the same 0, 2, 16.
4. **The 16 names (done):** `334500030 345500030 350500030 303345000 303455000 303505000 304444400 304550400 305504030 305040330 305050030 444440030 455040030 550403030 504033030 505003030`; all signature (9,1,1), all share the key of `P^(1,1)_(1,5,1)` (10,0,1). The n = 10 pair is `34504030`, `50505000`.
5. Status definition: `ct.lnaStatus` gives QUIPU if reached from a seed by known moves, NOT_QUIPU if no quipu has the LNA's Coxeter polynomial, UNPLACED otherwise. Title now names F-048 as the comparator. Not done: the n = 6..8 row stays merged (all zeros, same as before); Smith form of C+C^T for the F-010 pair vs E-152 data not compared (one line, unchanged).

## Claim

1. **Determined by the Cartan matrix: yes.** E-063's signature is that of `C + C^T`. The Euler form is `C^-1 + C^-T = C^-1 (C + C^T) C^-T`, congruent to it (Sylvester), so the signature is a function of `C`, and under derived equivalence the Cartan matrix changes by a congruence (E-063 cites Ladkani 3.15). It is therefore constant on every derived class, and so constant on the members of a class that H-017 counts (cords, relations) over.
2. **Consequence for H-017 itself.** The signature cannot see which member of a class is reached, so it can neither confirm nor refute "relations > cords". All it says about H-017 is E-063's weak bound (a quipu in a `pos <= n-2` class cannot exist; a relation moves `pos` by at most 1). H-017 stays OPEN and its status should not change on this evidence. What a "power check" of the signature can show is only its power as a *class separator*; it cannot show power against the search question.
3. **Power as a class separator (what it does separate).** Among LNAs whose Coxeter polynomial equals a quipu's, the signature differs from that quipu's for 0 LNAs at n = 6..9, **2 at n = 10 (`34504030`, `50505000`: (8,2,0) vs (10,0,0))** and **16 at n = 11** (all 16 UNPLACED LNAs; e.g. `334500030` (9,1,1) vs (10,0,1)). For these the polynomial is silent and the signature proves "not in that quipu's class" (compared against every quipu with that key; keys with several quipus all have one signature). The n = 10 pair is already certified by F-048 and the Brustle split (FINDINGS ~510): rediscovery. The 16 at n = 11 are not F-048 rows; their names are in the response above.
4. **Where it fails to separate (negative control).** The cospectral quipu pair of F-010 at n = 9 (`P^(1,4)_(1,0,1)`, `P^(1,2)_(1,1,2)`) has identical signature (8,1,0) and identical Smith form of `C+C^T`; they are non-equivalent per F-010, which finer congruence data separates (FINDINGS ~336), not the signature. At n = 6..9 no polynomial group of LNAs is split by the signature; at n = 10, 11 one group each.
5. **The iff of E-063 reproduces from the tables:** `pos <= n-2` holds exactly for status NOT_QUIPU and UNPLACED (n = 9: 9 of 9; n = 10: 260 + 2; n = 11: 2631 + 16) and for no QUIPU-status LNA; all hereditary quipus have `pos >= n-1` (n = 6..11). So at n = 9..11 the criterion is exactly "not provably quipu". Not new: E-063 records 262 and 2647; the split 260 + 2 and 2631 + 16 is bookkeeping.

Would refute: an LNA with `pos >= n-1` that is derived equivalent to no hereditary algebra of quipu type shows nothing against this; an LNA with `pos <= n-2` found in a quipu class (n <= 12) would refute the converse direction.

## Evidence

| n | LNAs | QUIPU, pos<=n-2 | NOT_QUIPU, pos<=n-2 | UNPLACED, pos<=n-2 | LNAs with quipu poly but different signature | poly groups split |
|---|---|---|---|---|---|---|
| 6-8 | 42/132/429 | 0 | 0 | 0 | 0 | 0 |
| 9 | 1430 | 0 | 9 of 9 | 0 | 0 | 0 |
| 10 | 4862 | 0 of 4600 | 260 of 260 | 2 of 2 | 2 | 1 |
| 11 | 16796 | 0 of 14149 | 2631 of 2631 | 16 of 16 | 16 | 1 |

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/050/maverick_sigpower.py 6 7 8 9 10   # about 1-2 min
timeout 10m .venv/bin/python workshop/rounds/050/maverick_sigpower.py 11           # under 10 min
timeout 10m .venv/bin/python workshop/rounds/050/maverick_sigpower.py 9            # also prints the F-010 pair
timeout 10m .venv/bin/python workshop/rounds/050/maverick_sixteen.py 9 10 11       # names, sharers, F-048 test (about 1.5 min)
```

## Prior record

E-063 (signature criterion, n = 8..11; counts 262, 2647), E-152 (one signature per class), F-048 (periodic + indefinite, 3 and 638 rows), E-069/E-072/E-073 (search controls), FINDINGS ~510 and ~1108 (the n = 10 pair and the polynomial/Euler-form separation), F-010 (cospectral pair). Items 1-2 are the answer to the assignment and I did not find the "congruent, hence constant on a class, hence blind to (cords, relations)" argument stated; it is short and elementary. Item 3 at n = 10 is known (F-048); at n = 11 the 16 are not F-048 rows and are plausibly new as a list. Items 1-2 are E-152/E-063 plus an explicit Sylvester step. Item 5 restates E-063.

## Code changed

None (new scripts `workshop/rounds/050/maverick_sigpower.py`, `maverick_sixteen.py` only; no tests touched).

## Next

- Ledger: H-017 stays OPEN; T6 note should say "signature is a class invariant, cannot probe (cords, relations); the open item is the search at n = 9 depth 7 (about 27 min per candidate, OVERNIGHT), not a better invariant." (chair)
- A positive control for the *search* (a monomial cord member at n = 8) is still missing (E-092); that is the actual power question and the signature cannot supply it.
- Whether the 16 n = 11 UNPLACED LNAs are the 2 n = 10 ones extended (vertex deletion/addition, S-1 kind): scholar or maverick.
