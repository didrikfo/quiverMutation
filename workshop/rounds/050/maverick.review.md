# Review of workshop/rounds/050/maverick.md

referee: scholar · round: 050
verdict: minor revision

## Reproduction

Re-ran `maverick_sigpower.py 6 7 8 9 10` (11 s), `... 11` (28 s) and `... 9`. All table cells match: n = 9 9/9 NOT_QUIPU with pos <= n-2; n = 10 260 + 2 UNPLACED, 2 quipu-polynomial LNAs with a different signature (`34504030`, `50505000`, (8,2,0) vs (10,0,0)); n = 11 2631 + 16, 16 different, e.g. `334500030` (9,1,1) vs (10,0,1); F-010 pair both (8,1,0), Smith (1^8, 12). n = 6..8: 0 differences, 0 pos <= n-2.

## True?

Yes. Items 1-2 are correct: C^-1 + C^-T = C^-1 (C + C^T) C^-T is a congruence, and a derived equivalence acts on Cartan matrices by C -> P C P^T (Ladkani 3.13/3.15), so the signature is constant on a derived class. The check is eigenvalue-based with tolerance 1e-9; the zero-count is stable in the output (zero eigenvalues occur only at n = 11, as (10,0,1)). I found no case the author missed within the stated range.

Two small overreaches:
- Item 3: "the signature proves 'in no quipu class'" for the 16 is right only given that the Coxeter polynomial is the quipu's; it holds, because the signature is a class invariant and differs from that quipu's. It does not show the LNA is outside every quipu class with that polynomial if several quipus share it (the script compares with one `qsig[key]`, last writer wins). At n = 9 the F-010 pair shares a signature, but the code does not check the case where cospectral quipus have different signatures. State that.
- Item 5 "the UNPLACED are all on the non-quipu side" is just a restatement of item 3 for pos <= n-2.

## New?

Mostly not new.
- Items 1-2: E-154 already records that "the class carries one signature" and that these invariants "cannot separate anything the key does not". E-065's own text says the signature criterion cannot give "relations > cords" (a relation is a rank-2 perturbation, lowers `pos` by at most 1). The explicit Sylvester argument for the Euler form is not stated there, but the conclusion is. Cite E-065 and E-154.
- Item 5: E-065 already records the counts n = 9: 9, 10: 262, 11: 2647 (= 260 + 2, 2631 + 16). The author's "UNPLACED part is new" is the split of E-065's 262/2647 into NOT_QUIPU and UNPLACED, which is a bookkeeping refinement.
- Item 3 at n = 10: FINDINGS ~510 (Brüstle split of `T^10+T^9+T+1`) and F-048, which certifies `34504030`, `50505000` (and `45050400`) as non-piecewise-hereditary by periodic Coxeter + indefinite Euler form. F-048 fires on 619 new LNAs at n = 11; the submission did not check whether the 16 lie among those 638. If they do, item 3 at n = 11 is known as exclusion, and the only new content is "the Coxeter polynomial coincides with a quipu's" for them.
- The n = 11 count of 16 as "LNAs with a quipu's polynomial but different signature": nothing found by grep for `334500030`, "different signature". Plausibly new, small.
- Item 4 (F-010 pair indistinguishable by signature and Smith form): E-154 and FINDINGS ~336 hold finer congruence data. Not checked against Smith form of C+C^T alone; the Smith result may be new but is one line.

## Evidenced?

Mostly. The table gives counts per n and the commands. Missing: a statement of how the UNPLACED/NOT_QUIPU status is defined (the author relies on `ct.lnaStatus`); the polynomial-group split is reported as counts only, with no list of the 16 names. The n = 6..8 row merges three lengths. The title says "beats the Coxeter polynomial", but polynomial vs signature is compared only for LNAs that share a polynomial with a quipu; F-048 (periodic + indefinite) is a stronger separator at n = 10, 11 and should be the named comparator.

## Scope

Title fits what was checked (n = 6..11, signature only). "H-017 stays OPEN" is a conclusion not drawn from computation but from the congruence argument, which is sound. Narrowed wording for the title: "...is a derived-class invariant (E-154, E-065), hence blind to (cords, relations); as a separator from a quipu's Coxeter polynomial it acts on 2 LNAs at n = 10 (known, F-048) and 16 at n = 11."

## Required for acceptance

1. Add E-065 and E-154 to the prior record for items 1-2 and 5; mark the counts 2647 = 2631 + 16 as already in E-065.
2. Check whether the 16 n = 11 LNAs are among F-048's 638 periodic-and-indefinite rows (one call to `periodic3`-style test on those 16) and say so; this decides whether item 3 is new at n = 11.
3. In item 3, either handle the case of several quipus sharing a polynomial (compare against every sharer) or state that it is assumed unique. At n = 9 the F-010 pair is exactly such a polynomial; at n = 10, 11 verify uniqueness.
4. List the 16 names (or a pointer to the script output) so the claim can be audited.
