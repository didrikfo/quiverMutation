# Review of workshop/rounds/007/experimentalist.md

referee: skeptic · round: 007
verdict: minor revision

## Reproduction

- `experimentalist_kd.py 14 34 3-9`: 55 s. Same classes, sizes, sums and k as the table (344..349 at n=14). Only the "closed" timings differ.
- Read the stored n=15 and n=17 outputs (`_n15.txt`, `_n17.txt`, `_n17_8/9`) and checked every singleton by hand. Each is the centre or an equal-size mirror pair: 344@15 (2,6 = 64,64; 4 centre), 344@17 (2,8 = 78; 4,6 = 86; 5 centre), 348@17 (2,4 = 9622; 3 centre), 349@17 (1,4 = 20769). The table matches the files.
- `experimentalist_core4046.py 13 4046`: {0,2}(63) {1}(104) {3}(2116), as stated. The E-060/STATE discrepancy is real; the author's 4046 result reproduces.
- Not re-run: n=16 and n=17 (6 to 8 min each), 45x, 4046 at n=14..16.

## True?

The data hold. Three wording problems.

1. **`k = x + 3` restates `s = hi`.** Since `hi = n - x - 3`, `k = n - s` equals `x + 3` exactly when `s = hi`. Every cell where 0 and hi share an orbit gives `s = hi` by construction. The n-independence of k therefore adds nothing beyond "0 and hi are in the same orbit". "Power: a wrong k would have to be wrong at each n" overstates this. The real content is that the orbit classes are the pairs `o <-> hi - o`, which is stronger, but the k table does not show it.
2. **`d = 0` is used two ways.** The headline says "d = 0 (no offset above s)", which is true by `s = hi`. The table reports raw `d` of 0 to 5, and the headline's `d_eff = 0` is a different quantity. The `d_eff` rule "equal-size singletons are a mirror pair" is a size coincidence, not an orbit identity (see Evidenced).
3. **The 4046 reflection claim rests on equal sizes.** `s = n - 11`, `d = 1` for 4046 at 12..16 pairs singletons by equal size (388/388, 149/149, 530/530, 534/534), which the author flags as unchecked. The "lone big orbit = hi" part is fine.

The `458@15` and `459@14,16` "k is an artefact" remark is a correct self-correction.

Not checked: x >= 10, n >= 18, and whether 346 being one orbit is a size effect or a join with `333@0`. The author says "the `333@0` orbit of E-065" without a check.

## New?

- `grep -E "34x|k\(34|x \+ 3|4046"` over FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS and literature found no record of `k(34x) = x + 3`.
- The E-056, E-060 and E-061 lines and STATE T2 (line 14) are the context. E-065, cited by the author, was not located in the files I grepped, only via the author's and STATE's citations.
- The 4046 discrepancy with STATE/E-060 is a new finding, and a useful one.
- I did not check the `literature/` hits for anything relevant to a reflection formula; the only `34x` hit there is a polynomial coefficient.

## Evidenced?

Mostly. Per-cell sizes are in the stored `.txt` files and the table is reproducible. Gaps:

- Mirror pairing of equal-size singletons is not verified as same-orbit. The author defers this to the skeptic (Next). It affects 344@15/17, 348@16/17, 349@16/17 and all of 4046. Equal size of 9622/9622 or 20769/20769 is what a true mirror symmetry would give, but it is also what two distinct isomorphic orbits would give, so the check is needed.
- No null was run. A baseline such as `33x` or random `s` with the same orbit sizes would show whether "some `s` pairs these classes" is easy to achieve. Only power-by-cell-count is given.
- Low-power cells (349@14, 348@14, 349@15) are flagged honestly.

## Required for acceptance

1. Reword the headline so `k = x + 3` is stated as equivalent to `s = hi`, and say the evidence is the pair structure `o <-> hi - o`. Drop "a wrong k would have to be wrong at each n".
2. Keep `d` (raw) and `d_eff` distinct in the headline.
3. Either run the mirror-join check for the equal-size singleton pairs (344, 348, 349 at n=15..17; 4046 at 14..16), or state in the Claim that "pairing" for singletons is by size only.
4. State the 4046 reflection claim as size-paired only, as in point 3.
5. Cite where the `333@0` join of 346 was verified, or mark it as unverified.
