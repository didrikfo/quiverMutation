# Conference round 020 position statement

author: theorist · round: 020 · kind: position

## Most promising question

**Why does the shuttle pattern (34 @k ↔ 403 @(k-1)) force exactly one placement per split word?** The placement's right gap g is n-independent and word-determined (E-098 table), but the mechanism is not yet derived: the boundary shapes (`4y` only at right gap 1, `334, 357` only at right gap 0) are read from tables (n = 12..16), and the proof that "R keeps g" is by formula and terminal check, not universal. The shuttle appears after R stops, so the question is why short shapes sit at the orbit boundary only — i.e., why does the closed `444` orbit under the reduced move set have its boundary at exactly g = 0 or 1?

## Weakest claim in the workshop

The claim that lemma R preserves the right gap g (E-098 claim 3). This is checked by the end formula (offset = o_max - g) and by the list of terminals, but not proved for all words. The exception `3344` shows the rule needs care: R is valid only when the shortened interval does not swallow a left relation; my first script error filtered R-inert words wrongly. Also: the boundary-slot table (`4y in S only at left gap 0 or right gap 1`) is verified at n = 12..16 for letters ≤ 8, never derived.

## What I need from another persona

- **Experimentalist:** Test words outside S (`3334, 2455, 3335`) for their R-terminal and right gap at n = 12..16, and n = 17; also 4-letter words with a 4 and 4+ placements. This will show whether "g word-only" holds beyond the split words, or whether the rule breaks.
- **Skeptic or I:** A selectivity null test of "g word-only, independent of n" on a random sample of 5-letter and 6-letter words.

