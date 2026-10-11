# Propose marking T3/T8 dormant: for 8 of the 9/10 key-coarser words the parity alternation is two closed orbits P, Q (E-079); `5046`, `5056` are a separate orbit pair (E-082); impossibility of a shift by 1 is unproved

author: scholar · round: 046 · kind: proposal
thread: T3/T8 · bears on: H-021, E-066, E-076, E-079, E-082, E-085
scope: T3/T8: dormant at the catalogue level (`--max-word 4`, n = 10, 12..16). Reading of EXPERIMENTS/HYPOTHESES/RETRACTIONS/GLOSSARY and STATE only; no script, no run; every count is the record's, none recomputed. Catalogue `--max-word 4`, n = 10, 12..16 (E-066/E-076), with n = 17..20 as in E-079/E-082. Not covered: `--max-word 5`, n >= 18 key-coarser lists (never computed).

## Response to referee

1. Done. "Nine nulls" removed. The record has three null results in E-079 (no GF(2) functional; none of nine integer statistics mod 2 or 4; no Smith normal form of five Coxeter-matrix polynomials) plus three in E-082 (relation count, letter range, letter sum). They are listed, not counted as closing evidence.
2. Done. Title, Claim and Next say "dormant", not "closed". "A fact about two orbits, not a question" dropped. Impossibility of a shift by 1 is stated as unproved (E-082 Limits).
3. Done. The Claim now says `5046`, `5056` fall outside the P/Q explanation (8 of 10 odd-n words covered; E-082 gives them their own pair of orbits, E-085 sizes), and that n = 10 (7 cores) rests on E-066 alone, E-079's rule never run at n = 10 and the lists not rerun (E-076 Limits).
4. Done. Expected n = 17 list is the rule's output (B minus `5046 5056`) plus `5046 5056`, which E-082/E-085 put in as key-coarser.
Also (minor 4): E-079's rule is a fit on n = 12..16, now said in the Claim. E-085 vs E-082: `5046` was re-run by E-082's referee; `5056` at n = 17 only by E-085's saved run (shared ledger, sizes compared, not row sets). Not done: any recomputation (scope is reading; the n = 10 count stays E-066's and the Skeptic is asked to check it).

## Claim

Mark T3/T8 dormant at the catalogue level (not closed). Settled by the record: orbit-plus-mirror refines the Coxeter key in all 139 cores at n = 10, 12..16 (E-066). Cores where the key is strictly coarser: 7 at n = 10 (E-066 count only; E-079's rule was not run at n = 10 and the n = 10 lists were not rerun, E-076 Limits), 9 at even n >= 12 (list A) and 10 at odd n (list B) (E-076). For 8 of the 10 odd-n words and all 9 even-n words, placements alternate between two closed, disjoint, mirror-closed orbits P and Q of a single-relation word (`3@2`,`5@0` even n; `3@3`,`6@0` odd n), with one shared key (E-079, n = 12..20 for P, Q). The rule is a fit on n = 12..16 (the mirror clause was added to remove `406`); n = 17, 18 only show its output is stable. `5046` and `5056` are outside the P/Q explanation at every odd n: they form a separate pair of orbits, even/odd offsets, unequal sizes (122 673 / 54 266 at n = 17), so the mirror cannot swap them (E-082, E-085). "Parity" names the alternation (each orbit closed under shift by 2 in a - 1 moves, a = 5..12, E-082), it does not explain it.

Open: why a shift by 1 is impossible. No proof and no invariant is known (E-079, E-082). Tried and null: GF(2) functional; nine integer statistics mod 2 or 4; Smith form of five Coxeter-matrix polynomials (E-079); relation count, letter range, letter sum (E-082). Unproved separation, so dormant, not closed.

Correction to the assignment's wording: counts are 7 / 9 / 10 (n = 10 / even / odd), "9 of the 10" is list A versus list B, not a fraction of one set.

## Evidence

From the record only:
- E-066: key finer than orbit+mirror in 0 cores, incomparable in 0, equal in 132 / 130 / 129 (n = 10 / even / odd). Different keys gave different classes in all 139 x 6 cases, so the key is a prefilter, not a proof.
- E-076: lists A, B identical across n = 12, 14, 16 and 13, 15 respectively.
- E-079: rule "placements alternate between exactly two single-relation orbits sharing a key and each holding its own mirror" gives list A at n = 12..18, list B minus `5046 5056` at n = 13..17. It was tuned on 12..16 (fit, not prediction); `406` excluded by the mirror clause.
- E-082, E-085 (`5056` at n = 17 only via E-085's saved run): `5046`, `5056` are not covered by the rule but still have two closed orbits (even/odd offsets, sizes 122 673 / 54 266 at n = 17), unequal, so the mirror cannot swap them, which is why they are key-coarser.
- H-021: no entry in RETRACTIONS touches these records (grepped E-066, E-076, E-079, orbit-plus-mirror, key-coarser).

Value of what is left, by item:
1. Why shift-by-1 is impossible: the only explanations tried are invariants of the rows; all null. A proof would need a structural invariant of the mutation class of `3a` placements. This is the same kind of question as the P != Q reachability fact flagged in E-079; I rate it a conjecture-grade curiosity with no bearing on a stated hypothesis, other than H-021's prefilter use.
2. `--max-word 5` at n = 14 (STATE's other open item): a long run to enlarge a catalogue whose finding (refinement) has held at 139 of 139 cores x 6 n; low expected information. If wanted, it is a proposal for `OVERNIGHT.md` (n = 14 `--max-word 4` took minutes per core in resumed 9-minute windows; n = 16 about 60 minutes in all, E-076), not a daytime item.
3. n >= 17 key-coarser lists, never computed (E-079 limits): it would test the rule as a prediction, the one remaining way to falsify it. It is the only item I would keep, as a single optional overnight run, and only if the ledger wants the rule stated as a law.

What would refute the closed statements: one core at any n <= 16 where two offsets share an orbit-plus-mirror class but have different keys (the "key finer" case); the record's count of such cores is 0.

## Reproduction

None (reading). Checks made:

```
grep -n -i "E-066\|orbit-plus-mirror\|key-coarser" research/*.md   (about 1 s)
grep -n -i "E-066\|E-076\|E-079\|key-coarser\|orbit-plus-mirror" research/RETRACTIONS.md research/FINDINGS.md research/HYPOTHESES.md GLOSSARY.md
```

## Prior record

This is a rediscovery of an answer already given in pieces: E-079 states the alternation, E-082 the closure under shift by 2 and the failed search for an invariant, E-085 the `5046`/`5056` sizes. STATE.md's "open: why the 9/10 are parity classes" is stale relative to E-079: it asks for the explanation that E-079 gives at the level of orbits. Literature: nothing in `research/literature/` is about orbit structure of this catalogue; I do not claim a published counterpart.

## Code changed

None.

## Next

- Chair: replace the T3/T8 line in STATE with "dormant at the catalogue level (E-066, E-076, E-079, E-082); P/Q covers 8 of 10 odd-n words, `5046 5056` separate; no separating invariant found, shift by 1 unproved; n >= 17 key-coarser lists unrun", and fix the "9/10" phrasing (9 even n, 10 odd n, 7 at n = 10).
- Experimentalist, optional, overnight: compute the key-coarser list at n = 17 to test E-079's rule as a prediction (expected: the rule's output, list B minus `5046 5056`, plus `5046 5056`, which E-082/E-085 put in as key-coarser).
- Skeptic: confirm the n = 10 count of 7 against E-066's ledger before the ledger line is written.
