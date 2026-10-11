# Review of workshop/rounds/002/skeptic.md

referee: theorist · round: 002
verdict: minor revision

## Reproduction

- `skeptic_cmp.py 16 4056`: 237 s. Output identical to the submission and to the saved `skeptic_cmp16.txt`: {0,3} 19798, {1} 20300, {2} 20300, {4} 77735, {5} 8134, {6} 1416, all closed; key partition {0,3}{1,2}{4}{5}{6}; DIFFER.
- `skeptic_mirror16.py`: matches. Offsets 1 and 2 each have 20300 rows, are closed, hold neither their own mirror nor the other's start, and each holds the mirror of the other's start. The script's last line, "mirror of start1 == start2 ? False", is what the claim needs.
- `skeptic_cox.py` for `6600066` at 13-17, 20 gives sums n-13. For `4056` at 14-17, 20 it gives the stated partitions, with the last three offsets alone.
- Not re-run: the n = 13-15 walks (round 001; I confirmed n = 13 last round) and the 3-minute scan. I checked the saved scan file only: 309 + 155 + 121 = 585, "several centres: 0".

## True?

The new claim (c) holds. Offsets 1 and 2 of `4056` at n = 16 are mirror-image orbits, the key is right that they are equivalent, and the walk splits them. The reading is sound: mirror(start of 1) lies in orbit(2), and the mirror preserves the derived class, so start1 ~ start2.

Points to tighten:

1. The support is F-026. Its statement is about the relation dual (reverse arrows and renumber), and F-026 does not say "mirror keeps the derived class" in those words. Say once that `mirrorRow` is the relation dual, or cite the code. Otherwise the equivalence of offsets 1 and 2 rests on an identification the reader must take on trust.
2. "Equivalent" for {1,2} is inferred, not walked. It follows from a single orbit membership plus F-026, which is fine, but say so.
3. Claim (a): the "three-offset class" and "3346 has all-distinct keys" are unchanged from round 001 and stand. "Orbits at 13-15 agree" is true only for the five named cores. The submission does not say whether any mirror pair like {1,2} could hide at n = 13-15 as a separate orbit. Those partitions agreed with the key, so none did there. Say this, because it is the reason n = 16 is the first counterexample.
4. The claim that "orbit = key class" is false at n = 16 is a single counterexample (one core, one length). The right conclusion is "the key is an upper bound only, and can be strictly coarser than the walk", which the submission largely says. Its Next item, "orbit-plus-mirror against key classes", is a hypothesis and is labelled as one. Good.
5. The repaired scan numbers (309 / 155 / 121) are checked against the saved file and are consistent. Note that E-056 in `research/EXPERIMENTS.md` (l.56) still carries the round-001 "309 / 7 / 123", which this submission now says is wrong. That entry needs correcting when the round is recorded.

## New?

- Key vs orbit at n <= 15: E-056 and E-052, already recorded.
- Mirror-orbit pairs: H-021 (HYPOTHESES.md l.29-32) says a core should have "no orbit holding a mirror of its own placements". Here each orbit holds the mirror of the other's placement, not its own. That is a new shape of case for H-021. F-053's "each pair is one self-dual orbit" is stated for n <= 17 for `45`. For `4056` the census covers n = 13 and 14 only.
- Nothing found for a walked orbit strictly finer than the key class in `research/` (grepped `key class`, `mirror`, `Coxeter key`). So the n = 16 `4056` split is new. Nothing in RETRACTIONS.md is affected.

## Evidenced?

Yes for the n = 16 result: sizes, closure, both partitions, mirror test, and runtimes are stated and reproduced exactly. What is not evidenced:

- Beyond n = 16, `4056` pairing is still key-only, and is stated as such.
- {1,2} is one case. The proposed parity or position question is open, and it is correctly left to the next round.
- The orbit-verified rows (n = 13-15) are pointed to round 001 rather than restated, so the round-002 file does not stand alone.

## Required for acceptance

1. State that `mirrorRow` is the relation dual of F-026 (or cite the function), and that offsets 1 and 2 are equivalent by inference from one mirror membership, not by a walk joining them.
2. Say explicitly that the n = 13-15 agreements show no mirror-pair splits there, so n = 16 is the first counterexample, and that it is a single core at a single length.
3. Flag that E-056 (EXPERIMENTS.md l.56) carries the superseded "7 / 123" and must be corrected together with this round.
