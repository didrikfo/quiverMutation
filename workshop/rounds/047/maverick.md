# Proposal (revised): T1/T2 (H-021' pairing cores, centre s = n - k(c)) have no open question with a cheap census; list the remaining items in the ledger

author: maverick · round: 047 · kind: proposal
thread: T1, T2 · bears on: H-021, H-020, F-053, F-032
scope: record read only (H-021 text and status line, E-055, E-058, E-061..E-064, E-066, E-067, E-069, E-070, E-072, E-073, plus E-076..E-100 headers via grep, round 2 of this note). No script was run this round and no new data; every count below is quoted from `research/`. The referee's re-run of E-067 (n = 14, closed, `held==pred`) is theirs, not mine.

## Response to referee

Verdict was minor revision; all six required items done, wording below.

1. **E-076..E-100 listed (done).** Relevant to T1/T2: E-076, E-079, E-082, E-085 (key-coarser lists A/B; `5046`/`5056` have two closed orbits of different sizes at n = 17; no invariant separating P from Q); E-077, E-081 (the `444` orbit; data cannot separate "has a 4" from "the `444 -> 34` route"); E-090 (lemma R; the end link `3x@0 -> 33(x-1)@0` seen on one labelled path, which is the part of E-067's undone "end link"); E-093, E-098, E-100 (split words with exactly one placement in the `444` orbit; `3344` exception; E-100: not enriched over a pooled-rate null). Residues they leave open for T1/T2: (a) the `444` orbit and what it contains; (b) `5046`/`5056`; (c) E-090's end link versus E-067's upper bound. I read headers and first lines only; I did not check these against E-067's table. In the referee's re-run, `33x` (x = 7, 8, high offsets) sits in an orbit of size 3767 = the n = 14 `444` orbit (E-081); not checked by me as the same orbit, but it bears on (c) and on `44x`.
2. **Reason 2 fixed (done).** The `k, d` fits are E-062/E-064 (12 cores, n = 13..16) and E-055/E-058 (109 pairing cores at n = 13). E-066 is the orbit+mirror vs key comparison, not a `k(c)` source; dropped. E-070 is a table of pair sums for `34x`, a restatement of the data (its own words), not a fit; dropped as a fit count. E-064 and E-079 also contain fits. The count "three times" is withdrawn; I do not have a clean count. E-069's chance-level result is a reason to distrust fits at small |R|, not a reason to stop; the interior formula and the n = 15/16 fits survive its null (stated in the table). Reason 2 now argues only for the held-out test below.
3. **E-072 and T3/T8 (done).** E-072 says "pair at even n, mirror-join at odd n" is not supported, and E-066 does not cover the 7 cores of E-061 (E-072: disjoint). So "odd-n loss of the 7 and the n = 16 loss are one phenomenon" is **inference**, marked as such; only the n = 16 `20300` pairs being orbit + mirror is recorded (E-066, E-076). "Close T3/T8 together" is dropped: STATE.md has T3/T8 done at the catalogue level with an open item (why no shift by 1; no invariant, E-082). T3/T8 stays separate.
4. **`44x` qualified (done).** The criterion (seed `aaa` collapses to `34` only for a = 4) has one positive datum (E-073 Limits) and a broader criterion failed (E-073 reproduction, `theorist_translator.py`; the referee counts 10 of 24 cells, which I did not recount). E-077/E-081 say the data cannot tell the `444 -> 34` route from "has a 4". "Matching a criterion" is now "consistent with a criterion that has one positive datum".
5. **Reason 1 (done, marked conjecture).** I have no citation that T10/T5 bear on the walk's pairing theory. STATE.md describes T10/T5 as key guard and J != 0 steps, not walk-orbit theory. Reason 1 is cut to: the pairing is a property of the forward reduced walk (E-067 Limits), not shown to be a property of the derived class. Dependency on T10/T5 is conjecture and removed from "new".
6. **Ledger wording (done).** STATE.md already lists T1 as dormant since round 021. This proposal changes only ledger wording (add the remaining-items line). The real change is H-021's status line (still OPEN, runs to round 021), the chair's call: iff REFUTED, H-021' SUPPORTED-AS-DESCRIPTION.

Scope sentence narrowed as the referee proposed ("not provable by more censuses" was an opinion; withdrawn).

## Claim

T1/T2: no open question with a cheap census. Settled: the literal H-021
("exactly when the orbit holds a mirror") is false on every reading (E-058);
the restated H-021' (for cores that pair, reduced-walk offset orbits are the
pairs `{o, s - o}`, `s = n - k(c)`, `n`-independent `k` and `d`) is a
description that holds on the cores tested (12 cores at n = 13..16, E-062/E-064;
the 7 cores at even n, E-061) and survives a null test for the interior formula
(E-069). Derived: for `33x` only, `k = 2x`, `d = x - 3`, lower bound, via the
double-mutation drift plus the self-dual seed `333` (E-067); the same argument
gives `aax` with `k = 2x + 3 - a` for a = 3, 5, 6 (E-073; `77x` at n = 15 in
the H-021 header). Remaining, to list in the ledger with owners: `k(c)` for
general c (`34x`: `k = x + 3` observed only, E-070); the upper bound of E-067
and the end link (E-090); the `44x`/`34x`/`45x` mechanism (E-073, E-077, E-081);
the odd-n and n = 16 fit losses (E-061, E-064; inference only for sameness);
`5046`/`5056` (E-082, E-085). Reopening test: a rule for `k(c)` that predicts a
held-out core. H-021 itself stays OPEN until the chair rewrites its header.
Speculation level: `supported` for the quoted counts, `idea` for the closure.

## Evidence (what the record says, by question)

| question | state | where |
|---|---|---|
| Is the mirror clause an iff? | no, on all three readings (loose: mirror => P fails 20; strict: P => strict fails 108 of 109, strict => P fails 7) | E-058, E-052 |
| Does the pairing exist for cores that pair? | yes: 12/12 at n = 13..15, 9/12 at 16 (3 lose the fit through an unmerged equal-size middle pair) | E-062, E-064 |
| Is the centre formula `first + last outside offset` real? | interior blocks 13/13 (n = 13), null-tested, chance 4.1/13; interior formula and n = 15/16 fits survive the null | E-062, E-069 |
| Is `k` n-independent? | yes for the 12 cores over 13..16; 7 cores flip with parity of n | E-064, E-061 |
| Is `k(33x) = 2x` explained? | lower bound derived; upper bound only computed, 135/135 placements at n = 14..16 (n = 14 re-run by referee) | E-067 |
| Does the drift extend? | `aax` for a = 3, 5, 6 (and 7 at n = 15) yes; `44x` joins the big orbit, criterion has one positive datum, and E-077/E-081 cannot separate routes; `34x` (`k = x + 3`, data only), `45x` no mechanism | E-073, E-077, E-081, E-070 |
| What did rounds 010..021 add? | `444` orbit; two-orbit `5046`/`5056`; end link; split words | E-076..E-100 |

## Why more rounds on T1/T2 look unattractive (argued, not evidenced)

1. The remaining questions are statements about the forward reduced-walk
   orbit (E-067 limits: "not derived equivalence"), a property of a rewrite
   system. Whether that bears on the derived class, and on T10/T5, is
   conjecture; no citation.
2. The `k, d` fits (E-055/E-058, E-062/E-064) gave a rule only for the drift
   families. E-069 warns that the most informative fits are chance-level at
   |R| <= 4 and 39 of 109 are vacuous; a further fit needs a held-out core.
3. The drift is exhausted by E-073 as far as it goes: it closes exactly when
   the seed `aaa` does not reach `34`; why `34` alone reaches the slider
   (`3333 <-> 3403`) is unexplained (E-073 Limits).

## Reproduction

Nothing run by me. Quoted rows: `grep -n "^## E-0\(53\|59\|60\|62\|65\|67\|68\|70\|71\|74\|75\|77\|79\|80\|83\|88\|91\|96\|98\)" research/EXPERIMENTS.md`.

## Prior record

Everything above is in `research/` except the argued reasons 1-3. New: only
the proposal to relist T1/T2 remaining items and the 10/24 and `444`-orbit
qualifications gathered from the record. I did not check the H-021 header
rewrite (chair, 048).

## Code changed

None.

## Next

- Chair, round 048: relist T1/T2 with the remaining items above, or name the
  held-out `k(c)` prediction (cores not used to fit it); rewrite H-021's header.
- Theorist: state E-067's upper bound and E-090's end link as one gap in
  H-021's text, so the ledger lists it; check whether `33x` at x = 7, 8 high
  offsets lies in the `444` orbit (size 3767, n = 14).
- Maverick: back to S-1 n = 15 K0 = 5 sizing.
