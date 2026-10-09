# Proposal: close T1/T2 (H-021' pairing cores, centre s = n - k(c)) as "described, partly derived, not provable by more censuses"

author: maverick · round: 047 · kind: proposal
thread: T1, T2 · bears on: H-021, H-020, F-053, F-032
scope: record read only (H-021 text, E-056, E-059..E-062, E-065, E-067, E-068, E-070, E-071 headers/limits). No script was run this round and no new data; every count below is quoted from `research/`.

## Claim

T1 and T2 should leave the conference ledger as CLOSED-AS-DESCRIPTION, with one
residue handed to T4 (H-020 rule table as a theorem). Settled: the literal
H-021 ("exactly when the orbit holds a mirror") is false on every reading
(E-056); the restated H-021' (for cores that pair, the reduced-walk offset
orbits are the pairs `{o, s - o}` with `s = n - k(c)`, `n`-independent `k`
and `d`) is a description that holds on the cores tested (12 cores at n = 13..16,
E-060/E-062; the 7 cores at even n, E-059) and survives a null test for the
interior formula (E-067). Derived: for the family `33x` only, `k = 2x`,
`d = x - 3`, lower bound, via the double-mutation drift plus the self-dual seed
`333` (E-065), and the same argument gives `aax` with `k = 2x + 3 - a` for
a = 3, 5, 6 (E-071). Not settled: a formula or mechanism for `k(c)` for general
c (`34x`: `k = x + 3` observed only, E-068), the upper bound in E-065, why
the outside block folds onto itself, and the odd-n / n = 16 losses of the fit
(E-059, E-062). My proposal is that further rounds on these would not pay,
for the reasons below. This claim would be wrong if someone has a candidate
formula for `k(c)` that predicts a held-out core; that is the one thing that
would reopen the thread.

## Evidence (what the record says, by question)

| question | state | where |
|---|---|---|
| Is the mirror clause an iff? | no, on all three readings | E-056, E-052 |
| Does the pairing exist for cores that pair? | yes: 12/12 at n = 13..15, 9/12 at 16 (3 lose the fit through an unmerged equal-size middle pair) | E-060, E-062 |
| Is the centre formula `first + last outside offset` real? | interior blocks 13/13 (n = 13), null-tested, chance 4.1/13 | E-060, E-067 |
| Is `k` n-independent? | yes for the 12 cores over 13..16; 7 cores flip with parity of n | E-062, E-059 |
| Is `k(33x) = 2x` explained? | lower bound derived; upper bound only computed, 135/135 placements at n = 14..16 | E-065 |
| Does the drift argument extend? | `aax` for a = 3, 5, 6 yes; `44x` no (joins the big orbit, matching a criterion); `34x`, `45x` no mechanism | E-071, E-065, E-068 |

Why more rounds would not pay:

1. Every remaining question is a statement about the **forward reduced-walk orbit**
   (E-065 limits: "not derived equivalence"). That is a property of a specific
   rewrite system, not of the derived category. The live work (T10, T5) is about
   whether such walks stay in a class; until that is understood, a closed form
   for the pairing centre of the walk is a theorem about a tool, and its payoff
   (halving censuses, E-051) is already collected by the committed reduction.
2. `k(c)` is the one open number, and the data already available (139 cores x
   n = 10, 12..16, E-064) has been fitted three times (E-060, E-061, E-068) without a
   rule beyond the drift families. A fourth fit is the pattern E-067 warned about:
   most informative fits are chance-level at |R| <= 4, and 39 of 109 are vacuous.
3. The only derivable piece, the drift, is exhausted by E-071: it closes exactly
   when the seed `aaa` does not collapse to `34`, i.e. it explains which families
   it explains.
4. The odd-n and n = 16 fit losses are one phenomenon (two singleton orbits of
   equal size, one orbit and its mirror, E-070, E-064). It belongs to T3/T8
   (orbit-plus-mirror vs key), which is dormant for the same reason and could be
   closed together with this.

## Reproduction

Nothing run. To check the quoted rows:
`grep -n "^## E-0\(59\|60\|61\|62\|65\|67\|71\)" research/EXPERIMENTS.md`.

## Prior record

Everything above is in `research/`. New: only the proposal to close, the grouping
of T1/T2 with T3/T8 for the ledger, and the reason 1 (the pairing is a statement
about the walk, so its theory waits on T10/T5). I did not check whether H-021's
header should be rewritten from OPEN to a status that records the iff as REFUTED
and H-021' as SUPPORTED-AS-DESCRIPTION; that is the chair's call at 048.

## Code changed

None.

## Next

- Chair, round 048: close T1/T2 as above, or name the held-out prediction that
  would keep them open (a rule for `k(c)` for cores outside the drift families,
  tested on cores not used to fit it).
- One cheap thing worth doing before closing, if the conference wants it:
  theorist, state E-065's upper bound as a stated gap in H-021's text rather
  than in an experiment's Limits line, so the closure ledger lists it.
- Maverick: back to S-1 n = 15 K0 = 5 sizing.
