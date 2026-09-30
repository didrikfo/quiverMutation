# Round 009 -- proceedings

Ordinary. Worked: experimentalist (T1/T3), theorist (T2/T4), maverick (T6). Referees: skeptic, experimentalist, toolsmith. All three: minor revision; I applied the wording fixes myself (chair notes at the end of each submission) and accepted all three. Step 0.5: the round-008 conference left no questions; the agenda proposed there stands.

## experimentalist -- accepted (with rescoping)
Claim: equal-size singleton pairs of `344 348 349` (n = 15..17) and `4046` (n = 14..16) are one orbit plus its mirror; orbit+mirror = key in all 12 cells; key-coarser cores at n = 12, 13 are disjoint from the 7 of E-059. Referee (skeptic) reproduced 348@16, 4046@15, 349@15 and n = 14; found it mostly in E-064 (only n = 17 new), claim 2 shown for n = 12, 13 only, and the "pair at even n, mirror-join at odd n" summary contradicted by its own table. Applied: rescoped, summary withdrawn. Promoted: **E-070**; H-021 status line.

## theorist -- accepted (computed criterion, not a theorem)
Claim: `aax` drift families close into pairs with `k = 2x + 3 - a` unless the seed `aaa` collapses to `34`; `33x`, `55x`, `66x` close, `44x` does not. Referee (experimentalist) reproduced the cheap scripts and ran `55x` n = 12, 13, `66x` n = 13, `44x` n = 12, `77x` n = 15: all agree; novelty confirmed. Applied: "exactly" softened, referee runs added. Promoted: **E-071**; H-021 status line. The reason `34` alone reaches the slider is unproved.

## maverick -- accepted (sizing corrected)
Claim: L = 5 control passes (42/42 at n = 6, 8/8 at n = 7, 0/50 at L-1); four n = 9 candidates reach nothing at depth 5. Referee (toolsmith) reproduced control and depth 4; found the depth 6 sizing wrong (about 30 min for 4 cells, 2 h for 16, one candidate per shard), "recorded paths shortest" tautological, n = 7 sample near-trivial. Applied in full. Promoted: **E-072**; H-017 status line.

## Questions for the steering committee
1. **H-017 depth 6 for all 16 n = 9 candidates** (about 2 h, 16 shards of one candidate; `maverick_verify.py` needs a candidate-index argument first, toolsmith). Recommend: toolsmith adds the argument and a `--budget-hours` wrapper in round 010; no overnight run yet -- depth 6 fits a chair's-slot run, and only depth 7 is overnight (already approved in Menu 4).
2. **Promote a criterion to H-021's text?** Recommend not yet: a skeptic's neighbour-aware null for "a = 4 is special" first.

## Decisions taken for the steering committee
None (no open questions from round 008).
