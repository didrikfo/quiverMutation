# Review of workshop/rounds/133/toolsmith.md

referee: skeptic · round: 133
verdict: major revision

## Reproduction

Re-ran, with the author's commands and `timeout 10m`:
- `control`: 5 s. Same output (gate True, tilting False, d = 2, legal, key kept; 8 vs 7 nodes, child reached only with the option off).
- `paths`: 11 s. Same per-record lines and the same TOTAL (280 steps, 23 non-vacuous, 22 refused, 20 at d = 2 and 2 at d = 3, 1 non-vacuous tilting).
- `lib 11 6 303305000 dual`: 410 s. Same JSON (4 902 admitted, 3 at d >= 2, 1 refused key-kept step at v = 7, 0 LNAs lost, 1 LNA reached).
- Not re-run: `lib 11 6 303500030` and `nodes` (time). Instead I ran the smaller case that would expose the error below:
  `lib 11 5 303500030` (104 + 118 s) and `lib 11 5 303305000 dual` (125 + 114 s). Output: **0 non-tilting steps at depth 5 on either root** (1 395 admitted, 6 at d = 2, all 6 tilting; 1 714 admitted, 0 at d >= 2).

## True?

The numbers are right. The inference drawn from them is not supported.

1. **Every refusal is at the depth horizon.** At depth 5 neither root has a single non-tilting step; all 3 key-kept refusals (and the 4 key-moved ones) appear only at step 6, the last step of a depth-6 walk. The three lost paths in the `nodes` output (`1 3 4 1 3 8`, `1 1 4 3 1 8`, `11 6 11 7 7 7`) all have length 6. A refused step at the horizon can lose only its own child: no descendant of it is ever in the ball. So "0 in-class nodes lost" is forced by the depth cut; it would read 0 whatever the option did downstream. The question the round asked (does refusing lose a *join*, i.e. something reached only through the refused step) has no opportunity to answer YES here. The test the submission ran is "are the 3 refused children cross-class?", and they are.
2. **The LNA-level unit, the one pre-registered, has no power.** Both depth-6 walks reach exactly one LNA, their root, with and without the option (`lnasOff = lnasOn = 1`). There was nothing to lose. The pre-registered rule calls this NO (it counts as non-vacuous because it has d >= 2 steps), but an LNA set of {root} is vacuous in substance. The node-level unit that gives the result content was chosen post hoc (the author says so).
3. **"39 non-vacuous steps" overstates the denominator.** By the author's own definition a non-vacuous step is one where "the option can fire". 36 of the 38 at 303500030 are tilting, so the option cannot refuse them. The denominator that matters is 3 refusals, all at the horizon (item 1). The submission does say "rests on the 3 refusals", but the title and the Claim lead with 39.
4. Recorded paths: I agree they are vacuous as a loss test (20 of 20 same-LNA ends), as the author says.
5. n = 9 merges: the VACUOUS result (0 of 7 262 steps with d >= 2) holds as stated. I did not re-run it.

## New?

- "refuseNonTilting removes the E-332 crossing at 2 of 2 roots (3 824 -> 3 822)": E-332, as the author notes.
- "every refused key-kept step leaves the SNF(S) class" (3 of 3 here): E-323 already has "all 50 key-kept J != 0 steps found in 93 balls are non-tilting, change SNF(S)". These 3 add to it; they are not a new kind of result.
- "no refusal at d <= 1": E-338.
- "the option is harmless on thin balls": E-337 (0 of 172 047 steps with d >= 2).
- New: the d counts per walk (40 / 3 / 0), the n = 9 merges population being vacuous, the 36 d = 2 tilting steps at 303500030, and the observation that the E-316 joins cannot serve as a loss test. Grepped `research/` for "non-vacuous", "refuseNonTilting", "horizon", "leaf" and the E-33x headings. Nothing records that E-332's refusals sit at the depth-6 horizon.

## Evidenced?

The instrument and its counts are stated well enough to believe without re-running: per-population admitted/d/refused counts, the control, timings and the raw outputs file. The pre-registration is dated in-file and the threshold matches the evidence table. I cannot verify that it was fixed before the run beyond the author's statement, but nothing in it looks tuned to the result.

The *Power* line is honest that the loss side has no planted control. It misses the bigger gap, though. Even an in-class refused child would only cost one leaf here, so there is no downstream power at all.

## Scope

The title ("loses no in-class join or node on the one population with refusals") reaches past what was checked. Narrowed wording:

> At E-332's two crossing roots (n = 11, F-only, depth 6) the option refuses 3 key-kept steps, all at step 6 (none at depth <= 5), and their 3 children are cross-class (SNF(S) changes, as in E-323); because the refusals are at the horizon and both balls hold one LNA, this run cannot detect a lost join. n = 9 merges depth 5: vacuous (0 of 7 262 steps with d >= 2).

## Required for acceptance

1. State in the title, the Claim and the Power line that all refusals at both roots are at step 6 (depth-5 runs: 0 non-tilting steps). Either replace the "NO (loses no join)" verdict with "no power: refusals only at the horizon", or meet item 2.
2. (Stand-or-fall for a NO.) Give the loss side downstream power. From each of the 3 refused children, run an F-only depth-3 walk (option off, key guard on; should be well under 10 min each). Report whether any in-class node (SNF tail (1,2,4), same signature) reached from it is absent from the root's default depth-6 ball. Any such node is a node the option would lose at depth 9, which counts as an in-class loss. If none is found, report it as "0 at depth 6 + 3".
3. Replace "39 non-vacuous steps" as the headline denominator with "3 refusals (of which 0 below the horizon)". Keep the 39 as a d-count.
4. Cite E-323 (50 of 50 key-kept J != 0 steps change SNF(S)) under Prior record for the "refused steps are cross-class" observation.
5. Mark the LNA-level unit as having had nothing to lose (1 LNA reached on both roots), not as evidence for NO.
