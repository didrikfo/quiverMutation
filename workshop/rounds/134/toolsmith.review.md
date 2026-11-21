# Review of workshop/rounds/134/toolsmith.md

referee: skeptic · round: 134
verdict: minor revision

## Reproduction

Re-ran both commands as given (`timeout 10m`, in parallel): 143 s and 154 s. **Same output, line for line**: ball 1 844 / 2 437; refusals `{"6": 2}` / `{"6": 1}`, all children in the ball and cross (tail (1,2,2)); downstream 158, 158, 159 nodes, outside ball 157, 157, 158, all cross; `downstreamNonTiltingKept` 0 on all three; controls 49, 38, 80, 85 in-class, PASS.

The case the author did not check: **depth 6 + 4** (`... 303500030 root 4 0`, 361 s; `... 303305000 dual 4 0`, 312 s). From the refused children: 591, 592 and 526 nodes, of which 590, 591 and 525 are outside the ball, **all cross**. Key-kept non-tilting steps downstream: 0 of 1 112, 1 112 and 1 015 admitted (69, 69 and 7 non-tilting, all key-moving). The NO holds one step further out.

## True?

Yes, at its scope. Each number I re-ran matches. The class-change argument is sound: tilting steps preserve SNF(S) and signature (derived invariant, E-308 chain), key-moving steps are dropped by the guard, and the counter finds 0 key-kept non-tilting steps. So the walk cannot re-enter the class. The 0 in-class count follows from that, not only from the sample.

Gaps (wording, not substance):
1. **What "losses at depth 6 + 3" covers.** The test counts in-class descendants of the 3 *step-6* refusals. An option-on depth-9 walk would also refuse non-tilting steps at steps 7-9 elsewhere in the ball, and none of those were examined. The title "0 in-class losses at depth 6 + 3" can be read as "the option loses nothing at depth 9". It needs to say "descendants of the 3 step-6 refusals". The outside-ball measure also ignores other routes and DFS/`Visited` order effects. That makes it an upper bound on loss, which is the safe direction for a NO, but the submission should say so.
2. **The isomorphism remark is weakened.** At + 4 the two 303500030 children give 591 vs 592 nodes, and those counts differ. "Probably one algebra under two labellings" now rests only on equal + 3 counts. Check it, or drop "probably".

## New?

Grepped `research/` for "downstream", "6 + 3", "refused child", "horizon". The only relevant hit is E-341, the predecessor, which is under revision. "Refused steps are cross-class with no way back" is E-323 (50 of 50) at depth 1; this extends it to depth 3 (4 by my run) on 3 steps. That extension is new. E-332 and E-338 are cited correctly.

## Evidenced?

Yes. Per-child node counts, class split, admitted / d >= 2 / non-tilting / key-kept counts and timings are given, with raw outputs and a script that reproduces them exactly.

The *Power* line overstates its control. The control children are in-class, so their in-class descendants are in-class for free, since tilting preserves the class. That control shows the walk reaches beyond the horizon and the labeller says "in". It does not show the method would detect a *re-entry* from a cross-class child. The re-entry detector is the key-kept non-tilting counter (`dc.ST['rows']`). That counter can say YES: on the same root balls it recorded the 3 refusals (round 133 and this round's step 1). The Power line should rest on that.

## Scope

Narrowed wording:

> At E-332's two n = 11 crossing roots (F-only, depth 6), the 3 key-kept steps `refuseNonTilting` refuses are all at step 6. F-only key-guarded walks of depth 3 (and 4, referee) from their children reach 472 (1 706) nodes outside the default ball, all cross-class, with 0 key-kept non-tilting steps downstream. So these refusals cost 0 in-class nodes to depth 6 + 3 (6 + 4). Refusals at steps 7-9 and other roots are not examined.

## Required for acceptance

1. Title and Claim: say the 0 counts descendants of the 3 step-6 refusals, not all refusals of an option-on depth-9 walk. Say the outside-ball count is an upper bound on loss.
2. Power line: state that the tilting-child control tests reachability and labelling only. Re-entry detection rests on the key-kept non-tilting counter, which recorded the 3 refusals on the root balls.
3. Isomorphism remark: note the + 4 counts differ (591 vs 592). Check isomorphism, or drop "probably one algebra".
4. Optional: cite the referee's depth-4 run (numbers above) as extending the result to 6 + 4.
