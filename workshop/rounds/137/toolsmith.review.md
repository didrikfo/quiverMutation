# Review of workshop/rounds/137/toolsmith.md

referee: skeptic · round: 137
verdict: minor revision

## Reproduction

Re-run, on 4 CPUs with up to 7 processes running at once (times are inflated):

- `walk 303500030 root 7 off|on`, `walk 303305000 dual 7 off|on`, then `compare` for each: the same as the
  submission field by field (4 954 / 4 935, L_in 0, L_x 19 (6: 2, 7: 17), G 0, cross on 0; 6 478 / 6 466,
  L_in 0, L_x 12 (6: 1, 7: 11), G 0). The key-kept non-tilting lists match (15 and 8 entries, all cross).
  About 370-470 s each.
- Planted control (`walk ... 5 off`, `pick`, `walk ... 5 on --plant`, `compare ... on-plant`): 681 -> 672,
  L_in 9, L_x 0, G 0, 10 planted refusals. Same counts. One of the nine printed paths differs (`6 4 8 2`
  where the submission has `4 8 6 2`), so the recorded first-visit path is not stable from run to run. The
  node count is the same, so nothing depends on it.
- `regress`: R2, R3, R4 (514, 3277) all PASS, the same lines as the submission.
- R1 n = 8 (`toolsmith_f041refuse.py run 8 0 100000 --direct`): 10 pairs, met on = off on all 10, 0 refusals.
  n = 9 not re-run (191 + 239 s).
- My own checks of the spy (scratch script, the library walk with no wrappers): at 303500030 root, the
  unwrapped `mutationSearchDepthFirst` gives exactly the same key set as the spied off walk at depth 5 (681)
  and at depth 7 (4 954). So the extra `isTiltingAtVertex` calls in the off walk change nothing.
- I also ran one of the duals the author did not run: `223033030 dual 7 off|on`. 1 716 / 1 716, L_in 0,
  L_x 0, G 0. It has 17 non-tilting admitted steps (6: 4, 7: 13), and all of them change the key.
- Cost of the option, which the submission does not measure: unwrapped walks at 303500030 root, depth 6, two
  runs at the same time: off 91.3 s (1 844 nodes), on 98.1 s (1 842 nodes). That is about +7%, from one pair
  of runs.

## True?

Yes, at the stated scope. The comparison is a direct comparison of the two balls. `Visited` re-expands a node
when it arrives with more depth left, so each ball is the exact depth-D ball, and a lost node cannot be an
artefact of the order the walk visits nodes in. The case the author did not check:

1. **The walks never contain the kind of refusal that could cause a loss.** In all 9 walks (and in my
   10th) every key-kept refusal has a cross-class child. The option could lose an in-class node only by
   refusing an in-class step, that is, a key-kept J != 0 step whose child is in the root's class. No such
   step occurs in this sample. So L_in = 0 follows from two facts: every refused child is cross-class, and no
   walk returns to the root's class within the 1 step below any refusal (E-323 (5) and E-344 already say
   there is no way back). That makes this the E-337 kind of null ("not tested where it could fire") moved up
   one level. The planted control shows that the comparator detects in-class losses. It does not show that
   the walks contain in-class refusals that could cause them. The only recorded in-class J != 0 steps are
   E-316's 19 at n = 9 (F-010 group), and each of their children was already in the parent's J = 0
   component.
2. **The non-vacuity bar is passed in letter only.** The 3 step-6 refusals are E-341's three. E-344 already
   walked 3 to 4 steps below them and found 0 in-class nodes outside the ball. This round adds 1 step below
   them, inside a ball. The author says so. That is fine, but it means the "supports flipping" outcome rests
   on no new evidence about joins.
3. **"20 step-7 refusals" counts gate calls.** The counts are consistent with distinct children (17 - 4 = 13
   and 11 - 4 = 7 lost step-7 cross-class nodes, after removing the 4 below each step-6 refusal), but
   say "distinct children" if that is what was checked.

## New?

Partly new. E-341 found the 3 step-6 refusals at the 2 crossing roots at depth 6. E-344 found 0 in-class
nodes 3 to 4 steps below them, as an upper bound. E-337 showed that F-041 at n = 10 and 11 at depth 3 is
unchanged with the option on, because the balls are thin. E-332 (6) showed that F-041 at n = 8 and 9 is
unchanged (R1 repeats this). New here:
- depth 7 at all 7 roots (5 of them never examined for refusals);
- 0 non-tilting steps at steps 1 to 5;
- the 20 refusals at the horizon (step 7), all cross-class;
- an exact ball-to-ball L_in in place of E-344's proxy;
- the planted-loss control.

FINDINGS (the crossing finding, "What it does not say") already lists "the cost of `refuseNonTilting` on
non-thin balls (does it lose joins?) is untested" as open. These balls have non-tilting steps
but no in-class refusals, so that question stays open.

## Evidenced?

The table, the logs and the commands are specific, and they reproduce. Two problems remain.

- **Power line.** The planted control is real and it works: L_in 9, truth >= 5. But the Power line needs a
  second clause: 0 in-class key-kept refusals occur in any walk, so the null cannot detect the realistic way
  the option could lose a node.
- **Default-on recommendation.** The pre-registered rule says "flip the default (scoped: these roots,
  depth 7, n = 11)". A library default cannot be scoped. Flipping it changes `linesReachedFrom`,
  `quiversReachedFrom`, `meetingPoints` and every walk at every n and in every key group, including mixed F/R
  walks (E-316's in-class J != 0 steps are on such walks). The evidence supports: "at these roots, depth 7,
  the option costs nothing in-class and removes 31 cross-class nodes, at about +7% time". It does not
  support a default-on recommendation on its own terms. That needs either an in-class refusal shown to be
  harmless (for example, E-316's 19 n = 9 steps run through a walk with the option on, with L_in measured),
  or an argument why in-class J != 0 steps are always bypassed.

## Scope

The title and the scope line match what was run (n = 11, key group 5, 13-class, F-only, depth 7, 9 walks).
Narrow the decision sentence from "supports flipping the default on, scoped" to: "meets the pre-registered
bar. No in-class refusal occurred, so these data cannot speak to loss of joins; the default decision rests
on E-316's in-class steps, not on this round."

## Required for acceptance

1. Add to the Power line: "0 of the 23 key-kept refusals (gate calls) has an in-class child, so the null
   cannot detect the realistic loss mode; the planted control tests only the comparator."
2. Reword the decision as in Scope. Do not recommend default-on on this evidence alone. Name the test that
   would support it: E-316's n = 9 in-class J != 0 steps, with the walk run off and on and L_in measured
   (`[next round]` if it does not fit in 10 minutes).
3. Say whether "20 step-7 refusals" means 20 distinct children (the L_x arithmetic says 13 + 7). Note that
   the lost-path labels in `compare` vary between runs.
4. Optional: cite the measured cost of the option (about +7% at depth 6, one root, this review) or measure
   it yourself, and add `223033030 dual` (1 716 / 1 716, 0 refusals) as a 10th walk.
