# At n = 10, 35 of 35 edges of the E-169 group-A paths pass the independent Hom(T,T[±1]) test, which also rejects 312 of 312 wrong-vertex/wrong-direction steps at n = 10 and 9 of 9 J != 0 steps at n = 7 (n = 10 pass and n = 7 control kept separate)

author: toolsmith · round: 061 · kind: result
thread: T10 (agenda item 1) · bears on: E-169, E-161, E-166, E-174, E-151
scope: n = 10, one start (LNA [0,5,0,4,0,3,3,0]), the 5 paths of E-169 (35 edges, no parallel-arrow quivers by construction); power control: (i) n = 10, 700 steps = every vertex of both a and dual(a) at the parent of each of the 35 edges; (ii) n = 7 class 2 only (9 failing steps of the rebuilt E-151 walk, 20 000 expansions). Cartan level + Hom dimensions; generation of K^b(proj) by T assumed (E-168); End(T) as an algebra is not retested here (E-174 did that). Class-1 failing steps (16) not rerun.

## Response to referee

1. Retitle: done; the n = 10 pass (35/35), the n = 10 power check (item 3) and the n = 7 control (9/9, n = 7 only) are separate in the title and table.
2. Cartan convention, once: "cart" means the labelled matrix H[i][j] = dim Hom(T_i,T_j) of End(T) equals the transpose of the child's Cartan matrix, H = Cartan(child)^T (the convention of E-161 / `skeptic_replay`); the claim's "equal" is read in that sense, nowhere else.
3. Power at n = 10: done. `toolsmith_homn10_power.py` runs `stepTest` at every vertex of a and of dual(a) (wrong vertex, reversed direction) at the parent of each of the 35 edges: 700 tests = 35 used steps (all pass) + 665 other. 312 of the 665 are rejected (Hom(T,T[-1]) = 1..6, Hom(T,T[1]) = 0 always); the other 353 pass and are themselves tilting by `tiltingPlus`. Agreement with `tiltingPlus`: 700 of 700 (every reject is a tiltingPlus fail, every pass a pass). So the test is not vacuous at n = 10 size. Caveat: the rejected steps are not E-169 edges and are not J != 0 steps of any walk; they show only that stepTest rejects non-tilting T at n = 10, with the same J-defect source (Hom(T,T[-1]) = sum J_i, E-128). Output `toolsmith_homn10_power_out.txt`; ~30 s.
4. Gap (2), "rejects 9 of 9 rests on the [-1] column": yes, Hom(T,T[1]) = 0 at all 9 (and at all 312 n = 10 rejections); stated in the table.

## Claim

On every one of the 35 edges of the five explicit 7-step tilting-only paths 05040330 -> 33460000 (E-169), the complex T replacing P_v has dim Hom(T,T[-1]) = dim Hom(T,T[1]) = 0, and the Cartan matrix of End(T) equals that of the child (labelled), as tested by the round-050 independent test (`skeptic_tilt.stepTest`, which does not use `tiltingPlus`, `perI` or the gate). So the J = 0 premise for the group-A merge is now checked at the Hom(T,T[m]) level, not only by the gate/`tiltingPlus` as in E-169 and not only by End(T) as in E-174 (which cannot discriminate). The same test returns Hom(T,T[-1]) = 1 and rejects 9 of 9 J != 0 failing steps of n = 7 class 2, with Cartan(End T) still equal to the child's (so Cartan alone would accept them). Power: at n = 10 the test rejects 312 of 312 deliberately wrong steps (not edges of the paths) and agrees with `tiltingPlus` on all 700 tested steps. It does not claim: derived equivalence of the end points (generation assumed), anything beyond these 5 paths and this start, or that a genuine J != 0 step at n = 10 would be caught by anything but the same mechanism (none is known at n = 10, E-169, E-175).

Would be refuted by: any edge of the 35 with a nonzero Hom(T,T[±1]) (none).

## Evidence

| set | edges | Hom(T,T[-1]) = 0 | Hom(T,T[1]) = 0 | Cartan(End T) = Cartan(child)^T | `tiltingPlus` passes |
|---|---|---|---|---|---|
| n = 10 path 0 (2 1 6 5 -4 -10 -8) | 7 | 7 | 7 | 7 | 7 |
| n = 10 path 1 (-10 -8 -4 -5 2 1 6) | 7 | 7 | 7 | 7 | 7 |
| n = 10 paths 2-4 (4 3 2 1 7 6 9; 4 3 7 2 1 6 9; 4 3 7 9 2 1 6) | 21 | 21 | 21 | 21 | 21 |
| total n = 10 | 35 | 35 | 35 | 35 | 35 |
| n = 10 wrong vertex/direction, 665 steps (power) | 665 | 353 pass, 312 rejected | 665 | - | 353 pass, 312 fail |
| n = 7 c2 failing steps (control, n = 7 only) | 9 | 0 (all Hom = 1) | 9 | 9 | 0 |

- Negative steps (s < 0) are taken on the dual algebra and carried back, as in E-169 `replay` and E-161.
- Control: the n = 7 class-2 failing steps have Hom(T,T[-1]) of dimension 1 at 9 of 9 (E-174 gave the same 9 and the 16 of class 1). The user's "one n = 7 J != 0 step" is thus 9 steps of c2; step 0 (v = 5) is the smallest to cite. Hom dimension 1 is the same value as the maximum seen on the n = 10 edges of E-174's End(T) check, but here all n = 10 Hom(T,T[±1]) are 0.
- Agreement with the gate/`tiltingPlus` on all 44 edges (35 pass, 9 fail): no disagreement. This is Ladkani's criterion restated (as in E-161), so it adds code independence at n = 10, not new information about J.
- Per-edge output: `workshop/rounds/061/toolsmith_homn10_out.txt` (35 lines + TOTAL), `toolsmith_homn10_ctrl_c2_out.txt`. Each edge takes < 0.2 s; the pickle rebuild for the control takes ~6.5 min (377 s, one run, DEADLINE 520 not hit).

## Reproduction

```
.venv/bin/python workshop/rounds/061/toolsmith_homn10.py plan          # 5 paths, 35 edges
timeout 10m .venv/bin/python workshop/rounds/061/toolsmith_homn10.py n10        # ~3 s
timeout 10m .venv/bin/python workshop/rounds/061/toolsmith_homn10_power.py      # ~30 s, wrong-step power at n = 10
timeout 10m .venv/bin/python workshop/rounds/050/toolsmith_collect.py 7 2 20000 /tmp/c2.pkl 20   # 377 s
timeout 10m .venv/bin/python workshop/rounds/061/toolsmith_homn10.py ctrl /tmp/c2.pkl            # < 1 s
```

## Prior record

E-161 / E-166 did this test for the n = 7 paths (40 and 33 edges) and E-166's scope says nothing at n = 10. E-169 tested the n = 10 edges by the gate and `tiltingPlus` only; E-174 did End(T) = child at all 35 (max Hom dimension 1) and says "only Hom(T,T[-1]) separates" but ran Hom(T,T[-1]) only for the 25 failing n = 7 steps (`toolsmith_endt2_hom.py`). So the n = 10 Hom(T,T[±1]) check on the 35 edges was not recorded (grep E-169, E-174, E-175: E-175 did Hom(T,T[±1]) on 1551 edges of the n = 9 F-010 balls, not on these paths). New: that check. It is expected given the gate (J = 0 at all 35), so it is an implementation cross-check, not a surprise.

## Code changed

New script `workshop/rounds/061/toolsmith_homn10.py` (reuses `skeptic_tilt.py`, `scholar_h015.py`). No library changes; no tests needed.

## Next

- skeptic: a genuine J != 0 step at n = 10 (none known) would be the only fully faithful control; the wrong-vertex power check is done.
- experimentalist: same Hom(T,T[±1]) replay for the 7 E-163 / the E-160 class-1 paths with the End(T)-checked dim (done n = 7 in E-166), and for any future n = 10 path with a J != 0 step.
- Remaining T10 gap is the premise itself (G2, E-181), not these edges.
