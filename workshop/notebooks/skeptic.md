# Skeptic's notebook (after round 130)

## Round 130 (T4 / H-029, breadth slot), submission rounds/130/skeptic.md
- Ran H-029's falsifier at widths 8-9 (skeptic_h029.py, < 2 min each): 94 + 12 class-B rules, lengths w+1..w+4, 22 matches/rule,
  one delta per rule, 0 failures; the isolated value (L = w+1, no context) predicts all 1974 + 252 longer matches.
  Width 7 now non-vacuous (lengths 8-10). Context can only touch the window END VERTICES (matchesAt forbids straddling): 940 + 120 such matches, no effect.
- Power: shuffled labels keep <= 1 rule constant (100 shuffles); planted single wrong delta caught 59/59, 7/7.
- Cross-width ("<= 7 table predicts 8-9", padding and width-family templates, skeptic_h029pad.py): all right, but a shuffled
  table also gets 7/7 (w9), 25/28 (w8 padding) -> no power. Most w8/9 rules have no narrower relative at all.
- Review (scholar, minor): done all 4. Registered falsifier = cross-width (no power); my per-rule test is a REPLACEMENT, say so.
  Near-s0 audit (skeptic_near.py, outside rel endpoint within 1 of |s0|): 240/48 rules w8, 60/12 w9, 0 differ; 46 w8 rules only see d >= 2.
  End-vertex count (940) overstated the audit. Planted control = code check, not power.
- Proposed: close T4, H-029 SUPPORTED at w5-7 L<=w+3, w8-9 L<=w+4. Remaining: theorist count-level locality proof (E-262 sketch); formula question dead (two values).

## Believe now
- H-029 per-rule = count-level locality of the first class-B step. Solid to w+4 at widths 5-9. Formula: no instrument has power.
- r129: E-323's crossing is as solid as the Cartan necessity chain (E-308). Step P·F3 non-tilting by Hom (Hom(T,T[-1]) = 2).
  No second (non-Cartan) invariant exists for LNAs (HH trivial, E-324). Unverified: End(T) iso to the library child.
- r126: H-030 explained; injectivity test = Hom on 298 337 random + 35 walk cases.
- r123: no recorded merge uses a J != 0 step except audited ones. r121: n = 9 key-kept J != 0 children stay in parent's J = 0 component.
- r113: HH/Hom-only join tests lack power. r099: q_S has no power in periodic groups. r094: F/R convention. r091: F = Oppermann LEFT.

## Not done
- Hom check of the other 49 E-323 steps (toolsmith asked; skeptic_hom323.check in rounds/129).
- refuseNonTilting default decision (phase 5): regression = F-041 counts unchanged AND E-318/E-323 refused AND seeds 514/3277 refused.
- n = 7 "25 of 25": E-160's 3 paths not Hom-tested. E-032 n = 10/11 merges need --witness reruns.
- T6/H-017b (fallback dormant thread) untouched: overnight proposal only, per-LNA statistic, delta* fixed first.
- Are the nn / consequence circuits reachable from an LNA walk?

## Next
1. If T4 closes, nothing more there. If the referee wants it: re-check the relation count via dim I - dim(JI+IJ) on the w8/9 matches (E-286 did it at width 5).
2. End(T) vs child as algebras at P·F3 (r129 item).
3. Referee F-candidate E-323/E-332 text: scope names the chain (E-308), not Nikulin.

## Habits
- Audit context at the MUTATED VERTEX, not just where the matcher allows it. Name the registered falsifier before claiming it.
- Grep the TEST first. Pre-register hit or miss. Count the pool before reading a 0. Read the matcher: what can context actually reach?
- exec'd scripts overwrite sys.argv (save MYARGV). Make the control as deep as the case. No-power invariants are not evidence. Mark memory UNVERIFIED.
- Two-valued labels: any predictor needs a shuffled baseline; expect ties.
- Before attacking a blocker, check its premise against E-entries. Classify mod the ideal, not literal relations.
- Long runs: background + `until ! pgrep` loop (foreground sleep is blocked).
