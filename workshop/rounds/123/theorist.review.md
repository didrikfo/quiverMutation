# Review of workshop/rounds/123/theorist.md

referee: scholar · round: 123
verdict: minor revision

## Reproduction

`timeout 600 .venv/bin/python -W ignore workshop/rounds/123/theorist_lnaparent.py`: 7.6 s wall (1.7 s CPU). The output matches the submission. I counted it independently from the printed lines:
- 82 instance lines. (LF) True on 82 (cycles False on all 82), H1w True on 82, ENG (0, 0, 0, ·, ·) on 82, and [MM] H = A YES on all 20 parents.
- Strict (H1) False on 26.
- Admitted is the F column on A lines and the R column on A^op lines, gate 1 and J0 1. That gives 48 admitted steps, and 0 of them fail strict (H1).
- Theta' is non-empty on 10 instances, 8 of them admitted.
- Under the swapped reading, 28 steps are admitted and exactly 2 fail strict (H1): A6 3300 at i = 4 (A) and at i = 2 (A^op).
- Controls as stated: [LF] -3 → NO; toy weight-0 cycle → NO; [MM] 1 failure; SHctrl NO on 3300 F1 and on the tree; ENGctrl bad cells 2 / 2 / 6.

Extra run, my own toy through the author's `lf.mk` / `lf_line`: kQ/(a1 c d - a2 d, a1 e) with a1: 0->1, c: 1->2, a2: 0->2, d: 2->3, e: 1->4. Homogeneity forces w(a2) = w(a1) + w(c), so A-weights cannot be equal.
- With w(e) = 1, (C) fails (w(q) = w(a2) = 2, letter q a2^-1 of weight 0), and (LF) still holds.
- With w(e) = 2, (C) holds.
- So (C) is sufficient but not necessary. On this parent it can be reached by reweighting.

## True?

The proof of the Proposition is correct. I checked each step against the letter list the engine produces:
- A-arrows are not kM letters. Only their stars are.
- Only stars end at i.
- Suffixed letters start at t(b), which is not i.
- The [LF] toy is a real parent, kQ/(a1 c). Its letter f a2^-1 is a loop at 2 of weight -3, so it is a negative cycle by itself.

The hypothesis "d p != 0" is load-bearing. A degree-1 letter g from i with d g = 0 gives g a2^-1 of weight 0 and kills (LF). I ran this. Minimal models exclude it, as the author says.

There is a shorter proof, and it gives a stronger result. Shift the weight of every letter by phi(t) - phi(s), with phi(i) = c, phi = 0 elsewhere, and max w(A) < c < min{w(p) : p a non-A letter from i}. (C) says exactly that such a c exists. After the shift:
- stars have weight c - w(b) > 0;
- unsuffixed letters from i have weight w(p) (+ w(a)) - c > 0;
- no other letter touches i, and loop letters at i are unchanged.

So under (C), every kM letter has positive weight after a vertex shift. That is a positive (Adams-connected) grading of kM, which is the standard hypothesis for Anick-type arguments, and it implies every cycle is positive. Steps (i)-(iii) and the pairing argument are a longer route to the same place. I recommend stating it this way, with (C) as the explicit condition for the shift.

Smaller points:
- The tree's docstring says "H = A not machine-checked here (d^2 = 0 is)". The submission does not say this. The model is the Anick resolution of a monomial algebra (2-chains y al, al z; 3-chain y al z → c), so it is right, but say so.
- (H1w) on LNAs is a sketch. "At most one Anick chain per cell and degree" holds for monomial LNAs, and the machine check agrees on all 82.

## New?

The record has no proof of (LF) for Construction-2.2 parents:
- E-319 says "(LF) checked per model, not proved for Construction-2.2 parents (T5 open)".
- H-031's status line says "(LF) unchecked for Construction-2.2 parents".
- So the Proposition is new to the record, and so is the (H1) → (H1w) correction: no hit for "H1w", and "pivot" hits only E-319 / H-031 / E-309-era entries.
- H(I) = 0 on no-loop steps is E-235, which the author concedes.
- No entry uses "Adams grading", "vertex potential" or "reweight" (grepped EXPERIMENTS, FINDINGS, HYPOTHESES, literature/).

Citation gap: rows 30, 3300 and 2400 have relations of length 3-4, which E-252 does not cover (its scope is length-2 zero relations). The F/R reading for them is E-258 (length 3-4 monomial, n = 5-7), which also covers the J/gate matching. Cite E-258 alongside E-252.

## Evidenced?

The machine evidence is stated specifically, per instance, and it reproduced exactly.

But "(LF) YES on 82 instances" has no power for the Proposition: every letter quiver on a line is acyclic, so the check is vacuous there. The Proposition rests on the proof plus one 2-arrow toy, and the toy tests necessity of the hypothesis, not the argument. That is acceptable for a proof, but the Power line should say it.

The [SH], [ENG] and [MM] controls bear on the per-instance checks, not on (LF). The author is right that strict (H1)'s 27 NOs are not power for cancellation.

## Scope

The title overreaches. "(LF) holds for every graded Construction-2.2 parent" is false on the author's own toy, which is graded, has w(a2) = 5, and has a weight -3 cycle. Narrowed wording:

"(LF) holds for every Construction-2.2 parent with a positive weight making d homogeneous and satisfying (C) (e.g. all A-arrows of equal weight); under (C) a vertex shift at i makes every kM letter positive."

"Loop parents covered" is by proof only. No finite-dimensional loop parent was run, and the submission says so.

## Required for acceptance

1. Narrow the title and the Answer line to the wording in Scope ("given (C)", not "every graded").
2. State that (C) is sufficient, not necessary (my toy: kQ/(a1 c d - a2 d, a1 e), w(e) = 1). Either leave open whether (C) is achievable by reweighting for every gradable parent, or prove it.
3. Add (or replace step (iii) with) the vertex-shift argument: under (C), all kM letters have positive weight after phi(i) = c.
4. Say the tree's H = A is not machine-checked, and that it is the monomial Anick resolution.
5. Cite E-258 for the F/R and gate reading on rows 30, 3300 and 2400.
6. Power line: say the 82-instance (LF) check is vacuous (acyclic letter quivers), so it is not evidence for the Proposition.
