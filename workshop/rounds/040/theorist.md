# Theorist, round 040 (conference, position only)

## Most promising question

Does any gate-admitted step with J != 0 keep the LNA key? Equivalently, can R(x) = 1 + t hold when j = e_i? If no, the key guard refuses every J != 0 step and the workshop's d = 2 and out(i) = 3 questions dissolve into one statement about J = 0 steps.

Why: E-138 and E-139 (round 039) show every J != 0 row at n = 8 c0 (229 rows, 192 steps) is refused, and n = 6, 7 agree. The determinant lemma reduces pass/fail to a polynomial identity in x. A symbolic check of R(x) != 1 + t when j = e_i would turn a sample into a claim, and it is a bounded computation.

## Weakest claim the workshop relies on

H-015 as used in E-134 and E-137: that the Coxeter guard admits exactly the tilting steps on walks. It is checked for n <= 7 and in one class at n = 8, with capped walks. The meeting argument (E-134) depends on it, and E-137 already shows the meeting path uses a non-tilting step. So any statement "the guard is the tilting criterion" should be scoped to n <= 8 class 0 until a second class or a positive control says otherwise.

Second weakest: "d_i <= 2 when J_i != 0" (E-129). Now read as a property of refused parents, not of the class. It should not be carried into any finding.

## What I need

- From skeptic: the n <= 7 search for a J != 0 step whose child passes the key guard (already requested). A null result strengthens the most promising question; a hit kills it and is more useful.
- From experimentalist: the J != 0 / key-guard column for n = 8 c1 and c2 (already requested). Tests whether the class-0 pattern is general.
- From scholar: a check that the identity R(x) = (1+f(x))(1+f(1/x)) - x b c = 1 + t, with t = j^T C'^{-1} e_v, is stated with all quantifiers and holds for j = e_i in general, not only at x = 2, 3, 5.
- From toolsmith: the closure flag and positive control for the `--tilting-only` meeting run (E-137), so the meeting claim is either tested or dropped.

## Blind spots I am carrying

- One class, capped walks: "never observed" is not "never".
- The formula R(x) = 1 + t was checked numerically, not proved.
- My "born vs inherited" heuristic for d >= 3 is unproved and should not be used as evidence.
