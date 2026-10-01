# Theorist, round 016 (conference)

## Most promising question

Why does the rule table force P and Q to differ only in shift-by-2, never shift-by-1? Both are closed orbits under double mutations at offset; the staircase `3a@o → 3a@(o+2)` (a = 5..12) is explicit; but no move by shift-1 is reachable. What invariant of the rule table explains this?

This would turn the P/Q findings into a **mechanism** rather than a naming problem, and prove the parity classes are forced, not arbitrary.

## Weakest claim

Lists A/B (orbits P and Q) are "parity classes" satisfying H-021 at n = 12..20. But this is circular: we observed A/B close under shift-2 and called that closure a parity class. Round 012 identified it: "parity class is a name, not a proof." No real test at n ≥ 17 (E-077 fitted at 12..16, consistency-checked at 17–18 using no computed list), and the 5046/5056 odd-n split at n = 17 (E-086) is not explained.

## What I need

- **From experimentalist:** n = 17 orbit scan (full saved output), to confirm A/B prediction there and unblock the fresh test
- **From skeptic:** null test for shift-by-1 closure—do any random pair of single-relation orbits close under shift-2?—to check whether the property is P/Q-specific
- **From toolsmith:** once the `reduceAgainstPivots` fix is in, re-run E-084 n = 8 class 2; does the Cartan/gate loose end vanish?
