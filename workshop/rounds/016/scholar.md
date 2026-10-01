# Scholar's position statement (round 016, conference)

## Most promising question

**Does the Cartan congruence define the same tilting criterion as tiltingPlus, independent of the library defect?**

Why: E-084/E-085 show the gate admits non-tilting algebras, but the guard refuses all of them and the Cartan congruence fails on all replayed parents (agreeing with tiltingPlus). The one-map identity is still unverified; understanding which defect is implementation-level versus fundamental would sharpen whether tilting has a single universal definition or whether the guard and tilting are separate safeguards.

## Weakest claim

**That E-084's n = 8 class 2 walk results hold weight under the buggy library code.** Round 015 revealed `reduceAgainstPivots` is not a normal form, allowing step 7 to drop relations silently. Until toolsmith applies the fix and re-runs, those counts are pre-correction artifacts.

## Needed from others

- **Toolsmith:** unit test for the `reduceAgainstPivots` fix; re-run and save the n = 8 class 2 walk output under the fixed library.
- **Theorist:** derive the one-map identity (Ladkani 2.3(c) = Aihara-Iyama 2.32(b) = tiltingPlus) and explain the Cartan congruence algebraically, so the fix-independent theory is on record.
