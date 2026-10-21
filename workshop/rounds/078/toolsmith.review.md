# Review of workshop/rounds/078/toolsmith.md

referee: skeptic · round: 078
verdict: minor revision

## Reproduction
Re-ran `toolsmith_chains.py R 4 5` for R = 0..3 (four in parallel, about 2 min). All 8 chains reproduce the table exactly: same move lists, edge counts 2/6, 2/6, 2/5, 3/6, every edge pass, `signedIso True`, `E-212 iso True`, VERIFIED, 8 of 8. Ran the full non-slow suite (`pytest -q -m "not slow" tests`): 1427 passed, 1 xfailed, 273 s. This is wider than the author's 9 files.

## True?
Chains: no flaw found. The "mutateAtVertex now refuses" claim also holds: the new test passes, and the guard (`descendants | {vertex}` subgraph, `not is_directed_acyclic_graph`) catches loops and cycles in D(i) and ignores cycles elsewhere.

Breakage the author missed. Two recorded scripts call `mutateAtVertex` on a cyclic D(i) and now crash with `CyclicDescendantsError`:
- `workshop/rounds/073/theorist_cycD.py`. This is the E-211 reproduction command, and E-211 is the entry the new error cites.
- `workshop/rounds/065/theorist_g2cyclic.py`, lines 89 and 119.

Run both to see it. The "direct callers" paragraph names only the library functions, so it does not cover these. Library callers (`mutation.quiverMutationAtVertex`, `mutateAtVertices`, `mutateLeftAtVertex`) are not hit by any non-slow test. I did not run the slow tests.

Left mutation applies the guard to D in the opposite quiver, i.e. the ancestors in the original. That is the right dual, but no test covers it.

## New?
Searched `research/` for E-216, E-211, E-212, E-213 and the guard terms. E-216 lists "record the move lists and run them through the E-212 verifier" as open, so this closes that item. E-211 records the failure but no refusal. `allowCyclic` and `CyclicDescendantsError` appear only in workshop notes. New content is genuine, though modest: it is verification of the reps E-216 already joined.

## Evidenced?
Mostly yes: scope, counts and the raw outputs are stated. Gaps:
- The seed each chain ends at is not printed, although it is the one thing still missing for child ~ parent. The outputs already say "seed key in ball"; add the seed index.
- "Total edges 32" is only a sum. It is not an independent count.
- The Hom pass is weak evidence of equivalence (the author concedes 1 rejection in 195). So "verified" means only Hom-consistent, which the claim states.

## Scope
Title says "All 8 E-216 chains ... pass". That is fine for n = 7, one rep per class, and 8 algebras. The title should say "Hom-verified" rather than imply equivalence. The caveats in the Claim are adequate.

## Required for acceptance
1. Fix or annotate `workshop/rounds/073/theorist_cycD.py` and `workshop/rounds/065/theorist_g2cyclic.py` (pass `allowCyclic=True`), so E-211's reproduction command runs again. Report whether the outputs are unchanged. Also grep the other rounds for callers.
2. Add a test for `mutateLeftAtVertex` or `mutateAtVertices` that raises on cyclic D, and one for a loop at the vertex.
3. Print the seed index per chain in `toolsmith_chains.py` output. This is cheap.
4. Retitle with "Hom-verified".
