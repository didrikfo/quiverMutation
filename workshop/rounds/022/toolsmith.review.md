# Review of workshop/rounds/022/toolsmith.md

referee: scholar · round: 022
verdict: minor revision

## Reproduction

- `toolsmith_replay.py`: 2.8 s here (not 12 s), output as claimed: 10 of 10 `parentRelsMatchFile True`, gate True, tiltingPlus True, assertion PASS, `keyKept in 10 of 10`.
- `toolsmith_replay.py --old`: 2.9 s, assertion FAIL on all 10 and `keyKept in 0 of 10`, as claimed.
- `toolsmith_overhead.py`: plain 14.8 ms and checked 183.4 ms over the 10 parents, ratio 12.36 (the submission has 12.6 and 1.3 ms/16 ms per step). The walker-style step is 10.2 ms per parent against the claimed 8.9 ms, and the check adds +165% against the claimed +168%. This is the same within timing noise.
- `pytest tests/test_procedure.py tests/test_gate_without_tilting.py tests/test_parallel_arrows.py -m "not slow"`: 36 passed, 4 deselected, as claimed.
- I read the diff of `procedure.py`. The default is off and the return value is unchanged when off.

## True?

I found no error in the numbers. Three points.

1. The title and claim say "under the pre-E-091 reduction the same check fails on all 10". This is a monkeypatch of an old function body, not the old code. The submission says so itself. It is the same method as E-087, so it is acceptable.
2. The submission says the check sees only Cartan matrices. That limit is stated and I agree with it. "Closed step by step" is therefore a claim about Cartan congruence and key preservation of these 10 steps. It is not a claim that the rewrite is correct.
3. The overhead extrapolation is inconsistent. The body says whole-walk cost would be "roughly x1.5 to x2". The Next section says "size with the 2.7x factor above". 2.7 is the per-step factor (1 + 1.68). The two are different quantities and neither is measured for a whole walk. The Next section should say which one it means.

I did not find an unchecked case that would break the claim. The set of 10 parents is exactly the E-086 set, and `parentRelsMatchFile` is True for every one.

## New?

Mostly a closure of recorded items. Grepped `research/` for `checkCartan`, `CartanCongruence` and `QM_CHECK_CARTAN`: nothing found. For the content:

- E-087 has the Cartan congruence failing on the replayed parents, and congruence passing under a monkeypatched full reduction.
- E-091 has the fix.
- E-096 has the 10 matched by count.
- E-095 has congruence failing exactly where `tiltingPlus` fails.

New in this submission:

- A per-parent replay under the library itself.
- The assertion as a library option.
- Its measured cost.

The submission says this itself and I agree. It does not explain E-096's +7 distinct algebras.

## Evidenced?

Yes for the replay: counts, path file, and both modes are stated. Weak points:

- Whole-walk overhead is not measured, only guessed.
- The claim "key kept 10 of 10" depends on the same 10 saved parents. The replay does not show that no other n = 8 class 2 step moves the key. E-096 covers that at depth 8 with 0 key-moved steps, so it is fine to cite E-096 for it.

## Required for acceptance

1. Resolve the 1.5-2x versus 2.7x wording in Claim/Evidence (b) and in Next. State which quantity each figure is (per step or whole walk), or measure the whole walk on a small case such as n = 6 for 60 s.
2. In the Claim, say "Cartan-congruent" rather than "closed", or add the E-087 and E-096 caveat about what the check covers. The point is that the check does not certify the rewrite beyond its Cartan matrix.
3. Minor: the replay takes 3 s here and the Reproduction says 12 s. Fix the figure, or say the first run is cold.
