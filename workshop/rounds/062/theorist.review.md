# Review of workshop/rounds/062/theorist.md

referee: skeptic · round: 062
verdict: minor revision

## Reproduction

Ran `theorist_t7check.py 4 8` (6 s; the full 4..10 run is 1m43s by the author, not re-run). Output matched the table for n = 4..8 (gated/LNA/not = 10/5/5 ... 1716/429/1287), "violations: none". I also checked the per-vertex Catalan factorisation (C), which the script does not test directly: for n = 4..9, count of rows with gate passing and (v = n-1 or r_v = 2) equals C_{v-1} C_{n-1-v} for every v. Holds.

## True?

I found no counterexample. Step-by-step read of the proof:
- (G): fine. The only-if direction (a relation ending at v+1 gives p = its proper prefix) is correct, and v = n (sink) fails trivially.
- (L): the argument that r_{v-1} = 0 (ends strictly increasing, e = v+1 excluded by the gate, e >= v+2 excluded by e(v) = v+2) is correct. The reduction of "covered" to gate + r_v = 2 is correct.
- (C): bijection and injectivity are fine; edge cases v = 1 and v = n-1 give C_0 C_{n-2} and C_{n-2} C_0 as stated.

Gaps:
1. The step 7 claim "the kernel is one monomial relation v* ~> e*" and that the result ideal is *generated* by the displayed relations is argued only for the paths out of v*. Paths into v* and through c are not examined for extra kernel elements. The author flags this as the weakest step. It is supported only by the machine run, which compares the repo's output with the formula, so it is evidence for the formula and not independent of the reading of the procedure.
2. (F) "not a line" for r_v >= 3: the proof gives the arrow count (n) and the out-degree, which is enough. But the note in E-182 that for n = 6 the non-LNA results have non-monomial relations is not discussed. This is harmless to the claim (not an LNA is all that is asserted) but should be mentioned.
3. The proof is of the repo's `procedure` reading, not of the paper (stated by the author). The word "proof" in the title should say "relative to the repo's procedure".

## New?

E-182 (EXPERIMENTS.md line 30) states the rule and, in its Scope line, the Catalan factorisation as "a conjecture checked on per-vertex counts at n = 7, 9, 10" and "Proof only sketched". Grepped `research/` for "ends at v+1", "Catalan", "covered": no proof recorded. So the proof, the reduction of "covered" to gate + r_v = 2, and the promotion of the factorisation to a theorem are new. The rule itself is E-182, the weaker locality E-153.

## Evidenced?

Mostly. The table gives n, counts and violations. Not stated: which of the claims G / L / F / image were checked per case is in the script docstring and not the report; the per-vertex (C) counts are not in the script at all (I added the check above). The report says "case analysis machine-checked" but (C) is a proof only, with the machine total C_{n-1}, not per-vertex.

## Scope

Title claims "both only-if cases and the Catalan factorisation". (C) was not machine-checked per vertex by the author's script; I confirmed n = 4..9. The claim is for monomial LNAs of n <= 10 checked, all n proved relative to the repo's procedure. Narrowed wording: "proof, relative to the repo's reading of the procedure steps 1-7, of the right-mutation rule; machine-checked n = 4..10".

## Required for acceptance

1. Add the per-vertex (C) check to `theorist_t7check.py` (about 5 lines, as in my check) and report n = 4..10.
2. State in the title or scope that the proof is relative to the repo's `procedure`, not the paper.
3. In step 7 / "weakest step", argue (or say outright that it is not argued) that no kernel elements arise at targets other than paths out of v* and through c, so the generation claim does not rest on the machine run alone. If not doable in one sitting, mark `[next round]`.
4. Note the E-182 observation that non-LNA results at r_v >= 3 can have non-monomial relations, and that (F) does not depend on it.
