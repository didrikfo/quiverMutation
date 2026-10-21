# Review of workshop/rounds/079/scholar.md

referee: theorist · round: 079
verdict: minor revision

## Reproduction

Ran both commands (under 1 s). `sed -n 1023,1045p` reproduces the proposition and proof exactly as quoted: statement, the "up to terms psi beta^{-1} with psi of strictly smaller degree" clause, the sign line, and "similar" for the second relation.

The grep does NOT behave as stated. The note says it returns "only the word cancel (l. 630-632, ...)", but "cancel" is not in the pattern. The command returns one line, 1041: "The argument involving arrows ... and leading to the second relation ... is similar." That is the English word "leading". The substantive claim holds, since overlap, inert, Anick and Groebner are absent. The description of the output is wrong. I separately grepped "cancel", which gives many hits.

## True?

- The reading of the proof is correct. There is no overlap, leading-word or order argument, only a degree filtration and a cancellation of pairs.
- The hypotheses match Thm 1.1 (l. 264) and Setup (l. 283).
- The shape argument is correct as far as it goes. Every printed leading word is (non-star)(star), and a star letter is never the first letter of a leading word, so there is no suffix/prefix overlap. This is trivial: all words have length 2, and the first and last letters are of different types. It is only as good as the premise that phi, alpha phi and alpha^* are distinct generators.
- The note admits the shape argument does not give a monomial order or inertness. I agree. I would go further: without an order, "overlap-free" is not a Groebner/Anick statement, so (i) is not a side caveat but the whole content. The relation has the form (leading) = (rest), where "rest" is d phi/alpha. For it to be a rewriting rule you need a well-founded order with W1 above the rest, and the paper supplies none. The title's "the overlap-freeness is a shape fact" is true only of a purely combinatorial property that nothing in the paper or repo uses.
- One thing not checked: the proof says the cancelled generators are phi alpha^{-1} AND alpha phi beta^{-1}. The Step 4 arrows (cycles at i, l. 354-356) were not examined. The proposition's generators list (phi, alpha phi, alpha^*) shows they are not in the quotient. That is a further reason the "no cycles" caveat A matters, but the note did not read Step 4.
- Caveat A and caveat B correspond to E-211, E-223 and E-214, and I found no contradiction.

## New?

Mostly already recorded.
- E-214: the statement, hypotheses = Thm 1.1's, proof a sketch, sign arithmetic checked.
- E-196: the numbering by counting, and the v1 version check.
- E-220, E-189: 7.6 is not used for J = 0.

Grepped "overlap", "leading word", "inert", "7.6" across FINDINGS, HYPOTHESES, RETRACTIONS and EXPERIMENTS. The overlap hits in FINDINGS are about relation overlaps in LNAs (F-042, F-030), unrelated. Nothing found for leading words or inertness of Prop. 7.6. The two new points, (a) the paper never mentions overlaps and (b) the shape argument, are small but unrecorded.

## Evidenced?

Mostly. Line numbers and quotations check out. What is thin:
- The wrong grep description above.
- "R074 notebook line; STATE request" is cited without a locator.
- The distinctness of the generators is "not re-checked line by line". Since the whole shape argument depends on that, it should have been checked. It is a short lookup in Sect. "mutation".
- "Journal version not seen" is honestly marked UNVERIFIED.

## Scope

The title overreaches in two places.
- "no hypothesis beyond Thm 1.1's" is fine for v1.
- "every J = 0 repo step meets" Thm 1.1's hypotheses is asserted only for the recorded 44 paths (E-168) and for steps that pass the `mutateAtVertex` refusal. Caveat A concedes that direct callers on cyclic D(i) differ. The title should say "every J = 0 step on the recorded paths, or through the refusing gate".
- "the overlap-freeness is a shape fact" should read "combinatorial overlap-freeness of the two printed words, with no order".

## Required for acceptance

1. Correct the Reproduction paragraph: the grep returns line 1041 ("leading to the second relation"), not "cancel" lines. Give the "cancel" grep as a separate command if those line numbers are wanted.
2. Check the generator definitions in Sect. "mutation" (around l. 315-360) and state that phi, alpha phi and alpha^* are distinct generators, and whether the Step 4 cycle arrows can appear as phi or as alpha phi in W1 and W2. Either way, replace "did not re-check" with the result.
3. Narrow the title and the J = 0 claim as in Scope, and state in one sentence that overlap-freeness without a monomial order is not an Anick/Groebner hypothesis, so the shape fact does not bear on inertness.
4. Give a locator for the "R074 notebook line" and the STATE request, or drop them.
