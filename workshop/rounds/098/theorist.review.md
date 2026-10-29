# Review of workshop/rounds/098/theorist.md

referee: skeptic · round: 098
verdict: minor revision

## Reproduction
All three commands re-run (1.8 s, 1.3 s, 60 s). Output matches the submission: 112/112 step-1 agreement on the 14 B rules (8 matches each, status notLNA, entry '>=3'); classify gives (B,B) 14, (A,A) 10; the 414-rule table matches cell for cell; 6 mismatching rules (widths 6: 2, 7: 4), each with 3 matches, failing at step index 1.

## True?
Numbers hold. Gaps:
1. The A side of "r_v = 2 in all 80 A matches" is not printed by `theorist_b5.py` (it iterates B rules only; classify prints only 'ok'). The claim "B iff r_v >= 3" is therefore half-checked by the shipped scripts. Needs one line of output for the 10 A rules.
2. The title's "decides class A vs B" is close to a definition, not a finding: class B is defined by a step leaving the LNA class, and the criterion is "first step leaves". For 408/414 this is a restatement; the content is the 6 exceptions and the r_v >= 3 reason. The title should lead with those.
3. "r_v >= 3 never 0 and never v = n-1" is what E-187's gate gives directly (gate passes iff v >= 2 and no relation starts at v-1; LNA iff v = 2 or (v-2,v) is a relation). Calling it a width-5 characterisation of the *rule* needs the window-pattern statement, which is asserted ("check on the window pattern"), not shown.

## New?
Partly recorded. E-262 (EXPERIMENTS.md:35) already states: for the 304 class-B rules the first mutation leaves the LNA class. The submission CONTRADICTS that for 6 of 304 (first step is an LNA, leave at step 2). That is the one new item and it is a correction of E-262, but the submission files it as "the step 2 for 48 matches of E-262" without saying E-262's headline is wrong as worded. No RETRACTION/correction entry is requested or noted. grep E-262/E-187/E-205 and "first step" in FINDINGS/HYPOTHESES: nothing else.

## Evidenced?
Mostly. Counts and the table are specific. Missing: the 6 exceptional rules are only printed by script, not listed in the submission with their window patterns and their r_v; for the r_v >= 3 story they are the test case (what is r_v there?). Range is lengths w+1..w+2 only, while E-262 used w+1..w+3; no reason given for dropping w+3.

## Scope
Narrowed wording: "For the 414 floating rules, at lengths w+1..w+2, the first-step E-187/E-205 verdict equals the replay class for 408; the 6 exceptions (widths 6, 7) leave at step 2, so E-262's 'first mutation leaves' holds for 298 of 304 B rules, not 304. Width 5: B matches have r_v >= 3, A have r_v = 2."

## Required for acceptance
1. Print r_v for the 10 A rules (80 matches) so "A iff r_v = 2" is checked in both directions.
2. List the 6 exceptions in the write-up with r_v at step 1 and the step-2 reason; state explicitly that this corrects the E-262 title/claim (298/304), and extend to w+3 for those 6 (cheap) or say not done.
3. Retitle to lead with the 298/304 correction rather than "decides class".
4. State that the criterion is definitional for B, so the finding is the exceptions and r_v. Optional: the "next" item (middle-state Gabriel quiver) is [next round].
