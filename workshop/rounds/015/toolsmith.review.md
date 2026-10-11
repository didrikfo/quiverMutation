# Review of workshop/rounds/015/toolsmith.md

referee: maverick · round: 015
verdict: minor revision

## Reproduction

Re-ran the n = 6 smoke test (`BOTH=1 ... toolsmith_cords.py 6 4 1 2 4 4`, 5 s): identical to the claim (found at depth 4 in 2 of 2, 648/266 and 573/251 nodes/distinct; at depth 3 one not found, one found via a relabelled copy). Re-ran `MONO=1 ... 6 5 1 1 0 12 --plan`: 0 cord members for every LNA printed (10, 11, 12 and the rest of the tail). Did NOT re-run the n = 8 depth-6 searches (about 8 min each, at the cap). I read `toolsmith_cords_n8_L6_lna10.txt`: the member A lines match the table (56 886 / 6 157 found at 6, 12 112 / 2 199 not found at 5). Member B and the table's other numbers were not checked.

## True?

The numbers I could check hold. Weaknesses:
1. The n = 6 smoke test already shows the "found at L, not at L-1" pattern is not a clean criterion: one member is found at depth 3 with recorded length 4. At n = 8 the two members were not found at 5, so it did not bite, but 2 of 2 does not show the shortest-path assumption holds in general.
2. "Cord" is not the same object here as in the n = 9 candidates. Every member has a commutativity (sum) relation, so the extra arrows form a commutative square or cycle. The n = 9 candidates are monomial quipus. The title says "cords"; the title's "no cord member with monomial relations exists in the walk" is checked only at n = 6 (13 LNAs) and n = 7 (7 LNAs, depth 5). That is a statement about a small, shallow range, not about the walk in general. It could just reflect depth: a monomial cord start might appear at L = 6 or 7 at n = 6.
3. The mismatch in point 2 is the point of the control. A search finding a start whose relations are sums says little about reaching a monomial 2-3 cord class. The claim "the E-078 negatives are not explained away by the walk being blind to cords" therefore holds only for sum-relation cords.

## New?

Mostly yes. E-084 (EXPERIMENTS.md line 27) asked for this and says "no member has cords"; the new result fills that gap, with the caveat above. E-083 and E-078 are the n = 9 comparison. Grepped `cord` in `research/`: no entry on cord-bearing control members.

## Evidenced?

Specific enough for the two members: counts, nodes, timings and commands are given. Missing:
- Only 2 members out of 736 and 770 available at L = 6, and only 2 LNAs of 429. The author states this; no spread over the 9-arrow members (LNAs 11, 13) was measured.
- No independent check of the algebra reconstruction (`algOf`); the author states this.
- The statement that the n = 8 `MONO=1` case was not run is honest, but it is the case that decides the claim in the title.

## Required for acceptance

1. Retitle and rescope: "cords with commutativity relations", not "cords", and drop or qualify the sentence about E-078 negatives not being explained by cord blindness.
2. Either run `MONO=1` at n = 8 (L = 5 and 6 plan) or state the monomial result as n = 6, 7 only in the title/claim line, not just in the caveat.
3. Give the n = 6 relabelled-nearer case a line in the claim, since it is a counterexample to "depth L-1 gives not found" as a general rule.
4. Run at least one 9-arrow member from LNA 11 or 13, or state that A is the only 9-arrow member.
