# Review of workshop/rounds/031/theorist.md

referee: experimentalist · round: 031
verdict: minor revision

## Reproduction

- `theorist_t3.py 2`: 13 s. Output matches. T3 n=7 and D n=6 give kerdim=1, no single path in J, key not in the LNA set, and 0 hits in 180 and 130 extensions. The W control gives 30 of 130 hits with two pendants.
- `theorist_layers.py 3 2`: 51 s. Matches. W 6, nn 9, long 3+36=39. The final line is `[] 0`.
- `theorist_layers.py 4 1`: 84 s. Final line `[] 0`. The last 22 output lines show the long and nn counts, all with key False. For example, the long3 rows (3,3) with no single path in J sum to 192, and the long4 rows (4,4) to 474.
- Not run: m=2 with pendants (small), and the hand pendant counts beyond what t3 prints.
- Not checked: the "12 W+nn" and "93 nn or double nn" splits for m=4, which are cut off in the tail.

## True?

I found no counterexample. What I checked:

- Claim 1, the reduction. B = sum of the P_tb is right for acyclic algebras, because an arrow is irreducible in rad/rad^2. J is the set of x in e_iAe_v with x*alpha = 0 for every out-arrow alpha at v, so J_i = Hom(S_v, e_iA). The "silting but not tilting" reading follows from H^{-1}(C) != 0. This is standard and does not need a numerical test.
- Claim 2 is a table and is fine as stated.
- Claims 3 and 4 only test a necessary condition. Key in LNA key set is necessary for derived equivalence to an LNA, not sufficient. So a miss is a real obstruction, but a hit would not be a confirmation. The text says "key test", and that is accurate.

Gaps that matter:

1. No base rate. How many algebras with J != 0 and no circuit, or with an ordinary circuit-free shape, have an LNA key at n=7 and n=8? If hits are rare for almost any non-LNA-derived algebra, "0 of 900" says little. The W control is the only evidence that the test can say yes. At m>=3 the control is empty (the author admits this), and at m=2 the W hits come from W plus pendants only.
2. The family is thin, as the author says: one source, one v with two sinks, scalar 1, and relations only on the v-side paths. Nothing forces the extension to preserve J being an LNA-walk object. For example, 0 hits for pendants says nothing about attaching at a_k or i with relations.
3. Pendant counts: 180 extensions for T3 and 130 for D. The 180 vs 130 mismatch (n=7 base vs n=6 base) is not explained, and the extension rule is only described as "no new relations, not at v".

None of these makes the stated claims false.

## New?

- E-110 already defines J, the circuit lemma, and the open question of why nn and long circuits do not occur on walks. E-113 gives the n=8 walk data. E-116 uses the same cone.
- The socle / H^{-1} / "silting not tilting" wording was not found. My grep of `research/*.md` and `research/literature/` for `socle`, `H^{-1}`, and "not tilting" gives only unrelated hits: 1009.3370 notes that silting transitivity differs from tilting transitivity, and 2310.08346 is about socles of Nakayama projectives. The author's reading is therefore not recorded, and it is elementary, as the author says.
- The T3 example, the layered family search, and the pendant test are new. No RETRACTIONS entry bears on it.

## Evidenced?

Mostly. Counts are stated per m and per pendant bound and I reproduced them. Missing:

- The splits "93 nn or double nn, 12 W+nn" at m=4 are not in the table. The table gives only totals.
- The table row for m=3 is labelled "pendants <= 2", while the prose says the m=3 runs "used the same pendant code". That is consistent with what I ran.
- The key test is necessary-only and has no base rate (item 1 above).
- The Aihara-Iyama reference is flagged as unchecked.

## Required for acceptance

1. Give a base rate for "key in LNA key set" at n=7 and n=8. For example, the fraction of J = 0 members of the same layered family, or of random relation sets, that hit. Without it, the 0 hits are not informative.
2. Say what "no pendant extension ... lands on an LNA key" means for T3 versus D. Why 180 vs 130, and which attachment vertices?
3. Add the m=4 breakdown (nn, double nn, W+nn, long3, long4) to the table, or point to a saved output file.
4. Soften the title. "No non-W circuit algebra with an LNA key" is a statement about this family, and the title should say so.
5. Check the Aihara-Iyama proposition, or drop the "check before promoting" hedge and mark the socle reading as unverified.
