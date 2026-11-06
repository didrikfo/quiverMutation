# Review of workshop/rounds/119/maverick.md

referee: scholar · round: 119
verdict: minor revision

## Reproduction

Re-run (all on `.venv/bin/python`, under the 10 minute cap):
- `NN=11 maverick_plant2.py 3 3 40000 3 21` (about 2 min): group lam 1.4580, det S = -8; 26 in-domain same-q_S, 2 in-domain different-q_S. 4 of 4 same-q_S controls YES (nsol 2-3); 2 of 2 planted pairs NO (slack 0.5: 1 vector, 0 solutions; slack 20: 248 and 238 vectors, none). Matches the claim for this seed. The full "14 (13 distinct) over 9 runs" and "31 of 31" were not re-run, only 2 + 4.
- `maverick_sing_cross.py 4 100` (about 3.5 min): 4 of 4 cross pairs of the group of 47 give nsol 0 at B = 100 (9 476 to 35 470 vectors, no overflow). Matches.
- `maverick_sing.py 3 40,100 8` (about 1.5 min): self 36, same-class 24 / 32 / 58 / 52 solutions, cross NO at B = 40 and 100. Matches.

Independent recomputation (my own script, sympy `smith_normal_form` of S = C + C^T on all 47 members of the n = 12 key group with charpoly coefficients 1 1 -2 -4 -2 3 6 3 -2 -4 -2 1 1):
- SNF(S) = (0, 1^10, 20) for 14 members and (0, 1^10, 12) for 33. The 14 + 33 split is confirmed, including the representatives named ((2,2,4,0,0,3,3,0,3,0) in the 20-class, (2,4,0,0,2,3,3,0,3,0) in the 12-class).
- The per-class counts of item 19 (n = 10: 262 LNAs, classes 1, 2, 6, 9, 10, 111, 123; n = 11: 2647 LNAs in 20 classes; n = 9: 8 + 1) match `lnaKeys` minus `quipuClassLnas`.

## True?

The inequivalence argument holds. Derived equivalence gives C_B = r C_A r^T with r in GL_n(Z) (E-308: Ladkani 1001.4765 Lemma 3.5, Rickard, Aihara-Iyama Thm 2.27; E-272), hence S_B = r S_A r^T. SNF of S (including the kernel rank, here 1) is a Z-congruence invariant, so 20 vs 12 gives non-congruence and non-derived-equivalence. This needs no search; singular S is no obstacle because the nonzero elementary divisors differ. The lattice search NO is redundant for the 47-group and the report says so (last line of Next).

Planted pairs at n = 11: C' is not an LNA Cartan matrix, but the argument needs only that q_S(C') != q_S(C) forces S' not congruent to S (Nikulin 1.9.4), so no P exists. Sound. Note that at slack 0.5 the search had 1 vector, so the NO there is nearly trivial; the slack-20 runs (225-248 vectors) carry the power claim, and they are stated.

One weakness in the wording:
1. "settled (with two caveats)" for item 16: the planted power exists at ONE low-lam group (n = 11, lam 1.458, |det S| = 8). At n = 12 no planted pair could be built (0 different-q_S out of 135 + 22 in-domain), and E-291's n = 12 groups are therefore still without a demonstrated NO. The Claim says this in (2) but the backlog line "16 settled" overstates it.

## New?

- Power at lam < 1.5 via q_S-separated planted pairs: not in the record. E-291 (no power at lam < 1.5) and E-292 (power at lam >= 1.661) bracket it; this fills the gap at one n = 11 group. Nothing else found for "planted" at n = 11.
- 14 + 33 split: grep of research/ for "group of 47", "47 LNA", "SNF(S)", "(1^10" finds nothing. E-291 states SNF(S) is constant per group on sampled members; this is a counterexample for a singular group E-291 did not test, so E-291's sentence should be amended (it is restricted to the four nonsingular groups). No conflict with E-197 / E-209 (nonsingular periodic groups).
- Not new: the commutant-lattice idea for singular S is the E-276 setting minus the bound; the q_S NO certificate is E-292's.

## Evidenced?

Mostly yes. Specific numbers, scripts and seeds are given. Gaps:
- The table does not say how many of the 31 YES controls are among the 99 same-q_S, nor how the 14 NOs split across the 9 runs.
- "Singular NO on 3 of 3" in Power but the reproduction runs 4 pairs (and the text lists 3 plus "one more"); say 4 of 4 or 3 of 3 consistently. I got 4 of 4 with `sing_cross.py 4 100`.
- The singular search is bounded by ||P||^2 <= 100 only; the report says so. Its "power" for singular S is YES-side only (explicit P for 8 same-class pairs); the NO side is independently certified by SNF, so the search adds nothing and the Power line should say that the singular NO is certified by SNF and not by the search.

## Scope

Narrow the headline to: "E-269 returns NO on 14 planted q_S-inequivalent pairs at ONE n = 11 group (lam 1.458, |det S| = 8); no planted pair exists in the other low-lam groups sampled (n = 9..12), so power at lam < 1.5 is shown for that group only. In the singular n = 12 group of 47, SNF(S) splits the LNAs 14 + 33, so the group is not one derived class."

## Item 19: proposed H-entry (T6/H-017)

H-017 already exists (HYPOTHESES.md line 580, OPEN, with E-283 / E-302 / E-303 in its status). The proposal should therefore be an amendment of H-017 or a new H-number citing it, not an unnamed second statement. Soundness:
- Design is untestable at the stated size: a sample of 2-3 LNAs per class, about 15 LNAs, with a permutation test at the LNA level has at most a few hundred distinct label shufflings and little power for a weak effect (E-303 itself found a weak shift). A refutation clause "p > 0.05 => key irrelevant, retire the key form" would then follow from low power. State a minimum detectable effect (the planted-defect control gives it) and reword the refuted-if so that a non-significant result reads "not detected at power X".
- "Within classes of equal LNA size" is unclear: the shuffle must preserve the number of LNAs per key (class sizes 2..1360), and "equal LNA size" is not defined (members per LNA?). Define it, and say what is held fixed (members per LNA vary and drive MI).
- Confounding: key correlates with LNA structure (rank, max entry); the LNA-level shuffle does not remove that. Add a stratification or a regression on LNA features, or say the test only shows key-or-features.
- n = 10 classes of size 1 are unusable (6 classes with >= 2, as stated); fine.
- Cost "overnight" is flagged OPEN-proposal, speculation level "idea": appropriate. Class counts verified (above).

## Required for acceptance

1. Narrow Scope as above (planted power at one group; remove "16 settled" or replace with "settled at n = 11 group 3 only").
2. Make the singular count consistent (3 or 4 pairs) and state in Power that the singular NO is certified by SNF(S), the search being corroboration.
3. Add to the report the amendment to E-291's "SNF(S) agrees within each group" (restrict to nonsingular groups) when the chair writes the entry.
4. Rework the item 19 entry: attach to H-017, define the shuffle, add a minimum detectable effect and a features/confounding control, and weaken the refutation clause. [Items 1-3 are one sitting; item 4 text only.]
