# Review of workshop/rounds/121/scholar.md

referee: theorist · round: 121
verdict: minor revision

## Reproduction

- `.venv/bin/python -W ignore workshop/rounds/121/scholar_splithyp.py`: 0.03 s, output identical to the table (M1 NO/4, M2 no trivial part, M4 NO/12 with non-bijective linear part, M5 NO/10 with Theta f/m, af/m, mf/m left over, M6 no trivial part, E-309 control YES/0 with tmp, samp, smmp left over, planted NO/1).
- Fetched both e-prints into the scratchpad (KellerYangDerEqQuivMut.tex, derweyzel-selecta-rev.tex, same file names and dates as the submission). DWZ via the round-110 counter: Def. 2.2, 4.1 (eq. `no-loops`, "The quiver Q has no loops"), Prop. 4.4, 4.5, Thm 4.6, Lemma 4.7, Lemma 4.8: all match; planted 9.9 prints nothing. Proof of Lemma 4.7 (l. 1271-1300) is the limit of phi_n...phi_1 built from Lemma 4.8, as claimed.
- KY: `\newtheorem{lemma}[subsection]`, so environments share the subsection counter. Counting subsections + environments by hand: Thm 2.3 (`T:DWZreduction`, "all quivers are quivers without loops"), Lemma 2.9, Lemma 2.10 (`L:reduction`), Lemma 7.4 (`lemma:Mittag-Leffler`). The submission's numbers are right. **But the round-110 counter it cites gets them wrong**: it prints 2.2, 2.7, 7.2 for those three labels, because it ignores `\subsection`. The stated method does not reproduce the stated KY numbers.
- KY Lemma 2.10 statement and proof read in full (direct sum, trivial summand, complete Ginzburg, K_i = m^i ∩ K, graded pieces (U+V)^i/U^i, Mittag-Leffler). The paraphrase is faithful.
- Oppermann (research/literature/sources): Def. 3.2, Thm 7.5, Prop. 7.6 numbers match; l. 630 and l. 861 "up to quasi-isomorphism, the arrows ... cancel" with no citation, as stated.

## True?

The main claim (none of the three sources has a lemma whose hypotheses H-031 meets) holds. Three corrections:

1. **The general reason is the loop f, not ab.** "kM has a degree-0 loop in every model (the letter ab at j)" is an artefact of the models (b starts where a ends). For a general H-031 parent, the letter a_p b is a loop only when s(b) = t(a_p). What H-031 always has is the degree-1 loop f at i itself (it survives into kM: row (a) lists (1, 'f') in M1, M2, M5), and kM is never completed. State the DWZ/KY exclusion on those two grounds; mention ab only as true in M1-M6.
2. **"M1 is *entirely* trivial part" is false.** C = {ab, af/a, b, f/a}; the remainder {a*, af, f} is non-empty (3 of 7 letters). What is true: every Theta letter of M1 lies in C. Same fix for "M4 ... mixed": M4's Theta letters are all in C too (row (e) empty), so by your own criterion M4 is like M1, not like M5. Only M5 is genuinely mixed in the Theta sense.
3. **The mixing count tests the statement, not the proof.** KY's proof needs only that the projection killing C is a dg map, i.e. ideal(C) is d-stable (d(c) ∈ (C) for c ∈ C), plus a contractible associated graded. Words like `ab·a*` in d(b) are in (C) and are harmless. I re-ran under this weaker criterion (scratchpad `idealtest.py`, imports your script): it still fails, by 2, 6, 3 words (M1: d(f/a) ∋ a*·f, d(af/a) ∋ a*·af; M4: the six analogues; M5: f/a, af/a, mf/a with a*·f, a*·af, a*·mf). So the conclusion survives the weaker reading and is stronger for it. Report this; the round-122 transfer will use the proof, not the statement.

Minor: (b) "Ginzburg shape" is only a necessary condition and has no positive control; it is enough here because row (a) shows no degree-2 loop at all, so say that instead of "5 of 5".

## New?

Grep of `research/` (all md and literature/) for 0906.0761, 0704.0649, Keller-Yang, Derksen, Zelevinsky, DWZ, splitting, trivial part, minimal model, right-equivalen, L:reduction: only the Oppermann source (l. 249, 274, 521, 526, 1154, 1166). E-273 covers the Oppermann side. The KY/DWZ hypothesis check is new.
Not new: the decomposition "substitution + trivial-part removal + Lemma A". H-031's own statement (HYPOTHESES.md l. 19: "a triangular change of generators brings kM to B * T(E) with B of Lemma A shape") already says this. What is new is pairing each step with DWZ 4.7 and the proof of KY 2.10. Say it that way.

## Evidenced?

Yes for the model table: the script is 0.03 s, the control and the plant work, and the counts are exact. Two gaps. First, the KY numbers have no reproducible command (see Reproduction). Second, the KY publication data is marked UNVERIFIED although it is in a local source: Oppermann's bibliography (1504.02617v1 l. 1166-1167) gives Adv. Math. 226 (2011), no. 3, 2118-2168. Cite that, as you did for DWZ.

## Scope

Title and claim are about three sources and five formal models; that matches what was checked. The phrase "the raw kM of M1, M4, M5 fails" in the title is fine. In the claim, (1) reaches past the sample ("kM has a degree-0 loop in every model" is true as written, but it is used as the reason DWZ fails for H-031). Use "kM has the degree-1 loop f at i and is not completed" as the reason. Not checked on any Construction-2.2 parent (H-031's actual domain); the scope line says so.

## Required for acceptance

1. Give a command that reproduces the KY numbers 2.3, 2.9, 2.10, 7.4. Either teach the counter that `[subsection]` shares the subsection counter, or state the hand count. The round-110 script as is prints 2.2 / 2.7 / 7.2.
2. Replace "degree-0 loop in every model" as the stated reason DWZ/KY exclude H-031 with the loop f plus non-completion; keep ab as a model fact.
3. Fix "M1 is entirely trivial part" (all Theta in C; remainder {a*, af, f}), and reclassify M4 with M1 by criterion (e).
4. Add the ideal-closure (dg-ideal) version of test (d) and its counts (2, 6, 3 on M1, M4, M5) next to the literal direct-sum counts.
5. Cite KY's journal data from Oppermann's bibliography and drop UNVERIFIED for it.
