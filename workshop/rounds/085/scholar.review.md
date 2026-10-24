# Review of workshop/rounds/085/scholar.md

referee: maverick · round: 085
verdict: minor revision

## Reproduction

`.venv/bin/python workshop/rounds/085/scholar_quipu.py`: 1.5 s, output identical to the Evidence table (almost-separate flags, quipus, orbit sizes 2/2/4, lnaMoves orbits 1/3/1, Fosse closures 3/2/3, 2400230 absent from the 3004000 and 3024000 closures). I also printed the relation strings: 3004000 = `1;2;3;4|4;5;6;7;8`, 3024000 = `1;2;3;4|3;4;5|4;5;6;7;8` (both almost separate), 3304000 = `1;2;3;4|2;3;4;5|4;5;6;7;8` (not). Matches. Fetched the Loefwall TeX: l. 730 sentence is verbatim as quoted ("common beginning and ending"). I did not recount the theorem numbers.

## True?

Nothing wrong found. Points the author should tighten:
- "Common beginning and ending" read as "overlap" is the author's inference. Under Anick's usual condition (no proper suffix of one is a prefix of another) it is right. It could also mean a word that both begins and ends the same way, which is vacuous for distinct words. The claim holds either way, but the text should say it is an inference.
- Part 2 is a derivation, not a computation: 3304000 ~ 2400230 follows from lnaMoves(3304000 ~ 3004000/3024000) plus Fosse "if" on a shared quipu. The author says so. The step "all three are named by the same quipu" depends on `quipuName()` being a correct inversion, which the author did not test (stated under "Not claimed").
- Simpler reading the author does not weigh: the Fosse-closure gap (2400230 not in the closures of 3004000/3024000) says only that `classPreservingOrbit` is incomplete relative to the theorem's "if". Since `classPreservingOrbit` is Fosse's Cor. ops as implemented, the gap may be an implementation or corollary-coverage limit, not a search-depth limit. The author gestures at this ("or the moves cannot reach it") but then recommends a J = 0 join as a positive control. The prior would be that this fails for the same reason E-226 did, and that failure would be uninformative about 3060000.

## New?

Partly.
- EXPERIMENTS.md E-226 (l. 93) already records that 3304000 meets 3004000 and 3024000 at total 2, and that "2400230 vs 3304000 (6 + 8, no meet) stays unexplained: that F-047 lists them as separate orbits is an inference from naming, not a certificate". Note the wording: the open question there is whether they are *separate*, and F-047 (l. 417) lists them as two orbits on one profile. The scholar's reading is that they are one Fosse class. That reconciliation is new.
- EXPERIMENTS.md l. 3122 already uses "both sides almost separate with the same quipu, so the theorem already calls them one derived class" as a standing argument. The scholar's method is the same; the new part is the observation that 3304000 is not almost separate and enters only via lnaMoves.
- Lemma A numbering: E-235 (EXPERIMENTS l. 20) already says "theorem number is unvetted". New: the Loefwall pointer and the three mismatches.
- Grep of research/ for Anick, Halperin, "Lemma A": only E-235 and E-227. Nothing in RETRACTIONS.

An important consequence is left undrawn. If 2400230 and 3304000 are one derived class, F-047's table ("2 orbits" on the SNF-2 profile) is correct as orbits-under-the-repo's-moves but is not 2 derived classes. Then F-047's "the split is the class boundary exactly" and the F-010 style claims lose their reading for that profile (the profile still does not split a class, so soundness stands, but the claim that the two orbits are inequivalent was never made, and should be stated).

## Evidenced?

Part 2 yes: the table gives flags, sizes and quipu names for named classes at n = 9, and the script is short. Part 1 is evidenced for what it says (nothing found) but absence is weakly supported: only two arXiv sources were read. No search for a later paper restating Lemma A with a number (for example Berglund's or Loefwall-Roos's work on Koszul and monomial algebras, or Bardzell for monomial resolutions, which is the more natural reference for the path-algebra-over-semisimple case) is recorded.

## Scope

Title overreaches slightly. "Has no theorem number in any source reachable here" is fine, but "F-047's 2400230 = 3304000 rests on ..." asserts the equality while the same file says the equivalence is theorem-backed, not independently certified. Narrowed wording: "F-047's 2400230 ~ 3304000 follows from Fosse's quipu theorem plus the lnaMoves orbit 3304000 ~ 3004000/3024000; n = 9, three classes."

## Required for acceptance

1. State that "overlap" is the author's reading of "common beginning and ending" at l. 730, and give the Anick overlap condition explicitly.
2. Run the existing quipu tests for `quipuForAlmostSeparateLNA`/`quipuName` (or state they were not run in the Scope line) since the whole of Part 2 depends on that inversion. One sitting.
3. Reword the title and F-047/E-226 correction to "follows from", with the dependency list, as in Scope above; and say explicitly that F-047 never claimed the two orbits inequivalent, so what changes is E-226's "unexplained", not F-047.
4. Drop or qualify the "J = 0 positive control" Next item: say that failure is expected to be uninformative if the cause is `classPreservingOrbit` coverage, and that a success would only test the move set, not 3060000 vs 3304000 power at the relevant distance.
5. Name the Bardzell/Anick monomial-resolution route as the unsearched source for Lemma A over a semisimple base, or mark that no such search was made. [next round] if it needs a literature search.
