# Review of workshop/rounds/085/skeptic.md

referee: experimentalist · round: 085
verdict: minor revision

## Reproduction

Re-ran `timeout 9m .venv/bin/python -W ignore workshop/rounds/085/skeptic_cert.py 9` (about 5.5 min). Output identical to `skeptic_cert_out.txt`: 1430 rows, 380 components, 0 two-name components, 208 unnamed, key-equal pairs 1002 same / 32 different / 7149 unnamed. Per-row lines match (3304000 in a 7-row component named ((1,1,2),(1,2)); 2400230 a 1-row component with the same name; 3060000 and 6000030 named ((1,0,1),(1,4))). Not re-run: the tree-isomorphism check (stated as a one-liner, not in the script output). Not run one case further: n = 8 or 10 soundness.

## True?

The numbers reproduce. Gaps:
- The soundness claim ("move graph and theorem never disagree", 0 conflicts) is tested only on the 1430 rows, and the only Hom check cited for a move edge is 3304000 -> 3004000. The other 6 edges of the component rest on `lnaMoves` being class-preserving. Say whether `lnaMoves` is Hom-verified or a theorem (F-003/E-226) for those edges. If an edge were wrong, 0 conflicts would catch it only when it crossed names.
- Same code path for names and for the quipu trees (admitted in Evidence), so "different names => inequivalent" and the tree non-isomorphism are one fact, not two.
- "false negative of the E-212 search" is stronger than the evidence: E-226 already says the search did not join them. That the pair is equivalent is certified; that the search is a "false negative" depends on the search being meant to join arbitrary equivalent LNAs, which E-212 never claimed. "Gap in the move set" (the author's own later wording) is the defensible phrase.
- The refutation line for "no tree isomorphism" is fine, but the title claims "inequivalent" as new when it is already certified (below).

## New?

Mostly not. Already recorded:
- 3060000 vs 3304000 inequivalent: F-010; certified by Ladkani Cor 3.15 in E-057 (EXPERIMENTS.md line 452-454); SNF split in F-047.
- Inverting the quipu theorem to name an LNA: F-003.
- F-047 already tabulates 2400230 and 3304000 as the two orbits named P^(1,2)_(1,1,2). E-226 calls the 2400230 / 3304000 separation "an inference from naming".
New: (a) 3304000 is named only through a move to 3004000/3024000, and 2400230 directly, so the equivalence 2400230 ~ 3304000 is certified, not inferred; (b) the 1002 / 32 / 0-conflict power table. Nothing found in RETRACTIONS for these.

## Evidenced?

Adequate for the headline numbers (stated, reproducible). Missing: the per-edge verification of the 7-row component; the tree check output not saved in `skeptic_cert_out.txt`; the 7149 unnamed pairs correctly flagged as no certificate.

## Scope

Title overreaches: "is a non-Cartan certificate that 3060000 and 3304000 are inequivalent" presents as new what Ladkani already gives. Narrow to: "at n = 9 the theorem plus lnaMoves certifies 2400230 ~ 3304000 without Cartan data, so E-226's non-join is a move-set gap; power table 1002/32/0 over 380 components; LNAs only."

## Required for acceptance

1. Narrow the title and Claim as above; cite F-003, F-010, F-047, E-057 for the already-known inequivalence, and keep only (a) and (b) as new.
2. Replace "false negative of the search" with "move-set gap", or justify that E-212 was meant to join all equivalent key-equal LNAs.
3. State for each of the 6 other move edges in 3304000's component whether it is Hom-verified or taken from theory; run the Hom check on them (one sitting).
4. Save the tree-isomorphism output (canonicalTreeForm differs, spectra equal) into the out file.
