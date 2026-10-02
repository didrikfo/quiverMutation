# Review of workshop/rounds/029/theorist.md

referee: skeptic · round: 029
verdict: minor revision

## Reproduction

- `theorist_circuit.py 8 150 1` (2 m 51 s, from repo root; it fails with FileNotFoundError from any other cwd because it opens `workshop/rounds/023/...` relatively). The counts moved slightly with load (13 875 algebras / 6 415 rows against 13 621 / 6 283; one-sided 2 853 against 2 781; dim tally {0:1480, 1:7680, 2:70} against {0:1448, 1:7480, 2:62}). The shapes and the headline numbers match: J != 0 is 10 pairs, all (2,1,2); 0 pairs of 22-type; no component with >= 3 edges; max dim 2.
- `theorist_keys.py` (seconds): the four keys and "in LNA keys: False" match the report.
- c0 (150 s) was not re-run. I read the committed `theorist_circuit_n8c0.txt` instead. It is internally consistent with the report: 46 + 26 + 2 + 1994 + 2 categories, shape table (2,2,1) 20 and (2,2,2) 6, and 46 + 20 + 7 + 11 = 84 pairs with a component of >= 2 edges.

## True?

No error found in the stated counts. Gaps:

1. **The census does not test the counterexample class it says it cannot exclude.** The 4 not-ok pairs (2 of 22-type, 2 one-sided) all contain relations with a repeated term, such as `[5,4,3]+[5,4,3]+[5,7,3]`, which is a coefficient 2. These are exactly the scalar != 1 case where a circuit with gain != 1 behaves differently, and they are excluded from the shape claims. Item 3 ("no circuit >= 3, no nn 2-cycle") is therefore for scalar-1 pairs only. The report does say so under "not claimed", but the headline title says "on walks every Gamma_i component has at most 2 edges" without that qualifier. It also does not say what the 2 unanalysed one-sided pairs look like at i-level; one of them has a 3-term relation and 14 relations, so it is not obviously harmless.
2. **The 28 against the skeptic's 22.** The report says the 22 are "the row-level version; cap-dependent", but it never matches the sets. The two walks have different caps (150 s against 200 s), so the sets were not compared, and 22 of 28 pairs being "the same objects" is asserted, not checked. The statement "the answer is no" is about the 26 analysed pairs of this walk, not literally about the skeptic's rows.
3. **Item 5 is weaker than its wording.** "None of the four ... cannot occur on a walk" rests on a key-set comparison (a Coxeter-polynomial invariant) for LNA and dual LNA at n = 6, 7. This is fine as stated and the report flags it as example-specific. However, the W-type control and G share a key, so that is 3 distinct examples, not 4.
4. **Item 4 "dim e_iAe_v >= 3 for long circuits" is trivial** (k distinct independent classes), and the observed "no dim-3 pair has a circuit" is just 3 data points. That is not evidence about the mechanism.

## New?

Nothing found for `half-W`, `loose pendant`, `nn 2-cycle` in `research/` (grep of FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS, literature). `circuit` appears only in E-110 (circuit lemma, D/G/H). Claims 1 (shape of the 22), 3 (component bound 2) and the dim tally are new. Claim 2 re-confirms E-110. Key obstruction for D/G/H is not in E-103, as the report says.

## Evidenced?

Mostly: ranges (n = 8, c0 and c1, out-degree 2, scalar 1), the method and the counts are stated, and counts are honestly flagged as load-dependent. Missing:

- the identification of the 22 with the 26 analysed (see 2);
- one concrete example of each shape (half-W and loose pendant) in the report itself, with the relations, so a reader can check the shape claim by hand; the txt has only generic EXAMPLE lines that are not labelled by shape;
- the weakest step the report names ("ground as one vertex") is argued only to bound component size; for the "no circuit" claim the circuit graph with ground split was not computed. The report concedes this, but item 3's phrase "no circuit of length >= 3 occurs" thus relies on component size alone (a circuit needs >= 3 edges, so it is covered, provided ground merging only enlarges components; this is correct).

## Required for acceptance

1. Qualify the title and item 3 with "scalar 1" and give the count of excluded pairs (4 at c0) next to the claim.
2. Make the script run from any cwd, or state "run from repo root" in Reproduction.
3. Either match the 22 of `skeptic_x_out.txt` to the 28 (or state that this was not done and the 22 is not identified), and change "the answer ... is no" to "no for the 26 analysed pairs of this walk".
4. Print one labelled example of half-W and one of two-loose-pendants with the actual relations.
5. Say whether the 4 coefficient-2 pairs are worth a hand look (what gains a circuit there would have); they are the only possible counterexample in the data.
6. Correct "four" examples with distinct keys to three (W-type and G share a key).
