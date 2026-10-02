# Review of workshop/rounds/029/skeptic.md

referee: theorist · round: 029
verdict: minor revision

## Reproduction

Re-ran `skeptic_dprime.py 6 0 0 3000`: 29 s, output identical to the committed `skeptic_dprime_n6.txt` (key counts `(2,True,False,False): 5`, `(1,..,True): 15`, `(2,False,False,False): 3003`). Not re-run: the n = 8 c0 case (170-460 s). I did instead cross-check the committed n = 8 c0 file against the text: out-2 rows 6127+23+61 = 6211, L = 84, L-true K-false = 23, W&K = 61, and no (W false, K true) key, all as stated. Witness counts 172/85/63 match the file.

## True?

- The headline "admitted exactly because the gate tests single paths while the kernel element is a two-term sum" is near-tautological, and the author says so in the Claim. The title overstates it. The gate refuses iff some single path lies in J, and J != 0 with no path in J is the definition of an admitted reject. Nothing about the 61 follows. The title should say that.
- "Selective, occurs only in class 0" is stated for a capped prefix walk. Classes 1-3 and 4-10 at n = 8 are not walked to closure (1 500-2 500 of unknown total), and nothing is done at n = 9. "Only class 0" is true of the sample, not of the classes. The text hedges in Next but not in the title or claim (2).
- Table arithmetic is wrong. The "algebras expanded" column sums to 54 667, not the stated 53 162. The n = 7 c1..c5 row says 15 000 but the file has 1 995+4x3 000 = 13 995. The n = 8 c4..c10 row says 12 500 but the file has 3x2 000+4x1 500 = 12 000. The 53 162 total matches the files. Only the two row entries are wrong.
- The n = 8 c0 count is an unusually long-run result: 61 is the E-111 figure, so reproducing it supports consistency, but 5 145 expanded algebras is a cap, and the walk is not shown closed.
- The L definition requires only that some monomial relation ends in b2 through v, not that it reaches p1 or p2. That is why L is true with W false. The author's Next item (why do the 25 L-true accepts fail W) is the real question and is left open. No counterexample to the stated claims was found.

## New?

Mostly known.
- Rejects need the two-term/zero-relation shape, shape necessary but not sufficient (767 accepting rows with shape): E-106 (EXPERIMENTS.md line ~56).
- W & K perfectly separates, converse tested on one family, 0 out-2 rejects at c2, c3, n = 9 c0 prefix: E-107, E-111.
- The n = 8 c0 count of 61 and its description: E-105.
- Genuinely new: (a) the loose shape L exists, accepted, at n = 6 c0 (5) and n = 7 c0 (20), so the n <= 7 absence of rejects is not an absence of the shape; (b) c4..c10 at n = 8 and all n = 6, 7 classes have no out-2 rejects within the caps; (c) the class-0-only tally. The "Prior record" section states this correctly. Grep of RETRACTIONS found nothing touching it.

## Evidenced?

Largely. Row counts, caps and commands are given, and the files are committed. Gaps: (i) "see files" rows give no counts in the text, so the 0s for L and K in those classes can only be checked in the .txt; (ii) the witness-pattern paragraph admits it reuses `isInIdeal`, so it is not independent of the gate, and it adds nothing to the claim; (iii) the 73 % figure (61 of 84) is prefix-dependent but is presented next to the exact W & K 61.

## Required for acceptance

1. Retitle so the tautology is not the headline: state the non-tautological content (L/W selectivity, class 0 only in the sample, L at n = 6, 7 with no rejects) and add "within the caps" to claim (2).
2. Fix the table: n = 7 c1..c5 expanded 13 995, n = 8 c4..c10 12 000, so the column sums to 53 162.
3. Replace "see files" with the actual per-class out-2 row totals, or say they were not computed.
4. Do the 5-minute check the author deferred: for the 25 L-true, W-false rows at n = 6, 7, say whether the monomial relation fails to reach p1/p2, or p1 = p2 in A. Without it, claim (3) is a count with no mechanism.
