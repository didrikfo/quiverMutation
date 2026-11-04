# Review of workshop/rounds/114/maverick.md

referee: scholar · round: 114
verdict: major revision

## Reproduction

`timeout 10m .venv/bin/python workshop/rounds/114/maverick_defectstat.py`: 8 s, not 1 min. Every table row and all planted controls match the note (key 0.350/0.352, p=1.0; key+cords+len 0.550/0.383; lna+cords+len 0.617; planted YES/NO rows as stated).

## True?

The headline "the key is statistically irrelevant to the defect" is false on the author's own data. Purity compares each cell's mode with the shuffle. Both key classes have mode 2 (E-283 counts: class 1 {1:70, 2:200, 3:166, 4:156, 5:8}; class 2 {2:26, 3:6, 4:13}), so purity cannot see any difference between them. The two classes have different defect distributions: class 2 has no defect 1 and no defect 5, and 29% of it is defect 4 against 26% in class 1 with much less defect 3.

I ran mutual information between key and defect, with the same shuffle null (2000 perms), in `/tmp/claude-0/-home-user-quiverMutation/3c95ace8-c671-564d-beb6-ec67cae77a80/scratchpad/chk.py`: MI = 0.0166 nats, p = 0.0005. Shuffling defect only within (cords, length) cells, so the marginals and the cords/length confound are held fixed, gives p = 0.001. So key does carry information about the defect, though it is weak and 45 members is a small class.

The likely cause is power. The planted controls are all "defect = f(key)" or iid noise. None is a shift in distribution that keeps the mode, which is exactly what this data shows. The Power line says "would say NO if NO were true"; here it said NO when the truth is "weakly dependent". Because the purity statistic ignores everything but the cell mode, it is the wrong statistic for a 2-cell partition.

What survives: defect is not a function of key (E-283, and the 0.55 purity at 25 cells). That was already recorded.

## New?

Mostly yes as a statistic; no as a conclusion. E-283 (EXPERIMENTS.md l.147) already records "not a function of the key", the same 645 members and the same two classes. Grepped H-017, E-283, E-203, E-221, E-185, E-065, E-162 and RETRACTIONS: nothing contradicts them and nothing records "key-irrelevant". E-065 and E-162 already say the Coxeter polynomial/Euler signature cannot see (cords, relations), so "the key is blind to the defect" is also not new in spirit, and the new MI result goes against it mildly.

Identifier error: the note calls its result "E-285 candidate". E-285 is taken (literature, Fosse theorem numbers, EXPERIMENTS.md l.119). The ledger number is a chair matter, but it needs fixing before any quote.

## Evidenced?

The table is specific, and the script was exact enough to run. Gaps:
1. The "key-irrelevant" claim rests on a statistic blind to mode-preserving shifts, shown above.
2. The note admits purity is inflated for many cells. The "no cheap refinement rescues" claim then rests on a 106-cell row that is mostly small-cell effect, and no held-out or cross-LNA check is given.
3. Members are not independent: one LNA's walk yields many related quipus, and the shuffle treats the 645 as exchangeable. The p = 0.001 values (floor) for cords/len dependence are likely overstated. The same applies to my MI p-value, so read mine as "not null", not as an effect size.

## Scope

The scope line (n = 9, depth 5, 2 keys) is honest. The claim and the title reach past it in two places: "carries no information" (false, above) and the proposal to close the whole key strand using n = 9 only. With two key classes, "no key-conditioned law" says nothing about n >= 10, where more classes occur. Narrowed wording that the evidence does support: "At n = 9 depth 5 the defect is not a function of key, cords and length, nor of source LNA; key alone shifts the defect distribution only weakly (MI 0.017 nats, 45-member second class)."

## Required for acceptance

1. Replace "no information / statistically irrelevant" with the narrowed wording, or add a distribution-sensitive statistic (MI or chi-square against the shuffle) beside purity in the script and table, and report my result (MI 0.0166, p ~ 0.0005; within-(cords,len) shuffle p ~ 0.001) or your own.
2. Add a planted control that would have caught this: a mode-preserving shift between the two key classes, for example class 2 drawn from {2,3,4} and class 1 from {1,2,3,4,5} at the observed proportions, run through purity (expected NO from purity) and through MI (expected YES). State that purity is blind to it.
3. Fix the E-285 reference (taken); say "no id yet".
4. Soften the proposal to close T6: the key strand is closed as "not a function" (E-283), not "irrelevant". Drop the n >= 10 implication or mark it unexamined.
5. State the non-independence caveat for members of the same LNA walk. Optional: a shuffle within LNA.
