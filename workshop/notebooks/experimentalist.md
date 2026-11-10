# Experimentalist notebook (rewritten round 123)

## What I now believe
- S-1/H-026 class level (round 123): the image key is NOT a class invariant at K = K0. n = 13 h = 2, 3 and n = 15 h = 3: the key class is ONE q_S class and 2-3 mutation orbits (lnaMoves rules + edgeMoves + doubleMutation, symmetric union-find over the whole class); the big orbit carries both image keys at every K <= K0, also on reduced rows. First counterexample n = 15 K = 3: 0003000000000 (i1) ~ 0003335000000 (i0), 26 moves (rules + double mutations at a repeated vertex). K > K0: one key, trivially constant.
- Rules-only symmetric orbits keep the image key at K = 2, 3 for n = 13 h = 3 and n = 15 h = 3 (E-307 extended), but NOT at n = 13 h = 2 (00300000000 ~ 00350000000, 13 rule moves via K = 0 rows). So constancy is a property of a move subset, not of the class.
- H-026 itself (key level, odd n, K0 = min(h, n-4-h)) stands; full-class K tables reproduce K0 = 3 at n = 15 h = 3. Even n: self-mirror class has no split (E-307). Two-3 classes: no law (115).
- T7 / H-010: no gap-2 pair lowered overlap (E-288/E-313, no power); k >= 5 hours per cell; liftPath ValueError on parallel arrows at k = 5, 6.
- Practical: full key-class scan at n = 15 = 2.67 M rows, ~110-165 s on 4 cores (numpy key 190 us/row); n = 13 = 7 s. Move neighbours 0.5 ms/row. The 200-shuffle control on ~7000 orbits costs ~7 min: cut to 50 shuffles next time.
- Lesson: "constant on components" depends entirely on which move families are in the union; always run the widest family (mut/der) before claiming an orbit invariant.

## What I tried
- 123: `rounds/123/experimentalist_classlevel.py` (plan / run N H [rules|mut|der]) + 7 *_out.txt.
- 119 gap2; 117 kpower (table/power/orbits); 115 kloc, nonlone; 114 tail; 113 kthresh; 110 refusal; 107 chain.

## What I would do next
1. Ask theorist for the right class-level statement (image key *set* per orbit? image derived class, E-236?) before more runs.
2. Why rules-orbits mix at n = 13 h = 2 but not n = 15 h = 3: run rules at n = 15 h = 2, n = 17 h = 2, 3 (n = 17 needs faster key / overnight).
3. Mutation-engine replay of the counterexample paths (toolsmith).
4. T7 overnight k = 5, 6 gap-2 pairs if approved.
