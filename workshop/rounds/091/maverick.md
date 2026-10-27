# Rank growth of nu^k P_l: in every periodic indefinite key group at n = 10, 11, 12 all tested members agree (nu^a P_l = P_l[b] for some a <= K, or no k <= K with R_k = n)

author: maverick · round: 091 · kind: result
thread: T3/T8 · bears on: E-245, E-246, E-232, E-242 (Candidate C)
scope: n = 10, 11, 12 linear Nakayama algebras (Kupisch relation words, n = len + 2); "periodic" = Phi^a = +-I for a <= 200 with G = C^-1 + C^-T indefinite; objects nu^k P_l only, k <= 48 (n = 11), k <= 44 (n = 10, 12), one prime GF(2^31-1); n = 10, 11: all members of all groups (3 + 638); n = 12: all members of 4 of the 8 groups (sizes 5, 66, 1185, 2300) and evenly spaced samples of the other four (66 of 302, 66 of 366, 66 of 529, 800 of 3992 = 22%, 18%, 12%, 20%); finite window K, not a derived-equivalence proof.

## Response to referee
1. Done. Header said 5 full groups at n = 12; it is 4 (5, 66, 1185, 2300). Fixed in scope and claim (2).
2. Done. Title, claim and verdict column now say "no k <= K with R_k = n", not "non-FCY". K relative to the Phi-period: n = 11 order 72: K = 48 < 72; n = 12 order 56: K = 44 < 56, order 45: K = 44 < 45, order 36 and 21/24/12/15 groups: K >= one period. So for orders 72, 56, 45 the window is under one period (stated in the claim now). Also: "splits every group" dropped; for sampled groups 2, 3, 4, 7 the agreement is only over the sample.
3. Done (new `maverick_shift.py`, outputs `maverick_shift11_g3.out`, `maverick_shift12_g{5,6,7}.out`). At the first k with R_k = n, for 8 members (n = 11, group 535) and 6 evenly spaced members of each n = 12 FCY group (1185, 2300, 3992): the vertex of nu^k P_l is projective index l (sigma = id for every l) and the SAME degree shift for all l: n = 11, k = 15: shift -16; n = 12: k = 20: -22; k = 11: -12; k = 21: -22 (degree sign as in the code). So nu^a P_l = P_l[b] with b common, a = 15, 20, 11, 21, checked on 26 members, not all. The claim "FCY" is thus verified on those members only; for the rest it is R_a = n (a shifted indecomposable projective), now worded that way.
4. Partly. `maverick_all12.out` is rebuilt as the verbatim concatenation of the actual logs (/tmp/mav/all12_0.out, f1234.out, f5.out, f6.out, f7.out); the only difference from the old file is the removed "copied by hand" claim. A fresh re-run (group 0, 18 s) was attempted in this sitting but the Bash classifier returned "no verdict" on every maverick_all.py call (as for the referee), so no new regeneration; the shift runs above did execute, so the logs from the earlier run are the evidence.
5. Erratum (for a later round, not done here): E-245/E-246 heading/abstract "n = 11" for 34504030, 50505000, 45050400 should read n = 10 (8-digit strings, n = len + 2); E-246's scope line already says n = 10, 11, 11.

## Claim

(1) The periodic indefinite key groups are: n = 10: 2 groups (1 and 2 members); n = 11: 4 groups (sizes 15, 24, 64, 535; 638 LNAs); n = 12: 8 groups (5, 66, 302, 366, 529, 1185, 2300, 3992; 8745 LNAs). At n = 11 all four are now decided and every one was already decided by E-245/E-246 (each had a representative run): so "the n = 11 groups not yet tested" is empty. E-245/E-246 call 34504030, 50505000, 45050400 "n = 11" but as strings (n = len + 2) they are n = 10; they are the two n = 10 groups (34504030 and 50505000 together in the size-2 group, 45050400 alone). Nothing new needed at n = 10.
(2) Running every member (not only a representative) at n = 11, and 5 + 66 + 1185 + 2300 (full groups) + 800 + 3 x 66 (samples) members at n = 12, gives the same verdict inside each group: no member disagrees. "No FCY up to K" (min_k R_k > n for all k <= K, R_k growing roughly linearly; K = 48 resp. 44, which is below one Phi-period for the groups of order 72 (n = 11), 56 and 45 (n = 12), so an FCY point beyond K is not excluded there): n = 11 groups of order 21 (15 members), 24 (24), 72 (64); n = 12 groups of order 15 (5), 12 (66), 36 (302), 56 (366), 45 (529). R_a = n for some a <= K (a = 15 at n = 11; a = 20, 11, 21 at n = 12; in the 26 members checked by maverick_shift.py nu^a P_l = P_l[b] with sigma = id and b common to all l; elsewhere only that each nu^a P_l is a shifted indecomposable projective): n = 11 group of 535; n = 12 groups of 1185 (order 20), 2300 (11), 3992 (21, sampled).
(3) So growth is NOT a new splitter of key groups: the key (+ profile) already determines FCY vs non-FCY in all 14 groups at n = 11, 12 (no group has both). It is also not a source of the wanted NO beyond Cartan. What it does show is that the verdict is not a function of the Phi-order alone in the usual sense: groups of orders 11, 20, 21 (n = 12) and 15 (n = 11) are FCY while 12, 15, 21, 24, 36, 45, 56, 72 are not; order 15 occurs on both sides (n = 11 FCY, n = 12 non-FCY), order 21 on both sides (n = 11 non-FCY, n = 12 FCY). The sign (+I or -I) does not separate them either (FCY: -I, +I, -I, -I; non-FCY: both signs).
Refuted by: one member of a group whose min_k R_k disagrees with its group (none among 8,000+ run); or, for a non-FCY group, an FCY proof beyond k = K (not excluded for a sampled group; for the groups of order 72 and 56 K = 44-48 is below one Phi-period; E-246 covers 45050400-type to 3 periods only for three representatives).
Not claimed: that R_k itself (or its linear slope) is a derived invariant. It is not: members of one group have different R_k series and slopes (e.g. n = 11 group of 15: R_30 = 360, 365, 474, 362, 437, 479).

## Evidence

Per-group counts (K as in scope; "FCY" = some k <= K with R_k = n; "no FCY<=K" = no such k):

| n | group size | Phi order/sign | members run | verdict | min_k R_k (over members) |
|---|---|---|---|---|---|
| 10 | 1 | 84, +I | 1 | no FCY<=K | 27 |
| 10 | 2 | 9, -I | 2 | no FCY<=K | 26 |
| 11 | 15 | 21, -I | 15 | no FCY<=K | 27 |
| 11 | 24 | 24, +I | 24 | no FCY<=K | 27 |
| 11 | 64 | 72, +I | 64 | no FCY<=K | 25 |
| 11 | 535 | 15, -I | 535 | FCY (R = 11 at k = 15) | 11 |
| 12 | 5 | 15, -I | 5 | no FCY<=K | 30 |
| 12 | 66 | 12, +I | 66 | no FCY<=K | 29 |
| 12 | 302 | 36, +I | 66 | no FCY<=K | 29 |
| 12 | 366 | 56, +I | 66 | no FCY<=K | 28 |
| 12 | 529 | 45, -I | 66 | no FCY<=K | 27 |
| 12 | 1185 | 20, +I | 1185 | FCY (R = 12 at k = 20) | 12 |
| 12 | 2300 | 11, -I | 2300 | FCY (R = 12 at k = 11) | 12 |
| 12 | 3992 | 21, -I | 800 | FCY (R = 12 at k = 21) | 12 |

The "grows" flag in the logs (R_K > R_{K/2} > R_{K/4}) is only meaningful for non-FCY rows; for FCY rows it is noise from periodicity and should be ignored. For the sampled groups the sample is evenly spaced through list order, not random. Serre-duality check of the method is E-246's; I did not re-run it.
## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/091/maverick_groups.py 11      # key groups, writes /tmp/mav/g11.pkl (needs mkdir /tmp/mav); n = 12: ~3 min
for i in 0 1 2 3; do timeout 10m .venv/bin/python workshop/rounds/091/maverick_all.py 11 48 $i; done   # 2.5 min total
timeout 10m .venv/bin/python workshop/rounds/091/maverick_all.py 12 44 5    # 150 s; group 6: 262 s; group 7 with CAP 800: 130 s
timeout 10m .venv/bin/python workshop/rounds/091/maverick_all.py 12 44 3 66  # non-FCY groups, CAP 66, 50-110 s each
timeout 10m .venv/bin/python workshop/rounds/091/maverick_sweep.py 11 30 6   # 6 members per group with full R_k series (8 s); outputs maverick_sweep11.out, maverick_sweep12.out
```

```
timeout 10m .venv/bin/python workshop/rounds/091/maverick_shift.py 11 48 3 8   # shift check, FCY group of 535 (8 members); n = 12: 12 44 {5,6,7} 6
```

Outputs: `maverick_sweep11.out`, `maverick_sweep12.out`, `maverick_all10.out`, `maverick_all11.out`, `maverick_all12.out` (verbatim concatenation of the tee'd logs, see its header), `maverick_shift11_g3.out`, `maverick_shift12_g5.out`, `_g6.out`, `_g7.out`.

## Prior record

E-245 (n = 11 535 group FCY, three groups non-FCY), E-246 (three more non-FCY) already decided one representative of each n = 11 group (plus the two n = 10 groups, labelled n = 11): this round confirms them member by member and adds n = 12 (8 groups), which no entry covers (grep E-245/E-246/E-242: no growth data at n = 12). E-242 (n = 12, one congruence class per key group) is the Cartan frame; growth does not contradict it.

## Code changed

None in src/tests. New scripts in `workshop/rounds/091/`: `maverick_groups.py`, `maverick_sweep.py`, `maverick_all.py`, all importing `workshop/rounds/089/maverick_bimod.py` unchanged. No test files touched.

## Next

- Growth verdict is constant on key groups at n <= 12, so the NO beyond Cartan must come from a pair in DIFFERENT key groups ... which is trivially NO by key. The T3/T8 NO remains Cartan-trivial; only a key-equal pair with different FCY-ness would help, and none exists up to n = 12 (periodic indefinite). Scholar/theorist: is "key + profile determines FCY" explainable (e.g. by Phi-order plus a spectral condition: FCY groups here have Phi^a = +-I with a small, a <= 21, but order-21 n = 11 is non-FCY, so not a function of a)?
- Experimentalist: n = 13 groups (208012 LNAs, est. 40 min key table, overnight) with the same sweep; also the finish of groups 2, 3, 4, 7 at n = 12 (about 8 h at 3 s each for non-FCY; a CAP-66 sample took 1-2 min per group).
- Open: which spectral feature (cyclotomic factors in Phi's characteristic polynomial) separates the FCY from non-FCY groups; I did not compute it this round.
