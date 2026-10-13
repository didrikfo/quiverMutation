# At n = 10 (one group), the Phi^18 group has two derived classes, not three: 90000000 = kA_10/rad^9 is of type D_10 and joins class 0 (a correction of E-178's "3 candidates" via F-045 and Happel-Seidel); Candidate C is moot at n = 10 but NOT at n = 11

author: maverick · round: 063 · kind: result
thread: T3/T8 · bears on: E-172, E-178, F-045, F-047, F-048, Candidate C
scope: n = 10, the one E-172 group (the "moot" claim is n = 10 only; n = 11 was scanned at the Cartan level and C is NOT moot there, see the response, item 3) (Coxeter polynomial (T+1)^2 (T^2-T+1)(T^6-T^3+1)); Cartan-level exact arithmetic (n = 10) and a Coxeter-key/signature census (n = 11), one representative per orbit-plus-mirror class (00000030, 34504030, 50505000, 90000000). No object-level computation of S^a on any complex. Level: tested on small cases (exact matrix identities; the D_10 part rests on cited theorems).

## Response to referee

Verdict was minor revision; all three required items done.

1. Retitle/scope: done. Title now carries "n = 10 (one group)" and says this is a correction of E-178/STATE's "3 candidates" by F-045 plus Happel-Seidel, not a new equivalence. Credit: FINDINGS l. 595 already counts the tame D_10 class as 321 = E-178's 320-row class 0 + 90000000, and F-045 already records "n = 10: 2 orbits into 1" (this merge). What I add is only: the wording correction for E-178/STATE, the Cartan-level table (C is blind on the group), the (a, b) = (9, even) statement, and 90000000 ~ 00000030 as an unjoined equivalent control. The narrowed claim, as the referee phrases it: at n = 10 the Phi^18 group has two derived classes (D_10 incl. 90000000; and W = {34504030, 50505000}), already separated by Euler-form signature; Cartan-level C adds nothing. The "Scope" line now says the same.
2. Runtime: corrected to about 2 s (referee measured 1.8 s).
3. n = 11 check: DONE, new script `maverick_n11periodic.py` (37 s, `timeout 10m` first; n = 10 control run reproduces F-048's 3 rows in 9 s). Result at n = 11 over all 16796 LNAs: 638 rows with periodic Phi and indefinite Euler form (= F-048's 638 firing; 619 new + 19 already certified). They fall into exactly FOUR Coxeter-polynomial groups, of sizes 535 (ord Phi 30, Phi^15 = -I), 64 (ord 72, no Phi^a = -I), 24 (ord 24, none), 15 (ord 42, Phi^21 = -I). ALL 638 have the same Euler signature (2 neg, 0 zero, 9 pos), and every LNA sharing a key with a row is itself a row (0 exceptions). So the signature separates nothing inside any of the four groups, and the Cartan-level data of each group is shared by all its members: Candidate C is NOT moot at n = 11 (it is live in every group), the opposite of n = 10 where the only group with two same-signature members had a D_10 class on the other side. I did not compute derived classes inside these groups (that needs the search), so I do not know whether any group has two classes; a group could be a single derived class. Consequence: I withdraw "close the Candidate C part of T3/T8"; the recommendation is now "closed at n = 10", open at n = 11 with four groups as the target list. The earlier Next-bullet question (a group with periodic Phi, >= 2 candidates, all of equal signature, not D_n type) is thereby answered "yes, candidates at n = 11", pending a class count.
4. Smaller points: "(b even)" is conditional (what any shift relation must satisfy at K_0; realised only on the D_10 class), now in the Claim as written; "refuted if ... Lambda(10,9) not D_10" is settled by two independent routes and the script's G-check is the F-045 computation re-run, not a new check (stated in Evidence). The orbit-plus-mirror classification and the singleton status of 90000000 are taken from E-172/E-178, not tested by my script.

## Claim

(1) Cartan level. For all four class representatives Phi has order exactly 18, Phi^a = -I only for a = 9 (and Phi^a = I only for a = 18). With [S] = -Phi on K_0 (so that kA_n gives S^(n+1) = [n-1]), the smallest consistent shift relation is S^9 ~ [b] with b EVEN (a = 9: [S]^9 = I); a = 18 allows b even too. The Cartan data cannot say more: b is not seen by K_0. This is the same for all four classes, so the Cartan level has no separating power here (as expected, this is the content of the Coxeter-polynomial key).
(2) The group's spectrum is exactly that of kD_10 (h = 18, Phi^9 = -I, matrix-checked). For D_10: S = tau[1], tau^(-9) = [1] (n even), so S^9 = [8]: (a, b) = (9, 8) is realised, and fits (1) (b even).
(3) 90000000 is kA_10/rad^9 = Lambda(10, 9), piecewise hereditary of type D_10 (Happel-Seidel, `research/literature/happel-seidel-piecewise-hereditary-nakayama.md` l. 64). Independently, F-045: its G = C^-1 + C^-T is positive definite with det 4 and SNF (1^8, 2, 2), the D_10 data, identical to 00000030 (class 0). So 90000000 is derived equivalent to class 0: the group has TWO derived classes, D_10 = {class 0 (320 rows), 90000000} and W = {34504030, 50505000}. W has G with 2 negative eigenvalues (not piecewise hereditary of Dynkin type, and not derived equivalent to D_10).
(4) Hence a shift relation exists exactly on the D_10 class (theorem, S^9 = [8]), and Candidate C cannot separate anything the signature has not: the only inequivalent pair is D_10 vs W, already certified by the Euler form signature (F-045/F-048 style: positive definite vs signature (8,2)). Whether W is fractionally CY is a genuine but separation-irrelevant question; I did not compute it. Not claimed: W has no shift relation; S^9 = [8] checked at object level for any LNA (it is quoted for kD_10).
Refuted if: a Z-congruence class of G or an F-045 classification puts 90000000 outside the D_10 class; or Lambda(10,9) is not of type D_10.

## Evidence

| rep | ord Phi | Phi^a = -I | [S]^a = I | C+C^T (neg/zero/pos) | G=C^-1+C^-T (neg) | det G | SNF G |
|---|---|---|---|---|---|---|---|
| 00000030 (class 0) | 18 | a = 9 | a = 9, 18 | 0/0/10 | 0 | 4 | 1^8,2,2 |
| 34504030 | 18 | a = 9 | 9, 18 | 2/0/8 | 2 | 4 | 1^8,2,2 |
| 50505000 | 18 | a = 9 | 9, 18 | 2/0/8 | 2 | 4 | 1^8,2,2 |
| 90000000 | 18 | a = 9 | 9, 18 | 0/0/10 | 0 | 4 | 1^8,2,2 |
| kD_10 (reference Cartan) | Phi^9 = -I, Phi^18 = I, char poly (l+1)^2 (l^2-l+1)(l^6-l^3+1) | | | | | | |

Reading: the four Cartan-level periodicity data are identical (so C at K_0 is blind, as round 055 predicted); the signature splits {00000030, 90000000} from {34504030, 50505000}. Round 058 counted 3 candidates / 2 certified inequivalences using the F-047 profile, which has 0 and 3 sharing a profile; the missing step was to read that as a derived equivalence (Happel-Seidel / F-045), which makes it 2 classes / 1 inequivalence. This also matches FINDINGS l. 595 ("321-member tame D_10 class": 320 + 1) and the "2 orbits into 1 at n = 10" merge of F-045.
Consequence for E-178: the pair 90000000 vs 50505000 is the pair D_10 vs W; its empty join to 6 + 6 is a non-join of an inequivalent pair, as stated, but 90000000 vs 00000030 is an EQUIVALENT pair that the F-047 moves (orbit plus mirror) do not join (singleton class 3). That pair is a long-path equivalence control for the n = 10 search that the workshop has not used.

## Reproduction

```
.venv/bin/python workshop/rounds/063/maverick_cartanC.py     # about 2 s
timeout 10m .venv/bin/python workshop/rounds/063/maverick_n11periodic.py   # n = 11 census, about 40 s
```

## Prior record

Known: S^a ~ [b] forces Phi periodic (0911.5137 Cor 1.9; 1310.1557 Thm 2a): used, not new. Known: the 321-member tame D_10 class and the n = 10 merge (F-045, FINDINGS l. 595); Happel-Seidel Lambda(n, n-1) of type D_n (literature note). E-172/E-178/STATE.md treat 90000000 as a third candidate class: not contradicted by data, but superseded by F-045 plus Happel-Seidel. The conversion "Phi^a = (-1)^b" in the literature note depends on the sign convention of Phi; with Phi = -C^-T C as in this repo, [S] = -Phi and the relation is Phi^a = (-1)^(a+b) (checked against kA_n and kD_10: a = 9, b = 8). New: the Cartan-level table, the (a, b) = (9, even) statement, and the 2-class reading of the group.

## Code changed

None (new scripts `workshop/rounds/063/maverick_cartanC.py` and `workshop/rounds/063/maverick_n11periodic.py` only).

## Next

- Recommendation (revised): Candidate C is closed at n = 10 only. At n = 11 it is live: the 638 periodic-indefinite LNAs form four Coxeter-key groups (535, 64, 24, 15), all of Euler signature (2,0,9). T3/T8 as a catalogue thread stays "done at catalogue level" (E-066), with the unproved "why no shift by 1" as the standing open item.
- Open question: do the groups of 535/64/24/15 contain more than one derived class? If a group is one class, C is moot there too; if two, signature cannot separate them and C (or a search) is the tool. Groups of ord Phi 30 and 42 admit S^a ~ [b] (Phi^15, Phi^21 = -I); those of ord 72 and 24 have no Phi^a = -I.
- Experimentalist: use 90000000 ~ 00000030 (equivalent by Happel-Seidel, not joined by orbit+mirror) as the long-path equivalence control for the J0 join search; correct the "3 candidates" wording in E-178 and STATE.
