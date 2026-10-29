# Nikulin's "class = genus" step for q_S: Prop. 1.13.3 not found in 28 arXiv e-prints (UNVERIFIED); Thm 1.14.2 and Cor. 1.13.4 are the statements cited by other authors, and its hypotheses hold for the 140 pair; the shift slope b/a is constant on the 26 logged members of each of the 4 FCY-type groups (E-253, re-read)

## Response to referee

Verdict minor revision; all required items done.
1. Done. Title and Claim (2) now say "not found in 28 arXiv e-prints" and "statements cited by other authors"; "verified" removed everywhere.
2. Done. For the 140 pair `workshop/rounds/097/maverick_pair140.out` gives det S = -40 for both members, rank 12. Negative det forces an odd number of negative eigenvalues; diag(S) = 2 > 0 forces a positive one; so S is indefinite. "Modulo S being indefinite" dropped. Hypothesis check of 1.14.2 is complete, up to the theorem's own wording being unread (only secondary citations seen).
3. Done. `scholar_slope.py` now parses the logged "order/sign" field of each SUMMARY line and compares it with (-1)^(a+b) and a = order; prints PASS/FAIL. Re-run: PASS in 4 of 4 groups (15/-1, 20/+1, 11/-1, 21/-1).
4. Done. Next no longer proposes citing 1.14.2 in E-265; it stays UNVERIFIED with e-print citations as pointers only.
5. Done. Claim (3) states the slope is of nu^a P_l on logged members, not a functor-level statement.
6. Done by note: none of the 28 fetched e-prints quotes the statement of Thm 1.14.2 or Cor. 1.13.3 (they only cite them); fetching further e-prints was not done this sitting.



author: scholar · round: 098 · kind: negative (verification) + rediscovery
thread: T3/T8 · bears on: E-265, E-259, E-253, E-261
scope: sources = arXiv e-prints only (not the original paper; "not found" means not found in these 28) (Nikulin's 1979 Izv. paper is not on arXiv); 28 e-prints fetched and grepped, none quotes the text of 1.13.3. Slope part: recorded logs of round 091 only (26 members; 535 group at n = 11, 1185/2300/3992 at n = 12).

## Claim

(1) I did not find "Nikulin Prop. 1.13.3" in 28 arXiv e-prints. Nikulin, *Integral symmetric bilinear forms and some of their applications*, Izv. Akad. Nauk SSSR 43 (1979), is not on arXiv, and no fetched e-print reproduces the statement of 1.13.3 (or 1.13.2). So E-265's "class = genus (Nikulin 1.13.3)" stays UNVERIFIED as a theorem number. (2) What the e-prints give, as statements cited by other authors (secondary; the statements are not quoted): the numbers 1.13.4 and 1.14.2 are used for "unique in its genus" of even indefinite lattices under a rank-versus-length condition (below). (3) Slope: nu^a P_l = P_l[b] (a statement about the object nu^a P_l on logged members, not a functor-level statement) has the same (a, b), hence the same b/a, on every logged member of each of the 4 FCY-type groups; the 4 non-FCY groups are untested; this is E-253 / E-261 (4) re-read, not new. It does not say that b/a is constant on unsampled members, nor that it is a derived invariant of the Cartan data.

## Evidence

Fetched with the standing arXiv permission (STEERING.md, arXiv fetches), why: to find the theorem text of Nikulin 1.13.3 in a paper that restates it. `export.arxiv.org` API searches (K3 / Nikulin / discriminant form / genus), then `curl -sL https://arxiv.org/e-print/ID`, unpacked into their own scratchpad directories, grepped (never executed): 1109.2879, alg-geom/9503003, 1707.05677, math/0309348, math/0307355, math/0609233, 0810.2945, 1902.00229, 1403.6061, math/0312396, 0907.2691, 1010.3904, 1208.2570, 1301.1118, 1607.01410, 1707.02563, 1803.10550, 1906.08081, 2102.01193, 2103.07324, 2103.13382, 2202.08030, 2402.07611, 2412.06921, 2502.02143, 2510.16279, math/9808129, math/9905193. The arXiv index is thin; searches returned few hits.

Source hits for the "1.13 / 1.14" family (file:line in those e-prints):
- 2202.08030 main.tex l.863: "N is the unique lattice of signature (2,10) whose quadratic form is isomorphic to q_N, see [Nikulin, Cor. 1.13.4]". Nearest to a class = genus use; the statement of 1.13.4 itself is not quoted. l.845: Thm 1.10.1 (existence of an even lattice with given signature and q) is quoted in full (conditions (1)-(4) on t+, t-, l(A_p)); that is existence, not uniqueness.
- 1902.00229 davide1.tex l.101-113 (Lemma lem:surjectivity): for an even lattice S, the bound l_p(A_S) <= rank S - 10 for odd p and an l_2 condition give "S is unique in its genus and O(S) -> O(q_S) is surjective", "by Theorem 1.14.2 in [Nikulin1979]".
- math/0312396 nspap1.tex l.2700, 2729 and 2103.07324 main.tex l.353: Thm 1.14.2 used for uniqueness of an indefinite even T when rk T > l(A_q), or rk = l with a splitting condition on q_2 ("T unique by [Nikulin, Thm. 1.14.2]").
- math/0609233 l.1308, 1010.3904 l.1180: Thm 1.14.4 gives uniqueness of a primitive embedding for rk >= l + 2 (a different statement).
- 1109.2879 and 1403.6061 quote Thm 1.12.2 in full (primitive embedding / discriminant gluing), not 1.13.3.

My memory of the original (Cor. 1.13.3: an even indefinite lattice with rank >= 3 and rank >= l(A_p) + 2 for odd p, plus a 2-adic condition, is determined by (t+, t-, q); Thm 1.13.2 the general form) is NOT checked and I do not rely on it. It is plausible that E-265's "1.13.3" is right, but the citation by other authors for "S unique in its genus" is Thm 1.14.2 (as used in 1902.00229) with the hypotheses restated there.

Do the hypotheses hold for the open pair (n = 12, key group 140, |det S| = 40)? S is indefinite (it is: det S = -40 < 0 and diag S = 2 > 0, see Response), |A| = 40 = 2^3 * 5 forces l_2 <= 3 and l_5 = 1, so rank 12 >= l_p + 2 for every p, comfortably. So the cited conditions would be met, and "same genus => Z-congruent S" would follow from the cited theorem. This is the Z-congruence of S = C + C^T only; it says nothing about Z-congruence of C itself or P Phi_a = Phi_b P, which is what E-265 already says is open. Also "same genus" requires (signature, q_S) equal; Nikulin's genus-by-(signature, q) statement (Cor. 1.9.4 in my memory) is likewise unmatched.

Slope table (from `workshop/rounds/091/maverick_shift*.out`, via `scholar_slope.py`):

| group (n, size) | Phi order/sign | members logged | (a, b) | b/a in absolute value | a + b parity |
|---|---|---|---|---|---|
| 11, 535 | 15, -1 | 8 | (15, -16) all | 16/15 | odd (Phi^a = -I) |
| 12, 1185 | 20, +1 | 6 | (20, -22) all | 11/10 | even (+I) |
| 12, 2300 | 11, -1 | 6 | (11, -12) all | 12/11 | odd (-I) |
| 12, 3992 | 21, -1 | 6 | (21, -22) all | 22/21 | odd (-I) |

Each group has one (a, b) and sigma = id on all logged members; the parity column is now checked by the script against the logged sign (PASS 4 of 4). The slope differs between groups (4 values), so it is not a global constant. Sign: b is negative in the repo's convention; whether FCY convention gives positive b was left unchecked by E-261 and I did not check it.

## Reproduction

```
.venv/bin/python workshop/rounds/098/scholar_slope.py     # under 1 s
```
The literature part has no script: the fetch/grep commands are as in the Evidence list (curl of the e-print, tar/gunzip into a fresh directory, grep for "1\.1[0-9]\.[0-9]" and "genus").

## Prior record

E-265 scope line and Consequence already mark Nikulin 1.13.3 "cited from memory, UNVERIFIED" and list its verification as open (STATE.md T3/T8). This submission leaves it open for the theorem number, narrows it to the 1.14.2 and 1.13.4 cited by other authors, and shows the hypotheses are met for the 140 group (S indefinite: det S = -40 in maverick_pair140.out, diag S = 2 > 0). The slope part is E-253 (one shift per key, 26 of 26) and E-261 Result (4) (16/15 on 8 of 8); nothing new beyond tabulating the other three groups' slopes (11/10, 12/11, 22/21) and the parity check on all four. The E-253 entry text gives them as shifts only (-16, -22, -12, -22), not as slopes. Not in RETRACTIONS.

## Code changed

Added `workshop/rounds/098/scholar_slope.py` (log parser, no tests). Nothing else.

## Next

- Chair: E-265 keeps "Nikulin 1.13.3" UNVERIFIED; pointers only (not citations to adopt): 1.14.2 as used in arXiv:1902.00229 (davide1.tex l.101-113), Cor. 1.13.4 in 2202.08030 l.863.
- Anyone with library access: the printed Nikulin 1979 (Math. USSR Izv. 14 (1980), 103-167), §1.13, to match the numbers; or Conway-Sloane SPLAG ch. 15 for the indefinite class = genus theorem as an independent source.
- Maverick/theorist: the 140 pair needs the Phi-equivariant question (Z[Phi]-module isomorphism); S-congruence is not the obstruction.
