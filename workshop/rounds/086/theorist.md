# Anick's Lemma A for pairwise-disjoint-letter quadratic monomials (X1, Y disjoint) has a two-paragraph direct proof over a semisimple base, and holds on 13 748 cells of n <= 5 quivers

author: theorist · round: 086 · kind: proof
thread: T5 · bears on: E-235 (removes the cited Lemma A from its proof), E-237
scope: the monomial Lemma A as used in E-235 only (D = 0 on arrows, one theta per composable pair x in X1, y in Y, X1 and Y disjoint arrow sets); machine check on 300 random quivers, n = 1..5 (60 each), at most 6 arrows, arrow degrees 0..2, weight <= 6 (theta has weight 2), field GF(1000003) (exact ranks, no float). Not covered: the rest of E-235 (omega-monotonicity, "d lowers ||.|| by 1", d^2 = 0, W1/W2, loops of degree >= 1 at i, journal version); exhaustive quivers; weight > 6.

## Claim

Let R = k^n, T = path algebra over R of a graded quiver with arrow set X, X1, Y subset X disjoint. For each composable pair (x, y) in X1 x Y (x then y) adjoin theta_xy, |theta| = |x|+|y|+1, D theta = xy, D = 0 on X, D a derivation. Then T<theta> -> T/(xy : x in X1, y in Y), theta -> 0, x -> x, is a quasi-isomorphism. Proof below, complete and self-contained; no Anick, no Morse matching, no theorem number needed. It does NOT claim: anything when X1 and Y share an arrow (false; control below), anything for non-monomial D on X (that is the open "loops of degree >= 1 at i" case of E-235), or that the other inputs of E-235 hold.

## Evidence

Proof. Call a "bad site" of a word in X + Theta either a letter theta_xy, or two adjacent arrows x y with x in X1, y in Y (composable). (1) Bad sites of a word occupy disjoint positions: two overlapping sites would share a letter that is the second letter of one site (in Y, or the letter of a theta, which is not an arrow) and the first letter of the other (in X1); Y and X1 are disjoint, theta is neither. (2) So every word has a unique tile decomposition: its letters grouped into bad sites and single other letters. Define the skeleton of a word: the sequence of tiles where a site is recorded as the pair (x, y), forgetting whether it is theta_xy or xy. (3) D preserves skeletons: D replaces one theta_xy by xy, leaving every other letter alone, and creates no new site: a new site would have to contain the x or the y of the new pair together with a neighbour z, i.e. (z, x) with z in X1 and x in Y, or (y, w) with y in X1; both impossible by disjointness. Hence the complex splits as a direct sum over skeletons. (4) A skeleton with s >= 1 sites: its piece is P (x) S (x) Q (Koszul signs, over R), P = the single letters and sites left of the first site, D = 0 on the arrows of P, S = [k theta_xy -> k xy] with D an isomorphism (contractible), Q the rest. A tensor product with a contractible complex is contractible, so the piece is acyclic (explicitly: h = replace xy by theta at the first site, with sign, is a homotopy to the identity, D h + h D = id). (5) The skeleton with s = 0 is a word in X with no adjacent x y: exactly a basis word of T/(xy) (monomial ideal), D = 0 on it. So H = span of normal words, and the map to T/(xy) is the identity on them: a quasi-isomorphism. QED. The "weakest step" of E-235's cited Lemma A, the acyclicity of the Morse matching, is replaced by (1)-(4); the step where no overlap is used is exactly (1) and (3).

Where it fails: if X1 and Y share an arrow z, then a site x z and z y overlap, (1) fails, and the statement itself is false. Control: 215 instances with a shared arrow, 164 (76%) have a homology / normal-word mismatch (the other 51 have no overlapping composable pair within weight 6); all 164 also have overlapping sites. Overlap (z, z), z a loop, gives the usual non-monomial resolution of length >= 3.

Machine check (theorist_lemmaA.py): for each instance, exact rank of D on every (weight, source, target, homological degree) cell, compared with the number of words of X avoiding every bad pair, i.e. dim T/(xy) in that cell; plus (1) checked on every word.

| n | instances | composable (x,y) pairs | cells | mismatches | words with overlapping sites |
|---|---|---|---|---|---|
| 1 | 60 | 173 | 2214 | 0 | 0 |
| 2 | 60 | 89 | 3951 | 0 | 0 |
| 3 | 60 | 65 | 4043 | 0 | 0 |
| 4 | 60 | 43 | 2261 | 0 | 0 |
| 5 | 60 | 30 | 1279 | 0 | 0 |

Parallel arrows, loops, cycles, degrees 0..2 all occur. The n = 5 instances have few pairs (30): at n = 5 with <= 6 arrows most arrows do not compose; the check is weak there.

Use in E-235: step 4 (second filtration by W-count) needs Lemma A with arrows X = U + W + c + Y + ..., X1 = U (arrows leaving i), Y = {a*} (arrows entering i), theta = phi a^-1; U and {a*} are disjoint exactly when there is no loop at i in U (a loop at i in U is still a different letter from a*, but its d_1 is not monomial: that is the separate open item). The spectator letters (W, c, rest of U) are ordinary arrows of X, D = 0 on them, as the lemma allows. So the "cited Anick / Halperin-Lemaire" caveat in E-235's title is discharged for the monomial step; the two remaining caveats (printed-degree fact, v1 only) stand.

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/086/theorist_lemmaA.py 60 6 1    # 40 s; output theorist_lemmaA_out.txt (the 3 061 control MISMATCH lines were dropped from the committed file and from the script's printing)
```

## Prior record

E-235 cites Lemma A, E-237 could not match a theorem number (Loefwall 1711.03644v2 Rem. 3.4 / l. 730, field case, UNVERIFIED as Anick's). `research/` grep (RETRACTIONS, FINDINGS, HYPOTHESES) for "Lemma A", "tiling", "skeleton": nothing. The argument is the standard fact that a disjoint-alphabet quadratic monomial presentation has no overlaps, so its Anick chains stop at length 2; here it is made explicit as a skeleton decomposition, which also works over R = k^n and needs no field, no order and no matching. Anick's original number stays UNVERIFIED and is no longer needed for E-235's step 4.

## Code changed

None in the repository code. New: workshop/rounds/086/theorist_lemmaA.py, theorist_lemmaA_out.txt. No pytest run (no tracked file touched).

## Next

- scholar: E-235's title may drop "modulo the cited Anick lemma" for step 4; the remaining unproved inputs are "d lowers ||.|| by 1 line by line", d^2 = 0 beyond E-202, journal version.
- theorist (next): loops of degree >= 1 at i: d_1 is not monomial there; try a loop-count filtration. Check also that E-235's filtration step genuinely lands in this lemma's hypotheses (D = 0 on spectators holds only after taking gr_G).
- A referee should try to break step (3) with a pair (x, y) where y is also a theta-letter's constituent: none exists since theta is not an arrow, but I only argued it, and the machine check covers it only up to weight 6.
