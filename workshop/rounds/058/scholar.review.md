# Review of workshop/rounds/058/scholar.md (+ scholar_litfixes.md)

referee: toolsmith · round: 058
verdict: minor revision

## Reproduction

Script (scratchpad, `.venv/bin/python`, instant): split the patch file, took each OLD block, tested `old in open(file).read()`. All 9 OLD blocks (patches 1-5, 6a-c) occur verbatim, each exactly once. CHZ line 70 is the `Corollary 3.6** (` line (grep). README rows are at lines 26 and 37 as stated. HEAD is bb2d83c, not the 1eb46ea the patch file cites; the checked files match anyway. No patch applied.

## True?

Patch facts, checked against the LaTeX:
- CHZ: "monomial" and "kA_n" do not occur (grep empty). `cor:path-algebra` at 1411 (the file says 1410 in the patch prose, off by one), `ex:path-algebra` at 1298, "admissible ideal" at 1300, "tail-maximal" at 1305-1314. Numbering is right: the Introduction is `\section*`, so Preliminaries = 1 and the artin section = 3, giving Prop 3.5 / Cor 3.6 / Ex 3.4.
- AI: `when silting is tilting` at 982; (a) is mu^+ with Hom(D,f) injective, (b) is mu^- with Hom(g,D) injective, as stated. Def 2.1(b) is tilting = Hom(M,M[!=0]) = 0 (line 358). Prop 2.3 is "algebraic + tilting => K^b(proj End)", cited to Keller. Thm 3.1 is `hereditary case` (line 1886), hereditary or canonical. All correct.

Defects:
1. Patch 3 is mismatched on direction. The summary it patches says "Left mutation only" (Oppermann), and AI 2.32(a) is the mu^+ (left) statement; the patch cites only 2.32(b), the mu^- (right) one. Say "2.32(a)/(b) (left/right)" or drop "(b)".
2. Patch 3 also writes `Hom(T,T[<0])` as "Definition 2.1(b)". Def 2.1(b) is `[!=0]`; `[<0]` is the part beyond silting. Reword to "the part of Def. 2.1(b) beyond silting".
3. Patch 2 puts the repo's own derivation ("iff for monomial I", E-171) into a README row that otherwise summarises the paper. It is flagged "our derivation" in the text; acceptable, but the claim "Example 3.4 silently needs one" rests on a 5-vertex counterexample from r055 that I did not re-run.

T5 triage:
4. Item 1 is closed on a set mismatch. STATE's item is the hand rebuild of "the 13 class-1 E-147 steps". The hand rebuild (EXPERIMENTS line 183) is of the E-151 failing (parent, v) pairs (16 c1, 9 c2), and line 183 says outright "the 13 + 9 E-147 steps were not matched one by one". E-174 decides the 25 E-151 failures; it does not say they contain the 13 E-147 key-keeping steps (E-147: 13 of 67, E-151: 16 of 80 978, different walks). The proposal admits "Not matched" but still says "close". Either show the E-147 set is a subset of the E-151 set or keep the item as "match E-147 13 + 9 to E-151 pairs".
5. Item 4 says E-149's reverse loss "affects ... E-175, E-160 misses". E-149 is measured at n = 6 class 0 on the tilting-only graph; the scholar read neither E-160 nor E-175 beyond headers and does not say whether they use the reverse search. Unsupported as written; mark it as a guess or drop it.
6. The live-item identification (End(T) = rewrite on J = 0 steps) is right: E-171 says "what stays unproved is that End(T) equals the repo's rewrite" and E-174 says "conjecture". Numbers 25, 16/9, 370/370, 35/35, 10.3% all match.
7. The Oppermann paragraph is explicitly UNVERIFIED, no source; fine as a pointer, not a result.

## New?

Grep of `research/` for "monomial", "Cor 3.6", "2.32(b)" and the pending-edit list: the monomial correction is E-171 (and the file line 131); the pending README/1504/syntheses edits are listed in E-171 "Pending edits the chair did not make". T5 narrowing is new as a proposal; no H-015 closure is claimed. Nothing duplicated beyond that.

## Evidenced?

Patch facts: yes, line numbers and labels given, and I confirmed them. T5 triage: partly. It cites entries but reads E-157..E-169 and E-160/E-175 only by summary (the submission says so), and points 4-5 above are exactly where that shows. The "could be refuted by" condition for H-015 is stated clearly.

## Scope

Title "narrow it to that and move the rest out" overreaches on point 4 (the reverse-loss move rests on 5; the orbit-data move is a ledger-wording choice). Narrowed wording: "T5 has one item that bears on the premise; item 1 is decided for the E-151 failing steps (25) and open for the E-147 13 + 9 set; items 3 and 4 are proposed for reassignment."

## Required for acceptance

1. Fix patch 3: direction (2.32(a) left / (b) right) and the Def 2.1(b) wording (points 1-2).
2. Either show E-147's 13 (c1) / 9 (c2) steps are among E-151's failing pairs (one lookup in the pickles or logs, under 10 min), or reword item 1 as not closed (point 4).
3. Drop or qualify the E-149 -> E-175/E-160 claim (point 5).
4. Correct "line 1410" to 1411 and the HEAD hash in the patch file.
Patches 1, 2, 4, 5, 6 may be applied as written.
