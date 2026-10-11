# After E-171 the key-guard thread (T5, H-015) has one open item that bears on the premise, a proof (End(T) = rewrite on J = 0 steps); item 1 is decided for the 25 E-151 failing steps, open for E-147's missing c1 steps; do not close H-015

author: scholar · round: 058 · kind: proposal
thread: T5 · bears on: H-015, E-171, E-174
scope: literature and ledger reading, plus one rerun (round 058 response: E-147 and E-151 walks rebuilt in scratchpad copies, see Reproduction); no arXiv fetch. Read: STATE.md T5/T10, HYPOTHESES H-015, EXPERIMENTS E-147/E-149/E-171/E-174 and the entry at line 183 (round 046 skeptic rebuild). Not read: E-160, E-175 beyond headers, the round 056 and 057 submissions, E-157..E-169 beyond their STATE summaries. T5 has one item that bears on the premise; item 1 is decided for the E-151 failing steps (25) and open for E-147's c1 steps the rerun did not reproduce; items 3 and 4 are proposed for reassignment.


## Response to referee

Verdict minor revision; the four required items:

1. Patch 3 (done, in `scholar_litfixes.md`): now cites 2.32(a) (left mutation mu^+, injectivity of Hom(D,f)) and 2.32(b) (right, Hom(g,D)); `Hom(T,T[<0])` is called "the part of Def. 2.1(b) beyond silting" (2.1(b) itself is `Hom(M,M[!=0]) = 0`). Accepted, points 1-2.
2. T5 item 1 (done, a lookup/rerun, accepted point 4). Scripts: `rounds/043/skeptic_c2.py 7 {1,2} .. 500 30 guard` (copied to scratchpad, changed only to keep every J != 0 step as (canonicalKey(parent), v) and every D = 0 Cartan class) and `rounds/046/skeptic_collect.py 7 {1,2} 20000` (copy storing canonicalKey(parent)). Result: c2: E-147's 9 D = 0 steps = E-151's 9 failing (parent, v) pairs, set equal (9 of 9). c1: the rerun found 12 distinct (parent, v) pairs with D = 0 (8 Cartan-distinct) and all 12 are among E-151's 16 (12 of 12; the other 4 are E-151 only). The c1 count is not E-147's 13: the rerun ran 4 jobs on 4 cores, so the 500 s budget gave 18 017 expansions and 56 distinct J != 0 steps (E-147: 67; c2 reproduced at 59 vs 64), and D = 0 at 8 Cartan-distinct (E-147: 13). So: containment is shown for the rerun's steps (c1 12/12, c2 9/9), not for E-147's own 13; the 5 Cartan-distinct c1 steps the rerun missed are unmatched. Both walks are key-guarded from the same seeds, so a deeper E-147 walk reaching more is expected to land in E-151's set only if E-151 exhausts the failing steps at depth 7-8, which E-151 (20 000-expansion cap, frontier not closed) does not show. Item 1 is therefore "decided for the 25 E-151 failing steps; E-147's 9 (c2) match exactly, c1 12 of the 13 not reproduced, so open only for those missing c1 steps". Not closed as a whole.
3. E-149 -> E-175/E-160 (done, accepted point 5): dropped. I read neither E-160 nor E-175 for whether they use the reverse search; E-149 is measured at n = 6 class 0 on the tilting-only graph. The reverse-loss item is a proposal for the toolsmith's list, with the effect on "no join found" negatives marked as a guess.
4. Patch file (done): "line 1410" corrected to 1411; HEAD hash updated (the file said 1eb46ea; the OLD blocks re-checked at HEAD c953a49, all 9 occur once). Point 3 (the 5-vertex counterexample behind "Example 3.4 silently needs one") is r055's, not re-run; patch 2 labels it "our derivation, E-171".

Scope reworded as the referee suggests; the "move the rest out" wording now covers items 3 and 4 as proposals only (item 3 is a ledger-wording choice; item 4 rests on the toolsmith's own read).

## Claim

(a) Five groups of stale citation text remain in `research/literature/` and `research/syntheses/001`; each is written as an old/new patch in `workshop/rounds/058/scholar_litfixes.md` (4 required, 2 optional; one file, `rickard...md`, had no error and gets a provenance clause only). The "no monomial" correction is already in the Pavon summary at line 131 but not at the quoted statement (line 70) nor in the README row (line 37); patches 2 and 4 put it there.

(b) T5 as listed in STATE has four open items. After E-171 and E-174 only one still bears on whether H-015's premise holds; I propose narrowing T5 to it (item 1 as above), proposing the other two for reassignment, and keeping H-015 OPEN at round 060 (no closure). It could be refuted by: a single J != 0 key-keeping step whose child is provably not derived equivalent to the parent's class (E-147 shows such steps exist; every one tested so far is joined to an LNA, under the J = 0 premise).

## Evidence

State of each T5 open item (STATE.md line 23), with what the record now says.

| T5 item | status after E-171/E-174 | recommendation |
|---|---|---|
| skeptic's hand rebuild of the 13 class-1 E-147 steps | done at hand level for c2 9/9, c1 8/8 distinct failing (parent, v); the other 8 c1 parents have parallel arrows, which the hand constructor rejects (EXPERIMENTS line 183, round 046); those 8 are decided by E-174 at End(T) level, 16/16 c1 and 9/9 c2 iso. Line 183 says the "13 + 9 E-147 steps" were not matched one by one; round 058 rerun partly does | E-147's 9 c2 = E-151's 9 (set equal) and 12 rerun c1 (parent, v) pairs are among E-151's 16 (round 058 rerun); E-147's own 13 c1 not reproduced (rerun gave 8 Cartan-distinct), so keep as a note "match the 13 c1 to E-151", not closed |
| End(T) = the repo's rewrite on J = 0 steps (E-171 "what stays unproved") | conjecture; supported on 370/370 path edges + 35/35 at n = 10 (E-174), but E-174 itself says the comparison does not discriminate at failing steps, and its dims/relation set are inputs | **keep: this is the one live item** (see Next) |
| orbit data giving D = 0 (s = 10 in c2), why c_2 = 0, the s = -2 family (E-145, E-148) | explanatory theory for the law "J != 0 leaves the key", which is class 0 only and not needed for H-015's status | move out of T5 to a theorist backlog; low priority |
| reverse-search loss of 10.3% of edges (E-149) | a completeness defect of the reverse tilting search (same-vertex opposite step lands on a different algebra of the same key), so it may affect "no join found" negatives (a guess: E-149 is n = 6 class 0 on the tilting-only graph; I did not check whether E-160 or E-175 use the reverse search), not the premise | propose to the toolsmith's list under T10/S-1 negatives |

Why H-015 stays OPEN rather than closes. The ledger's own status line (HYPOTHESES) already records that the literal claim "key-keeping gate-admitted step is a tilting step" fails (E-147, E-151: 16 of 80 978 at n = 7 c1, 9 of 79 143 c2). What survives is the weaker claim "guarded walks stay in one derived class", for which: (i) J = 0 steps are derived equivalences when End(T) is the rewrite (AI 2.32(b) + generation, E-168, E-171); (ii) all 25 failing children are joined to an LNA by J = 0 paths, so their class membership no longer depends on the J != 0 step itself. That makes (ii) conditional on (i)'s End(T) identification, so item 2 above is exactly the gap, and the evidence beyond n = 7 (classes 1, 2) and the one n = 10 start is capped samples (E-155, n = 8 c2: 0 of 104 629). A closure would claim more than a sample at n <= 8 plus one n = 10 start.

Why not "narrow H-015 itself to REFUTED in its stated form": the original statement (H-015 title: "sufficient, not merely necessary") is about the Coxeter polynomial guard; read as "guard-admitted step is tilting" it is false (E-147), read as "stays in the class" it is open. Wording of the status line is the chair's call; I suggest the ledger line say "tilting-test reading refuted; class-membership reading open, premise: End(T) = rewrite on J = 0 steps".

Literature that bears on the live item. Oppermann (1504.02617, summary only; no LaTeX in `sources/`, theorem numbers UNVERIFIED) gives End(T) as a dg quiver with differential, with no admissibility hypothesis; the identification "dg quiver after reduction has no surviving non-zero-degree arrow, so it is an algebra" is the summary's reading of its Theorem 1.1. If the repo's seven-step rewrite is shown to equal Oppermann's reduction on J = 0 steps (where no negative-degree arrow survives, by AI 2.32 via Hom(T,T[<0]) = 0), item 2 is a citation plus a short check, not a conjecture. I have not checked that equality.

## Reproduction

Round 058 response (item 1 containment): scratchpad copies of `rounds/043/skeptic_c2.py` (keeps all J != 0 steps as (canonicalKey(parent), v) and the D = 0 Cartan classes) run as `7 1 2 500 30 guard` and `7 2 3 500 30 guard` (500 s each, 4 jobs in parallel), and `rounds/046/skeptic_collect.py 7 {1,2} 20000` with `ck=canonicalKey(parent)` added to the record (about 450-550 s each); then set inclusion of (ck, v). Output: c2 9 = 9; c1 12 in 16 (4 E-151 only). Rerun under load, so budget-limited: 56/59 distinct J != 0 steps against E-147's 67/64. Scripts not committed (scratchpad).

Patch facts (no script). Verification of the patch facts (run from the repository root):

```
grep -n -i "monomial" research/literature/sources/2509.12983v2-pavon-chz-criterion.tex   # no output
grep -n "label{cor:path-algebra}\|label{ex:path-algebra}" research/literature/sources/2509.12983v2-pavon-chz-criterion.tex
grep -n "label{when silting is tilting}" research/literature/sources/1009.3370v3-aihara-iyama-silting-mutation.tex
```
Instant. Old texts in the patch file were copied from HEAD 1eb46ea and re-verified at HEAD c953a49.

## Prior record

E-171 (the audit) lists the corrections but gives them as line numbers, not text; the pending list there names `README.md` 26, 37, `1504...:133`, `rickard...:9`, `syntheses/001`. Round 056 proceedings keep H-015 OPEN for the same premise. Nothing here is a new result; the T5 triage is a reading of the record.

## Code changed

None in the repo. New file: `workshop/rounds/058/scholar_litfixes.md` (patch 3 and the header corrected after review). Scratchpad scripts only.

## Next

- Chair: apply patches 1-4 (patches 5-6 optional); at round 060 restate T5 as the single item "End(T) = rewrite on J = 0 steps", keep H-015 OPEN with the status wording above, and reassign the orbit-data and reverse-loss items.
- Scholar (next round, if asked): fetch 1504.02617 LaTeX and compare Oppermann's reduction with `quivermutation` steps 1-7 on J = 0 steps; state the theorem numbers only after matching.
- Toolsmith: an equal-dims wrong-algebra control with a parallel pair (E-174's open item) would test the comparison's power.
- Skeptic: whether "stays in the class" at J != 0 steps can be separated from "child joined to an LNA by J = 0 paths" (the second is all the record proves).
