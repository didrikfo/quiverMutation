# After E-171 the key-guard thread (T5, H-015) has one open item that is a proof, not a computation; narrow it to that and move the rest out; do not close H-015

author: scholar · round: 058 · kind: proposal
thread: T5 · bears on: H-015, E-171, E-174
scope: literature and ledger reading only; no computation, no arXiv fetch (the two LaTeX sources in `research/literature/sources/` sufficed). Read: STATE.md T5/T10, HYPOTHESES H-015, EXPERIMENTS E-147/E-149/E-171/E-174 and the entry at line 183 (round 046 skeptic rebuild). Not read: the round 056 and 057 submissions, E-157..E-169 beyond their STATE summaries.

## Claim

(a) Five groups of stale citation text remain in `research/literature/` and `research/syntheses/001`; each is written as an old/new patch in `workshop/rounds/058/scholar_litfixes.md` (4 required, 2 optional; one file, `rickard...md`, had no error and gets a provenance clause only). The "no monomial" correction is already in the Pavon summary at line 131 but not at the quoted statement (line 70) nor in the README row (line 37); patches 2 and 4 put it there.

(b) T5 as listed in STATE has four open items. After E-171 and E-174 only one still bears on whether H-015's premise holds; I propose narrowing T5 to it, handing the other three to their natural owners, and keeping H-015 OPEN at round 060 (no closure). It could be refuted by: a single J != 0 key-keeping step whose child is provably not derived equivalent to the parent's class (E-147 shows such steps exist; every one tested so far is joined to an LNA, under the J = 0 premise).

## Evidence

State of each T5 open item (STATE.md line 23), with what the record now says.

| T5 item | status after E-171/E-174 | recommendation |
|---|---|---|
| skeptic's hand rebuild of the 13 class-1 E-147 steps | done at hand level for c2 9/9, c1 8/8 distinct failing (parent, v); the other 8 c1 parents have parallel arrows, which the hand constructor rejects (EXPERIMENTS line 183, round 046); those 8 are decided by E-174 at End(T) level, 16/16 c1 and 9/9 c2 iso. Not matched: n = 8, "13 + 9 E-147 steps one by one" (line 183 says so) | close as a T5 item; the residue (n = 8 matching) is a note, not a thread |
| End(T) = the repo's rewrite on J = 0 steps (E-171 "what stays unproved") | conjecture; supported on 370/370 path edges + 35/35 at n = 10 (E-174), but E-174 itself says the comparison does not discriminate at failing steps, and its dims/relation set are inputs | **keep: this is the one live item** (see Next) |
| orbit data giving D = 0 (s = 10 in c2), why c_2 = 0, the s = -2 family (E-145, E-148) | explanatory theory for the law "J != 0 leaves the key", which is class 0 only and not needed for H-015's status | move out of T5 to a theorist backlog; low priority |
| reverse-search loss of 10.3% of edges (E-149) | a completeness defect of the reverse tilting search (same-vertex opposite step lands on a different algebra of the same key), so it affects "no join found" negatives (E-175, E-160 misses), not the premise | move to the toolsmith's list under T10/S-1 negatives |

Why H-015 stays OPEN rather than closes. The ledger's own status line (HYPOTHESES) already records that the literal claim "key-keeping gate-admitted step is a tilting step" fails (E-147, E-151: 16 of 80 978 at n = 7 c1, 9 of 79 143 c2). What survives is the weaker claim "guarded walks stay in one derived class", for which: (i) J = 0 steps are derived equivalences when End(T) is the rewrite (AI 2.32(b) + generation, E-168, E-171); (ii) all 25 failing children are joined to an LNA by J = 0 paths, so their class membership no longer depends on the J != 0 step itself. That makes (ii) conditional on (i)'s End(T) identification, so item 2 above is exactly the gap, and the evidence beyond n = 7 (classes 1, 2) and the one n = 10 start is capped samples (E-155, n = 8 c2: 0 of 104 629). A closure would claim more than a sample at n <= 8 plus one n = 10 start.

Why not "narrow H-015 itself to REFUTED in its stated form": the original statement (H-015 title: "sufficient, not merely necessary") is about the Coxeter polynomial guard; read as "guard-admitted step is tilting" it is false (E-147), read as "stays in the class" it is open. Wording of the status line is the chair's call; I suggest the ledger line say "tilting-test reading refuted; class-membership reading open, premise: End(T) = rewrite on J = 0 steps".

Literature that bears on the live item. Oppermann (1504.02617, summary only; no LaTeX in `sources/`, theorem numbers UNVERIFIED) gives End(T) as a dg quiver with differential, with no admissibility hypothesis; the identification "dg quiver after reduction has no surviving non-zero-degree arrow, so it is an algebra" is the summary's reading of its Theorem 1.1. If the repo's seven-step rewrite is shown to equal Oppermann's reduction on J = 0 steps (where no negative-degree arrow survives, by AI 2.32 via Hom(T,T[<0]) = 0), item 2 is a citation plus a short check, not a conjecture. I have not checked that equality.

## Reproduction

No script. Verification of the patch facts (run from the repository root):

```
grep -n -i "monomial" research/literature/sources/2509.12983v2-pavon-chz-criterion.tex   # no output
grep -n "label{cor:path-algebra}\|label{ex:path-algebra}" research/literature/sources/2509.12983v2-pavon-chz-criterion.tex
grep -n "label{when silting is tilting}" research/literature/sources/1009.3370v3-aihara-iyama-silting-mutation.tex
```
Instant. Old texts in the patch file were copied from HEAD 1eb46ea.

## Prior record

E-171 (the audit) lists the corrections but gives them as line numbers, not text; the pending list there names `README.md` 26, 37, `1504...:133`, `rickard...:9`, `syntheses/001`. Round 056 proceedings keep H-015 OPEN for the same premise. Nothing here is a new result; the T5 triage is a reading of the record.

## Code changed

None. New file: `workshop/rounds/058/scholar_litfixes.md`.

## Next

- Chair: apply patches 1-4 (patches 5-6 optional); at round 060 restate T5 as the single item "End(T) = rewrite on J = 0 steps", keep H-015 OPEN with the status wording above, and reassign the orbit-data and reverse-loss items.
- Scholar (next round, if asked): fetch 1504.02617 LaTeX and compare Oppermann's reduction with `quivermutation` steps 1-7 on J = 0 steps; state the theorem numbers only after matching.
- Toolsmith: an equal-dims wrong-algebra control with a parallel pair (E-174's open item) would test the comparison's power.
- Skeptic: whether "stays in the class" at J != 0 steps can be separated from "child joined to an LNA by J = 0 paths" (the second is all the record proves).
