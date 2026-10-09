# Digest

What happened, one entry per round, newest first. Written by the chair for the
human coming back after a while: what was claimed, what survived review, what
was promoted into `research/`, what the chair needs from you. Each entry at
most `max_digest_entry_lines` lines, linking to `rounds/NNN/` for the rest.

---

## Round 053 -- 2026-10-09 -- ordinary
Chair settled round 052's questions (agenda approved; PDFs wanted; `canonicalKey` cap not yet). Skeptic, toolsmith, scholar on T10; all three minor revision, responses done.
- **Skeptic (E-164, accepted):** the Hom(T,T[±1]) test accepts all 33 edges of the 3 E-158 paths (Cartan level, J = 0 premise); 25 of 25 failing children = 23 by an edge-tested path + 2 (c1 13, 15) by key equality.
- **Toolsmith (E-165, accepted narrowed):** End(T) as a quiver with relations is iso to the next algebra on 13/13 E-161 edges, arrows generate End(T) on each; the same holds at the 8 decided of 16 J != 0 failing steps, so this check does not see the J = 0 premise. 8 parallel-arrow cases undecided.
- **Scholar (note):** arXiv still 403; from memory the J = 0 step needs Hom in K^b(proj A) only; the live gap is generation of K^b(proj) by T and coverage. Provenance lines of the two literature summaries are inconsistent.
- Consequences: H-015 stays OPEN; "25 of 25" must keep "under the J = 0 premise (generation assumed)".
- **Question for the human:** can you supply PDFs of arXiv:1009.3370 and 2509.12983? Recommend yes; nothing else needs you.

## Round 052 -- 2026-10-09 -- conference
All six personas gave position statements (no runs). Details in `rounds/052/`.
- **Common view:** the weakest claim is the J = 0 premise (a `tiltingPlus` step is a derived equivalence; Hom test Cartan-level, generation assumed, AI 2.31/2.32 and CHZ 3.6 unread). The key-guard law is class 0 only. "25 of 25 children joined" must carry the premise.
- **Ledger:** no status changes. H-015 OPEN (open point: vertex-level vs whole-T tilting), H-010 SUPPORTED, H-017 OPEN, H-020 SUPPORTED (theorist: is the rule table needed for any verdict?), H-021 OPEN. No F- or R- entries, no thread closed.
- **Proposed agenda:** (1) non-Cartan test of the premise (End(T), E-158 Hom replay); (2) group-A witness path; (3) T3/T8 P vs Q finer invariant with power control (outside T10; four personas named it); (4) H-020 rules = [] ablation; (5) T7/H-010 k = 4 breadth.
- Step 0.5: round-051 questions unanswered; decided no overnight, PDFs wanted but literature parked.
- **Questions for the human:** approve or change the agenda in `STEERING.md` (recommend approve); **PDFs of arXiv:1009.3370, 2509.12983** if you can supply them; `canonicalKey` cap 5040 (recommend not yet, with the docstring rewords in a toolsmith round).

---

## Round 051 -- 2026-10-09 -- ordinary
Called: toolsmith (T10 i), experimentalist (T10 ii), theorist (T4). Referees: skeptic (x2), scholar. All reviews accept or minor revision; responses done. Details in `rounds/051/`.
- **toolsmith**: under the J = 0 premise, c1 children 14/15 (key b32eca) are joined to an LNA/dual by a replayed path of total length 13 (7 + 6), so **all 25** failing n = 7 children are joined. The skeptic replayed all 13 edges with the independent Hom test. Hom-tested only for the E-155 paths and this one. Promoted **E-161**.
- **experimentalist**: `merges.py 10 --depths 3 4 5` found no link (repeats E-032); F-037 (one of two n = 10 merges) uses no J != 0 step on 5 of 19 paths. The group-A merge is unreplayed; `meetingPoints` compares labelled keys. Promoted **E-162** (narrowed).
- **theorist**: floating rules of width 6..8 hold at one length each (0 failures in 43 116 applications; gaps listed); the H6 ablation changes nothing in the reduced walk, which is already F-032, but does in a rules-only walk. Promoted **E-163** (narrowed).
- Consequences: H-015 stays OPEN, status line extended with a pointer (the failing children look like class members, conditional on the premise). Docstring rewords (guard, `canonicalKey`) still pending.
- Step 0.5: round-050 questions unanswered; decided no overnight, PDFs wanted but literature parked.
- **Questions for the human:** (1) overnight `merges.py 10 --depths 5 6 7 --witness` (recommend no; a group-A witness path first by a relabelling-aware search); (2) **PDFs of arXiv:1009.3370 and 2509.12983** if you can supply them.
- Round 052 is a conference.

---

## Round 050 -- 2026-10-09 -- ordinary
Called: toolsmith (T10 i, depth-7 ball for the misses), skeptic (independent test of the E-155 premise), maverick (T6/H-017 signature, breadth slot). Referees: theorist, experimentalist, scholar; all minor revision, all answered.
- **toolsmith**: under the J = 0 premise child 12 and 3 of the 5 misses are joined (total <= 12): 23 of 25 children joined; c1 14/15 (one key) still open. `canonicalKey` returns None on parallel-arrow bundles, blinding the child side; its docstring is wrong. Promoted E-158 (narrowed: no-key control not built).
- **skeptic**: an independent Hom(T,T[m]) test accepts all 324 edges of the 40 printed E-155 paths and rejects the 25 failing J != 0 steps (Ladkani's criterion restated; generation assumed). Promoted E-159.
- **maverick**: the H-017 Euler signature is Cartan-determined, so a class invariant blind to (cords, relations); separates 16 n = 11 LNAs from quipu polynomials (new), 2 at n = 10 (known). H-017 stays OPEN. Promoted E-160.
- Consequences: T10 premise now independently checked on the E-155 paths; guard-docstring and `canonicalKey`-docstring rewords wait for a toolsmith round; no hypothesis or finding changed.
- Decided for the committee: no overnight (sharded instead); PDFs cannot be arranged from here.
- **Questions for the human**: (1) overnight c1 14/15 at horizon 13, recommend no (toolsmith tries depth-6 target ball first); (2) PDFs of arXiv:1009.3370, 2509.12983 if you can supply them.

## Round 049 -- 2026-10-09 -- ordinary

Worked: toolsmith (T10 i), experimentalist (T10 ii), theorist (T4, breadth slot). Referees: skeptic (x2), scholar. All three: minor revision, answered, accepted narrowed. Details in `rounds/049/`.
- **Toolsmith:** a tilting-only meet-in-the-middle search joins **19 of the 25** key-keeping, `tiltingPlus`-failing n = 7 children (c2 8/9, c1 11/16) to an LNA of their class by printed paths of length 6-11, and all 25 parents too; 5 miss at bound 11, 1 undecided. Positive controls 12/12 (length 9) and 3/3 (length 11); skeptic replayed 7 paths. Conditional on the untested premise that J = 0 + `tiltingPlus` steps are derived equivalences. Promoted **E-155**.
- **Experimentalist:** the 25 failing edges lie on no seed-to-seed walk shorter than 15 (shortest merge: 4); 11 are leaves of the capped graph. Partly a restatement of E-149's depths. Promoted **E-156**.
- **Theorist:** six falsifiable hypotheses behind the H-020 rule table; H1 (length independence) holds for width <= 5 to length 11 (0 failures in 11 170); orbit sizes of core `45` at n = 13; "outside the derived class" defined (needs a derived invariant with a power control; none exists). Promoted **E-157**.
- **Consequences:** H-015 stays OPEN with a pointer: the key guard's failures look like class members, not leaks (conditional). Docstring not yet reworded.
- Step 0.5: round-048 questions unanswered; agenda approved, literature parked, no overnight.
- **Questions for the human:** (1) overnight depth-7 child ball for the 5 misses (recommend no; toolsmith shards it first); (2) **PDFs of arXiv:1009.3370 and 2509.12983** if you can supply them.
- Next round 050 is ordinary; 052 is the next conference.

---

## Round 048 -- 2026-10-09 -- conference

All six personas wrote position statements (`rounds/048/`); nothing run, nothing promoted to E-/F-/R-.

- **Convergence:** all six name one question: are the 25 key-keeping, `tiltingPlus`-failing n = 7 children (E-145, E-149) in the derived class? Preferred test: a tilting-only path search back to an LNA with a positive control (a hit settles it; a miss is only bounded). Weakest claim named by all: "the key guard keeps a walk in one derived class" (class 0 only; E-153 one class; AI/CHZ citations unverified).
- **Disagreement:** whether H-017 (T6) can be revived as a second invariant; scholar and maverick warn its Euler signature is Cartan-determined and may fail as E-152 did. First step: check what the invariant is.
- **Ledger:** H-015 stays OPEN (no R-entry: it is about derived equivalence); H-021 stays OPEN, status line rewritten short, history moved into the body; no other status change, no F-entry, no thread closed.
- **Proposed agenda:** (1) T10 tilting path back + skeptic hand-check + theorist's exact "outside the class"; (2) n = 7 c1, c2 edge tally, `merges.py --witness` on a real link; (3) H-017 power check; (4) T4 falsifiable statement; (5) S-1 n = 15 sizing.
- Step 0.5: round-047 questions unanswered; decided no depth-9 overnight, literature parked.
- **Questions for the human:** (1) **approve or change the proposed agenda in `STEERING.md`** (recommend approve); (2) **PDFs of arXiv:1009.3370 and 2509.12983** if you can supply them; (3) overnight: none.
- Next round 049 is ordinary; 052 is the next conference.

---

## Round 047 -- 2026-10-09 -- ordinary

Worked: experimentalist (T10 ii), toolsmith (T10 iii), maverick (T1/T2 breadth slot). Referees: skeptic, theorist, scholar. Details in `rounds/047/`.
- **Experimentalist:** E-094's n = 8 class 2 depth-8 guarded walk reproduces; 0 of 104 629 key-kept edges fail `tiltingPlus` (the 2 failures are key-refused); depth 9 not run; slice 2 not independently re-run. Promoted **E-153**.
- **Toolsmith:** the 10 F-041 n = 8 merges have shortest witness 3 + 3 = 6 (inside the key-guarded gate graph), 60/60 edges J = 0; `merges.py --witness` added (opt-in, tested). Little beyond E-148. Promoted **E-154**.
- **Maverick:** note, nothing run: T1/T2 have no cheap census left; reopen only for a rule for k(c) predicting a held-out core. Referee corrected citations. Nothing promoted; H-021 header rewrite goes to the ledger.
- Consequences: none; H-015 stays OPEN (these point the E-148 way, E-145/E-149 stand at n = 7). H-015 entry gains a pointer.
- Next: round 048 is a **conference** (ledger: H-021 status, H-015, T10).
- **Questions:** (1) overnight depth-9 n = 8 c2 walk? Recommend no; do the n = 7 c1, c2 tally first. (2) PDFs of 1009.3370 and 2509.12983 if you can supply them.

## Round 046 -- 2026-10-08 -- ordinary

Worked: theorist (T7 revision), skeptic (T10 i), scholar (T3/T8 breadth slot). Referees: maverick, experimentalist, toolsmith. Details in `rounds/046/`.

- **theorist (E-150, E-151, accepted narrowed):** for an interior pair `(8:3)(9:3)` only a bystander sharing >= 2 arrows (a run of three) lowers the overlap (21 + 14 placements inert at k <= 2; 5 of 5 runs of three lower by k = 3, `(10:4)` only at k = 3). Refines F-022. One-step locality: for all LNAs n <= 10 exactly 2 of n mutations give an LNA, changed starts within {v-2, v-1, v}. H-010 still unproved.
- **skeptic (E-152, accepted as a null):** the 25 key-keeping, `tiltingPlus`-failing children (n = 7 c1, c2) have equal Cartan invariants and an integral P, but every pair of LNAs in those key classes does, so the test says nothing about class membership. Whether they leave the class is still open.
- **scholar (note, T3/T8):** key-coarser counts are 7 / 9 / 10 (n = 10 / even / odd); the "parity classes" are two mirror-closed orbits P, Q with one key (E-077, E-080); 5046/5056 outside. Marked dormant, not closed; no invariant separates P from Q.
- **Consequences:** none new. The `mutationSearchDepthFirst` docstring stays unreworded (unsupported as tilting, untested as class); H-015 stays OPEN.
- Step 0.5: round-045 questions unanswered; I kept `tiltingPlus` unpromoted and literature parked (arxiv.org still 403).
- **Questions for the human:** (1) promote `tiltingPlus` as opt-in `tiltingGuard`? Recommend not yet. (2) **arxiv.org access or PDFs of 1009.3370, 2509.12983?** Recommend yes. (3) Overnight: none; the deep E-094 replay is sized first.
- Next round 047 is ordinary; 048 is the next conference.

---

## Round 045 -- 2026-10-08 -- ordinary

Worked: experimentalist, toolsmith (guard audit T10), theorist (T7 breadth slot). Referees: skeptic (x2), maverick. Details in `rounds/045/`. Also done per the special request: STATE.md rewritten from scratch; H-015 ledger pass; thread T10 opened.

- **experimentalist (E-148, accepted narrowed):** the 10 F-041 n = 8 merges also meet in a tilting-only walk (gate + `tiltingPlus`, key guard off); all 2396 edges are J = 0 and key-keeping; guarded depth-4 edges from LNAs at n = 6..8 (n = 8 partial) are all J = 0. So no merge examined depends on a J != 0 step. Deep merges (distance >= 5), where E-145's steps live, are not covered. Library `meetingPoints` agrees on 10 of 10.
- **toolsmith (E-149, accepted narrowed):** an added J = 0 / `tiltingPlus` check costs 8-11% of a step. 16 of 80 978 (n = 7 c1) and 9 of 79 143 (c2) key-keeping steps fail `tiltingPlus`, first at parent depth 7-8; none at n = 8 samples. Reconciled with E-084 (its walks stopped at the first rejecting level). Whether those children leave the derived class is untested.
- **theorist (T7, major revision, nothing promoted):** no proof of H-010; the proposed "run-of-three" lemma mostly restates F-022 and is tested only to k = 2; closure of T7 not accepted. Revision due next round.
- **Consequences:** **H-015 moved SUPPORTED -> OPEN** (sufficiency is refuted for tilting at n = 7 c1, c2; class membership untested). The docstring of `search.mutationSearchDepthFirst` ("coxeterGuard keeps the walk in one derived class") is unsupported as a tilting claim; rewording waits on the skeptic's out-of-class test (T10).
- Step 0.5: round-044 questions unanswered; I approved the proposed agenda with T10 first, no overnight, literature parked.
- **Questions for the human:** (1) **promote `tiltingPlus` into the library as an opt-in `tiltingGuard` keyword?** I recommend yes after the out-of-class test; (2) **allow arxiv.org or supply PDFs of 1009.3370 and 2509.12983** (recommend yes).
- Next round (046) is ordinary; 048 is the next conference.

---

## Round 044 -- 2026-10-08 -- conference

All six personas wrote position statements (`rounds/044/`); nothing run, nothing promoted.

- **Convergence:** five of six name the same live question: E-145 found gate-admitted J != 0 steps that keep the LNA key at n = 7 classes 1 and 2, so the "key guard" law holds for class 0 only and may be no evidence for H-015. Four ask the skeptic for an independent hand rebuild of the class-1 steps (E-145 rests on one script).
- **Weakest claims named:** the key-guard law (the guard is also the BFS filter); E-145 class-1 count; H-015 cited from memory (AI 2.31/2.32, CHZ 3.6 unread); the S-1 K-threshold law (two lengths).
- **Proposed agenda (round 044):** (1) skeptic rebuild + guard-off n = 7 census with `tiltingPlus` and Cartan per step; (2) theorist: orbit data for D = 0, why c_2 = 0; (3) toolsmith: can a key-preserving step fail `tiltingPlus`, reverse loss by depth; (4) scholar retries arXiv; (5) S-1 n = 15 K0 = 5 sizing. Maverick alone favoured S-1 first.
- Step 0.5: round-043 questions unanswered; I decided agenda replaced by the proposed one, no overnight, literature parked with an arXiv retry.
- **Questions for the human:** (1) **approve or change the proposed agenda in `STEERING.md`** (recommend approve); (2) overnight: none; (3) **PDFs of arXiv:1009.3370 and 2509.12983, if you can supply them.**

---

## Round 043 -- 2026-10-08 -- ordinary

Worked: skeptic, experimentalist, toolsmith. Referees: theorist, skeptic, maverick (all minor revision; all accepted with qualifications). Details in `rounds/043/`.

- **skeptic (E-145):** at n = 7 classes 1 and 2 there are gate-admitted J != 0 steps that keep the class key (13 of 67 and 9 of 64 distinct). They are walk descendants, fail tiltingPlus, and are Cartan-incongruent. So E-140's "none keeps the key" and the E-138/E-141 law are class-0 statements, and the key guard is no evidence for H-015 off J = 0 steps. The referee reproduced the counts and rebuilt the class-2 steps by hand; the class-1 steps rest on one script. The orbit relation e_i = F^s e_w fails on most H1/H2 steps in c1, c2.
- **experimentalist (E-146, null):** at n = 8 none of 27 J != 0 steps (25 distinct) has the E-143 shape (|out v| = 2 or |supp J| = 2); Q's lowest term is x^3 for some class-1 steps. Shallow, time-capped sample. Corrects E-143: "orbit relation never absent" holds only at n = 6, 7 c0.
- **toolsmith (E-147):** the reverse search passes a positive control at depth 2-4 (12/12). 154 of 1500 edges (10.3%) are lost in reverse, never to a filter: the opposite step lands on a different same-key algebra. The control is shallow, so deep reach is not shown.
- Step 0.5: the round-042 questions were unanswered; I decided agenda kept, literature item parked, no overnight.
- Next round (044) is a conference.
- **Questions for the human (recommendations in `rounds/043/proceedings.md`):** (1) keep the round-040 agenda, item 1 now "why key-preserving J != 0 steps exist at n = 7 c1, c2 and what that means for H-015" (recommend yes); (2) overnight: none, the 3 h reverse job waits for a deep control (recommend no); (3) **PDFs of arXiv:1009.3370 and 2509.12983 in `research/literature/`, or item 3 stays parked**.

---

## Round 042 -- 2026-10-08 -- ordinary

Worked: scholar, theorist, maverick. Referees: skeptic, experimentalist, toolsmith (all minor revision). Details in `rounds/042/`.

- **scholar (note, nothing promoted):** arXiv is still blocked by the proxy (403), so Aihara-Iyama 2.31/2.32 and CHZ 3.6 stay unverified. By hand, "monomial" is needed for the path-wise Cor 3.6; the flag was already in `research/literature/2509.12983` and E-066/E-122, so it is not new.
- **theorist (E-143, accepted with qualifications):** under shape hypotheses on the v-row, Q(x) = x(adj S_ii - adj S_wi - adj S_iw) is proved. With F e_w = e_i, Q's x^2 coefficient is 1 iff c_2 = 0, also proved. The orbit relation e_i = F^s e_w and c_2 = 0 are observed only. Reproduced at n = 6, 7 and guard off; "exactly when" is proved for s = 1 only.
- **maverick (E-144, accepted with corrections):** at n = 12 the S-1 K-threshold law holds. The lone-3 class (2746 LNAs, one orbit) sends K >= 4 ends to one n = 11 key class and K >= 3 ends to two. The single-3 table gives first failures at n = 11, 13, 15 for K0 = 3, 4, 5; n = 15 is a prediction only. Reproduced by the referee.
- Step 0.5: the round-041 questions were unanswered. I decided: agenda kept, no overnight, scholar tries the arXiv fetch (blocked).
- **Questions for the human (recommendations in `rounds/042/proceedings.md`):** (1) keep the round-040 agenda (recommend yes); (2) **can you place PDFs of arXiv:1009.3370 and 2509.12983 in `research/literature/`? Without them agenda item 3 stays parked**; (3) overnight: none.

---

## Round 041 -- 2026-10-08 -- ordinary

Worked: experimentalist, theorist, toolsmith. Referees: skeptic, scholar, maverick (all minor revision; all accepted with qualifications). Details in `rounds/041/`.

- **experimentalist (E-140):** with the key guard off, none of 278 gate-admitted J != 0 steps at n = 6, 7, 8 (classes 0-2, 5 of 9 cells non-vacuous) has a child with the class key; a random-algebra control keeps the key in 55 of 1264 steps, so the test can say yes. No LNA-parent control. The "282" in the draft was a miscount (278).
- **theorist (E-141):** "J != 0 forces a key change" is false in general (4-vertex counterexample, but its parent is not an LNA key, so the walk law is untouched). On LNA walks the Coxeter polynomials differ by a term starting at x^2 with coefficient 1; no proof. The trace-coefficient route to a proof is closed.
- **toolsmith (E-142):** `--tilting-only` meet now has a positive control and closure flags. 16 of 16 hits share 0 keys with the tilting-only LNA side, forward and reverse, but no search closed: a bounded miss. E-137's small hit sides were a 12 s cap. The reverse search loses 12% of edges and has no control.
- **Not worked:** agenda items 3 (literature) and 4 (S-1 n = 12).
- Round 040's questions were decided by the chair: agenda approved, literature fetch allowed if the network does, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** keep the round-040 agenda with item 1 reshaped to "why does Q(x) start at x^2 on walks" (recommend keep); the toolsmith's 3 h reverse-search overnight job (recommend no until a reverse positive control exists); scholar to fetch arXiv:1009.3370 and arXiv:2509.12983 (recommend yes if the network allows).

---

## Round 040 -- 2026-10-08 -- conference

All six personas wrote position statements (`rounds/040/`); no new work, no referees, no `research/` changes.

- **Convergence:** five of six name the same most promising question: does any gate-admitted J_i != 0 step keep the LNA key (E-137/E-138 say no on walks, but the key guard is also the BFS filter, so the refusal is partly built in)? Maverick alone presses S-1 (the K-threshold law at n = 11, 13, 15).
- **Weakest claims named:** E-134's membership claim (unsupported, not refuted); "d_i = 2 whenever J_i != 0" (capped, E-131 counterexamples); H-015 as the tilting test; the AI Thm 2.31 citation (from memory, PDF never read); the S-1 threshold law (two lengths).
- **Proposed agenda (round 040), ranked:** 1 key-guard-off tabulation of J != 0 children at n = 6-8 with positive control; 2 `--tilting-only` meet with closure flag; 3 AI/Ladkani vs `tiltingPlus`; 4 S-1 n = 12 after `--plan`.
- Round 039's questions were decided by the chair: keep the agenda, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** approve the round-040 agenda in `STEERING.md` (recommend approve); may the scholar fetch arXiv:1009.3370 and arXiv:2509.12983 into `research/literature/` (recommend yes if the network allows); overnight: none (recommend none).

---

## Round 039 -- 2026-10-07 -- ordinary

Worked: skeptic, theorist, maverick. Referees: experimentalist (x2), skeptic (all minor revision; all accepted with qualifications). Details in `rounds/039/`.

- **skeptic (E-137):** E-134's meeting path is not a tilting path: the LNA -> M side passes `tiltingPlus` at every step, but the first hit -> M step is gate-admitted, key-preserving and fails `tiltingPlus`/Cartan congruence. With non-tilting steps removed, no hit meets the LNA side (capped). So "the 16 fans lie in the LNA class" is unsupported, not refuted. No positive control yet.
- **theorist (E-138):** on the n = 8 c0 walk all 229 J != 0 rows (192 steps, incl. all 9 with d >= 3) are refused by the key guard; same at n = 6, 7. So E-129's d = 2 and E-135's out(i) = 3 describe parents of refused steps. out(i) = 3 is not forced by the gate (hand algebra with out(i) = 2). Class 0 only.
- **maverick (E-139):** S-1 at n = 13: the lone-3 key class (5023 LNAs) is two orbits (4349 + 674); the 4349-orbit's K >= 4 ends map to two different n = 12 key classes, so transport by deletion fails here, as E-133 predicted.
- **Promoted:** E-137, E-138, E-139 (E-134 annotated). Library: no changes.
- Round 038's questions were decided by the chair: keep the agenda, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** keep the round-036 agenda, item 1 now "does any J != 0 step keep the key; tilting-only meet with a positive control" (recommend keep); overnight: none (recommend none).

---

## Round 038 -- 2026-10-07 -- ordinary

Worked: experimentalist, toolsmith, scholar. Referees: theorist, skeptic, experimentalist (all minor revision; all accepted with qualifications). Details in `rounds/038/`.

- **toolsmith (E-134):** all 16 gate-admitted n = 6 fans of E-132 (d >= 3, J != 0, LNA key) meet an LNA in the key-preserving forward graph (25-27 shared algebras), so E-132's "no LNA reached" was a bounded miss. Equivalence rests on the Coxeter guard (H-015), not the gate; no step checked with `tiltingPlus`. E-129's d = 2 is therefore a forward-walk statement, not a class one. n = 6 classes do not close in 240 s. `maverick_single` crash fixed (test added).
- **experimentalist (E-135):** on the n = 8 c0 capped walk, J != 0 with d >= 3 occurs in 9 rows / 5 algebras, all with out-degree(i) = 3, including two (3,1) rows with no parallel arrows; necessary in the sample, not sufficient (9 of 17-19). Referee reproduced every d >= 3 cell.
- **scholar (E-136):** `C_B = r C_A r^T + H`, `H_{vi} = dim J_i` (28 000 walk steps, 0 failures); d_i is a Cartan entry, not a class invariant, and no published bound is known. Largely E-093/E-127 as observation; the closed formula is new.
- **Promoted:** E-134, E-135, E-136 (E-132 annotated). Library: only `tests/test_toolsmith_single.py`.
- Round 037's questions were decided by the chair: keep the agenda, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** keep the round-036 agenda, item 1 now "check one meeting path with `tiltingPlus`; why out(i) = 3" (recommend keep); overnight: none yet (recommend none; class-0 growth run proposed by toolsmith can wait).

---

## Round 037 -- 2026-10-06 -- ordinary

Worked: theorist, skeptic, maverick. Referees: skeptic, theorist, experimentalist (all minor revision; all accepted with qualifications). Details in `rounds/037/`.

- **skeptic (E-131):** the 5 parallel rows of E-127 are 3 algebras; J_i != 0 comes from one relation on the single non-parallel out-arrow, not from the parallel multiplicity. Three rows past E-129's cap have J_i != 0 with d_i = 4, 4, 5, so "d_i = 2 at J_i != 0" holds only for the capped sample. Referee re-ran the 8.5-minute walk.
- **theorist (E-132):** the gate does not force d_i = 2: a hand algebra with (d, dim J) = (3, 1) is gate-admitted but not on a walk (key not LNA); 183 such n = 6 fans, 167 off the LNA keys, 16 on them, no LNA found by a bounded BFS (not a proof).
- **maverick (E-133):** with `mirrorRow` the S-1 image is still a function of (K, word); the I1/I2 split fits "lone 3 or 7 with K_eff = 3". Key-level prediction that "K >= K0 holds" fails first at n = 11, 13, 15. Also corrects E-125 (94 ends, not 82). Not run at n = 12/13.
- **Promoted:** E-131, E-132, E-133. No library change. Open: `maverick_single.py` label block crashes; d_i tables unsaved.
- Round 036's questions were decided by the chair: approve agenda, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** keep the round-036 agenda, item 1 reshaped to "which derived-class invariant governs d_i" (recommend keep); overnight: none (recommend none).

---

## Round 036 -- 2026-10-06 -- conference

All six personas wrote position statements; no new work, nothing promoted. Details in `rounds/036/`.

- **Convergence:** four of six name the same question: why d_i = 2 whenever J_i != 0 on walks (E-129), with the 28 rows at d_i >= 3 (all J_i = 0) unexplained.
- **Weakest claims named:** the d_i <= 2 bound (empirical, capped data); E-123's "no LNA keys on circuit members" (no base rate, arguably circular); E-130's "0/44 reached" (keys of the 44 unchecked); "room to move" in S-1 (hypothesis).
- **Proposed agenda (replaces round 032):** (1) why d_i = 2 at J_i != 0; (2) audit the 28 d_i >= 3 rows and walk past the cap; (3) n = 7 closure: key check on the 44, then reverse search; (4) the 5 parallel rows of E-127; (5) S-1: the last-letter-3 I1/I2 split, n = 12 after `--plan`.
- Round 035's questions were decided by the chair: keep the agenda, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** please approve or change the round-036 agenda in `STEERING.md` (recommend approve); overnight: none (recommend none).

---

## Round 035 -- 2026-10-05 -- ordinary

Worked: toolsmith, experimentalist, scholar. Referees: skeptic (x2), theorist (scholar accepted; the others minor revision, accepted with qualifications). Details in `rounds/035/`.

- **scholar (E-128):** Hom(N,N[-1]) = {(y_b): y_b in J_t(b), sum b y_b = 0}, so "silting-not-tilting iff some J_i != 0" holds for any A and loopless v, cyclic or not (given AI 2.31); only the dimension count changes on a cyclic quiver. Computed on six cases, referee re-ran. The repo gate sees simple paths only, so E-126's L1 is unguaranteed on cyclic quivers.
- **experimentalist (E-129):** on capped walks (n = 8 c0/c1/c2, n = 9 c0) all 285 rows with J_i != 0 have d_i = 2, dim J_i = 1; the 28 rows with d_i >= 3 have J_i = 0. Thin at the cap's edge; no distinct-algebra counts.
- **toolsmith (E-130):** the n = 7 BFS closure of the both-die classes does not fit one command (frontier ratio about 2.4 and 1.9); 0 of 44 targets reached in 240 s per class ("not found").
- **Promoted:** E-128, E-129, E-130. No library change.
- Round 034's questions were decided by the chair: keep the agenda, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-032 agenda until the round-036 conference (recommend yes); (2) overnight `toolsmith_closure.py run 1 --budget-hours 8`: recommend no, the reverse-direction search is sized first.

---

## Round 034 -- 2026-10-05 -- ordinary

Worked: theorist, skeptic, maverick. Referees: experimentalist, scholar, skeptic (all: minor revision; all accepted with qualifications). Details in `rounds/034/`.

- **theorist (E-126):** on a gate-admitted v, dim J_i <= d_i - 1, so J_i != 0 needs d_i >= 2 and d_i = 2 gives dim J_i = 1; for acyclic algebras Hom(T,T[-1]) = sum J_i (silting-not-tilting iff J != 0, conditional on AI 2.31). "dim J_i = 1 on walks" still needs d_i <= 2 there: a layered algebra has dim J = 2. No cyclic test of the Hom(N,N[-1]) step.
- **skeptic (E-127):** the 5 out-degree 3/4 parallel rows with J != 0 (n = 8 c0) are real Cartan failures; the defect sits at the J support (socle reading not refuted). Control is a 200-row sample; the 5 may be 3 orbits.
- **maverick (E-125):** in the failing n = 11 class the image at n = 10 is a function of (K, core word); 12 core words go I1 at K = 3, I2 at K = 4. "Room to move" is a hypothesis; E-118 already holds most of the counts.
- **Promoted:** E-125, E-126, E-127. No library change.
- Round 033's questions were decided by the chair: keep the agenda, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-032 agenda (recommend yes); (2) overnight: none; the toolsmith sizes the n = 7 BFS closure first (recommend none).

---

## Round 033 -- 2026-10-04 -- ordinary

Worked: scholar, toolsmith, experimentalist. Referees: skeptic (x2), theorist (all: minor revision; all accepted with qualifications). Details in `rounds/033/`.

- **scholar (E-122):** AI 2.32(b) at a vertex reduces exactly to J_i = Hom(S_v, e_iA), with no monomial hypothesis; matches step 7 and E-078. Closes E-066's "no derivation". Caveats: two parents only (coker 0, dim J = 1); Hom(N,N[-1]) not treated. No new obstruction to circuits.
- **toolsmith (E-123):** the LNA-key test has no base rate: walk parents carry the key by construction, and in E-121's layered family 0 of 2 704 members have one (even circuit-free, pure-W). E-121's key absence is therefore not evidence.
- **experimentalist (E-124):** a both-die square exists at n = 6 (gate-admitted, J != 0) but none of 42 enumerated has an LNA key; 48 at n = 7, 1 408 at n = 8; capped n = 7 walks found none ("not found"). dim J_i = 1 on every walk row. Unreconciled: out-degree 3 parallel rows with J != 0 (not shown to be real failures).
- **Promoted:** E-122, E-123, E-124. No library change.
- Round 032's questions were decided by the chair: keep the agenda, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-032 agenda (recommend yes); (2) overnight: none; the toolsmith sizes the n = 7 BFS closure first (recommend none).

---

## Round 032 -- 2026-10-04 -- conference

All six personas filed position statements; no new research, nothing promoted. Details in `rounds/032/`.

- **Convergence:** four of six land on T5: why LNA-derived walks have no long circuits. The socle reading (J_i = H^{-1}(cone), Aihara-Iyama 2.31/2.32) must be checked against step 7, and an invariant (Coxeter polynomial / Euler form) must separate the 962 non-W circuit members from LNA keys.
- **Weakest claims named:** no base rate for LNA keys among circuit members; E-120's "no both-die row at n = 6, 7" is unexplained; "derived equivalence is the obstruction" is circular without a mechanism; K >= 4 deletion rule untested past n = 11.
- **Proposed agenda (round 032):** (1) circuits/invariant/socle check, (2) LNA-key base rate at n = 8 c0, (3) why both-die rows start at n = 8, (4) Gamma_i two-edge bound, (5) S-1 threshold K.
- Round 031's unanswered questions were decided by the chair: keep the agenda, no overnight.

**Questions for you (the chair takes the recommended option if unanswered):** (1) approve the proposed agenda above in STEERING.md (recommend yes); (2) overnight: none (recommend none). Round 033 is ordinary.

---

## Round 031 -- 2026-10-03 -- ordinary

Worked: theorist, skeptic, toolsmith. Referees: experimentalist, scholar, skeptic (all: minor revision; all accepted with qualifications). Details in `rounds/031/`.

- **theorist (accepted, E-121):** J_i is the socle part Hom(S_v, e_iA) = H^{-1} of the mutation cone, so a circuit means "silting but not tilting". The two-term kernel structure alone does not forbid nn 2-cycles or circuits >= 3 (layered hand algebras, n = 6..8), and 962 non-W circuit members plus pendants have no LNA Coxeter key. Necessary test only, no base rate, thin family; no proof.
- **skeptic (accepted, E-120):** the 25 loose W-false rows at n = 6, 7 are all "half-W" (the monomial kills one term only), so J = 0; same pattern as the 23 loose accepts at n = 8; the 61 rejects are the both-terms-die rows, which the capped n = 6, 7 walks never contain.
- **toolsmith (accepted, E-119):** doubled-arrow controls (W-type, cancels, parallel out-arrows, G) are gate-admitted and fail the Cartan check; the code's W is True on "cancels" too and misses a tripled-arrow H chain. The `dim e_iAe_v` count is exact on 33 580 random pairs, lifting E-117's parallel-arrow caveat for n <= 6.
- **Promoted:** E-119, E-120, E-121. No library change.

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-028 agenda (recommend yes); (2) overnight runs: recommend none. Round 032 is a conference.

---

## Round 030 -- 2026-10-03 -- ordinary

Worked: experimentalist, scholar, maverick. Referees: theorist, experimentalist, skeptic (all: minor revision; all accepted with qualifications). Details in `rounds/030/`.

- **experimentalist (accepted, E-117):** at n = 8 classes 0, 1, max `dim e_iAe_v` is 1 to BFS depth 3, 2 first at depth 4, then grows (8 at depth 8, c0). Caveats: acyclic-only depth; dim count unvalidated for parallel arrows.
- **scholar (accepted, E-116):** cone estimate `dim(child) <= max(d, M d)`, M = out-degree at the mutated vertex; referee found 0 violations on 46 586 mutation edges. Thin Hom is forced only at depth 1, so a dimension bound cannot explain "no long circuits".
- **maverick (accepted, E-118):** with the E-115 labels the free-end K >= 3 deletion rule survives at n = 8, 9, 10 but fails at n = 11 (one class splits inside one orbit); K >= 4 holds on the 10 resolved n = 11 classes (418 of 442 sources unresolved).
- **Promoted:** E-116, E-117, E-118. No library change.

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-028 agenda (recommend yes); (2) overnight runs: recommend none. Decided for you (round 029): agenda kept, no overnight.

---

## Round 029 -- 2026-10-02 -- ordinary

Worked: theorist, skeptic, toolsmith. Referees: skeptic, theorist, experimentalist (all: minor revision; all accepted with qualifications). Details in `rounds/029/`.

- **theorist (accepted, E-113):** the skeptic's 22 odd rows are not length-2 ground paths but "half-W" (20) or two loose pendants (6). At n = 8 classes 0, 1 no circuit-graph component has more than 2 edges (no nn 2-cycle, no circuit >= 3), for scalar-1 relations; hand-built D, W-type/G, H cannot occur on walks (their Coxeter keys are in no LNA key set). The general obstruction is not proved.
- **skeptic (accepted, E-114):** the 61 D' rejects are admitted because the gate checks single paths (near-tautology). The loose D' shape also occurs at n = 6, 7 class 0 (5, 20 rows) and is never rejected there; no other out-degree 2 reject in capped walks.
- **toolsmith (accepted, E-115):** the F-047 Smith profile places all 16 unresolved n = 9 and 176 unresolved n = 10 LNAs (n = 9 was already F-047's table), so the S-1 class labels there are independently supported; at n = 11 it places only 24 of 442.
- **Promoted:** E-113, E-114, E-115; glossary: half-W / loose pendants. No library change.

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-028 agenda (recommend yes); (2) overnight runs: recommend none. Decided for you (round 028): agenda approved, no overnight.

---

## Round 028 -- 2026-10-02 -- conference

All six personas wrote position statements (`rounds/028/`); no new work, no referees, nothing promoted.

- **Convergence:** four of six name the same weak claim: the converse of rule W (reject => W) rests on 70 same-shape rejects, the "cancels" branch and parallel-arrow rejects are untested. Theorist and scholar both ask why LNA-derived algebras have no circuit of length > 2; toolsmith and maverick agree the 16 n = 9 and 176 n = 10 unresolved cospectral LNAs block validation of S-1.
- **Proposed agenda (ranked):** (1) no long circuits / converse of W (theorist, scholar); (2) hand-built parallel-arrow and "cancels" controls for W (experimentalist, toolsmith; skeptic checks gate admission of the 61 D' rejects); (3) resolve the cospectral LNAs (toolsmith, experimentalist); (4) S-1: derive the K >= 3 rule (maverick, theorist); (5) carried: mirror chain, `k(33x) = 2x`, n = 15..17. No disagreement between personas.

**Questions for you (the chair takes the recommended option if unanswered):** (1) approve or change the proposed agenda in `STEERING.md` (recommend yes); (2) overnight runs: recommend none. Decided for you (round 027): no overnight, agenda kept.

---

## Round 027 -- 2026-10-02 -- ordinary

Worked: scholar, experimentalist, maverick. Referees: theorist, skeptic, toolsmith (all: minor revision). Details in `rounds/027/`.

- **scholar (accepted, partial proof):** "reject => W" cannot come from step 7 (J is defined by the parent alone) and is false for hand-built algebras (D = E-103's "cancels", G, H). For monomial + two-term relations (scalar 1) J != 0 iff a circuit graph has a circuit; W is its length-2 ground-path case. On walks every J != 0 is a length-2 ground path. Referee: lemma needs "balanced" for scalars != 1.
- **experimentalist (accepted, null extension):** W has 0 mismatches on 32 132 fresh out-degree 2 rows (n = 8 classes 0, 2, 3; n = 9 class 0 prefix), but only the 61 n = 8 class 0 rows are positives; no out-degree >= 3 or parallel-arrow reject; no positive control for parallel arrows.
- **maverick (accepted with caveats; suggested question S-1):** no vertex-deletion rule by position (same-class pairs kept 0.43 vs 0.19 chance at n = 9). Deleting a free end vertex with >= 3 free vertices gives an image class that depends only on the source class at n = 8, 9, 10 (3, 5, 10 resolved classes); all numbers re-run by the referee. Not derived; labelling not independently validated.
- **Promoted:** E-110, E-111, E-112; glossary: circuit graph, free-end strip. No library change.

**Questions for you (the chair takes the recommended option if unanswered):** (1) overnight n = 9 class 0 walk: recommend no; (2) keep the round-024 agenda (recommend yes). Round 028 is a conference. Decided for you (round 026): no overnight, agenda kept.

---

## Round 026 -- 2026-10-02 -- ordinary

Worked: skeptic (revision), theorist, toolsmith. Referees: theorist, skeptic, experimentalist. Details in `rounds/026/`.

- **skeptic (accepted):** the 42 n = 8 class-0 rejects each have a kernel element x that is a two-term sum, killed termwise by one out-arrow and only as a sum by the other; 15 of 42 survive the shorter presentation. Referee: reproduces, marginally new over E-105/E-106; wording fixes applied.
- **theorist (accepted as a conjecture):** rule W (a two-term relation through v whose difference x is killed by the other out-arrow) matches kerdim > 0 on 17 802 out-degree 2 rows with 0 mismatches. Referee re-ran n = 8 class 1 only; the converse rests on 70 same-shape rejects and the "cancels" branch is untested.
- **toolsmith (accepted):** `longSquare` now handles parallel arrows (it missed 4 of 16 steps at n = 8 class 1); a resumable reject walk with `--budget-hours` exists, resume verified at n = 7. n = 9 class 0 would be multi-night and may never close.
- **Promoted:** E-107 (conjecture W), E-108, E-109. No library change.

**Questions for you (the chair takes the recommended option if unanswered):** (1) overnight n = 9 class 0 reject walk: recommend no, test W on a short n = 9 prefix first; (2) keep the round-024 agenda (recommend yes). Decided for you (round 025): agenda kept; no overnight.

---

## Round 025 -- 2026-10-02 -- ordinary

Worked: scholar (revision), skeptic, experimentalist. Referees: experimentalist, theorist, skeptic. Details in `rounds/025/`.

- **scholar (accepted):** at n = 8 class 0 the walk-level "reject iff long square" fails: 61 distinct out-degree 2 rejects, each commuting into one outgoing arrow and killed into the other by a zero relation (called D'), all failing the Cartan test through the real rewrite. The referee checked all 61 independently. n = 9 class 0 (500 s prefix) shows none, so presence is coverage-dependent.
- **skeptic (revise):** the 42 reproduce and are BFS-reached (vacuous: the key test is the BFS filter). "Same mechanism as E-103" is only a shape match: a control has 754 same-shape steps that accept vs 42 that reject. Title said "minimal"; the author's caveat says otherwise. Due: extract the kernel element, test the shorter presentation.
- **experimentalist (runs kept, claim not):** out-degree >= 2 rejects recur at n = 8 class 1. "A second way to break the iff" was wrong (parallel-arrow long squares the test misses); counts shift with machine load (15 vs 38); the zeros at n = 8 class 3 and n = 9 class 0 support nothing.
- **Promoted:** E-105 (scholar), E-106 (the capped runs, with the corrections), glossary "D' reject", a range note on E-100.

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-024 agenda (recommend yes); (2) overnight: none until a checkpoint exists for an n = 9 class 0 walk (recommend none). Decided for you (round 024): agenda kept; no overnight.

---

## Round 024 -- 2026-10-02 -- conference

All six personas wrote position statements (`rounds/024/`); no new work, nothing promoted. Four of six (experimentalist, scholar, skeptic, toolsmith) name the same question: the 42 n = 8 class-0 gate-admitted rejections with out-degree 2 and no long square. Theorist and maverick name the mirror-chain rule (E-104: a fit with no mechanism). Weakest claim named: "reject iff long square" (true at n <= 7 only). Disagreement: the skeptic suspects the 42 are unreachable by key-preserving walks (so the iff may survive on walks); the experimentalist and scholar read them as a mechanism change at n = 8.

Proposed agenda (round 025 onwards): (1) classify the 42, test their reachability, run them through the Cartan rewrite, capped walks at n = 8 c1/c3 and n = 9 c0; (2) mirror chain: sort failing LNAs by big blocker, deeper n = 11/12 shapes; (3) split words/shuttle/`k(33x) = 2x`; (4) cord = commutativity cycle.

**Questions for you (the chair takes the recommended option if unanswered):** (1) approve this agenda in STEERING.md (recommend yes); (2) overnight: none (recommend). Decided for you (round 023 questions): agenda kept; no overnight.

---

## Round 023 -- 2026-10-02 -- ordinary

Worked: skeptic, scholar, maverick. Referees: theorist, experimentalist, skeptic. Details in `rounds/023/`.

- **skeptic:** off the walks, a genuine long relation always makes the step reject (362 apparent counterexamples were redundant presentations). E-100's `hasLongSquare` is only a one-out-arrow test: 26 hand-built gate-admitted rejections have no such square (19 with two out-arrows); two examples are unreachable by the key-preserving walk (n = 6). Random hand-built parents, so not a refutation of E-100's counts. Accepted, scoped.
- **scholar:** derivation sketch: square => reject needs a minimal relation; reject => square also needs out-degree 1; checked on walks at n = 5..7. The referee ran n = 8 class 0 and found 42 out-degree-2 rejections with no square, so the walk-level iff holds at n <= 7 only. **Revise**, not promoted.
- **maverick:** E-101's peeling formula fails on 1/6/24 LNAs with several big relations at n = 8/9/10; a fitted mirror-chain rule has 0 mismatches there and 13 of 13 out-of-sample depth-3/4 predictions up to n = 12 (referee ran n = 12). A fit, not a proof.

Promoted: E-103, E-104.

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-020 agenda (recommend yes); (2) overnight: none (recommend). Decided for you (round 022 questions): agenda kept; no overnight. Round 024 is a conference.

---

## Round 022 -- 2026-10-02 -- ordinary

Worked: experimentalist, theorist, toolsmith. Referees: theorist, skeptic, scholar. All minor revision; I scoped the claims and accepted. Details in `rounds/022/`.

- **experimentalist:** control for the long-sided square: 0 of 479 761 tilting steps (n = 5..7 guarded walks) have it, while all 2 104 rejecting steps (n = 6, 7 class 0) do; strict A5 is not selective. Capped walks, lower bounds; rejecting side from two runs only.
- **theorist:** the cord criterion splits into a depth-1 rule (a relation of >= 3 arrows that is not blocked; fitted among 4 variants, 0 mismatches n = 6..10 and an n = 11 sample) and a peeling depth 1 + min(a, b) for single-big-relation LNAs. E-099's "within 3 steps" is true only for n <= 9 (n = 10 needs 4). "Cord" here means a cycle member. A sketch for why 3 arrows, no proof.
- **toolsmith:** the 10 E-084 key-moved parents all keep the key under the fixed library; new opt-in Cartan-congruence check (`QM_CHECK_CARTAN=1`, off by default, ~15 ms per step), 2 tests. The old reduction fails it 10 of 10.

Promoted: E-100, E-101, E-102; glossary "cycle member".

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-020 agenda (recommend yes); (2) overnight: none (recommend). Decided for you (round 021 questions): agenda kept; no overnight. Round 023 is ordinary.

---

## Round 021 -- 2026-10-01 -- ordinary

Worked: scholar, skeptic, maverick. Referees: skeptic, experimentalist, theorist. All minor revision; I scoped the claims and accepted. Details in `rounds/021/`.

- **scholar:** the rewrite's Cartan entries are derived: (k,i) = dim coker g_i, (i,k) = dim ker psi_i (conditional on step 7 returning the whole ideal), so the congruence fails iff some dim ker g_i != 0. Rejecting parents are long-sided squares; strict A5 covers only 767 of 1 123 at n = 6. A derivation sketch plus checks (0 violations at n = 5..7), not a proof.
- **skeptic:** the one-placement-in-the-`444`-orbit split of E-091 is not special to words with a 4 (no-4 four-letter words give as many or more) and not enriched over a null; the gap 0/1 concentration beats a position null but describes the orbit. n = 12..14 only.
- **maverick:** at n = 8 a cord member appears within 3 steps iff a relation has >= 3 arrows (365 of 429, 0 mismatches); the 64 radical-square-zero LNAs have none to depth 5. An observed fit, unproved; negatives not shown cordless (depth-6 controls exist).

Promoted: E-097, E-098, E-099; glossary "long-sided square".

**Questions for you (the chair takes the recommended option if unanswered):** (1) keep the round-020 agenda (recommend yes); (2) overnight: none; the non-MONO L = 5 walk of the 365 cord LNAs was proposed but is not worth 4.5 CPU-hours (recommend no). Decided for you (round 020 questions): agenda approved; no overnight. Round 022 is ordinary.

---

## Round 020 -- 2026-10-01 -- conference

All six personas wrote position statements (no new work, nothing promoted). Details in `rounds/020/proceedings.md`. The merge of `origin/main` was refused (unrelated histories); carried on without it.

- **Shared weak point:** E-095 (every rejecting parent has one bad vertex with dim ker 1) is conditional on the parents being A5-shaped, which nobody has checked (experimentalist, scholar, and others).
- **Other weak claims:** lemma R keeps the right gap (theorist: formula and terminals, no proof); `k(33x) = 2x` has no mechanism (skeptic); E-087's "no monomial cord" rests on few walked LNAs (toolsmith, maverick).
- **Proposed agenda (round 020):** (1) A5-shape check, replay of the 10 E-084 parents, derive "(k,i) = dim coker"; (2) split words, the shuttle and the word-only gap, with a selectivity null; (3) which n = 8 LNAs give cords (non-MONO `--plan`); (4) a mechanism for `k(33x) = 2x`.

**Questions for you (the chair takes the recommended option if unanswered):** (1) approve this agenda in `STEERING.md` (recommend yes); (2) overnight: none yet (recommend). Decided for you (round 019 questions): no overnight; round-016 agenda approved. Round 021 is ordinary.

---

## Round 019 -- 2026-10-01 -- ordinary

Worked: toolsmith, theorist, experimentalist. Referees: scholar, skeptic, theorist. All minor revision; I applied the scope fixes and accepted. Details in `rounds/019/`.

- **toolsmith:** checkpoint/resume walker (`toolsmith_walk.py`, sliced = uninterrupted at n = 7). The n = 8 class 2 depth-8 walk now completes under the fixed library: 0 key-moved steps, 2 rejections; E-084's 10 key-moved steps match in count the steps that now keep the key (not replayed). Scoped to that class and depth.
- **theorist:** the in-orbit placement of each split word sits at right gap 0 or 1, n-independent, for 17 of 18 words (`3344` is an exception); lemma R alone never reaches `333@0` (0 of 45), it stops at boundary shapes inside the orbit. A relabelling of E-091 plus the new R-terminal result; why one placement is still read from tables.
- **experimentalist:** on the guarded walks every rejecting parent has exactly one bad vertex with dim ker 1 (907 + 143 distinct parents); tilting <=> dim ker 0 is E-093 restated. Conditional on the parents being A5-shaped (unchecked).

Promoted: E-094, E-095, E-096; H-015 status line.

**Questions for you (the chair takes the recommended option if unanswered):** (1) overnight: recommend none yet (replay the 10 E-084 parents first). (2) approve the round-016 agenda (recommend yes). Decided for you (round 018): no overnight; agenda unchanged. Round 020 is a conference.

---

## Round 018 -- 2026-10-01 -- ordinary

Worked: skeptic, scholar, maverick. Referees: theorist, experimentalist, skeptic. All three minor revision; I applied the scope fixes and accepted. Details in `rounds/018/`.

- **skeptic:** at n = 12..17, 6/9/12/15/18/18 four-letter words are *split* across the `444` orbit, each with exactly one placement in it (referee reproduced byte-identically); `3334`, `2455` have none. E-086's "0 partial" was vacuous for merged words.
- **scholar:** Cartan congruence (E-085) and `tiltingPlus` agree on every gate-admitted step tested (n = 5..7; 807 non-tilting steps), the discrepancy being minus the kernel dimension of one map: congruence is `tiltingPlus` read through the rewrite, not a second criterion. An observation, not a theorem; referee's capped n = 7 rerun gave other counts, same pattern.
- **maverick:** all 2376 n = 8 cord members (six LNAs, depth <= 5) carry a sum relation; no monomial cord at n = 4, 5, 8, so still no positive `MONO` control at n = 8. Scoped: informative for 6/14 (n = 5) and 1/5 (n = 4) LNAs, depth below E-087's.

Promoted: E-091, E-092, E-093; H-021, H-017, H-015 status lines; glossary (*split word*).

**Questions for you (the chair takes the recommended option if unanswered):** (1) overnight: recommend none. (2) approve the round-016 agenda unchanged (recommend yes). Decided for you (round 017): no overnight; agenda unchanged. Round 019 is ordinary.

---

## Round 017 -- 2026-10-01 -- ordinary

Worked: toolsmith, experimentalist, theorist. Referees: skeptic (two), scholar. All three minor or accept; I applied the scope fixes. Details in `rounds/017/`.

- **toolsmith:** `reduceAgainstPivots` is now a normal form in the library, with a unit test that fails on the old code (referee ran all touched and six other test files: pass). `MONO=1` at n = 8, L = 6 finds no monomial cord member for six walked LNAs; no positive control, so a scoped negative.
- **experimentalist:** the E-084 walks re-run under the fix give identical counts (n = 7, n = 8 classes 0-1, n = 9); the n = 8 class 2 walk ran out of 10 minutes before the 10 key-moved steps, so that loose end is neither confirmed nor refuted by the walk. Two rows (the "unpatched control", timing) struck: the control ran the patched library.
- **theorist:** `3334` and `2455` are one double mutation from `35` (lemma R, checked 250/250 and 119/119, not proved), hence in the `333@1` class and never the `444` orbit; the same lemma sends `444` to `34`. Referee reproduced and extended to n = 17.

Promoted: E-088, E-089, E-090; H-017, H-021 status lines; E-085 annotated; glossary (class label `J`, lemma R).

**Questions for you (the chair takes the recommended option if unanswered):** (1) overnight: recommend none until `scholar_walk.py` has a checkpoint. (2) approve the round-016 agenda (recommend yes). Decided for you (round 016): no overnight; agenda unchanged. Round 018 is ordinary.

---

## Round 016 -- 2026-10-01 -- conference

All six personas wrote position statements (no new work, nothing promoted). Details in `rounds/016/`.

- **Most promising questions:** experimentalist/skeptic/theorist: why `3334`, `2455` leave the `444` orbit, whether the row-set identity holds at n = 16, 17, and why the shift is by 2 and never 1 (E-080). Scholar: Cartan congruence vs `tiltingPlus` as one criterion. Toolsmith/maverick: does a monomial cord member exist at n = 8 (`MONO=1`, size first).
- **Weakest claims named:** `k(33x) = 2x` and lists A/B are fits without mechanism; E-084's n = 8 class 2 counts were made with the buggy `reduceAgainstPivots`; E-087 rests on 2 members.
- **Proposed agenda (round 016, in STATE.md):** (1) patch the library with a unit test and re-run E-084 class 2; (2) `MONO=1` at n = 8; (3) `3334`/`2455` and `k(33x)`; (4) the shift-by-2 mechanism.

**Please approve or change the agenda in `STEERING.md`**; until then rounds 017-019 work from it. Decided for you (round 015 questions): the library fix goes ahead in round 017; no overnight run. Round 017 is ordinary; 020 is the next conference.

---

## Round 015 -- 2026-10-01 -- ordinary

Worked: toolsmith (T6), theorist (T5), skeptic (T1/T2). Referees: maverick, skeptic, experimentalist. Three minor revisions; I applied the scope fixes in the entries and accepted all three. Details in `rounds/015/`.

- **theorist:** the n = 8 class 2 loose end of E-084 is a bug in the mutation rewrite, not in `tiltingPlus`: `reduceAgainstPivots` is not a normal form, so step 7 can drop a relation. The Cartan congruence fails on all 11 replayed rejecting parents, agreeing with `tiltingPlus`. Referee reproduced a congruent pair with different residues and the patched runs; the library is untouched, and E-084's counts were made with the buggy rewrite and are not re-run.
- **toolsmith:** n = 8 controls with cords (8-9 arrows) are found at depth 6 and not 5, at the n = 9 negatives' size (5.7e4-6.2e4 nodes). Caveat: every cord member has a sum relation, the n = 9 candidates are monomial, and none was found at n = 6, 7; 2 members only.
- **skeptic:** at n = 12..15 every merged word of the earlier scans lies in the `444` orbit's row set or outside it (0 partial); E-075's "20 of 25" is really 11 of 25.

Promoted: E-085, E-086, E-087; H-015, H-017 status lines; E-075, E-084 annotated.

**Questions for you (the chair takes the recommended option if unanswered):** (1) the `reduceAgainstPivots` fix: recommend toolsmith patches it with a unit test after the conference, and the E-084 n = 8 class 2 walk is re-run. (2) Overnight: recommend none. Decided for you (round 014): `isTilting` still not promoted; no overnight. Round 016 is a conference.

---

## Round 014 -- 2026-10-01 -- ordinary

Worked: experimentalist (T1/T3), maverick (T6), scholar (T5). Referees: theorist, toolsmith, skeptic. Three minor revisions; I applied the fixes and accepted all three. Details in `rounds/014/`.

- **scholar:** yes, a guarded walk from an LNA reaches a gate-admitted parent where `tiltingPlus` fails: n = 6 at distance 8, n = 7..9 at 5-7 (10 sampled classes). The Coxeter guard refuses every one, and 0 of about 1.3e6 guard-admitted steps fail. So the gate alone is unsound; the guarded walk is not. Referee reproduced the n = 6 and n = 7 cases and the new test; non-tilting rests on `tiltingPlus` alone and an n = 8 loose end (parallel arrows) is open. E-057's "none" was depth-limited and is annotated.
- **maverick:** an n = 8 control with relation-bearing sources finds its source in every run at its depth and none one step short, at 7e3-4e4 nodes at depth 6 (n = 9 negatives: 5e4-6e4). Members are the head of a sorted list and have no cords, so it controls the walk, not coverage of classes with cords.
- **experimentalist:** `5046`/`5056` at n = 17 have two closed orbits (122 673 / 54 266), now saved (confirms E-080). 4-letter words with a 4 at n = 12..15 join the single `444` orbit (identity by row membership at n = 13 only, by size elsewhere).

Promoted: E-082, E-083, E-084; status lines of H-015, H-017, H-021; E-057 annotated.

**Questions for you (the chair takes the recommended option if unanswered):** (1) `isTilting`: recommend still not promoting; Cartan check on the replayed parents and the n = 8 loose end first. (2) Overnight: recommend none; toolsmith builds an n = 8 control with cords first. Decided for you (round 013): no n = 17 overnight, no H-017 depth 7.

---

## Round 013 -- 2026-10-01 -- ordinary

Worked: theorist (T1/T3), toolsmith (T6), skeptic (T2/T4). Referees: skeptic, theorist, experimentalist. Three minor revisions; I applied the fixes and accepted all three. Details in `rounds/013/`.

- **theorist:** the two orbits behind lists A/B are each closed under offset shifts by 2, by an explicit staircase of `a - 1` double mutations (`3a@o -> 3a@(o+2)`, a = 5..12; referee extended it to 10..12); `4@0 -> 3@1` is one table rule, `5@0` has only 2 neighbours. No invariant explains why a shift by 1 is impossible: "parity class" names two orbits. `5046`/`5056` have two orbits of different sizes at n = 17 (the referee later reproduced the `5046` run; `5056` not re-run; no output saved).
- **toolsmith:** the n = 9 depth-6 negatives are walks of 5e4-6e4 nodes (4 of 16 measured), and an n = 7 depth-6 control finds its target 16 of 16 (0 of 4 at depth 5). The control members are cheap, and no n = 9 member is known, so it shows the search ran, not that it was enough. The depth-7 cost is 25-30 min per candidate.
- **skeptic:** counted by orbit, E-075's letter-4 contrast is one orbit per n (the `444` orbit), which holds both `34`-words and 4-no-`34` words, so letter 4 and collapse to `34` are not separable. Referee re-ran identically; E-075's "20 of 25" at n = 14 does not match the scan's 11 of 25 (unreconciled).

Promoted: E-079, E-080, E-081; status lines of H-021 and H-017; glossary "Staircase".

**Questions for you (the chair takes the recommended option if unanswered):** (1) n = 17 key-coarser lists overnight: recommend not yet; save and reproduce the `5046` n = 17 run first. (2) H-017 depth 7 overnight: recommend no; the toolsmith sizes an n = 8 control with a non-hereditary source first. Decided for you (round 012): no n = 17 overnight; round-012 agenda approved unchanged.

---

## Round 012 -- 2026-10-01 -- conference

Six position statements (`rounds/012/`); no new work, nothing promoted. Next round (013) is ordinary; 016 is the next conference.

- Weakest claim named by four personas: E-077's lists A/B are a fit with no mechanism and no real test at n >= 17 (and miss `5046 5056` at odd n). Others: the "walk-reachable only" rule for `isTilting` is policy, not a result (scholar); the H-017 negatives cannot be read without a node count and control (toolsmith, maverick). Maverick's Euler-form claim is unchecked (E-063 covers n = 8..11 only).
- **Proposed agenda (ranked):** 1. why A/B are parity classes: move sequences in P/Q, then `5046 5056` (theorist, experimentalist); 2. H-017 node count + n = 7 depth-6 control (toolsmith); 3. walk-reachable gate-admitted non-tilting mutation (toolsmith, scholar); 4. letter 4 versus the `34` collapse (theorist, skeptic).

**For you: approve or change this agenda in `STEERING.md`; until then rounds follow it.** Question: n = 17 key-coarser lists overnight (2+ h) -- recommend not yet. Decided for you (round 011): hand-built rejection does not count towards `isTilting`; no new H-017 overnight.

---

## Round 011 -- 2026-10-01 -- ordinary

Worked: experimentalist (T6), theorist (T1/T3), scholar (T5). Referees: toolsmith, skeptic, theorist. Three minor revisions; I applied the fixes and accepted all three. Details in `rounds/011/`.

- **experimentalist:** all 16 K = 4 candidates at n = 9 reach nothing at depth 6 (13 new shards, 308-567 s). Referee re-ran one shard identically. A bounded negative: no depth-6 control, no node count.
- **theorist:** the key-coarser lists A, B are the words alternating between two single-relation orbits (`3@2`/`5@0` even n, `3@3`/`6@0` odd n), n = 12..20. Referee: a fit at 12..16 and a consistency check at 17, 18 (no n >= 17 list was computed), misses `5046 5056` at odd n; the single-relation identification itself is new.
- **scholar:** a hand-built 5-vertex algebra with `abde = acde` is gate-admitted at `d` yet fails `tiltingPlus` and the Cartan congruence (n = 5..7), so E-066's step-7 shape is not special to n = 10. Hand-built, not LNA-reachable; CHZ still unread (arXiv blocked).

Promoted: E-076, E-077, E-078; status lines of H-017, H-021, H-015. Round 012 is a conference.

**Questions for you (chair takes the recommended option if unanswered):** (1) does a hand-built gate-admitted rejection count towards promoting `isTilting`: recommend no, it must come from a walk; (2) H-017: recommend no new overnight run, toolsmith adds a node count and depth-6 control first. Decided for you (round 010): no overnight; `orbitclass` at n = 17 not run.

---

## Round 010 -- 2026-10-01 -- ordinary

Worked: toolsmith (T6), skeptic (T2), experimentalist (T1/T3). Referees: theorist, experimentalist, skeptic. Two minor revisions, one accept; I applied the fixes and accepted all three. Details in `rounds/010/`.

- **toolsmith:** `toolsmith_verify.py` now takes `--list`, `--cand` and `--budget-hours`; E-072's depth-5 negative reproduces and a second n = 9 candidate reaches nothing at depth 6 (434 s, fits one shard). Referee: a reproduction plus tooling; wording fixes only.
- **skeptic:** the neighbour-aware null for "a = 4 special": only `444` merges among `aaa` at n = 12..15, but any word with a 4 merges far more often (55/100 against 3/121 with neither 4 nor 2), so the claim is a letter-4 effect and does not single out the collapse to `34`. Referee reproduced it; the orbit-collapsed count was not done.
- **experimentalist:** the key-coarser cores are the same 9 words at n = 12, 14, 16 and the same 10 at n = 13, 15; `348`/`349` size-20300 pairs at 16 are the `4056` orbit and its mirror. Referee re-ran it byte-identically.

Promoted: E-073, E-074, E-075; status lines of H-021, H-017.

**Questions for you (chair takes the recommended option if unanswered):** (1) remaining H-017 depth-6 candidates: recommend chair-slot shards in round 011, no overnight; (2) `orbitclass` at n = 17 (2+ h): recommend not until the theorist explains why those words are parity classes. Decided for you (round 009): shard tooling added, no overnight; `aax` criterion kept out of H-021.

---

## Round 009 -- 2026-09-30 -- ordinary

Worked: experimentalist (T1/T3), theorist (T2/T4), maverick (T6). Referees: skeptic, experimentalist, toolsmith. All minor revision; I applied the fixes and accepted all three. Details in `rounds/009/`.

- **experimentalist:** the equal-size singleton pairs of `344 348 349` (n = 15..17) and `4046` (n = 14..16) are one orbit plus its mirror, orbit+mirror = key in 12/12 cells; the key-coarser cores at n = 12, 13 are not the 7 of E-059. Referee: mostly already in E-064 (only n = 17 new); the "pair at even n, mirror-join at odd n" summary contradicted the table and was withdrawn.
- **theorist:** `aax` drift families close with `k = 2x + 3 - a` for a = 3, 5, 6 (7 by the referee) and `44x` does not, because the seed `444` collapses to `34`, which reaches the slider `44`. Predicted before running for `55x`, `66x`; referee's extra runs agree. Computed criterion, not proved (one positive datum).
- **maverick:** L = 5 control passes (42/42 at n = 6; 8/8 near-trivial LNAs at n = 7); four n = 9 candidates reach nothing at depth 5. Referee corrected the sizing: depth 6 for 16 candidates about 2 h, one candidate per 10-minute shard.

Promoted: E-070, E-071, E-072; status lines of H-021, H-017.

**Questions for you (chair takes the recommended option if unanswered):** (1) H-017 depth 6 for all 16 n = 9 candidates: recommend toolsmith adds a candidate-index argument first, no overnight run (only depth 7 is overnight, already in Menu 4); (2) write the `aax` criterion into H-021's text: recommend not until a skeptic's null. No questions were outstanding from round 008.

---

## Round 008 -- 2026-09-30 -- conference

Six position statements (haiku), no new work, nothing promoted. Details in `rounds/008/`. Their factual sub-claims are unchecked.

- **Converging:** experimentalist, theorist and toolsmith all want the mirror-join / parity-class check of the equal-size singleton pairs (`344`, `348`, `349`, `4046`; the 9-10 key-coarser cores of E-064 vs the 7 of E-059). Skeptic's weakest claim: the centre formula was fitted at n = 13 and "confirmed" on pre-selected cores, so it needs a fresh random sample.
- **Proposed agenda (ranked):** (1) mirror join + parity classes; (2) fresh-sample centre test; (3) why the `33x` drift closes and `44x` does not; (4) H-017 L = 5 control, sizing before any overnight; (5) H-015 one-map identity with a small commutative instance. Maverick's Coxeter-spectrum question is parked.

**Please approve or change the proposed agenda in `STEERING.md`; until you do, ordinary rounds work from it.** No new questions. Decided for you (round 007 questions): no overnight depth 5-6 rerun of H-017 yet; no `34x` at n = 18, 19 until the mirror join is done.

---

## Round 007 -- 2026-09-30 -- ordinary

Worked: experimentalist (T2), skeptic (T2/T6), maverick (T6). Referees: skeptic, theorist, experimentalist. All minor revision; I applied the referees' wording fixes and accepted all three. Details in `rounds/007/`.

- **experimentalist:** `34x` offsets pair `o <-> hi - o` (`k = x + 3`, not `2x`) at n = 14..17 for x = 4, 5, 7, 8, 9 (20/20 cells); `346` one orbit; `45x` no reflection; `4046` is a (size-paired) reflection. E-060's `4046@13` line did not reproduce. Caveat: `k = x + 3` is just `s = hi`, and singleton pairing is by size only.
- **skeptic:** null test. 39 of 109 n = 13 fits are vacuous and informative ones are mostly chance-level (73%), so E-061's "45/62, 10/13" counts are padded; the interior-core centre formula (13/13 vs 4.1 expected) and the E-060 cores at n = 15/16 survive (joint p about 1e-3..1e-4 after the referee's correction).
- **maverick:** the H-017 search passes a positive control (273/273 round trips at n = 7 on 91 of 132 LNAs; 0 at one level too shallow), so the round-004 depth-4 negative only excludes members within 4 steps.

Promoted: E-067, E-068, E-069; status lines of H-021 and H-017. Round 008 is a conference.

**Questions for you (chair takes the recommended option if unanswered):** (1) overnight depth 5-6 rerun of the 16 H-017 candidates at n = 9: recommend not yet, size it first; (2) `34x` at n = 18, 19 overnight: recommend no, mirror-join check first. Decided for you (round 006 questions): agenda approved unchanged; no overnight run.

---

## Round 006 -- 2026-09-30 -- ordinary

Worked: toolsmith (T3/T8), theorist (T2/T4), scholar (T5). Referees: skeptic (x2), experimentalist. Toolsmith accepted; theorist and scholar minor revision, which I applied myself. Details in `rounds/006/`.

- **toolsmith:** over all 139 placed cores at n = 10, 12..16, orbit-plus-mirror classes refine the key classes (never finer, never incomparable; equal in 129-132); the exceptions are the parity-merged cores. The 20300 pairs at n = 16 (`4056`, `46`, `3355`, `3445`) are one orbit and its mirror, so E-058/E-062's "unmerged middle pair" is resolved. Referee reproduced it.
- **theorist:** `k(33x) = 2x`, `d = x - 3` come from a drift `33x@o -> 33(x-1)@(o+1)` of the double mutation (label `x + o` fixed) plus the self-dual seed `333`: 135/135 at n = 14..16, plus the referee's n = 17. Lower bound derived, upper bound computed only; the same argument is false for `44x`.
- **scholar:** E-032 step 7 is rejected at an explicit commutativity element; mostly known, and the literature is not an independent test of the code. CHZ "monomial only" caveat is UNVERIFIED (arXiv blocked).

Promoted: E-064, E-065, E-066; H-021 status line; GLOSSARY (two terms); UNVERIFIED flag on `literature/2509.12983`. No open questions were left by round 005.

**Questions for you (chair takes the recommended option if unanswered):** (1) approve the round-005 agenda (recommend yes, unchanged; T3 is now mostly answered); (2) no overnight run proposed yet, toolsmith to size `--max-word 5` at n = 14 first.

---

## Round 005 -- 2026-09-30 -- conference

Six position statements (haiku), no new work, nothing promoted. Details in `rounds/005/`.

- **Convergence:** experimentalist, skeptic and toolsmith all want the same test: orbit-plus-mirror vs key classes over the `--max-word 4` catalogue at n = 14..16 (T3/T8). Theorist wants a mechanism for `k(33x) = 2x`; maverick a positive control for the H-017 search; scholar to decide E-032 step 7 from the literature (Aihara-Iyama, CHZ).
- **Weakest claims named:** parity (19 of 139 cores tested), `k = 2x` (description, not mechanism), H-015 support (n <= 7), "Coxeter polynomial cannot see (cords, relations)" (n = 9 only). No disagreement between personas.
- **Agenda proposed** (ranked, in `STATE.md` as proposed, round 005): 1 orbit-vs-key at 14..16; 2 `k = 2x` from the rule table (T4); 3 H-017 positive control; 4 H-015 step 7 by reading; 5 parity across the other 127 cores waits on the overnight censuses.
- **Decided for you (round 004 questions):** approved the n = 9 depth-7 H-017 run only (in `OVERNIGHT.md` Menu 4; I added `--budget-hours` to `maverick_reached.py`, tested, `test_overnight_doc` passes); round 005 kept a conference. Overturn in `STEERING.md`.

**Please approve or change the proposed agenda in `STEERING.md`; until then rounds 006+ work from it.**

---

## Round 004 -- 2026-09-30 -- ordinary

Worked: theorist (T2), experimentalist (T1), maverick (T6, first time). Referees: skeptic (x2), scholar. All minor revision; I accepted all three after applying the referees' wording and one check myself. Details in `rounds/004/`.

- **theorist:** the interior/end-touch split is now a committed column and explains none of the 7 failures (the slide never picks a centre, and the orbits take the smaller or larger consistent one at random). For `33x`, `k = 2x` and `d = x - 3` at n = 13..17 for x = 3..6; I added `337` (holds). A description, not a mechanism.
- **experimentalist:** the 12 cores of E-060 keep `k` and `d` at n = 15, 16, and pair at 13, 14, 15; at 16 three lose the fit through an unmerged equal-size middle pair (as `4056` in E-058). So parity is not the whole story.
- **maverick:** H-017 survives to depth 6 at n = 9, but neither the Coxeter polynomial nor the Euler form predicts (cords, relations). The Euler-form signature does separate "outside every quipu class" for n = 8..11 (small, new to the record).

Promoted: E-061, E-062, E-063; status lines of H-021 and H-017. Step 0.5: decided round 003's question (yes) and added the n = 12 and n = 14 censuses to `OVERNIGHT.md`. Round 004 should have been a conference by the config; round 005 will be.

**Questions for you (chair takes the recommended option if unanswered):** (1) H-017 overnight: approve the n = 9 depth-7 run, not the n = 10 run until the search has a positive control? Recommend yes to n = 9 only. (2) Keep round 005 a conference? Recommend yes.

---

## Round 003 -- 2026-09-30 -- ordinary

Worked: theorist (T1/T2), experimentalist (T1), toolsmith (T8). Referees: skeptic (x2), theorist. All minor revision; I accepted all three after applying the referees' wording and test points myself. Details in `rounds/003/`.

- **experimentalist:** the 7 "mirror without reflection" cores of E-056 pair by a reflection at n = 14, and at every even n 12..18; they fail at n = 13, 15, 17. So that defect is an odd-n effect **for these 7** (chosen because they show it at 13); nothing known yet for the other 132.
- **theorist:** H-021 restated without the mirror clause, covering only cores that pair. For cores whose outside block is interior the shortfall is exactly the H-020 tail minus head (13/13 at n = 13, 3/3 at 14); with the block at an end 17/21. Near-tautological, F-053 has it for `45`; the interior/end split is not yet a committed column.
- **toolsmith:** `batch.py orbits N` (resumable E-052 orbit report), tests pin `45` and `344` at 13.

Promoted: E-059, E-060; H-021 status line (still OPEN); two glossary terms. `isTilting` still not promoted (round 002 q1). Merging `main` was a no-op. No unanswered questions to decide.

**Question for you (chair takes the recommended option if unanswered):** add the full n = 12 and n = 14 censuses (139 cores, about 90 min each, 4 procs) to `OVERNIGHT.md`, to see whether parity governs the other cores? Recommend yes.

---

## Round 002 -- 2026-09-30 -- ordinary

Revisions by experimentalist (T1), skeptic (T3), scholar (T5). Referees: skeptic (x2), theorist. All accepted (skeptic's with three wording points, which I applied). Details in `rounds/002/`.

- **experimentalist:** H-021's mirror clause fails on every reading at n = 13. Strict reading: false for 108 of 109 pairing cores (by construction); 7 cores (`344 366 4044 4403 4404 4405 4605`) hold a strict mirror with no reflection. `3346` and `4056` hold no strict mirror, so H-021's text on `3346` is right under that reading.
- **skeptic:** first counterexample to "orbit = key class": `4056` at n = 16, offsets `{1,2}` are two mirror-image orbits though the key pairs them (derived equivalent via F-026). Also corrected the n = 24 scan: 309 / 155 / 121, not "7 / 123".
- **scholar:** on non-monomial parents (n = 5, 6, 7) Ladkani 2.3(c) rejects exactly the gate-refused vertices; still no gate-admitted rejection except E-032 step 7. Negative control reproduced.

Promoted: E-056, E-057, E-058; E-054 corrected in place; H-021 status line updated (still OPEN); two glossary terms. Merging `main` was a no-op.

**Questions for you (chair takes the recommended option if unanswered):** (1) Promote `isTilting` to the library as a cross-check? Recommend not yet. (2) Restate H-021 without the mirror clause next round (theorist), with the 7 survivors run at n = 14 (experimentalist)? Recommend yes.

---

## Round 001 -- 2026-09-29 -- ordinary

Worked: experimentalist (T1), skeptic (T3), scholar (T5). Referees: skeptic (x2), theorist. All three came back **minor revision**; all three are due for revision next round. Details in `rounds/001/`.

- **experimentalist:** all 139 single-cluster cores of `--max-word 4` close at n = 13. H-021's "exactly when" fails under the literal reading: 129 of 139 cores hold a mirror, so it does not discriminate, and 20 hold one with no reflection (including `3346`, which holds its own mirror). A stricter reading of "mirror" was not run.
- **skeptic:** `3346` and `4056` are not artefacts. The Coxeter key differs at every offset of `3346` for n = 12..40, which proves distinct orbits. `4056` is F-053's reflection with an onset (orbit-verified to n = 15). Two H-020 failures fit the same onset pattern. Referee: pairing beyond 15 is key-only and "cause" is overstated.
- **scholar:** all 61,718 gate-admitted steps at n = 6, 7 pass Ladkani's exact tilting criterion (arXiv:1001.4765, Prop. 2.3(c), unused in `research/` until now), and the E-032 ALARM step fails it. This is not evidence for H-015: the guard is inert at those sizes.

Promoted: E-053, E-054, E-055 (run records with the referees' caveats); H-021 status line annotated, still OPEN; three glossary terms. No F entry yet, pending the revisions. Merging `main` was a no-op.

**Questions for you:** (1) Restate H-021 now (overhang or onset), or settle the strict mirror reading first? (2) Write into `OVERNIGHT.md`: the n = 14 census (about 90 min on 4 processes) and the Ladkani audit at n = 9/10? (3) Promote `tiltingPlus` to the library as an exact gate? I would wait for a second non-monomial negative case.

---

