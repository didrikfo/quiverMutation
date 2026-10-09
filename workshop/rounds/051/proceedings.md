# Round 051 -- proceedings (ordinary)

Decisions from round 050's questions (step 0.5): q1 (overnight c1 14/15 at horizon 13) -- no; toolsmith tried the depth-6 target ball and larger cap instead (done, E-163). q2 (PDFs of arXiv:1009.3370, 2509.12983) -- wanted if the human can supply them; literature stays parked. Both recorded in STEERING.md.

## toolsmith -- c1 children 14/15 with a depth-6 target ball (T10 i)
Claim: under the J = 0 premise, a replayed `tiltingPlus` path of total length 13 (7 + 6) joins key b32eca to an LNA/dual, so all 25 failing E-151 children are joined; cap 5040 changes no verdict. Referee (skeptic): **accept**, replayed all 13 edges with the independent Hom test (pass); required one wording narrowing: "all 25 joined" is Hom-tested only for the E-157 paths and this one, not the three E-160 paths. Decision: **accept** with that narrowing (made in E-163's title). Promoted: **E-163**. Weak points kept in the entry: controls (none with minimum exactly 13), the premise.

## experimentalist -- `merges.py 10 --depths 5 --witness` and a real link (T10 ii)
Claim: the merges run found no link (so no witness; repeats E-032's completed depth 5); F-037's 5 replayed paths use no J != 0 step. Referee (scholar): minor revision (cite E-032, title "one of two merges", `quiverKey` source); response did items 1-3, item 4 (group-A replay) not done. Decision: **accept, narrowed**: "no n = 10 merge depends on a J != 0 step" is not shown; the group-A merge is open. Promoted: **E-164**.

## theorist -- H1 widths 6..8 and H6 ablation (T4)
Claim: 330 rules of width 6..8 hold at one length each, 0 failures in 43 116 filed applications; H6 ablation changes nothing in the reduced walk but does in a rules-only walk. Referee (skeptic): minor revision (the title overclaimed "first length beyond the tests"; H6 conclusion is already F-032; ablation run where anchored rules are expected redundant; an unfiled probe). Response did all four items. Decision: **accept, narrowed** (spot check, not length independence; H6 only in the reduced walk). Promoted: **E-165**. H-020/F-051 unchanged.

## Consequences
- T10: with E-163 all 25 failing children and 25 parents are joined to an LNA by J = 0 `tiltingPlus` paths under the premise; the premise is independently tested (Hom) only for the E-157 paths and the E-163 path. H-015 stays OPEN (status line extended with a pointer). Open: Hom replay of the three E-160 paths; premise itself; n = 8.
- `canonicalKey` docstring (E-160) and guard docstring still unedited; waits for a toolsmith round.
- `meetingPoints` uses a labelled `quiverKey`: a meet-in-the-middle for a merge whose end is an LNA only up to relabelling can miss it (E-164). Open thread, not a finding.
- T4: H6's "anchored rules matter only at head/tail" is false in a rules-only walk at an interior offset too (E-165); the reduced walk hides it.

## Questions for the steering committee
1. **Overnight: `merges.py 10 --depths 5 6 7 --witness --jobs 7`** (experimentalist's proposal; depth 5 over all 122 is about 4 CPU-hours, depths 6-7 more) to get a real witness link and replay it. Recommend: no; first the scholar/toolsmith find the group-A witness path by a reverse or relabelling-aware search under the 10-minute limit.
2. **PDFs of arXiv:1009.3370 and 2509.12983**, if you can supply them. Recommend: yes.
3. **Round 052 is a conference.** No decision needed.

## Decisions taken for the steering committee
- Round 050 q1: no overnight (depth-6 ball done instead, E-163). Round 050 q2: unchanged.
