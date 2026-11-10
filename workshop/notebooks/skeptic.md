# Skeptic's notebook (after round 123)

## Round 123 (T10 replay of recorded merges, H-015 / H-032)
- Two kinds of "merge" on record, and the call conflated them:
  (a) n = 9 F-010 JOINS built around an audited J != 0 step (E-233/234/239/244/254/255): root -path-> P -STEP-> X -path-> M ~iso~ M' <-path- root.
      Mixed F/R; only the root halves are mostly F-only (16 of 19). Every one contains exactly one J != 0 step BY DESIGN (the audited step).
      Replayed 19 (13 printed + 5 E-234 terminal-only regenerated, totals 8 8 7 7 7 reproduce + E-233 X~P) and 8427 via iso to 8063 (True):
      258 steps, 19 J != 0 (the audited ones), 239 J = 0 / tiltingPlus / key kept / gate; INV (SNF, signature, q_S) constant at EVERY node.
      SNF(S) and signature are equal for the two F-010 halves, so only q_S has power at n = 9 (planted cross merges: meet iso False, q_S differ, 2 of 2).
  (b) The library's own LNA-to-LNA merges (no J filter): F-037 (5 paths, 35 edges), F-054 meet + mixed certificate, C(2,4,5) atlas meet,
      F-041 meetingPoints n = 8 (10), n = 9 (52; 13 none), n = 10 (206; 98 none) = F-041 counts; 276 merges. 0 J != 0 steps, 0 INV changes.
- Control: E-318 walk flagged (step 4 F3 J != 0, NOT-tp, SNF tail (1,2,4)->(1,2,2)).
- Scripts: rounds/123/skeptic_mergereplay.py (records|extra|f037|control|mkpkl), skeptic_f041.py N [lo hi], skeptic_jtilt.py (copy of E-317 classifier).
- Mostly confirmation: join steps were already Hom-verified in their records (skeptic_join replay), F-037 edges tiltingPlus (E-164), F-041 n = 8 (E-156).
  New: perI + INV at every node; F-041 n = 9/10 witness edges; C(2,4,5) meet path replayed (E-164 said unreplayed).

## Believe now
- No recorded merge (of either kind) uses a J != 0 step except the audited steps of (a), and those stay in q_S class. E-318's walk (F3 R8 F3, then J != 0 F3) is
  outside every recorded path; no recorded LNA merge passes a J != 0 step -- but that is FORCED (r123 response, toolsmith): the depth-3 key-guarded
  balls of all F-041 pairs n = 8, 9, 10 have 0 J != 0 walk edges (2 396 / 18 518 / 117 642 distinct; skeptic_ballJ.py). Exactly 1 meeting per pair.
  Lesson: before reading a 0 count on paths, count the pool the paths are drawn from. C(2,4,5) merge was already E-169.
- r121: n = 9 key-kept J != 0 children land in parent's J = 0 component. r113: HH/Hom-only tests lack power. r099 q_S no power in periodic groups.
  r094 F/R convention; r091 F = Oppermann LEFT. Unit of evidence = orbit.

## Not done
- n = 7 "25 of 25": E-160's 3 paths still not Hom-tested (needs E-151 c1 rebuild > 10 min); INV has no power at n = 7 anyway.
- E-032 n = 10/11 overnight merges and n = 11 merges.py runs: no paths recorded; would need --witness reruns.
- H-030 seeds rebuild (871 514 1629 1806 1944 2257 2794 3277); n = 12 47-group second instrument; H-029 widths 8-9.

## Next
1. Library default: refuse J != 0 steps in mutationSearchDepthFirst (E-151 cost 8-11%) -- prunes 0 edges in F-041 balls; regression = F-041 counts unchanged AND E-318 walk refused.
2. H-030 seeds. 3. Referee any n = 12 T10 run: pre-registered YES, q_S/SNF per kept child.

## Habits
- Pre-register hit/miss; grep the step class, not my E-list; check what a "merge" is in each record before counting; exec'd scripts overwrite sys.argv/ARGV
  (save mine as MYARGV); control as deep as the case; no-power invariants (SNF at n = 9 F-010) are not evidence; mark memory UNVERIFIED.
