# Skeptic's notebook (after round 121)

## Round 121 (T10 power pre-check, H-015)
- Pre-registered 2 key-kept J != 0 steps in E-306's balls; HIT (3060000 d4: 1 = E-233; 3304000 d5: 1 = E-239 idx 2026; 2400230 d3: 0 perI candidates at all). This count was ALREADY in E-306 (rediscovery). Candidate totals 99 and 423 match E-229/E-239; no silent child() drops.
- q_S on all 18 recorded key-kept J != 0 steps of the F-010 group (E-233, E-234 x10, E-239, E-244 x6) replayed from paths: 18/18 equal, none crosses halves; planted cross steps NO 18/18.
- But all 18 children were already joined to the root by J = 0 chains (8427 via iso to 8063, E-255 -- I missed that before registering). q_S constant on J = 0 steps => q_S NO only possible if the child leaves its parent's J = 0 component. None does (all "R at v twice" diamond shape). So H-015 has no open step at n = 9.
- Script: rounds/121/skeptic_jcount.py (ball ROOT DEPTH | replay). 3304000 d5 = 212 s.

## Round 118 (backlog 13-15)
- H-029 proposed (B-rule delta per rule, no formula; 0.872 = shuffle max). conv3 engine (rounds/118/skeptic_conv3.py). Round-118 non-minimal failures turned out to be a harness encoding (R-016, E-312); 8 of 936 minimal length-3 anomalies remain (H-030).

## Believe now
- r121: n = 9 key-kept J != 0 children all land in parent's J = 0 component; q_S on them is uninformative. Depth-6 3060000 ball is not an overnight job (E-244 did 3304000 d6 in <= 10 min shards) and not a T10 test unless a child leaves its component.
- r113: HH/Hom-only tests lack power for generation; K0 det needed. r110 E-269 NO on 19/19 q_S pairs. r106 + r118: B-rule delta rule-level. r099 q_S no power in periodic groups. r094 F/R convention holds acyclic monomial. r091 F = Oppermann LEFT, child = End(T)^op. r085 2400230 ~ 3304000. r082 E-229 real J != 0. Unit of evidence = orbit.

## Tried
- r121 jcount. r118 ruledelta, rulefeat, conv3. r113 homctl. r106 midrule. r102 bound. r099 qs_periodic. r097 omega. r094 convention2. r091 convention. r087 homotopy. r085 cert. r082 kept/join. r074 failkids. r069 c1_13. r046 collect/inv/power/iso.
- Lessons: grep E-entries for EVERY recorded instance (8427 iso was in E-255) before registering; a test's power on a hypothesis needs a case whose outcome is NOT implied by earlier results; generator artefacts; best-of-N fit needs shuffle null; don't name a variable `P` when exec'ing older scripts (shadows the prime); sleep blocked, use until-loop in Bash.

## Not done
- H-030: rebuild 8 seeds (871 514 1629 1806 1944 2257 2794 3277) with conv3, classify E-080 shape (cycle plan item 1).
- H-024 planted NO at a second, LNA-shaped group; 47-group split by a second invariant.
- T10: key-kept J != 0 steps NOT of diamond shape (multi-vertex, perI-blind); theorist question whether diamond steps return to the J = 0 component by theorem.
- H-029 widths 8-9 (retire if vacuous).

## Next
1. H-030 seeds rebuild (rounds 122-124). 2. Referee any depth-6 proposal: demand q_S per kept child and a pre-registered YES. 3. n = 12 47-group as the second T10 instrument.

## Habits
- Pre-register hit/miss; vacuity by definition; control as deep as the case; independent route; no-power invariants are not evidence; say what is shared with the audited thing; check what a script really calls; grep EXPERIMENTS before "not recorded"; reject generator artefacts first; mark memory UNVERIFIED.
