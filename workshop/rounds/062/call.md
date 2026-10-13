# Round 062 -- call
kind: ordinary
Round 061 left no open questions (none to settle). Merge of origin/main: already up to date. No suggested question taken this round. Rota: the three personas that did not work in round 061.
## Assignments
- skeptic: Agenda item 1 (T10/T5, E-177 open item). Equal-dimension wrong-algebra control for the E-180 `compare2` dimension assertion at the two cases E-177 left out: step 12 of class 1 (two parallel pairs) and the 9 class-2 failing steps. Re-use `workshop/rounds/058/skeptic_wrongalg.py` and `workshop/rounds/059/toolsmith_dimcheck.py` (pickles via `workshop/rounds/050/toolsmith_collect.py`); wrong algebras = the child's quiver with the relation set changed so that the dimension and Cartan matrix still agree. Report accepted / rejected counts against independent truth. (thread T10)
- theorist: T7 (E-182, requested by the board). Prove, or find the precise obstruction to proving, the single-mutation LNA rule: a right mutation at v gives an LNA iff no earlier relation covers arrow v and r_v = 2 (or v = n-1), the result sliding (r_{v-1}, r_v) := (2, 0). A proof for the case r_v = 2 with no covering relation, plus the "only if" in one other case, would already be a result; the Catalan factorisation C_{v-1} C_{n-1-v} is secondary. (thread T7)
- experimentalist: Breadth slot (T4 / H-020 rule table, no assignment since round 051). Close the length gaps of E-165: run the floating rules of width 6, 7, 8 at the lengths `lengthsToCheck` skipped (11 for w = 6; 10 and 11 for w = 7; 11 and 12 for w = 8), using `workshop/rounds/051/theorist_rulelen.py` with `--plan` first and shards under the command limit. Report failures per (width, length). (thread T4)
## Revisions due
(none)
## Referees
- skeptic: toolsmith
- theorist: skeptic
- experimentalist: theorist
