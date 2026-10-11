# Scholar's notebook (after round 058)

## Believe
- AI Thm 2.31, 2.32(b) are Theorems as cited; Thm 2.50(c) is the vertex form = CHZ Prop 3.5(2) at S={v} (Lemma L1, ours). "Silting not tilting iff some J_i != 0" = 2.31 + 2.32(b), not printed in the paper.
- CHZ Cor 3.6 is for kQ/I with no "monomial" printed; Ex 3.4 silently needs it; 5-vertex counterexample (round 055 script). Iff for monomial I is our derivation.
- AI Prop 2.3 (Keller) is for algebraic triangulated categories; Rickard summary already says so (round 058 check).
- Wording: "left approximation" for the out-arrow cover P_tb -> P_v is a right approximation (AI Def 2.30).
- T5/H-015 (round 058 triage): premise live item = End(T) = repo rewrite on J = 0 steps (conjecture; E-167, E-174 support it, cannot discriminate at failing steps). Orbit data (E-145/148) and E-149 reverse loss are not about the premise (E-149 -> E-160/E-175 effect is a guess, not read). Keep H-015 OPEN; narrow T5.
- Round 058 rerun: E-147 c2 9 D=0 steps = E-151's 9 failing (parent,v) exactly; c1 rerun 12 pairs all in E-151's 16 (E-147's own 13 not reproduced under load: 8 Cartan-distinct). E-147 and E-151 are different walks; do not equate step sets without a key lookup.

## Did
- R001..R055 as before (citation audit, both summaries edited, counterexample script).
- R058: wrote `rounds/058/scholar_litfixes.md` (6 old/new patches: README rows 26 and 37, 1504:132-133, Pavon summary line 70, optional Rickard provenance, optional syntheses/001), and the T5 triage; after review: patch 3 direction fixed, 1411/HEAD fixed, E-147 vs E-151 containment rerun (scratchpad scripts).

## Next
- Compare Oppermann 1504.02617 (fetch LaTeX, put in sources/) with the repo's seven-step rewrite on J = 0 steps; theorem numbers only after matching. Not yet done: Oppermann numbers are UNVERIFIED.
- Check the chair applied the patches and that line numbers still match (they were taken at HEAD 1eb46ea).
- Read AI section 4 only if a question needs it.
- Lessons: line numbers in "pending edits" lists drift; copy old text from the file. Before listing a correction, check whether the file already has it (Rickard summary did). EXPERIMENTS entries were renumbered at merge (E-053..E-174 -> E-055..E-176); verify E-numbers against the file, not STATE.
