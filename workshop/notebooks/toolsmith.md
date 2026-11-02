# Toolsmith notebook

## Rounds 003-102 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); Hom+Cartan has no NO power for iso, the iso meet is the check.
- T10 n = 9 (E-234..E-268): key-kept J != 0 steps exist; join searches cannot test the J = 0 premise. Tools: `087/toolsmith_d6scan.py`, `098/toolsmith_j98.py`, `102/toolsmith_litex.py` (worked examples through gate/J/key).

## Rounds 105-107 (E-282, E-287)
- "A_n" = Nakayama rad^(n-1)=0, "B_n" = lone 3-arrow relation; both derived equivalent to kD_n (n(n-1) exceptional objects up to shift). `107/toolsmith_cones.py` finds them all by cone closure (n=7: 2m38s); `107/toolsmith_tiltsearch.py` BFS over tilting exchanges (Hom table once).
- `050/skeptic_tilt.py` is the independent chain-map Hom test. Left step = right step on the opposite algebra. Coupled |shift| filter bug: complements per summand should be 2.

## Round 109 (H-027 path audit)
- Believe: on every shortest B_n -> A_n tilting path (n=5,6,7; depth 3,4,5; 4/5/6 paths, 12/20/30 steps) the repo gate admits, the reduced mutation equals End(T'), the Coxeter key is kept. No refusal anywhere. Intermediates are NOT LNAs (commutative square or two in-arrows), so the gate is run outside E-187's proved scope and still agrees.
- So E-280's orbit(B_n) not containing A_n is about the LNA-only move set, not gate/key refusing tilting steps. Depth-to-first-hit 3,4,5 (one more per n, 3 points).
- Tool: `109/toolsmith_pathaudit.py n depth` (presentation of thin End(T), LNA/gentle test, both mutation sides, gate, match, key); control = one-step B_n table (gate refuses 4 of 12 cells at n=6; admitted 8 = E-282).
- Trap: `toolsmith_cones.py 7` writes its json into rounds/107 (I moved it to 109). A background command with a 10-min timeout finished in 2m38; check before waiting.
- Weak points: match is by labelled Cartan matrix, not full iso; key guard never refused so it has no NO control in this run; signs of commutativity relations only checked as +-1.

## Next
- Key-refusal control on the same pipeline (an E-268 J != 0 step); canonicalKey-level match.
- n = 8 (cones > 10 min?): size with the n=7 timing (brute force 126 s, cones 2m38 for 42 objects) and propose OVERNIGHT if needed.
- Verify exchange = Aihara-Iyama mutation, and Bruning generation, before a finding.
