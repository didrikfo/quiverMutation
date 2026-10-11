# E-136's meeting path is not a derived-equivalence path: its first hit-side step is gate-admitted, key-preserving and fails Ladkani 2.3(c), and with non-tilting steps removed no hit meets an LNA

author: skeptic · round: 039 · kind: negative (to E-136's reading, and a boundary of H-015)
thread: T5 · bears on: H-015, E-086, E-095, E-134, E-136, E-131

## Claim

E-136 gave 25 algebras shared by the key-preserving forward BFS from the two class-0 LNAs and from each of E-134's 16 hits (n = 6, d >= 3, J != 0, LNA key), and read it as "all 16 lie in the LNA derived class if the Coxeter guard holds". Replayed step by step: (1) an LNA -> M path (M = {1->5->2 = 1->6->2 commutative square, one relation}, 11 steps) passes all four tests at every step (gate, `tiltingPlus`, Cartan congruence C' = R C R^T, recomputed child equals stored child); (2) the hit -> M path (2 steps, identical for hits 0, 5, 11) has a FIRST step at v = 2 that is gate-admitted AND key-preserving but fails `tiltingPlus` (kernel dim 2 at i = 1) and fails the Cartan congruence; only its second step is tilting. So the meeting runs through a step the guard admits and Ladkani 2.3(c) rejects. With every non-`tiltingPlus` step removed from both walks, none of the 16 hits meets the LNA side (0 shared, hit side 323-501 algebras in 12 s, LNA side 27 518 in 150 s). Consequence: E-136's "derived equivalent to an LNA" is NOT supported by this meeting; H-015's guard is not sufficient at these parents (it was tested only on parents reached by guarded walks from LNAs, where E-086 found 0 of about 1.3e6 failures). Not claimed: that the 16 hits are outside the LNA derived class (a bounded miss again, 12 s per hit); nor that H-015 fails on LNA-reachable parents.

## Evidence

The hits are exactly the algebras with J_i != 0 at their mutation vertex v = 2, which by AI 2.32(b)/E-123 is what makes the mutation silting-not-tilting; the key-preserving BFS of E-136 expanded them anyway (it applies gate + key only). So the first step was non-tilting by construction of the hit set; the replay makes it explicit (E-130/E-123 already say J != 0 iff not tilting, so the algebra is old; the application to E-136's meeting is new).

Replay, hit 0 (`skeptic_n6path.txt`; hits 5 and 11 identical):
| path | step | v | gate | tiltingPlus | Cartan congruence | key kept | ker (i: dim) |
|---|---|---|---|---|---|---|---|
| LNA -> M | 1-11 | 1,2,2,2,3,4,1,3,5,1,5 | all True | all True | all True | all True | all {} |
| hit -> M | 1 | 2 | True | **False** | **False** | True | {1: 2} |
| hit -> M | 2 | 1 | True | True | True | True | {} |
Shared keys printed (first three by depth sum, full tuples in the txt): hit-depth 2 / LNA-depth 11, hit-depth 3 / 11, hit-depth 5 / 10; 16 shared in this run (LNA side 23 526 in 120 s; E-136 had 25, 27 for the referee: wall clock). `None` keys are dropped on both sides (69 on the hit side), so no None-equality.

Tilting-only BFS (`skeptic_n6tilt.txt`): both sides keep a step only if `tiltingPlus` holds, plus gate and key. LNA side 27 518 (150 s). Hits 0-15: tilting-only seen 323-501, shared with LNA side 0 for all 16. The unfiltered runs of E-136 (835-866 seen, 25 shared) therefore owe every meeting to the non-tilting step(s) out of the hit (it is the only non-tilting step on the shortest path; I did not enumerate all meeting paths).

A step that is not tilting can still land in the same derived class, so the meeting is not shown false; it is unsupported. Also note that the unfiltered LNA side contains J != 0 rows ((2,1) x 631) whose key-preserving children are likewise unvetted; I did not test whether LNA-forward steps ever fail `tiltingPlus` here (E-086 says the guarded walk from LNAs does not).

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/039/skeptic_n6path.py 120 20 0,5,11 > workshop/rounds/039/skeptic_n6path.txt   # about 4 min
timeout 10m .venv/bin/python -u workshop/rounds/039/skeptic_n6tilt.py 150 12 0 > workshop/rounds/039/skeptic_n6tilt.txt     # about 9 min (16 hits)
```
Both execute the E-136 preamble (`toolsmith_n6close.py`) and `rounds/001/scholar_h015.py` (`tiltingPlus`, `rplus`, `cartan`); counts drift with the wall-clock caps.

## Prior record

H-015 (SUPPORTED; E-086 guard refuses every gate-admitted non-tilting parent reachable from LNAs; E-095 `tiltingPlus` == Cartan congruence), E-123/E-130 (J_i != 0 iff silting-not-tilting), E-134 (hits, gate-admitted, off every walk), E-136 (the meeting). Grepped HYPOTHESES/EXPERIMENTS/RETRACTIONS for "key-preserving ... not tilting": no entry records a key-preserving, gate-admitted, non-tilting step. This is new as an observation about H-015's domain: the guard is not a substitute for `tiltingPlus` off LNA-reachable parents (the E-134 hits are the examples), so any meet-in-the-middle must filter by `tiltingPlus`. It also gives a concrete reason E-134's forward walks from the hits showed d >= 3: nothing else.

## Code changed

New only: `workshop/rounds/039/skeptic_n6path.py`, `skeptic_n6tilt.py` (+ outputs). No library changes, no tests needed.

## Next

- toolsmith: put a `tiltingPlus` filter in `toolsmith_n6meet.py` (`--tilting-only`) and retest at longer caps and in the reverse direction (mutate each hit at every vertex with J = 0 only, look for an LNA); overnight candidate: hit side to closure (it is small, 323-501 at 12 s).
- theorist: the 16 hits share one key with the 2-LNA class; is there a tilting path at all? A proof either way would decide whether (3,1),(3,2) algebras are in the class.
- referee note: statements "E-131's d = 2 is only a walk statement" stands; E-136's "lie in the LNA class" should read "reached through a non-tilting step".
