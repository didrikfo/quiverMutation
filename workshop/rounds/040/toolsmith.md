# Toolsmith, round 040 (conference)

No runs made this round. Read: persona, notebook, STATE, DIGEST (last four entries), STEERING.

## Most promising question (next few rounds)

Does any gate-admitted, non-tilting (J != 0) step keep the LNA Coxeter key when it is taken from an LNA-side walk? In other words, does the key-preserving graph ever carry a non-tilting step, or does the key guard always refuse it?

Why: E-136 put the 16 n = 6 fans into the LNA derived class through a meeting path, and E-139 showed that the path's hit -> M step fails `tiltingPlus`. E-140 showed that all 229 J != 0 rows on the n = 8 c0 walk are refused by the key guard. So the open question is whether a tilting-only meet finds anything at all. That is a yes/no question with a computable answer at n = 6, 7, and it decides whether the meet-in-the-middle idea (the toolsmith's own lesson from round 038) is usable here.

Next concrete step (toolsmith, if the chair assigns it): `--tilting-only` in `toolsmith_n6meet.py`, with a positive control (a known tilting-only meet at n = 6) and a closure flag for the LNA side. No overnight run.

## Weakest claim the workshop relies on

"A key-preserving walk stays in the same derived class." E-136's "the 16 fans lie in the LNA class" rests on the Coxeter guard (H-015), not on the gate and not on `tiltingPlus` (E-139 found a gate-admitted, key-preserving step that fails `tiltingPlus`). The key is an exact presentation invariant, but preserving it along a step is not shown to preserve the derived class. This weakness also sits under E-131 (d = 2 at J != 0 on walks), which is a walk statement only.

Secondary: "0 hits reached" in bounded BFS (E-134, E-132) is a bounded miss, as E-136 already showed; any new "no hit" needs a positive control in the same family.

## What I need

- skeptic: a J != 0 step whose child passes the key guard, at n <= 7 (the request already open from round 039); or confirmation that none exists in the tilting-only sample.
- theorist: R(x) = 1 + t for j = e_i, or a statement of which relation makes the hit -> M step non-tilting. This tells me whether a tilting-only meet can ever reach the hits.
- experimentalist: the n = 8 c1 and c2 column of J != 0 rows with the key-guard result, so the E-140 pattern can be checked beyond class 0.
- scholar: whether `tiltingPlus` and Ladkani 2.3(c) are the same test on gate-admitted steps (E-068 open item), since that decides which check the meet should use.

## Notes for the chair

- The S-1 request (deletion transport) needs a key at n = 12 that I have not computed; I can build it after the T5 meet is settled, not now.
- No change to the library or the notebook is proposed this round.
