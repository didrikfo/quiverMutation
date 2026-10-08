# Round 044 (conference): Toolsmith position

**Most promising question.** Is the Coxeter-key guard sound on walk-reachable steps: does a key-preserving step ever fail `tiltingPlus`? If yes, the BFS filter that every closure and meet search (E-134..E-142) relies on is measuring a larger set than tilting, and the T5 meet has to be redone with `tiltingPlus` as the filter.

Why this one: E-145 already shows gate-admitted, key-preserving, non-tilting walk descendants at n = 7 classes 1 and 2 (13 of 67 and 9 of 64 distinct). That breaks the E-138/E-141 law outside class 0. It is cheap to test for every class at n = 6..8 with the key guard off and `tiltingPlus` on each step. It also decides whether the T5 search results mean anything.

**Weakest claim the workshop relies on.** That "key preserved" implies "tilting" for walk steps, so the LNA-key BFS is a valid stand-in for the tilting class. E-145 contradicts it for J != 0 steps; E-138 and E-141 only hold for class 0. Every E-134..E-142 meet result inherits this. My own reverse search adds a second gap: 154 of 1500 edges (10.3%) are lost in reverse with no filter behind it (E-147), and the control is shallow (depth 2-4).

**Why not the other candidates.** The orbit-plus-mirror question (T3/T8) is already answered at catalogue level (E-064). The H-017 mutation search needs a node count and depth before its negatives mean anything, but that is a sizing task, not a question. S-1 is parked behind literature I cannot read.

**What I need.**
- From the experimentalist: a key-guard-off census at n = 7 (classes 0-2) with per-step `tiltingPlus` and Cartan congruence, sized with `--plan` first. Counts, not one-off examples.
- From the skeptic: an independent rebuild of one class-1 E-145 step by hand, since the class-1 counts rest on one script.
- From the theorist: a precise statement of the property the meet needs (tilting, silting, or key) and which steps must be tilting for H-015 to hold. Without it I cannot say which failure matters.
- From the scholar: the statement of the `tiltingPlus` criterion (AI 2.32(b) / Ladkani 2.3(c)) once the PDFs are readable. The library criterion is only as good as its match to that statement.
- From the chair: no overnight run. The census is expected to finish within the 10-minute command limit at n = 6..7.

**My commitments this round.** Name the reverse-search loss by depth as the first task of the next toolsmith slot (E-147's open item). Add a `--key-off` flag to `toolsmith_n6meet.py` only if the experimentalist's census needs it. No library changes.
