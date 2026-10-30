# <Title: the claim, as a sentence>

author: <persona id> · round: <NNN> · kind: <result | negative | conjecture | proof | tool | proposal | rediscovery | literature | backlog>
thread: <thread id from STATE.md, or "new"> · bears on: <H-nnn, F-nnn, ...>
scope: <n range · classes · caps on depth, time or sample · what was not covered>
power: <for a test that can say NO (no join, no separation, no failure): the control on which this method does say NO
       where NO is known to be true, or does catch a planted error, with its outcome; or "none" -- then the title
       says "(no power)". For a proof or literature submission: "n/a">
backlog: <phase 5 only: the backlog items taken, each marked settled | retired | converted | narrowed only>

## Claim

One paragraph. What is true, stated so that it could be wrong. Say exactly
what it does *not* claim if a reader would be tempted to read more into it.
A pattern seen only on a capped sample, one class, or one n is stated as
that, in the title as well ("at n = 7, class 0"), not as a law. Say what
result would refute it.

## Evidence

What was checked, over what range, how many cases, how many failures. Tables
of data are welcome and do not count towards the length limit. For a proof,
the proof, with the step you are least sure of marked.

## Reproduction

The exact command(s), each with how long it took:

```
timeout 10m .venv/bin/python ...
```

## Prior record

What `research/` already says about this (identifiers), and why this is new --
or that it is not, and what that means.

## Code changed

Files touched and why; which tests were run. "None" if none.

## Next

What should be asked next, and by which kind of researcher. Requests to other
personas go here ("toolsmith: a way to ...").
