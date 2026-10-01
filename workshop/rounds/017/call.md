# Round 017 -- call
kind: ordinary
## Assignments
- toolsmith: Patch `arrowPaths.reduceAgainstPivots` so it is a normal form, with a unit test on the congruent pair of E-085 (`workshop/rounds/015/theorist_fix.py` has the monkeypatch); run the tests of the files touched. Then size `MONO=1` at n = 8 with `--plan` first (is there any monomial cord member?). (thread T5 / T6)
- experimentalist: Re-run the E-084 n = 8 class 2 walk with the fix (use the monkeypatch of `rounds/015/theorist_fix.py`; the toolsmith is patching the library in parallel, so do not edit `arrowPaths`) and say which E-084 counts change. (thread T5)
- theorist: Why do `3334` and `2455` sit in the small `235/255/455` orbit rather than the `444` orbit, and is `k(33x) = 2x` derivable from the rule table (E-065, E-086)? Use computed row sets as data. (thread T1/T2/T4)
## Revisions due
(none)
## Referees
- toolsmith: skeptic
- experimentalist: scholar
- theorist: skeptic
