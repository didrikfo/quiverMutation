"""The double mutation of arXiv:2310.08346, `proposition:doubleMutation`.

The paper's statement, verbatim in substance:

    Let r be a relation from vertex s to t, such that there is another relation
    starting in vertex s - 1, but no relation starting in t - 1.  Then Lambda is
    derived equivalent to the Nakayama algebra obtained from it by

    * if t is not the final vertex, adding a new relation, of the same length as
      r, from s + 1 to t + 1;
    * shortening at the start every relation starting properly within r;
    * lengthening toward its end any relation ending properly within r.

The proof is two left tilting mutations at t, so this is a mutation equivalence
and carries the sequence `[-t, -t]`.  The paper calls the result `L_t(Lambda)`
and names the dual `R_s(Lambda)` without stating it; here it is taken through
the relation dual of F-026 rather than written out by hand.

**Why `lnaMoves` could never hold it.**  By minimality a relation cannot lie
inside another, so every relation the proposition changes -- one starting or
ending properly within `r` -- runs out past `r`'s far end.  Every one of them
straddles any window around `r`, and `matchesAt` refuses a straddled window.
So no width of table rule states this, and no discovery run could have found it:
the rules that do exist are its instances with nothing crossing `r` (the pair
slide, the end collapse), which is why F-023 kept finding the rules blocked by a
bystander.

**Two things read off the statement, both left to the engine to decide.**

* The hypothesis "a relation starts at s - 1" is used in the proof only to make
  the new vertex slot in between `P_s` and `P_{s-1}`, so the quiver is again a
  line.  At `s = 1` there is no `P_0` and the quiver is a line regardless, so
  `allowSource` admits `s = 1` with no companion.  F-029's source doubling is
  exactly that case with nothing crossing `r`.
* Lengthening a relation can make it contain another -- the companion at
  `s - 1` ending at `t - 1` lengthens onto `r` itself -- and then it is no longer
  minimal and goes.  The paper says "a minimal set of relations" throughout, so
  the result is reduced to one.
"""

from . import lnaMoves


def intervalsOf(relLengths):
    """(start, end) vertices of each relation."""
    return [(start, start + arrows) for start, arrows in lnaMoves.relationsOf(list(relLengths))]


def relLengthsOf(length, intervals):
    """Back to the per-vertex row, keeping only a minimal set; None if not an LNA."""
    intervals = set(intervals)
    minimal = [(a, b) for a, b in intervals
               if not any((c, d) != (a, b) and a <= c and d <= b for c, d in intervals)]
    row = [0] * (length - 2)
    for a, b in minimal:
        if not (1 <= a <= length - 2) or row[a - 1]:
            return None
        row[a - 1] = b - a
    if not lnaMoves.isAdmissible(length, row):
        return None
    return tuple(row)


def dualIntervals(length, intervals):
    return [(length + 1 - b, length + 1 - a) for a, b in intervals]


def leftDoubleMutation(length, relLengths, s, allowSource = True):
    """`L_t` for the relation starting at `s`: (relation lengths, sequence), or None."""
    intervals = intervalsOf(relLengths)
    ends = dict(intervals)
    if s not in ends:
        return None
    t = ends[s]
    starts = set(ends)
    if (t - 1) in starts:
        return None
    if (s - 1) not in starts and not (allowSource and s == 1):
        return None
    result = []
    for a, b in intervals:
        if (a, b) == (s, t):
            result.append((a, b))
        elif s < a < t:
            result.append((a + 1, b))
        elif s < b < t:
            result.append((a, b + 1))
        else:
            result.append((a, b))
    if t < length:
        result.append((s + 1, t + 1))
    row = relLengthsOf(length, result)
    if row is None or row == tuple(relLengths):
        return None
    return row, [-t, -t]


def rightDoubleMutation(length, relLengths, t, allowSource = True):
    """`R_s` for the relation ending at `t`, as the dual of `L`; sequence `[s, s]`."""
    dual = relLengthsOf(length, dualIntervals(length, intervalsOf(relLengths)))
    answer = leftDoubleMutation(length, dual, length + 1 - t, allowSource)
    if answer is None:
        return None
    row, sequence = answer
    back = relLengthsOf(length, dualIntervals(length, intervalsOf(row)))
    # Reversing the line sends vertex v to n + 1 - v and left to right (F-026).
    return back, [length + 1 - abs(v) for v in sequence]


def rewritesOf(length, relLengths, allowSource = True):
    """Every `L_t` and `R_s` out of an LNA, as (relation lengths, sequence) pairs."""
    reached = []
    for a, b in intervalsOf(relLengths):
        for answer in (leftDoubleMutation(length, relLengths, a, allowSource),
                       rightDoubleMutation(length, relLengths, b, allowSource)):
            if answer is not None:
                reached.append(answer)
    return reached
