"""The Coxeter polynomials of a whole length, as a table to look a match up in.

The question both of the searches in `treeSearch` and `quipuRelations` ask is the
same one: *is this algebra's Coxeter polynomial carried by a linear Nakayama
algebra of the same length, and is that LNA one the quipu theorem already
places?*  The polynomial is a derived invariant and derived equivalence forces
the number of simples, so a candidate can only ever be equivalent to an LNA of
its own length, and a difference in the polynomial settles it outright.  A
*match* settles nothing on its own -- F-010 is where the polynomial stops
separating classes -- so a match is a lead, and a non-match is a result.

Three tables, all cached, all keyed by the integer coefficient tuples of
`invariants.coxeterCoefficients`:

* `lnaKeyIndex(length)` -- every LNA of the length, grouped by polynomial;
* `quipuKeySet(order)` -- the polynomials the quipus of the order carry, which
  are exactly the polynomials of the classes the quipu theorem names;
* `quipuClassLnas(length)` -- the LNAs that are *provably in* a quipu class,
  by seeding the theorem and closing under the free move and the double
  mutation (research F-032), which is what `freeMoves.coverage` computes.

`lnaStatus` puts the last two together and is the one a search calls: for each
LNA it says QUIPU (in a quipu class, by the moves), NOT_QUIPU (provably in no
quipu class, since no quipu of the order carries its polynomial), or UNPLACED
(the moves do not reach it and its polynomial is a quipu's, so the question is
open and a search would have to settle it).
"""

import functools

from . import freeMoves
from . import invariants
from . import nakayama
from . import quipuForms


QUIPU = "quipu"
NOT_QUIPU = "not quipu"
UNPLACED = "unplaced"


def lnaCartanMatrix(length, relLengths):
    """The Cartan matrix of an LNA, straight from its relation lengths.

    Building the path algebra and counting paths gives the same matrix and costs
    a hundred times as much; at `n = 11` that is the difference between seconds
    and half an hour over the whole length.  The path `i -> j` is nonzero exactly
    when no relation starts at or after `i` and ends at or before `j`, so the
    first relation starting at or after `i` bounds the row.
    """
    relations = [(start + 1, arrows) for start, arrows in enumerate(relLengths) if arrows]
    matrix = [[0] * length for _ in range(length)]
    for source in range(1, length + 1):
        stop = min((start + arrows for start, arrows in relations if start >= source),
                   default = length + 1)
        for target in range(source, min(length, stop - 1) + 1):
            matrix[target - 1][source - 1] = 1
    return matrix


def lnaCoxeterKey(length, relLengths):
    """The Coxeter polynomial of one LNA, as an integer coefficient tuple."""
    return invariants.coxeterCoefficients(lnaCartanMatrix(length, tuple(relLengths)))


@functools.lru_cache(maxsize = None)
def lnaKeys(length):
    """Every LNA of the length, as a dict from relation lengths to polynomial."""
    return {relLengths: lnaCoxeterKey(length, relLengths)
            for relLengths in nakayama.allRelationLengths(length)}


@functools.lru_cache(maxsize = None)
def lnaKeyIndex(length):
    """The same table read backwards: polynomial -> the LNAs carrying it.

    The LNAs are named the way the rest of the repo names them, by the digits of
    their relation lengths -- `3033030` -- and sorted, so the result is stable.
    """
    index = {}
    for relLengths, key in lnaKeys(length).items():
        index.setdefault(key, []).append(relLengths)
    return {key: tuple(sorted(members)) for key, members in index.items()}


@functools.lru_cache(maxsize = None)
def quipuKeyIndex(order):
    """Polynomial -> the quipus of the order carrying it, in the paper's notation.

    A quipu quiver has no relations, so its Cartan matrix is the reachability
    matrix of the quiver and any orientation gives the same polynomial -- the
    orientations of a tree are all derived equivalent by BGP reflection.  More
    than one quipu under one polynomial is F-010's cospectral pair, first seen at
    order 9.
    """
    index = {}
    for parameters in quipuForms.allQuipusOfOrder(order):
        algebra = nakayama.QuipuAlgebra(*parameters)
        key = invariants.coxeterKey(algebra)
        index.setdefault(key, []).append(quipuForms.formatQuipu(parameters))
    return {key: tuple(sorted(names)) for key, names in index.items()}


@functools.lru_cache(maxsize = None)
def quipuKeySet(order):
    """The set of polynomials carried by some quipu of the order."""
    return frozenset(quipuKeyIndex(order))


@functools.lru_cache(maxsize = None)
def quipuClassLnas(length):
    """The LNAs provably in a quipu class, as a frozenset of relation lengths.

    Seeding every LNA the quipu theorem names outright and closing under the
    free move and the double mutation, which is research F-032's route and
    places every LNA in a quipu class at `n <= 11`.  The rule table is left out:
    F-032 measured that it adds nothing on top of these two.
    """
    result = freeMoves.coverage(length, rules = [], free = True, edges = True, doubles = True)
    return frozenset(result['covered'])


@functools.lru_cache(maxsize = None)
def lnaStatus(length):
    """Relation lengths -> QUIPU, NOT_QUIPU or UNPLACED, for every LNA.

    NOT_QUIPU is a *proof*: the Coxeter polynomial is a derived invariant, so an
    LNA whose polynomial no quipu of the order carries is derived equivalent to
    no hereditary algebra of quipu type.  UNPLACED is an admission: the known
    moves do not reach the LNA from a seed, but its polynomial is a quipu's, so
    only a search or a finer invariant can say which it is.
    """
    covered = quipuClassLnas(length)
    keys = lnaKeys(length)
    quipuPolynomials = quipuKeySet(length)
    status = {}
    for relLengths, key in keys.items():
        if relLengths in covered:
            status[relLengths] = QUIPU
        elif key not in quipuPolynomials:
            status[relLengths] = NOT_QUIPU
        else:
            status[relLengths] = UNPLACED
    return status


def className(relLengths):
    """The repo's name for an LNA: the digits of its relation lengths."""
    return "".join(str(arrows) for arrows in relLengths)


def summary(length):
    """Counts of each status, for a report header."""
    counts = {QUIPU: 0, NOT_QUIPU: 0, UNPLACED: 0}
    for value in lnaStatus(length).values():
        counts[value] += 1
    return counts
