"""Asking questions of a finished classification.

`mutationClassTable` is the write side: one row per LNA, keyed by its relation
string, which is what the search needs.  Reading a classification wants the
other shape entirely -- one row per *class* -- and it wants to be able to leave
most of the table on disk.  The row count grows like the Catalan numbers (1430
LNAs at n = 9, 16796 at n = 11, 58786 at n = 12) while the class count does not
(20 at n = 9, 127 quipus at order 12), so the summary is small even where the
table is not.  Opening the CSV in a spreadsheet stopped being a way to look at
the answer somewhere around n = 10.

So: a polars scan over the parquet, grouped by class, with the filters the
questions actually take -- by Coxeter polynomial, by quipu, by class size, by
what kind of class it is.  `classes()` is one query; `members()` is the only
call that touches the per-LNA rows, and only for the one class asked for.

What a class' *kind* means, read off its hereditary form:

| form | kind | what it says |
|---|---|---|
| `P^(m)_(k)` | quipu | the quipu it is derived equivalent to. Complete: same quipu, same class |
| `C(w,...)` | canonical type | the weight type of the canonical algebra with its Coxeter polynomial |
| `not piecewise hereditary` | not piecewise hereditary | derived equivalent to no hereditary algebra, so to no quipu. A *negative* statement -- it separates, it never merges |
| an AHU string | tree, not a quipu | a relation-free quiver was reached and its tree is not a quipu (H-006 would be broken by one) |
| `wl:...` | not a tree | the underlying graph has a cycle, which no LNA search has produced |
| empty | unnamed | no relation-free quiver reached and the theorem does not cover it |

A form of `a|b` means a search reached two non-isomorphic hereditary algebras
in one class, which cannot happen and would be a bug; it is reported as its own
kind rather than hidden.
"""

import os

import polars as pl

from . import mutationClassTable as table

QUIPU = "quipu"
CANONICAL = "canonical type"
NOT_PIECEWISE_HEREDITARY = "not piecewise hereditary"
NON_QUIPU_TREE = "tree, not a quipu"
NOT_A_TREE = "not a tree"
UNNAMED = "unnamed"
CONTRADICTORY = "contradictory"

KINDS = [QUIPU, CANONICAL, NOT_PIECEWISE_HEREDITARY, NON_QUIPU_TREE, NOT_A_TREE,
         UNNAMED, CONTRADICTORY]


def kindOfForm(form):
    """Which kind a hereditary-form value names.  See the module docstring."""
    form = (form or "").strip()
    if not form:
        return UNNAMED
    if "|" in form:
        return CONTRADICTORY
    if form == table.NOT_PIECEWISE_HEREDITARY:
        return NOT_PIECEWISE_HEREDITARY
    if form.startswith("P^"):
        return QUIPU
    if form.startswith("C("):
        return CANONICAL
    if form.startswith("wl:"):
        return NOT_A_TREE
    return NON_QUIPU_TREE


class MutationClass:
    """One derived equivalence class, as the summary wants it."""

    def __init__(self, name, size, coxeterPolynomial, hereditaryForm, representative):
        self.name = name
        self.size = size
        self.coxeterPolynomial = coxeterPolynomial
        self.hereditaryForm = hereditaryForm
        self.representative = representative

    @property
    def kind(self):
        return kindOfForm(self.hereditaryForm)

    @property
    def isQuipuClass(self):
        return self.kind == QUIPU

    def __repr__(self):
        return "MutationClass({0!r}, size={1})".format(self.name, self.size)

    def __eq__(self, other):
        return isinstance(other, MutationClass) and vars(self) == vars(other)


class Classification:
    """A finished classification of one length, read rather than computed.

    Built over a polars `LazyFrame`, so the per-LNA rows are only materialised
    for the queries that need them.  `classes()` reads six columns and groups;
    `members()` reads the rows of one class.
    """

    def __init__(self, frame, length = None):
        self._frame = frame.lazy() if isinstance(frame, pl.DataFrame) else frame
        self.length = length

    # -- opening one ------------------------------------------------------

    @classmethod
    def forLength(cls, length, directory = "."):
        """The classification `classify.py` writes for a length.

        Prefers the parquet, which can be scanned, and falls back to the CSV,
        which cannot -- so a table only written as CSV is read whole.
        """
        parquet = os.path.join(directory, "A_{0}_mutation_classes.parquet".format(length))
        if os.path.exists(parquet):
            return cls(pl.scan_parquet(parquet), length)
        csv = os.path.join(directory, "A_{0}_mutation_classes.csv".format(length))
        if os.path.exists(csv):
            return cls.fromCSV(csv, length)
        raise FileNotFoundError(
            "no classification of length {0} in {1!r}; run classify.py {0}".format(
                length, directory))

    @classmethod
    def fromCSV(cls, fileName, length = None):
        return cls(table.MutationClassTable.fromCSV(fileName, length).toDataFrame(), length)

    @classmethod
    def fromTable(cls, mutationClassTable):
        """Straight from a table in memory, as `classifyLength` returns one."""
        return cls(mutationClassTable.toDataFrame(), mutationClassTable.lineLength)

    # -- the summary ------------------------------------------------------

    def classes(self, kind = None, minSize = None, maxSize = None,
                coxeterPolynomial = None, quipu = None):
        """The classes, largest first, filtered by whatever was asked for.

        `kind` is one of the constants above; `quipu` matches the hereditary
        form exactly, so it selects the one class a quipu name identifies;
        `coxeterPolynomial` matches the printed polynomial exactly. Size bounds
        are inclusive.
        """
        assigned = self._frame.filter(pl.col(table.CLASS) != "")
        grouped = (
            assigned.group_by(table.CLASS)
            .agg([
                pl.len().alias("size"),
                pl.col(table.COXETER).first().alias("coxeter"),
                # The form is written onto every row of a class, but a class
                # settled by a merge can have it on only some of them.
                pl.col(table.HEREDITARY).filter(pl.col(table.HEREDITARY) != "")
                  .first().alias("form"),
                pl.col(table.RELATIONS).sort_by(
                    [pl.col(table.PATH).str.len_chars(),
                     pl.col(table.RELATIONS).str.len_chars(),
                     pl.col(table.RELATIONS)]
                ).first().alias("representative"),
            ])
            .collect()
        )
        found = [
            MutationClass(row[table.CLASS], row["size"], row["coxeter"],
                          row["form"] or "", row["representative"])
            for row in grouped.iter_rows(named = True)
        ]
        if kind is not None:
            found = [c for c in found if c.kind == kind]
        if quipu is not None:
            found = [c for c in found if c.hereditaryForm == quipu]
        if coxeterPolynomial is not None:
            found = [c for c in found if c.coxeterPolynomial == coxeterPolynomial]
        if minSize is not None:
            found = [c for c in found if c.size >= minSize]
        if maxSize is not None:
            found = [c for c in found if c.size <= maxSize]
        return sorted(found, key = lambda c: (-c.size, c.name))

    def classNamed(self, name):
        """One class by name, or None."""
        return next((c for c in self.classes() if c.name == name), None)

    def members(self, className):
        """The LNAs of one class: (relations, mutation path, numbering).

        The only query that touches the per-LNA rows, and it reads the rows of
        one class.  Ordered by mutation path length, so the class
        representative -- which needs no mutations -- comes first.
        """
        rows = (
            self._frame.filter(pl.col(table.CLASS) == className)
            .select([table.RELATIONS, table.PATH, table.NUMBERING])
            .collect()
        )
        members = [
            (row[table.RELATIONS], row[table.PATH], row[table.NUMBERING])
            for row in rows.iter_rows(named = True)
        ]
        return sorted(members, key = lambda m: (len(m[1]), len(m[0]), m[0]))

    # -- the questions worth having a name ---------------------------------

    def coxeterCollisions(self):
        """Coxeter polynomial -> the classes sharing it, where more than one does.

        The polynomial is a derived invariant but not a complete one, so this is
        where it stops separating classes.  For the quipu classes that is exactly
        the cospectral quipus, and there are none below order 9 -- research
        F-010, which computes the same thing from the trees without running a
        classification at all.
        """
        grouped = {}
        for mutationClass in self.classes():
            grouped.setdefault(mutationClass.coxeterPolynomial, []).append(mutationClass)
        return {poly: found for poly, found in grouped.items() if len(found) > 1}

    def kindCounts(self):
        """How many classes of each kind, and how many LNAs they hold."""
        counts = {}
        for mutationClass in self.classes():
            classes, members = counts.get(mutationClass.kind, (0, 0))
            counts[mutationClass.kind] = (classes + 1, members + mutationClass.size)
        return counts

    def unsettled(self):
        """The LNAs no search has placed.  Empty for a finished classification."""
        rows = (
            self._frame.filter(pl.col(table.CLASS) == "")
            .select(table.RELATIONS)
            .collect()
        )
        return [row[0] for row in rows.iter_rows()]

    def rowCount(self):
        return self._frame.select(pl.len()).collect().item()

    def __repr__(self):
        return "Classification(length={0})".format(self.length)
