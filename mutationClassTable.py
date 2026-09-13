"""The table of mutation classes for the linear Nakayama algebras of one length.

One row per LNA, keyed by its relation string.  The table was a list of lists
read straight out of a CSV, scanned linearly on every lookup; it is now a class
with a polars DataFrame behind it and a dict index over the key column, so a
lookup is O(1) instead of O(rows).  That matters because the class search looks
each of the LNAs it reaches up in the table, so the old cost was
O(reached * rows) per search and O(reached * rows^2) per length, against a row
count that grows like the Catalan numbers.

The CSV layout is unchanged, so a table written by this class is readable by the
old code and vice versa.
"""

import csv

import polars as pl

RELATIONS = "Relations"
CLASS = "Mutation class"
PATH = "Mutation path from class representative"
COXETER = "Coxeter polynomial"
NUMBERING = "Numbering"
HEREDITARY = "Hereditary form"

COLUMNS = [RELATIONS, CLASS, PATH, COXETER, NUMBERING, HEREDITARY]

# Tables written before the hereditary column existed have five columns; they
# are read back with that column empty.
LEGACY_COLUMNS = [RELATIONS, CLASS, PATH, COXETER, NUMBERING]


class MutationClassTable:
    """One row per LNA of a fixed length.

    Every column holds a string, as it did when this was a CSV read with the
    csv module: the Coxeter polynomial is the printed form of the sympy
    expression, the mutation path and the vertex numbering are ';'-joined
    integers, and an empty class means the LNA has not been reached yet.
    """

    def __init__(self, rows, lineLength = None):
        self._rows = [list(row) + [""] * (len(COLUMNS) - len(row)) for row in rows]
        self.lineLength = lineLength
        self._reindex()

    def _reindex(self):
        self._indexByRelations = {}
        for i, row in enumerate(self._rows):
            # A relation string is unique per LNA, so the first occurrence wins
            # and a duplicate would be a bug in the table's construction.
            self._indexByRelations.setdefault(row[0], i)

    # -- construction ------------------------------------------------------

    @classmethod
    def forLength(cls, lineLength, relationStrings):
        """An empty table with one unassigned row per relation string."""
        rows = [[relations] + [""] * (len(COLUMNS) - 1) for relations in relationStrings]
        return cls(rows, lineLength)

    @classmethod
    def fromCSV(cls, fileName, lineLength = None):
        """Read a table written by this class or by the original code.

        The original createMutationClassCSV wrote a header row and
        saveLineRelationsAndMutationsToCSV did not, so a file may or may not
        have one.  Both are accepted.
        """
        with open(fileName, newline = "") as f:
            rows = list(csv.reader(f))
        if rows and rows[0][0] == RELATIONS:
            rows = rows[1:]
        return cls(rows, lineLength)

    @classmethod
    def fromDataFrame(cls, frame, lineLength = None):
        rows = [[row[column] for column in COLUMNS] for row in frame.iter_rows(named = True)]
        return cls(rows, lineLength)

    # -- reading -----------------------------------------------------------

    def __len__(self):
        return len(self._rows)

    def rows(self):
        """The rows in file order.  Each is the live list, so edits show up."""
        return self._rows

    def rowFor(self, relationString):
        """The row for one LNA, or None if the table has no such LNA."""
        index = self._indexByRelations.get(relationString)
        return None if index is None else self._rows[index]

    def indexFor(self, relationString):
        return self._indexByRelations.get(relationString)

    def unassignedRelationStrings(self):
        """The LNAs no search has reached yet, in file order."""
        return [row[0] for row in self._rows if not row[1]]

    def classNames(self):
        return {row[1] for row in self._rows if row[1]}

    def membersOfClass(self, className):
        return [row[0] for row in self._rows if row[1] == className]

    def classesByCoxeterPolynomial(self):
        """Coxeter polynomial -> the set of class names carrying it.

        A polynomial with more than one class name is a merge candidate: the
        polynomial is invariant under derived equivalence, so those classes are
        either the same class with a mutation path the search did not find, or
        genuinely different classes that happen to share a polynomial.
        """
        grouped = {}
        for row in self._rows:
            if row[1]:
                grouped.setdefault(row[3], set()).add(row[1])
        return grouped

    # -- writing -----------------------------------------------------------

    def assign(self, relationString, className, mutationPath, coxeterPolynomial,
               numbering, hereditaryForm = ""):
        """Fill in the result columns for one LNA."""
        index = self._indexByRelations[relationString]
        self._rows[index] = [
            relationString, className, mutationPath, coxeterPolynomial, numbering,
            hereditaryForm,
        ]

    def setHereditaryFormForClass(self, className, hereditaryForm):
        """Record the hereditary algebras a class' search reached, on every row.

        A relation-free quiver reached from an LNA pins down its derived
        equivalence class completely, so this is the sharpest invariant the
        search produces.  Two classes with different non-empty values here are
        certainly distinct, whatever their Coxeter polynomials.
        """
        for row in self._rows:
            if row[1] == className:
                row[5] = hereditaryForm

    def classesByHereditaryForm(self):
        """Hereditary form -> the set of class names that reached it.

        Classes sharing a form are certainly the same derived equivalence class,
        so this is a merge certificate rather than a merge candidate.  The empty
        form means the search reached no relation-free quiver and says nothing.
        """
        grouped = {}
        for row in self._rows:
            if row[1] and row[5]:
                grouped.setdefault(row[5], set()).add(row[1])
        return grouped

    def renameClass(self, oldName, newName):
        """Point every row of one class at another class' name."""
        for row in self._rows:
            if row[1] == oldName:
                row[1] = newName

    # -- interchange -------------------------------------------------------

    def toDataFrame(self):
        return pl.DataFrame(
            {column: [row[i] for row in self._rows] for i, column in enumerate(COLUMNS)},
            schema = {column: pl.String for column in COLUMNS},
        )

    def writeCSV(self, fileName, header = False):
        """Write the table as CSV.

        header defaults to False to match saveLineRelationsAndMutationsToCSV,
        which rewrote the file without one after every class.
        """
        with open(fileName, "w", newline = "") as f:
            writer = csv.writer(f)
            if header:
                writer.writerow(COLUMNS)
            writer.writerows(self._rows)

    def writeParquet(self, fileName):
        """Write the table as parquet, for a length too large to reread as CSV."""
        self.toDataFrame().write_parquet(fileName)

    @classmethod
    def fromParquet(cls, fileName, lineLength = None):
        return cls.fromDataFrame(pl.read_parquet(fileName), lineLength)
