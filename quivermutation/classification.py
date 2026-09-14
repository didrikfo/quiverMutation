"""Classifying every LNA of one length, end to end.

`classifyLength` is the entry point and runs four steps, each cheaper than the
one after it:

1. `seedTableFromQuipuTheorem` -- name every class the quipu theorem of
   arXiv:2305.06642 covers, with no searching at all;
2. `mutationSearch` -- search by mutation for the rows that are left;
3. `annotateHereditaryForms` -- give each class the hereditary algebra it is
   derived equivalent to, from the theorem where it applies and from a search
   where it does not;
4. `resolveMergeCandidates` -- settle whatever `mergeReport` could not, by
   searching harder for a mutation path between two candidate classes.

`mergeReport` is the certificate: `certain` for classes a shared hereditary form
proves equal, `separated` for ones differing forms prove distinct, `candidate`
for the rest.
"""

import math
import os

import numpy as np

from . import invariants
from . import lines
from . import lnaMoves
from . import mutation
from . import nakayama
from . import mutationClassTable
from . import pathAlgebra
from . import piecewiseHereditary
from . import search



def assignMutationClassInTable(table, mutationList, mutationClassName, printOutput = False):
    """Record a completed search in the table, and return the class it landed in.

    mutationList is the cleaned list of (path algebra, mutation path, vertex
    numbering) triples the search reached, all of them LNAs of the table's
    length and all derived equivalent to each other.

    Two passes, as the CSV version had.  The first looks for an LNA in the list
    that the table has already classified, and takes the shortest such link: if
    it finds one, this whole list belongs to that existing class rather than to
    a new one, and the mutation path from that class' representative has to be
    prefixed to every path recorded below.  The second pass fills in every LNA
    in the list the table has not reached yet.  If the first pass moved the
    class name, every row already carrying the old name is renamed at the end.

    This replaces saveLineRelationsAndMutationsToCSV.  Two things change:

    * the linear scan over every table row, per mutation, per pass, is now a
      dict lookup, since a relation string identifies at most one row;
    * the Coxeter polynomial is computed once rather than once per mutation.
      The original recomputed it inside the first pass and used whichever value
      the loop happened to leave behind, which is the one for the last entry, so
      that entry is the one used here.  Every entry has the same polynomial --
      they are mutations of each other -- but taking the last one keeps the
      output identical rather than merely equivalent.
    """
    inputMutationClassName = mutationClassName
    minMutationLength = np.inf
    baseMutationVertexString = ''
    baseVertexNumbering = {}
    if not bool(mutationList):
        return mutationClassName
    for n in range(1, len(mutationList[0][2]) + 1):
        baseVertexNumbering[n] = n
    coxPoly = invariants.coxeterPoly(mutationList[-1][0])

    for mut in mutationList:
        vertexNumbering = mut[2]
        mutationVertices = mut[1]
        row = table.rowFor(lines.relSetToString(mut[0].rels))
        if row is None or not bool(row[1]):
            continue
        mutationLength = len(row[2]) + len(mut[1])
        if mutationLength >= minMutationLength:
            continue
        mutationClassName = row[1]
        oldNumberingString = row[4]
        if bool(oldNumberingString):
            oldNumberingList = [int(v) for v in oldNumberingString.split(';')]
        else:
            oldNumberingList = range(1, len(mut[0].vertices()) + 1)
        renumberedReverseMutationVertexString = ''
        if bool(mut[1]):
            reverseMutationVertices = mutation.reverseMutationSequence(mutationVertices, vertexNumbering)
            renumbered = []
            for n in range(len(reverseMutationVertices)):
                if reverseMutationVertices[n] > 0:
                    renumberedVertex = oldNumberingList[reverseMutationVertices[n] - 1]
                else:
                    renumberedVertex = -oldNumberingList[-reverseMutationVertices[n] - 1]
                renumbered.append(str(renumberedVertex))
            renumberedReverseMutationVertexString = ';'.join(renumbered)
        if bool(row[2]) and bool(renumberedReverseMutationVertexString):
            baseMutationVertexString = ';'.join([row[2], renumberedReverseMutationVertexString])
        elif bool(renumberedReverseMutationVertexString):
            baseMutationVertexString = renumberedReverseMutationVertexString
        else:
            baseMutationVertexString = row[2]
        baseVertexNumbering = {}
        for n in range(1, len(vertexNumbering) + 1):
            baseVertexNumbering[n] = oldNumberingList[mutation.getVertexNumberingKeyFromValue(vertexNumbering, n) - 1]
        minMutationLength = mutationLength

    for mut in mutationList:
        localVertexNumbering = mut[2]
        localVertexNumberingList = [localVertexNumbering[n] for n in range(1, len(localVertexNumbering) + 1)]
        vertexNumberingList = [baseVertexNumbering[localVertexNumbering[n]] for n in range(1, len(baseVertexNumbering) + 1)]
        vertexNumberingString = ';'.join([str(n) for n in vertexNumberingList])
        relSetString = lines.relSetToString(mut[0].rels)
        row = table.rowFor(relSetString)
        if row is None or bool(row[1]):
            continue
        if bool(mut[1]):
            localMutationVertexString = ';'.join(
                str(vertexNumberingList[localVertexNumberingList.index(v)]) for v in mut[1]
            )
        else:
            localMutationVertexString = ''
        if bool(baseMutationVertexString) and bool(localMutationVertexString):
            totalMutationVertexString = ';'.join([baseMutationVertexString, localMutationVertexString])
        elif bool(localMutationVertexString):
            totalMutationVertexString = localMutationVertexString
        else:
            totalMutationVertexString = baseMutationVertexString
        table.assign(relSetString, mutationClassName, totalMutationVertexString,
                     str(coxPoly.as_expr()), vertexNumberingString)

    if printOutput:
        print('as part of the class ', mutationClassName)
    if mutationClassName != inputMutationClassName:
        table.renameClass(inputMutationClassName, mutationClassName)
    if printOutput:
        for row in table.rows():
            print(row)
    return mutationClassName


def expandClassByMoves(table, lineLength, relationString, className, coxeterPolynomial,
                       hereditaryForm = ''):
    """Fill in every LNA reachable from one by the verified moves of lnaMoves.

    Each move is a rewrite on the relation lengths with a mutation sequence that
    realises it, checked exhaustively against the mutation engine, so the whole
    orbit belongs to one derived equivalence class with a path to prove it.  This
    costs a table lookup per move where the depth-first search costs a subtree.

    Returns the number of rows filled in.
    """
    relLengths = lines.relationStringToLineRelLengths(lineLength, relationString)
    orbit = lnaMoves.closureUnderMoves(lineLength, relLengths)
    filled = 0
    for name, (sequence, numbering) in orbit.items():
        memberString = nakayama.LinearNakayamaAlgebra(lineLength, name).relationString()
        row = table.rowFor(memberString)
        if row is None or bool(row[1]):
            continue
        table.assign(memberString, className,
                     ';'.join(str(v) for v in sequence),
                     coxeterPolynomial,
                     ';'.join(str(numbering[p]) for p in range(1, lineLength + 1)),
                     hereditaryForm)
        filled += 1
    return filled


def adoptClassesByMoves(table, lineLength, maxRounds = 10):
    """Give every unclassified LNA the class of any classified LNA in its orbit.

    The outward direction -- expanding a known class along its move orbit --
    only reaches what the orbit of a *classified* LNA contains.  Running it
    inward as well catches the LNAs whose own orbit happens to touch something
    already classified, which is the same relation read the other way and costs
    the same lookup.

    Returns the number of rows placed.
    """
    placed = 0
    for _round in range(maxRounds):
        changed = 0
        for relationString in list(table.unassignedRelationStrings()):
            relLengths = lines.relationStringToLineRelLengths(lineLength, relationString)
            orbit = lnaMoves.closureUnderMoves(lineLength, relLengths)
            for name, (sequence, numbering) in orbit.items():
                memberString = nakayama.LinearNakayamaAlgebra(lineLength, name).relationString()
                row = table.rowFor(memberString)
                if row is None or not row[1]:
                    continue
                table.assign(relationString, row[1], '', row[3], '', row[5])
                changed += 1
                break
        placed += changed
        if not changed:
            break
    return placed


def seedTableFromQuipuTheorem(table, lineLength, printOutput = True, expandByMoves = True):
    """Assign every LNA the quipu theorem covers, before any searching.

    Theorem `thm:QuipuToAn` of arXiv:2305.06642 names the derived equivalence
    class of any LNA with almost separate relations outright, in O(1), so every
    such row can be filled in before a single mutation is computed.  The class is
    named by its quipu rather than by an LNA, which is both a better name and
    makes two seeded classes with the same quipu literally the same class.

    The depth-first search then only has to place the rows the theorem misses,
    and those inherit a seeded class as soon as the search reaches any seeded LNA
    -- which assignMutationClassInTable already does, since a seeded row is an
    already-classified row like any other.

    Coverage falls as the length grows (100% at n = 4, 54% at n = 8, 19% at
    n = 12) but the quipus it names do not: it already finds every class of
    every length checked so far.

    With expandByMoves, each seeded LNA also drags in its whole orbit under the
    verified moves of lnaMoves, which reaches LNAs the theorem does not cover at
    all -- the ones whose relations overlap too much -- without any searching.

    Returns the number of rows filled in.
    """
    seeded = expanded = 0
    for row in table.rows():
        if row[1]:
            continue
        form = search.hereditaryFormFromTheorem(lineLength, row[0])
        if not form:
            continue
        pathAlg = nakayama.LinearNakayamaAlgebra.fromRelationString(lineLength, row[0])
        polynomial = str(invariants.coxeterPoly(pathAlg).as_expr())
        table.assign(
            row[0], form, '', polynomial,
            ';'.join(str(v) for v in range(1, lineLength + 1)), form)
        seeded += 1
        if expandByMoves:
            expanded += expandClassByMoves(table, lineLength, row[0], form, polynomial, form)
    adopted = 0
    if expandByMoves:
        adopted = adoptClassesByMoves(table, lineLength)
    if printOutput:
        print('Seeded {0} rows from the quipu theorem, {1} more by expanding those '
              'classes along move orbits, {2} more by adopting a class through a '
              'move orbit: {3} of {4} rows in {5} classes, with no search'.format(
                  seeded, expanded, adopted, seeded + expanded + adopted, len(table),
                  len(table.classNames())))
    return seeded + expanded + adopted


def annotateHereditaryForms(table, lineLength, maxDepth = 0, printOutput = True):
    """Give every class in the table the hereditary algebra it is equivalent to.

    Three passes, cheapest first.

    1. The quipu theorem, for every LNA in the table with almost separate
       relations.  A class picks up the form of any such member it contains.
       Two members disagreeing would mean either the search merged two classes
       that are not equal or the inversion of the theorem is wrong, so that is
       reported rather than silently resolved.
    2. Whatever the class searches already found, which is kept where pass 1
       says nothing.
    3. Only if maxDepth > 0: an iterative-deepening search from the members of
       each class still without a form.  This is the expensive one -- a class
       that reaches no relation-free quiver at all makes it explore the whole
       tree from every member -- so it is off by default.

    Returns the table.
    """
    formsByClass = {}
    conflicts = {}
    # Shared across classes: many algebras delete down to the same smaller one.
    deletionCache = {}
    for row in table.rows():
        if not row[1]:
            continue
        form = search.hereditaryFormFromTheorem(lineLength, row[0])
        if not form:
            continue
        known = formsByClass.setdefault(row[1], form)
        if known != form:
            conflicts.setdefault(row[1], {known}).add(form)

    for className, forms in conflicts.items():
        print('WARNING: class {0} contains LNAs equivalent to different quipus: {1}. '
              'Either the search merged two distinct classes, or the theorem was '
              'inverted wrongly.'.format(className, sorted(forms)))

    for className in sorted(table.classNames()):
        rows = [row for row in table.rows() if row[1] == className]
        fromTheorem = formsByClass.get(className, '')
        fromSearch = next((row[5] for row in rows
                           if mutationClassTable.isIdentifyingForm(row[5])), '')
        if fromTheorem and fromSearch and fromTheorem != fromSearch:
            print('WARNING: class {0} is {1} by the theorem but the search reached '
                  '{2}'.format(className, fromTheorem, fromSearch))
        form = fromTheorem or fromSearch
        source = 'theorem' if fromTheorem else ('search' if form else '')
        if not form and maxDepth > 0:
            form = search.findHereditaryFormForClass(table, lineLength, className, maxDepth, printOutput)
            source = 'search' if form else ''
        if not form:
            # No quipu.  Either the class is not piecewise hereditary at all, in
            # which case it is in no quipu class and the certificate says so, or
            # it is of canonical type and the Coxeter polynomial names which.
            certified = None
            for row in rows:
                chain = piecewiseHereditary.notPiecewiseHereditaryByDeletion(
                    lineLength, lines.relationStringToLineRelLengths(lineLength, row[0]),
                    deletionCache)
                if chain is not None:
                    certified = (row[0], chain)
                    break
            if certified is not None:
                form = mutationClassTable.NOT_PIECEWISE_HEREDITARY
                witness, chain = certified
                source = 'not piecewise hereditary, witness {0!r} via {1}'.format(
                    witness, ' -> '.join(step[2] for step in chain))
            elif rows and rows[0][3]:
                weights = piecewiseHereditary.canonicalWeightType(lineLength, rows[0][3])
                if weights is not None:
                    form = 'C({0})'.format(','.join(str(w) for w in weights))
                    source = 'canonical algebra' + (
                        ', tubular' if piecewiseHereditary.isTubular(weights) else '')
        table.setHereditaryFormForClass(className, form)
        if printOutput:
            print('class {0}: {1} ({2})'.format(className, form or '-', source or 'nothing'))
    return table


def _memberAndItsDual(lineLength, relationString):
    """An LNA and its relation dual, both as path algebras.

    Reversing every arrow of an LNA keeps it in the same derived equivalence
    class -- one of the three class-preserving operations of arXiv:2305.06642 --
    and turns right mutations into left ones, so searching from both covers both
    directions of a reachability that is otherwise one-way.
    """
    algebra = nakayama.LinearNakayamaAlgebra.fromRelationString(lineLength, relationString)
    dual = algebra.relationDual()
    return [algebra] if dual == algebra else [algebra, dual]


def resolveMergeCandidates(table, lineLength, depth = 8, printOutput = True):
    """Settle the classes mergeReport could not, by searching harder for a link.

    A candidate is a group of classes sharing a Coxeter polynomial that the
    hereditary form does not separate, because at least one of them has no form.
    Such a class is one the quipu theorem does not cover and whose search never
    reached an LNA that it does, so the only way to place it is to find a
    mutation path from it to a classified LNA.

    For each such class this runs a deeper search from each of its members in
    turn and merges as soon as one reaches a row belonging to another class.
    Since the reached row already carries a hereditary form, the merge inherits
    it, and the group stops being a candidate.

    The search is run from each member *and from its relation dual*, because
    mutationSearchDepthFirst only walks right mutations, which makes reachability
    directional: A can reach B at depth d while B reaches nothing at that depth.
    Since rightMutate(dual(P)) = dual(leftMutate(P)), a right-mutation path out
    of dual(X) is a left-mutation path out of X, and the relation dual of an LNA
    is derived equivalent to it, so everything reached either way is in X's
    class.

    Returns the list of (class merged away, class merged into) pairs.
    """
    merges = []
    report = mergeReport(table)
    unresolved = set()
    for classNames in report['candidate'].values():
        for className in classNames:
            rows = [row for row in table.rows() if row[1] == className]
            if not any(row[5] for row in rows):
                unresolved.add(className)

    for className in sorted(unresolved):
        members = sorted(table.membersOfClass(className), key = lambda r: (len(r), r))
        if not members:
            continue
        merged = False
        for relationString in members:
            for startPoint in _memberAndItsDual(lineLength, relationString):
                reached = []
                search.mutationSearchDepthFirst(startPoint, depth, [], 'resolve', printOutput = False,
                                         collected = reached)
                for mut in lines.mutationListLineCleanup(reached, printOutput = False):
                    row = table.rowFor(lines.relSetToString(mut[0].rels))
                    if row is None or not row[1] or row[1] == className:
                        continue
                    if printOutput:
                        print('class {0} reaches {1} (class {2}) at depth {3} from {4!r}'.format(
                            className, row[0], row[1], depth, relationString))
                    target = row[1]
                    form = row[5]
                    table.renameClass(className, target)
                    if form:
                        table.setHereditaryFormForClass(target, form)
                    merges.append((className, target))
                    merged = True
                    break
                if merged:
                    break
            if merged:
                break
        if not merged and printOutput:
            print('class {0} still unresolved at depth {1}'.format(className, depth))
    return merges


def mergeReport(table):
    """What the table says about which classes should be merged.

    Returns a dict with three keys:

    * 'certain'   -- hereditary form -> class names that all reached it.  More
                     than one name means those classes are provably the same
                     class and the search simply missed the mutation path.
    * 'candidate' -- Coxeter polynomial -> class names sharing it that the
                     hereditary form does not settle, because at least one of
                     them has no form recorded.  These need a deeper search.
    * 'separated' -- Coxeter polynomial -> class names sharing it that the
                     hereditary form proves to be distinct classes.  These must
                     not be merged, and are the cases where the Coxeter
                     polynomial is not a complete invariant.
    """
    byForm = table.classesByHereditaryForm()
    formOfClass = table.formOfEachClass()

    certain = {form: names for form, names in byForm.items() if len(names) > 1}
    candidate = {}
    separated = {}
    for polynomial, classNames in table.classesByCoxeterPolynomial().items():
        if len(classNames) < 2:
            continue
        forms = {formOfClass.get(name, '') for name in classNames}
        if '' in forms:
            candidate[polynomial] = classNames
        elif len(forms) > 1:
            # Different forms, so distinct classes.  This includes one class
            # certified not piecewise hereditary against one with a quipu: the
            # certificate cannot merge classes but it does separate them from
            # every quipu class.
            separated[polynomial] = classNames
        elif forms == {mutationClassTable.NOT_PIECEWISE_HEREDITARY}:
            # Both only carry the negative certificate, which says nothing about
            # whether they are the same class.
            candidate[polynomial] = classNames
    return {'certain': certain, 'candidate': candidate, 'separated': separated}


def mutationSearch(lineLength, mutationDepthStart, startRow = 0, createNewCSVfile = False,
                   printMutations = False, fileName = None, table = None, writeEveryClass = True,
                   collectHereditary = False, seedFromQuipuTheorem = False,
                   printProgress = True):
    """Classify every LNA of the given length by depth-first tilting mutation.

    Walks the table of all Catalan(lineLength - 1) LNAs.  For each one that no
    earlier search has reached, runs a depth-first search of mutations out of it
    and records every LNA that search reaches as belonging to the same class.

    The depth decays as max(mutationDepthStart - floor(log10(row)), 2), which is
    what keeps later rows affordable, and is also why the result is a lower
    bound on each class rather than the classification itself: two LNAs in the
    same class end up in different classes here if no mutation path between them
    fits in the depth.  Merging those is a separate step.

    Returns the MutationClassTable.  Pass `table` to continue an existing one,
    and startRow to begin partway through it.

    seedFromQuipuTheorem fills in every row the quipu theorem covers before the
    search starts, so the search only has to place the rest.  The classes are
    then named by their quipu rather than by an LNA.

    collectHereditary makes the search also record every relation-free quiver it
    passes through, which identifies the class completely.  It is off by default
    because annotateHereditaryForms gets the same answer from the quipu theorem
    in O(1), while collecting during the search costs a canonical form at every
    relation-free node and roughly doubles the run time.
    """
    if fileName is None:
        fileName = 'A_{0}_mutation_classes.csv'.format(lineLength)
    if table is None:
        if createNewCSVfile:
            table = mutationClassTable.MutationClassTable.forLength(
                lineLength,
                [lines.relSetToString(relSet) for relSet in lines.generateAllPossibleLineRelations(lineLength)],
            )
            table.writeCSV(fileName, header=True)
        else:
            table = mutationClassTable.MutationClassTable.fromCSV(fileName, lineLength)

    if seedFromQuipuTheorem:
        seedTableFromQuipuTheorem(table, lineLength, printOutput=printProgress)
        table.writeCSV(fileName, header=True)

    mutationDepth = mutationDepthStart
    numberOfRows = len(table)
    for i in range(numberOfRows):
        row = table.rows()[(startRow + i) % numberOfRows]
        if not bool(row[1]):
            algebra = nakayama.LinearNakayamaAlgebra.fromRelationString(lineLength, row[0])
            lineNumberString = algebra.className()
            if printProgress:
                print('row {0}/{1}: searching from {2} at depth {3}'.format(
                    i + 1, numberOfRows, lineNumberString, mutationDepth))
            pathAlg = algebra
            if printMutations:
                pathAlgebra.printPathAlgebra(pathAlg)
            mutList = []
            hereditaryFound = [] if collectHereditary else None
            search.mutationSearchDepthFirst(pathAlg, mutationDepth, [],
                                     'A{0}_{1}'.format(lineLength, lineNumberString),
                                     printOutput=printMutations, collected=mutList,
                                     collectedHereditary=hereditaryFound)
            cleanMutList = lines.mutationListLineCleanup(mutList, printOutput=printMutations)
            className = assignMutationClassInTable(table, cleanMutList, lineNumberString,
                                                   printOutput=printMutations)
            if hereditaryFound:
                forms = {}
                for canonical, quipu, path in hereditaryFound:
                    if canonical not in forms or len(path) < len(forms[canonical][1]):
                        forms[canonical] = (quipu, path)
                table.setHereditaryFormForClass(className, search.formatHereditaryForms(forms))
            if writeEveryClass:
                table.writeCSV(fileName, header=True)
        mutationDepth = int(max(mutationDepthStart - math.floor(math.log10(i + 1)), 2))
    table.writeCSV(fileName, header=True)
    return table


def classifyLength(lineLength, mutationDepthStart = 6, resolveDepth = 6, fileName = None,
                   printOutput = True, resume = False):
    """The whole classification of one length, end to end.

    1. Seed every LNA the quipu theorem covers, naming each class by its quipu.
    2. Depth-first search from each remaining unclassified LNA, which inherits a
       seeded class as soon as it reaches a seeded LNA.
    3. Name any class the search created but the theorem did not cover.
    4. Resolve what is left: any group of classes sharing a Coxeter polynomial
       that the hereditary form does not settle gets a deeper search, from each
       member and from its relation dual.

    `resume` continues from an existing CSV rather than starting over.  The table
    is written after every class searched, so an interrupted run -- and a long one
    will be interrupted, since a length-10 classification takes hours and does not
    survive the machine going away -- picks up where it stopped.

    Returns (table, report).  A report with empty 'candidate' means every class
    is settled: 'certain' entries are classes proved equal, 'separated' entries
    are classes proved distinct despite sharing a Coxeter polynomial.

    This replaces the hand-merge step.  For n <= 8 it reproduces the published
    classification with nothing left over.
    """
    if fileName is None:
        fileName = 'A_{0}_mutation_classes.csv'.format(lineLength)
    existing = None
    if resume and os.path.exists(fileName):
        existing = mutationClassTable.MutationClassTable.fromCSV(fileName, lineLength)
        if printOutput:
            print('Resuming from {0}: {1} of {2} rows already placed'.format(
                fileName, len(existing) - len(existing.unassignedRelationStrings()),
                len(existing)))
    table = mutationSearch(lineLength, mutationDepthStart, 0,
                           createNewCSVfile = existing is None,
                           fileName = fileName, table = existing,
                           seedFromQuipuTheorem = True,
                           printProgress = printOutput)
    annotateHereditaryForms(table, lineLength, printOutput = printOutput)
    merges = resolveMergeCandidates(table, lineLength, resolveDepth, printOutput = printOutput)
    for merged, into in merges:
        if printOutput:
            print('merged class {0} into {1}'.format(merged, into))
    report = mergeReport(table)
    for polynomial, classNames in report['certain'].items():
        target = sorted(classNames)[0]
        for className in classNames:
            if className != target:
                table.renameClass(className, target)
        if printOutput:
            print('merged {0} into {1} (same hereditary form)'.format(
                sorted(classNames), target))
    table.writeCSV(fileName, header = True)
    table.writeParquet(fileName.replace('.csv', '.parquet'))
    report = mergeReport(table)
    if printOutput:
        print('{0} LNAs, {1} classes, {2} still candidates, {3} separated'.format(
            len(table), len(table.classNames()), len(report['candidate']),
            len(report['separated'])))
    return table, report
