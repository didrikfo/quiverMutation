#!/usr/bin/env python
"""Other families of quivers that could carry the classes the quipu theorem misses.

    python families.py trees 9 10 11        every tree, against every LNA
    python families.py quipus 9             every quipu with relations, against the
                                            LNAs that lie in no quipu class
    python families.py quipus 9 --verify 4  and search each lead for a mutation path
    python families.py free 8               are two-arrow relations free on a quipu?

The quipu theorem of arXiv:2305.06642 says an LNA whose relations are almost
separate is derived equivalent to the path algebra of a quipu -- a tree, with no
relations.  Every LNA it does not cover has to be classified some other way, and
the question both subcommands ask is whether some *other* family of quivers plays
the same role for those:

* `trees` -- every tree of the order, not just the quipus, hereditary as they all
  are.  F-031 did this for the trees of maximum degree three; this does the rest.
* `quipus` -- quipu quivers that do carry relations, in every orientation.  The
  theorem is about the relation-free ones only, and a walk from an LNA to its
  quipu passes through quivers that have relations, so there is no reason the
  shape should stop mattering when they do.

Both work by the Coxeter polynomial, which is a derived invariant: a candidate
whose polynomial no LNA of the length carries is ruled out outright, and one that
matches is a lead to be settled by a mutation search.  `--verify` runs that
search on the leads it finds.
"""

import argparse
import collections
import sys
import time

from quivermutation import coxeterTables as ct
from quivermutation import invariants as inv
from quivermutation import lines
from quivermutation import nakayama as nk
from quivermutation import quipuRelations as qr
from quivermutation import search as se
from quivermutation import treeSearch as ts


def relationsAsVertices(match):
    """The relations of a match, numbered from 1 as the rest of the repo does."""
    return ["-".join(str(vertex + 1) for vertex in path) for path in match['relations']]


def arrowsOf(match):
    """The quiver's arrows, numbered from 1."""
    arrows = []
    for (first, second), bit in zip(match['edges'], match['orientation']):
        tail, head = (first, second) if bit else (second, first)
        arrows.append("{0}->{1}".format(tail + 1, head + 1))
    return arrows


def reportTrees(orders):
    """Every tree of each order, matched against the LNAs of that length."""
    for order in orders:
        rows = ts.treeReport(order)
        nonQuipu = [row for row in rows if not row['isQuipu']]
        cospectral = [row for row in nonQuipu if row['cospectralQuipus']]
        leads = ts.leads(order)
        counts = ct.summary(order)
        print()
        print("order {0}: {1} trees, {2} not quipus; {3} LNAs, of which {4} in a quipu "
              "class, {5} in none, {6} unplaced".format(
                  order, len(rows), len(nonQuipu),
                  sum(counts.values()), counts[ct.QUIPU], counts[ct.NOT_QUIPU],
                  counts[ct.UNPLACED]))
        print("  {0} non-quipu trees share a polynomial with a quipu of the order "
              "(cospectral, so not derived equivalent to it)".format(len(cospectral)))
        for row in cospectral:
            print("     {0}  max degree {1}  cospectral with {2}".format(
                row['name'], row['maxDegree'], ", ".join(row['cospectralQuipus'])))
        if not leads:
            print("  no non-quipu tree shares a polynomial with an LNA outside a quipu "
                  "class -- nothing to follow up")
            continue
        print("  {0} LEADS:".format(len(leads)))
        for row in leads:
            print("     {0}  max degree {1}  {2}".format(
                row['name'], row['maxDegree'], inv.formatCoefficients(row['key'])))
            for statusName, members in sorted(row['byStatus'].items()):
                print("        {0}: {1}".format(
                    statusName, ", ".join(ct.className(m) for m in members)))
    return 0


def summariseLnas(lnas, limit = 8):
    """The LNAs under one polynomial, the ones outside a quipu class first.

    A polynomial shared with a quipu class can be carried by hundreds of LNAs,
    and printing them all buries the two or three the run is actually about.
    """
    ordered = sorted(lnas, key = lambda pair: (pair[1] == ct.QUIPU, pair[0]))
    shown = ["{0} ({1})".format(name, status) for name, status in ordered[:limit]]
    if len(ordered) > limit:
        shown.append("and {0} more".format(len(ordered) - limit))
    return ", ".join(shown)


def targetGroups(order, everyLna):
    """The polynomial groups a quipu-with-relations run is trying to hit.

    One entry per Coxeter polynomial, listing the LNAs that carry it and their
    status.  With `everyLna` off -- the default -- only the polynomials of LNAs
    in no quipu class are included, which are the classes the quipu theorem
    leaves unclassified and the ones a new family would be for.
    """
    status = ct.lnaStatus(order)
    groups = collections.defaultdict(list)
    for relLengths, key in ct.lnaKeys(order).items():
        if everyLna or status[relLengths] != ct.QUIPU:
            groups[key].append(relLengths)
    return groups


def reportQuipus(orders, minArrows, verifyDepth, showLimit, everyLna, keepLines,
                 dedupe = False):
    """Quipu quivers with relations, matched against the LNAs of the length."""
    statuses = None if everyLna else (ct.NOT_QUIPU, ct.UNPLACED)
    for order in orders:
        counts = ct.summary(order)
        groups = targetGroups(order, everyLna)
        print()
        print("order {0}, relations of {1} arrows or more".format(order, minArrows))
        print("  {0} LNAs: {1} in a quipu class, {2} in none, {3} unplaced".format(
            sum(counts.values()), counts[ct.QUIPU], counts[ct.NOT_QUIPU],
            counts[ct.UNPLACED]))
        print("  matching against {0} LNAs under {1} Coxeter polynomials".format(
            sum(len(members) for members in groups.values()), len(groups)))
        started = time.time()
        result = qr.search(order, minArrows = minArrows, statuses = statuses,
                           includeLines = keepLines,
                           keepPerKey = showLimit * 4, dedupe = dedupe)
        print("  {0} algebras walked, {1} through the sieve, {2} matching{3}, "
              "in {4:.0f}s".format(
                  result['walked'], result['sieved'], result['matched'],
                  " up to isomorphism" if dedupe else "", time.time() - started))
        byKey = result['examples']
        unmatched = [key for key in groups if key not in byKey]
        print("  {0} of the {1} polynomials are carried by a quipu with relations, "
              "{2} by none".format(len(byKey), len(groups), len(unmatched)))
        for key in sorted(unmatched):
            print("    no quipu algebra: {0}   {1}".format(
                ", ".join(ct.className(member) for member in groups[key]),
                inv.formatCoefficients(key)))
        for key, group in sorted(byKey.items()):
            group.sort(key = lambda match: (len(match['relations']),
                                            sum(len(path) for path in match['relations'])))
            lnaNames = summariseLnas(group[0]['lnas'])
            shapes = result['shapes'][key]
            print()
            print("  {0}".format(inv.formatCoefficients(key)))
            print("    {0} quipu algebras, over {1} quipu shapes; LNAs: {2}".format(
                result['counts'][key], len(shapes), lnaNames))
            print("    commonest shapes: {0}".format(", ".join(
                "{0} ({1})".format(name, count) for name, count
                in sorted(shapes.items(), key = lambda item: -item[1])[:5])))
            for match in group[:showLimit]:
                print("      {0}  arrows {1}  relations {2}".format(
                    match['quipu'], " ".join(arrowsOf(match)),
                    " ".join(relationsAsVertices(match))))
                if verifyDepth:
                    reached = verify(match, order, verifyDepth)
                    if reached:
                        print("        REACHES {0} at depth {1}".format(
                            ", ".join(sorted(reached)), verifyDepth))
                    else:
                        print("        reaches no LNA of the group at depth {0}".format(
                            verifyDepth))
    return 0


def verify(match, order, depth):
    """The LNAs of the match's group that a mutation search out of it reaches."""
    wanted = {name for name, _status in match['lnas']}
    algebra = qr.algebraFromMatch(match)
    reached = se.linesReachedFrom(algebra, depth)
    names = {"".join(str(arrows) for arrows
                     in lines.relationStringToLineRelLengths(order, relationString))
             for relationString in reached}
    return names & wanted


def reportMembers(orders, depth, showLimit):
    """The quipu algebras a mutation walk proves are in each unclassified class.

    The other side of `quipus`: instead of enumerating the family and asking
    which members *could* be in a class, this walks out of the LNAs themselves,
    so everything it prints is in the class it is printed under, with a mutation
    path behind it.  It is what a normal form would have to be read off (H-014).
    """
    for order in orders:
        status = ct.lnaStatus(order)
        outside = sorted(relLengths for relLengths, value in status.items()
                         if value != ct.QUIPU)
        print()
        print("order {0}: {1} LNAs in no quipu class, walked to depth {2}".format(
            order, len(outside), depth))
        reachingNone = []
        everything = set()
        for relLengths in outside:
            algebra = nk.LinearNakayamaAlgebra(order, list(relLengths))
            reached = qr.reachedQuipuAlgebras(algebra, depth)
            confirmed = {certificate: path for certificate, path in reached.items()
                         if not isLinearlyOrientedCertificate(certificate)}
            everything.update(confirmed)
            if not confirmed:
                reachingNone.append(ct.className(relLengths))
                continue
            simplest = sorted(confirmed, key = lambda certificate: (
                len(certificate[1]), sum(len(path) for path in certificate[1])))
            print("  {0}: {1} quipu algebras with relations".format(
                ct.className(relLengths), len(confirmed)))
            for certificate in simplest[:showLimit]:
                name, arrows, relations = qr.describeCertificate(certificate)
                print("      {0}  arrows {1}  relations {2}  via {3}".format(
                    name, " ".join(arrows), " ".join(relations), confirmed[certificate]))
        print("  {0} distinct quipu algebras confirmed; {1} LNAs reached none".format(
            len(everything), len(reachingNone)))
        for className in reachingNone:
            print("    reached none: {0}".format(className))
    return 0


def isLinearlyOrientedCertificate(certificate):
    """Whether a certificate is the line itself, which is an LNA and not news."""
    return all(head == tail + 1 for tail, head in certificate[0])


def reportFree(orders, minArrows, limit):
    """Whether deleting a two-arrow relation keeps the polynomial, on quipus."""
    for order in orders:
        checked, agreeing, failures = qr.freeRelationCheck(order, minArrows, limit)
        print("order {0}: {1} ideals with a two-arrow relation, {2} keep their Coxeter "
              "polynomial when it is deleted, {3} do not".format(
                  order, checked, agreeing, len(failures)))
        for failure in failures[:10]:
            print("    {0}  {1}".format(
                failure['quipu'],
                " ".join("-".join(str(v + 1) for v in path)
                         for path in failure['relations'])))
            print("      before {0}".format(inv.formatCoefficients(failure['before'])))
            print("      after  {0}".format(inv.formatCoefficients(failure['after'])))
    return 0


def main(argv = None):
    parser = argparse.ArgumentParser(
        description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest = "command", required = True)

    trees = subparsers.add_parser("trees", help = "every tree, against every LNA")
    trees.add_argument("orders", type = int, nargs = "+")

    quipus = subparsers.add_parser("quipus", help = "quipus with relations, against the LNAs")
    quipus.add_argument("orders", type = int, nargs = "+")
    quipus.add_argument("--min-arrows", type = int, default = 3, dest = "minArrows",
                        help = "the shortest relation to allow (default 3; 2 is the "
                               "whole enumeration and costs several times as much)")
    quipus.add_argument("--verify", type = int, default = 0, dest = "verifyDepth",
                        help = "search each lead printed for a mutation path to the "
                               "LNA, to this depth (default 0, off)")
    quipus.add_argument("--show", type = int, default = 5, dest = "showLimit",
                        help = "how many algebras to print per polynomial (default 5)")
    quipus.add_argument("--every-lna", action = "store_true", dest = "everyLna",
                        help = "match against every LNA of the length, not only the "
                               "ones that lie in no quipu class")
    quipus.add_argument("--dedupe", action = "store_true",
                        help = "count matches up to isomorphism of the algebra, which "
                               "needs every one of them in memory and is affordable "
                               "up to order 9")
    quipus.add_argument("--keep-lines", action = "store_true", dest = "keepLines",
                        help = "keep the linearly oriented line, whose algebras are "
                               "the LNAs themselves and so match trivially")

    members = subparsers.add_parser(
        "members", help = "the quipu algebras a walk proves are in each unclassified class")
    members.add_argument("orders", type = int, nargs = "+")
    members.add_argument("--depth", type = int, default = 3,
                         help = "how far to walk out of each LNA (default 3)")
    members.add_argument("--show", type = int, default = 3, dest = "showLimit",
                         help = "how many to print per class (default 3)")

    free = subparsers.add_parser("free", help = "are two-arrow relations free on a quipu?")
    free.add_argument("orders", type = int, nargs = "+")
    free.add_argument("--min-arrows", type = int, default = 2, dest = "minArrows")
    free.add_argument("--limit", type = int, default = None,
                      help = "stop after this many ideals per order")

    args = parser.parse_args(argv)
    if args.command == "trees":
        return reportTrees(args.orders)
    if args.command == "quipus":
        return reportQuipus(args.orders, args.minArrows, args.verifyDepth,
                            args.showLimit, args.everyLna, args.keepLines, args.dedupe)
    if args.command == "members":
        return reportMembers(args.orders, args.depth, args.showLimit)
    return reportFree(args.orders, args.minArrows, args.limit)


if __name__ == "__main__":
    sys.exit(main())
