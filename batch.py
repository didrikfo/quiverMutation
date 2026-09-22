#!/usr/bin/env python
"""Long runs, named and resumable, so a machine can be left working on them.

    python batch.py --list                        what there is to run
    python batch.py sample 14 --count 2000        probe 2000 LNAs of length 14
    python batch.py sample 14 --count 2000 --jobs 7 --budget-hours 9
    python batch.py sample 14 --summary           what the ledger says, no work
    python batch.py cores 15 --jobs 7             slide every core along a length

Every task writes an append-only ledger under `logs/`, one line per finished
unit, and every run resumes from it: the same command again does what is left
and nothing that is done.  `--budget-hours` stops cleanly between units and
exits 2, which is what `overnight.py` restarts on.

The tasks here are the ones that are *only* worth running long.  The three
existing long jobs keep their own front doors, because each has parameters that
do not fit a common shape, and `--list` says what they are.
"""

import argparse
import itertools
import sys

from quivermutation import freeMoves as fm
from quivermutation import jobs
from quivermutation import sampling


class SampleTask(jobs.Task):
    """Draw LNAs of an intractable length and see what the cheap pipeline says.

    Above about n = 13 there are too many LNAs to classify them all -- 208012 at
    n = 13 and 1767263190 at n = 20 -- and the short lengths that *can* be done
    completely are unrepresentative: in a quiver of length 8 every vertex is
    within three arrows of an end.  So the question a long length can be asked
    is a statistical one, and this asks it: what fraction of LNAs of length n
    does the quipu theorem name, what fraction do the moves carry to one it
    names, and what is left over.

    The leftovers are the point.  At n = 9, 10 and 11 they are a handful of
    orbits whose merging is research H-013, and nobody has seen whether their
    rate falls, holds or rises with the length, or whether at n = 16 they take
    shapes that n = 9 has no room for.  Each drawn leftover is recorded with its
    overlap profile, so the shapes can be counted afterwards.

    With `--depth` above 0 each leftover also gets a deduplicated mutation
    search, and what it reaches is recorded.  That is much the most expensive
    part, which is why it is off unless asked for.
    """

    name = 'sample'
    help = "draw LNAs of a long length and probe them"

    def addArguments(self, parser):
        parser.add_argument("length", type = int,
                            help = "the line length to draw from")
        parser.add_argument("--count", type = int, default = 1000,
                            help = "how many LNAs to draw (default 1000)")
        parser.add_argument("--seed", type = int, default = 0,
                            help = "the run's seed; a different seed is a "
                                   "different sample and a different ledger")
        parser.add_argument("--depth", type = int, default = 0,
                            help = "also search this deep from every leftover "
                                   "(default 0, meaning do not search)")
        parser.add_argument("--orbit-limit", type = int, default = 20000,
                            dest = "orbitLimit",
                            help = "how far to walk a move orbit before calling "
                                   "the row a leftover (default 20000)")
        parser.add_argument("--walk", choices = WALKS, default = "plain",
                            help = "how the free move is walked: 'plain' only "
                                   "deletes a two-arrow relation, as every run "
                                   "before 2026-09-22 did; 'reduced' treats an LNA "
                                   "and its stripped form as one state and may add "
                                   "one as well -- it places more, and costs about "
                                   "four times as much at n = 12 (E-049) "
                                   "(default plain)")

    def ledgerPath(self, args):
        # The depth is in the name because a run with a search and a run without
        # are different work on the same draws, and a ledger must not claim a
        # unit was done to a depth it was not.  The orbit limit is in it for the
        # stronger reason: it decides whether a draw is recorded as a leftover,
        # so two runs under different limits disagree about the answer and not
        # merely about the effort, and one must not resume from the other.
        # The walk is in it for the same reason: the reduced walk places rows the
        # plain one leaves over (E-049), so the two disagree about the answer.
        return "logs/sample-n{0}-s{1}-d{2}-o{3}{4}.jsonl".format(
            args.length, args.seed, args.depth, args.orbitLimit, _walkSuffix(args))

    def units(self, args):
        return ["{0}/{1}/{2}".format(args.length, args.seed, index)
                for index in range(args.count)]

    def run(self, unit, args):
        index = int(unit.rsplit("/", 1)[1])
        relLengths = sampling.drawFor(args.length, args.seed, index)
        record = sampling.probe(args.length, relLengths,
                                orbitLimit = args.orbitLimit,
                                free = _freeFor(args))
        record['index'] = index
        if args.depth > 0 and record['settledBy'] == 'leftover':
            record['search'] = _searchFrom(args.length, relLengths, args.depth)
        return record

    def summarise(self, records, args, out = sys.stdout):
        import collections
        results = [record['result'] for record in records]
        if not results:
            print("nothing in the ledger yet", file = out)
            return
        tally = collections.Counter(result['settledBy'] for result in results)
        total = len(results)
        print("n = {0}, seed {1}: {2} of {3} drawn".format(
            args.length, args.seed, total, args.count), file = out)
        print("  out of {0} LNAs of this length\n".format(
            sampling.countLNAs(args.length)), file = out)
        for settledBy in ('theorem', 'moves', 'leftover'):
            count = tally.get(settledBy, 0)
            print("  {0:<9} {1:6d}  {2:5.1f}%  +- {3:.1f}".format(
                settledBy, count, 100.0 * count / total,
                100.0 * _standardError(count, total)), file = out)

        leftovers = [result for result in results if result['settledBy'] == 'leftover']
        if not leftovers:
            print("\n  no leftovers drawn, so nothing to say about their shape",
                  file = out)
            return

        # Two different facts wear the one name.  A leftover whose orbit
        # closed is a statement about the moves; one whose walk hit the cap
        # is a statement about the budget, and E-037 is what it costs to
        # read the second as the first.  Rows written before the split
        # existed carry no `orbitClosed` and are counted as neither.
        closed = sum(1 for result in leftovers
                     if result.get("orbitClosed") is True)
        capped = sum(1 for result in leftovers
                     if result.get("orbitClosed") is False)
        print("\n  of the {0} leftovers:".format(len(leftovers)), file = out)
        print("    orbit closed  {0:6d}   the whole forward orbit held none"
              .format(closed), file = out)
        print("    orbit capped  {0:6d}   the walk ran out at {1} rows"
              .format(capped, args.orbitLimit), file = out)
        if len(leftovers) - closed - capped:
            print("    not recorded  {0:6d}   written before the split existed"
                  .format(len(leftovers) - closed - capped), file = out)
        if capped:
            print("    a rate read off the total above is part rate and part"
                  " cap;" + chr(10) + "    a larger --orbit-limit moves the line",
                  file = out)
        print("\n  leftovers by maximum relation overlap:", file = out)
        shapes = collections.Counter(result['maxOverlap'] for result in leftovers)
        for overlap in sorted(shapes):
            print("    overlap {0}: {1}".format(overlap, shapes[overlap]), file = out)
        print("\n  leftovers by number of relations:", file = out)
        counts = collections.Counter(result['relations'] for result in leftovers)
        for relations in sorted(counts):
            print("    {0} relations: {1}".format(relations, counts[relations]),
                  file = out)
        print("\n  the first few, by name:", file = out)
        for result in leftovers[:12]:
            print("    {0}  overlap {1}, orbit {2}".format(
                result['name'], result['maxOverlap'], result['orbit']), file = out)


class CoresTask(jobs.Task):
    """Slide an overlapping core along a long line, and ask where it is placeable.

    The sampler above asks what a *typical* LNA of a long length does.  This asks
    the opposite question, and it is the one the short lengths cannot be asked at
    all: take a small configuration of overlapping relations -- a **core** -- put
    it at every position of a long line with nothing else in the quiver, and see
    whether the moves carry it to an almost separate LNA from there.

    That is how F-042 was found.  Sliding the `45` core along the line at
    `n <= 13` gave a clean law -- the moves place it when it sits against the
    source or within one arrow of the sink, and nowhere else -- and the law is
    invisible at `n = 9`, where the core fits at three offsets and all three are
    inside.  H-018 says outright that what would settle its phrasing is the same
    census at `n = 14` and `n = 15`, over *every* core rather than the handful
    that were slid by hand.

    Two shapes have never been testable at any length the repo has run.  F-040
    counted them: at `n <= 11` not one LNA outside a quipu class has two heavy
    clusters with a relation-free stretch between them, and a **barricade** --
    two clusters walling in a two-arrow relation -- needs 4 + 1 + 2 + 1 + 4
    arrows and so first fits at `n = 13`.  A length of 15 has room for both, and
    `--gaps` is what puts them in the catalogue.

    **Why a failed walk is not recorded as a failure.**  `freeMoves.orbitOf`
    walks forwards only and stops at a row cap, and E-037 is the experiment that
    burned on exactly that: 49 barricades at `n = 15` and `n = 16` went down as
    outside, having merely run out of budget, and `movesJoin` then joined them in
    45 seconds walking from both ends.  So a placement here gets one of three
    verdicts and never two:

    * `inside`    -- the walk reached an almost separate LNA, or a two-ended walk
                     met one.  A certificate; the class is a quipu class.
    * `outside`   -- the forward orbit **closed** under the moves without holding
                     one.  A real statement about this move set, and still not a
                     statement about derived equivalence.
    * `undecided` -- the orbit hit the cap and the joins did not meet.  The
                     budget was measured, not the moves.  These are what to spend
                     a second night on, not what to draw a conclusion from.
    """

    name = 'cores'
    help = "slide overlapping cores along a long line, at every offset"

    def addArguments(self, parser):
        parser.add_argument("length", type = int,
                            help = "the line length to place the cores in")
        parser.add_argument("--max-word", type = int, default = 3,
                            dest = "maxWord",
                            help = "how many vertices a single core may span "
                                   "(default 3)")
        parser.add_argument("--max-arrows", type = int, default = 6,
                            dest = "maxArrows",
                            help = "the longest relation a core may hold "
                                   "(default 6)")
        parser.add_argument("--gaps", default = "1,2,3",
                            help = "the free-arrow gaps to try between two cores, "
                                   "comma separated; the empty string leaves "
                                   "pairs out (default 1,2,3)")
        parser.add_argument("--pair-word", type = int, default = 2,
                            dest = "pairWord",
                            help = "how many vertices each half of a pair may "
                                   "span (default 2)")
        parser.add_argument("--orbit-limit", type = int, default = 20000,
                            dest = "orbitLimit",
                            help = "how far to walk the forward orbit before "
                                   "calling the placement undecided (default 20000)")
        parser.add_argument("--join-limit", type = int, default = 6000,
                            dest = "joinLimit",
                            help = "rows per side for each two-ended join "
                                   "(default 6000)")
        parser.add_argument("--walk", choices = WALKS, default = "plain",
                            help = "how the free move is walked: 'plain' only "
                                   "deletes a two-arrow relation, as every run "
                                   "before 2026-09-22 did; 'reduced' treats an LNA "
                                   "and its stripped form as one state and may add "
                                   "one as well -- it places more, and costs about "
                                   "four times as much at n = 12 (E-049) "
                                   "(default plain)")
        parser.add_argument("--cores", default = "",
                            help = "run only these core words, comma separated; "
                                   "the ledger is the same one, so this is a "
                                   "way to do part of a census first")
        parser.add_argument("--core-limit", type = int, default = 0,
                            dest = "coreLimit",
                            help = "run only the first this many core words of "
                                   "the catalogue (default 0, meaning all)")
        parser.add_argument("--no-mirror", dest = "mirror", action = "store_false",
                            help = "ask both a row and its relation dual, rather "
                                   "than one and reading the other in the mirror")

    def ledgerPath(self, args):
        # Every parameter that changes what a unit *means* is in the name: the
        # two limits decide a verdict, so a ledger written under one must not be
        # read as though it answered under another.
        # The walk is a third: the reduced walk places rows the plain one calls
        # outside (E-049).  A plain ledger keeps the name it always had.
        return "logs/cores-n{0}-w{1}p{2}a{3}g{4}-o{5}j{6}{7}.jsonl".format(
            args.length, args.maxWord, args.pairWord, args.maxArrows,
            args.gaps.replace(",", "") or "none", args.orbitLimit, args.joinLimit,
            _walkSuffix(args))

    def units(self, args):
        return ["{0}@{1}".format(word, offset)
                for word, offset in _placements(args)]

    def run(self, unit, args):
        word, offset = unit.split("@")
        return _verdictFor(args.length, word, int(offset),
                           orbitLimit = args.orbitLimit,
                           joinLimit = args.joinLimit,
                           free = _freeFor(args))

    def summarise(self, records, args, out = sys.stdout):
        import collections
        results = [record['result'] for record in records]
        if not results:
            print("nothing in the ledger yet", file = out)
            return
        byWord = collections.defaultdict(dict)
        for result in results:
            byWord[result['core']][result['offset']] = result['verdict']
        # A mirrored census asks one row of each dual pair; the other half of
        # every slide is the same verdict read in the mirror.  Only filled in
        # where the ledger has nothing of its own for that placement.
        if getattr(args, 'mirror', False):
            for result in list(results):
                row = tuple(int(letter) for letter in result['name'])
                word, offset = _placementOf(_mirror(args.length, row))
                byWord[word].setdefault(offset, result['verdict'])
        tally = collections.Counter(result['verdict'] for result in results)
        print("n = {0}: {1} placements of {2} cores".format(
            args.length, len(results), len(byWord)), file = out)
        for verdict in ('inside', 'outside', 'undecided'):
            print("  {0:<10} {1}".format(verdict, tally.get(verdict, 0)), file = out)
        print("\n  offsets read left to right from the source; "
              "i inside, o outside, ? undecided", file = out)

        # The interesting cores are the ones whose verdict *moves* with the
        # offset: a core that is inside everywhere says nothing the theorem does
        # not, and one that is outside everywhere is a property of the core.  A
        # core that changes is a property of the placement, which is H-018.
        moving = {word: offsets for word, offsets in byWord.items()
                  if len(set(offsets.values()) - {'undecided'}) > 1}
        settled = {word: offsets for word, offsets in byWord.items()
                   if set(offsets.values()) == {'outside'}}
        print("\n  cores whose verdict changes with where they sit "
              "({0}):".format(len(moving)), file = out)
        print("    core           slide          head tail  interior",
              file = out)
        for word in sorted(moving, key = _coreSortKey):
            slide = _slideOf(moving[word])
            head, tail, interior = _headAndTail(slide)
            print("    {0:<14} {1:<18} {2:>4} {3:>4}  {4}".format(
                word, slide, head, tail, interior), file = out)
        print("\n    head is how many offsets from the source are inside, tail how many"
              "\n    from the sink.  Both are read off this run and neither is a"
              "\n    claim; what is worth comparing between two lengths is whether"
              "\n    they change.  A core whose interior holds an inside is the"
              "\n    one to look at first.", file = out)
        print("\n  cores outside at every offset tried ({0}):".format(len(settled)),
              file = out)
        for word in sorted(settled, key = _coreSortKey)[:40]:
            print("    {0}".format(_offsetLine(settled[word], word)), file = out)

        # F-040's count, asked at a length with room for the shape.  At n <= 11
        # not one LNA outside a quipu class had two heavy clusters with a
        # relation-free stretch between them.  Anything printed here is a row
        # that does, and that the moves did not place.
        separated = [result for result in results
                     if result.get('heavyClusters', 0) >= 2 and result.get('freeGap')]
        placed = collections.Counter(result['verdict'] for result in separated)
        print("\n  rows with two heavy clusters and a free arrow between them: "
              "{0}".format(len(separated)), file = out)
        for verdict in ('inside', 'outside', 'undecided'):
            print("    {0:<10} {1}".format(verdict, placed.get(verdict, 0)), file = out)
        if separated and not placed.get('outside'):
            print("    none outside -- F-040's count holding at a length "
                  "with room for the shape", file = out)
        for result in sorted((result for result in separated
                              if result['verdict'] == 'outside'),
                             key = lambda result: result['name'])[:40]:
            print("    {0:<16} core {1} at {2}, orbit {3} closed".format(
                result['name'], result['core'], result['offset'],
                result.get('orbit')), file = out)


WALKS = ("reduced", "plain")


def _freeFor(args):
    """The `free` argument a walk takes, from `--walk`."""
    return fm.REDUCED if getattr(args, 'walk', 'plain') == 'reduced' else True


def _walkSuffix(args):
    """What the walk adds to a ledger's name: nothing for the plain walk, so the
    ledgers written before the reduced one existed keep their names."""
    return "-reduced" if _freeFor(args) is fm.REDUCED else ""


def _slideOf(offsets):
    """One core's verdicts across its offsets, as `i` / `o` / `?` / `.`.

    Indexed by the offset itself and not by how many are in the ledger, so
    that a partial run reads as a gap rather than as a shifted table --
    which is the one way this display could quietly lie, since the whole
    finding it is meant to show is *which* offsets are which.
    """
    letters = {"inside": "i", "outside": "o", "undecided": "?"}
    return "".join(letters.get(offsets.get(offset), ".")
                   for offset in range(max(offsets) + 1))


def _offsetLine(offsets, word):
    """The core's name and its slide, for the lists that print nothing else."""
    return "{0:<14} {1}".format(word, _slideOf(offsets))


def _headAndTail(slide):
    """How many inside offsets at each end of a slide, and whether the rest are not.

    Returned for reading and not for concluding: at one length these two
    numbers describe one table and nothing more.  They are printed because
    the thing worth comparing between two lengths is whether they *change*,
    and that comparison is unreadable off two rows of letters of different
    widths.  An undecided counts as neither inside nor outside, so it ends a
    head and does not clear an interior, which is the cautious way round.
    """
    head = 0
    while head < len(slide) and slide[head] == "i":
        head += 1
    tail = 0
    while tail < len(slide) - head and slide[len(slide) - 1 - tail] == "i":
        tail += 1
    interior = slide[head:len(slide) - tail]
    if "i" in interior:
        reading = "HOLDS AN INSIDE"
    elif set(interior) - {"o"}:
        reading = "not finished"
    else:
        reading = "all outside"
    return head, tail, reading


def _coreSortKey(word):
    return (len(word), word)


def _placements(args):
    """Every (core, offset) this run covers, in a stable order.

    Cheap and deterministic: `jobs.runTask` calls this before doing any work and
    diffs it against the ledger, so it must not walk any orbits and must give the
    same list on a resumed run as on the first one.
    """
    # Under the reduced walk a word holding a two-arrow relation is not a core
    # of its own: it is the word without it, placed one or more vertices along,
    # and that placement is already in the catalogue -- `245` at 0 is `45` at 1.
    # They are one state of the walk, so asking both is the same work twice
    # (E-049: 489 of the 3034 placements at n = 16).
    withTwos = _freeFor(args) is not fm.REDUCED
    cores = _singleCores(args.maxWord, args.maxArrows, withTwos)
    gaps = [int(piece) for piece in args.gaps.split(",") if piece.strip()]
    if gaps:
        atoms = _singleCores(args.pairWord, args.maxArrows, withTwos)
        for gap in gaps:
            for first in atoms:
                for second in atoms:
                    cores.append(first + "0" * gap + second)
    cores = _selected(cores, args)
    placements = []
    for word in cores:
        for offset in range(0, max(0, args.length - 2)):
            if _rowFor(args.length, word, offset) is not None:
                placements.append((word, offset))
    if getattr(args, 'mirror', False):
        # The relation dual keeps the class and the move set is closed under
        # it (F-026), so a row and its mirror get one verdict: `45` at `o` and
        # `504` at `n - 7 - o` are one question.  Ask the smaller row of each
        # pair; `summarise` writes the answer in under both.  1370 pairs at
        # n = 11 and 12, under both walks, agreed without exception (E-049).
        # Of the two, the one earlier in catalogue order is asked, so that a
        # `--core-limit` prefix asks exactly what the whole census would.
        order = {_rowFor(args.length, word, offset): rank
                 for rank, (word, offset) in enumerate(placements)}
        placements = [(word, offset) for rank, (word, offset) in enumerate(placements)
                      if order.get(_mirror(args.length,
                                           _rowFor(args.length, word, offset)),
                                   rank) >= rank]
    return placements


def _mirror(length, row):
    """The relation dual of a row: vertex `v` goes to `n + 1 - v`."""
    mirrored = [0] * len(row)
    for position, arrows in enumerate(row):
        if arrows:
            start = position + 1
            mirrored[length - start - arrows] = arrows
    return tuple(mirrored)


def _placementOf(row):
    """(word, offset) for a row: the word is the row without its outer zeros."""
    nonzero = [index for index, value in enumerate(row) if value]
    first, last = nonzero[0], nonzero[-1]
    return "".join(str(value) for value in row[first:last + 1]), first


def _selected(cores, args):
    """The catalogue narrowed by `--cores` and `--core-limit`, in catalogue order.

    Both are filters and neither is in the ledger's name, which is deliberate:
    a placement means the same thing however the run was narrowed, so a night
    that does the first sixty cores and a night that does the rest write the
    same ledger and the second does not redo the first.  Anything that changes
    what a *verdict* means -- the two limits -- is in the name instead.

    The catalogue is built the same way at every length, so `--core-limit 60` at
    n = 12 and at n = 16 ask about the same sixty cores.  That is what makes two
    part-finished censuses comparable, and the second overnight run is why it is
    here: it left n = 13 at 62 cores and n = 15 at 125, and only the overlap
    could be read.
    """
    wanted = [piece.strip() for piece in getattr(args, 'cores', '').split(",")
              if piece.strip()]
    if wanted:
        known = set(cores)
        missing = [word for word in wanted if word not in known]
        aliases = [word for word in missing if "2" in word]
        if aliases and _freeFor(args) is fm.REDUCED:
            raise ValueError(
                "{0}: a relation of two arrows is free, so under the reduced "
                "walk this is the same state as the word without it, one or "
                "more offsets along -- ask for that word, or pass "
                "--walk plain".format(", ".join(aliases)))
        if missing:
            raise ValueError(
                "not in this catalogue: {0} -- widen --max-word/--max-arrows/"
                "--gaps, or check the spelling".format(", ".join(missing)))
        chosen = set(wanted)
        cores = [word for word in cores if word in chosen]
    limit = getattr(args, 'coreLimit', 0)
    if limit:
        cores = cores[:limit]
    return cores


def _singleCores(maxWord, maxArrows, withTwos = True):
    """The core words: short, heavily overlapping, no relation shorter than two.

    `withTwos = False` leaves out every word holding a relation of two arrows,
    which is the catalogue the reduced walk needs: such a word is its stripped
    form at another offset (see `_placements`).

    A core is a relation-length word with no leading or trailing zero whose
    relations overlap somewhere in **two or more arrows** -- which is exactly the
    condition that puts an LNA outside the quipu theorem's reach, so a word
    without it has nothing to ask.  Sorted, so the catalogue is stable and a
    resumed run asks for the same units in the same order.
    """
    from quivermutation import overlap as ov

    found = []
    for size in range(2, maxWord + 1):
        for code in itertools.product(range(0, maxArrows + 1), repeat = size):
            if code[0] == 0 or code[-1] == 0 or 1 in code:
                continue
            if not withTwos and 2 in code:
                continue
            if ov.maxOverlap(list(code)) < 2:
                continue
            found.append("".join(str(value) for value in code))
    return sorted(set(found), key = _coreSortKey)


def _rowFor(length, word, offset):
    """The LNA of that length holding this core at that offset, or `None`.

    `None` where the placement is not an LNA at all: a relation has to fit inside
    the line, and the enumeration's own condition is that the relations' starts
    **and** their ends both strictly increase.  A word that overruns the sink, or
    whose second relation ends no later than its first, is not a row, and is left
    out of the catalogue rather than failing in a worker.
    """
    row = [0] * max(0, length - 2)
    if offset < 0 or offset + len(word) > len(row):
        return None
    ends = []
    for position, letter in enumerate(word):
        arrows = int(letter)
        if arrows == 0:
            continue
        if arrows < 2:
            return None
        start = offset + position + 1
        if start + arrows > length:
            return None
        if ends and start + arrows <= ends[-1]:
            return None
        ends.append(start + arrows)
        row[offset + position] = arrows
    return tuple(row)


def _verdictFor(length, word, offset, orbitLimit, joinLimit, free = True):
    """Place the core, walk out of it, and say which of the three verdicts holds.

    Cheapest first, and the second step exists only because the first one's
    failure is ambiguous -- see the class docstring and E-037.
    """
    from quivermutation import freeMoves as fm, overlap as ov

    row = _rowFor(length, word, offset)
    if row is None:
        raise ValueError("{0} at offset {1} is not an LNA of length {2}".format(
            word, offset, length))
    clusters, freeGap = _clusterShape(row)
    record = {
        'length': length,
        'core': word,
        'offset': offset,
        'name': "".join(str(value) for value in row),
        'maxOverlap': ov.maxOverlap(list(row)),
        'relations': sum(1 for value in row if value),
        # The word is what was *asked for*; these are what the row actually is.
        # They differ, and the difference is the point: `5600055` looks like two
        # clusters with three free vertices between them and is one cluster,
        # because the six-arrow relation reaches across the gap.  F-040 counts
        # clusters, so a census that counted zeros in the word would not be
        # answering F-040's question.
        'heavyClusters': clusters,
        'freeGap': freeGap,
        'orbitLimit': orbitLimit,
        'joinLimit': joinLimit,
        'walk': 'reduced' if free == fm.REDUCED else 'plain',
    }
    if ov.isAlmostSeparate(length, row):
        # Not reachable from `_singleCores`, which only builds heavy words, but a
        # caller with its own catalogue would hit it and the answer is free.
        record.update(verdict = 'inside', by = 'theorem', orbit = 1,
                      certificate = record['name'])
        return record

    # Stop the walk at the first almost separate row instead of enumerating the
    # orbit to the cap and then looking through it.  The second overnight run
    # measured what the difference is worth: of its 932 placements the 666
    # `inside` ones cost 49 of the 54 core-hours, and every one of a re-run
    # sample of twelve found its certificate inside the first 351 rows of an
    # orbit the walk had taken to 20000.  The verdicts are unchanged; what
    # changes is that a census of a length now fits in a night instead of four.
    walk = fm.orbitReport(length, row, free = free, edges = True, doubles = True,
                          limit = orbitLimit,
                          stopWhen = lambda member: ov.isAlmostSeparate(length, member))
    # `orbit` is rows *walked*, which is the orbit's size only when it closed.
    # `stopped` says which, and is absent from the ledger rows written before
    # the walk stopped early -- where `orbit` meant something else.
    record['orbit'] = len(walk.rows)
    record['stopped'] = walk.stoppedBy
    if walk.found is not None:
        record.update(verdict = 'inside', by = 'orbit',
                      certificate = "".join(str(value) for value in walk.found))
        return record

    if walk.closed:
        # The frontier emptied, so this is the whole forward orbit and the
        # absence is the move set's doing and not the budget's.
        record.update(verdict = 'outside', by = 'closed orbit', certificate = None)
        return record

    targets = _separatedTargets(length, row)
    record['targetsTried'] = ["".join(str(value) for value in target)
                              for target in targets]
    for target in targets:
        meeting = fm.movesJoin(length, row, target, free = free, edges = True,
                               doubles = True, limit = joinLimit)
        if meeting is not None:
            record.update(verdict = 'inside', by = 'join',
                          certificate = "".join(str(value) for value in target),
                          meetsAt = "".join(str(value) for value in meeting))
            return record
    record.update(verdict = 'undecided', by = 'orbit capped', certificate = None)
    return record


def _clusterShape(row):
    """How many heavy clusters the row has, and whether a free arrow separates two.

    A *heavy cluster* is F-040's: a maximal run of relations linked by overlaps
    of two arrows or more, with at least two relations in it, which is exactly
    what puts an LNA outside the quipu theorem's reach.  F-040 counted them over
    every LNA of lengths 8 to 11 and found that **not one** LNA outside a quipu
    class has two of them with a relation-free stretch between them -- and noted
    that the shape barely fits at those lengths.  `freeGap` is the flag that
    makes the same count readable here, at a length with room.

    A relation covers arrows `start .. start + arrows - 1`, so two clusters have
    a free arrow between them when the second starts at least two arrows past
    where the first stops covering.
    """
    from quivermutation import overlap as ov

    runs = [run for run in ov.overlapRuns(list(row)) if len(run) >= 2]
    if len(runs) < 2:
        return len(runs), False
    freeGap = False
    for earlier, later in zip(runs, runs[1:]):
        stops = max(start + arrows - 1 for start, arrows in earlier)
        begins = min(start for start, _arrows in later)
        if begins - stops >= 2:
            freeGap = True
    return len(runs), freeGap


def _separatedTargets(length, row):
    """Almost separate rows near this one, to walk towards from both ends.

    Nothing here claims a target is equivalent to the row -- `movesJoin` is what
    would prove that, and failing to meet one proves nothing either way.  They
    are candidates, and the four are the four smallest edits that take a heavy
    overlap away: shorten the relation on either side of it until it is gone, or
    delete the relation on either side outright.  A relation shortened below two
    arrows is deleted instead, there being no such relation.
    """
    from quivermutation import overlap as ov

    ordered = []
    for shrinkLeft in (True, False):
        for delete in (True, False):
            candidate = _separate(length, row, shrinkLeft, delete)
            if candidate is None or candidate == tuple(row):
                continue
            if candidate in ordered:
                continue
            if ov.isAlmostSeparate(length, candidate):
                ordered.append(candidate)
    return ordered


def _separate(length, row, shrinkLeft, delete):
    """One edit of `_separatedTargets`, applied until no heavy overlap is left."""
    from quivermutation import overlap as ov

    current = list(row)
    # Every pass takes at least one arrow off some relation, so the total number
    # of arrows bounds the loop; the bound is there to stop a future change to
    # the edits turning a non-convergence into a hang in a worker.
    for _attempt in range(sum(row) + len(row) + 1):
        profile = ov.overlapProfile(current)
        if not profile or max(profile) < 2:
            return tuple(current)
        heavy = next(index for index, value in enumerate(profile) if value >= 2)
        positions = [index for index, value in enumerate(current) if value]
        target = positions[heavy] if shrinkLeft else positions[heavy + 1]
        if delete or current[target] <= 2:
            current[target] = 0
        else:
            current[target] -= 1
    return None


def _searchFrom(length, relLengths, depth):
    """A deduplicated mutation search from one LNA, and what it reached.

    Deduplicated because without it this is unaffordable: the walk revisits the
    same algebra along many routes, 11.4x at n = 9 and depth 6, and the factor
    roughly doubles per level (research E-042).
    """
    from quivermutation import fingerprint, lnaMoves as lm, nakayama as nk, search

    algebra = nk.LinearNakayamaAlgebra(length, list(relLengths))
    visited = fingerprint.Visited()
    collected, hereditary = [], []
    search.mutationSearchDepthFirst(algebra, depth, [], 'sample',
                                    printOutput = False, collected = collected,
                                    collectedHereditary = hereditary,
                                    visited = visited)
    reached = set()
    for pathAlg, _path, _numbering in collected:
        row = lm.asRelLengths(lm._copy(pathAlg), length)
        if row is not None:
            reached.add(''.join(str(value) for value in row))
    return {
        'depth': depth,
        'reached': sorted(reached),
        'hereditary': sorted({quipu for _form, quipu, _path in hereditary if quipu}),
        'walk': visited.summarise(),
    }


def _standardError(count, total):
    """The standard error of a proportion, for reading a sample honestly.

    A run that draws 2000 of 1767263190 LNAs and reports "3.1% leftovers" is
    reporting 3.1% plus or minus 0.4, and the difference between that and a
    later run's 2.6% is not a trend.  Printing the error is cheaper than
    explaining afterwards that it was never there.
    """
    if total <= 0:
        return 0.0
    proportion = count / total
    return (proportion * (1 - proportion) / total) ** 0.5


TASKS = {task.name: task for task in [SampleTask(), CoresTask()]}


#: The long jobs that keep their own front door, with the command that runs
#: them.  They are listed here so that `--list` is the whole inventory and a
#: session does not have to remember which script does what.
ELSEWHERE = [
    ("classify", "python classify.py 10 --resume --budget-hours 9",
     "classify a whole length up to derived equivalence"),
    ("merges", "python merges.py 10 --depths 5 6 7 8 --jobs 7 --budget-hours 9",
     "search for mutation paths between the orbits the moves leave over"),
    ("overlaps", "python overlaps.py 10 --free --doubles --no-rules",
     "which relation-overlap configurations a length leaves unplaced"),
    ("discover", "python discover.py --max-arrows 7 --max-width 8",
     "search for new rewrite rules"),
]


def main(argv = None):
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--list", action = "store_true", dest = "listTasks",
                        help = "print the tasks and exit")
    subparsers = parser.add_subparsers(dest = "task")
    for task in TASKS.values():
        sub = subparsers.add_parser(task.name, help = task.help,
                                    description = task.__doc__,
                                    formatter_class = argparse.RawDescriptionHelpFormatter)
        task.addArguments(sub)
        sub.add_argument("--jobs", type = int, default = 1,
                         help = "worker processes (default 1)")
        sub.add_argument("--budget-hours", type = float, default = None,
                         dest = "budgetHours",
                         help = "stop cleanly once this long has passed, "
                                "exiting 2 so a wrapper can restart it")
        sub.add_argument("--summary", action = "store_true",
                         help = "print what the ledger establishes and stop")
        sub.add_argument("--plan", action = "store_true",
                         help = "print the ledger, the size of the run "
                                "and how much of it is left, and stop")

    args = parser.parse_args(argv)
    if args.listTasks or args.task is None:
        _printList()
        return 0

    task = TASKS[args.task]
    if args.plan:
        # Sizing a night is the one thing that has to be answerable
        # without doing any of it, and `units` is cheap by contract.
        ledger = jobs.Ledger(task.ledgerPath(args))
        units = task.units(args)
        done = ledger.done()
        left = [unit for unit in units if unit not in done]
        print("ledger  {0}".format(ledger.path))
        print("units   {0}".format(len(units)))
        print("done    {0}".format(len(units) - len(left)))
        print("left    {0}".format(len(left)))
        return 0
    if args.summary:
        ledger = jobs.Ledger(task.ledgerPath(args))
        task.summarise(ledger.records(), args)
        return 0
    return jobs.runTask(task, args, budgetHours = args.budgetHours, jobs = args.jobs)


def _printList():
    print("Tasks that run here:\n")
    for task in TASKS.values():
        print("  {0:<10} {1}".format(task.name, task.help))
        print("             python batch.py {0} --help".format(task.name))
    print("\nLong jobs with their own front door:\n")
    for name, command, description in ELSEWHERE:
        print("  {0:<10} {1}".format(name, description))
        print("             {0}".format(command))
    print("\n`overnight.py` runs a set of these for a fixed number of hours and "
          "restarts\nwhatever dies.  Every one of them resumes from its own "
          "checkpoint.")


if __name__ == "__main__":
    sys.exit(main())
