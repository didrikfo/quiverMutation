# Glossary

The words this project uses, what each one means *here*, and where it is
defined or measured. Many are ordinary mathematical terms used in a narrower
sense; several are this project's own coinages, and a few mean something
different in the code from what they mean in the literature. When a term is
backed by a finding or a function, the pointer is given so the definition can
be checked rather than trusted.

Identifiers like `F-028` are entries in [`research/`](research/): `F` findings,
`H` hypotheses, `R` retractions, `E` experiments. Module names are under
[`quivermutation/`](quivermutation/).

**[Equivalences that save work](#equivalences-that-save-work)** at the end lists
every symmetry the project knows about, which parts of the code already use
each one, and which do not yet. Read it before you design a run.

---

## Contents

1. [Notation](#notation)
2. [Quivers, relations and algebras](#quivers-relations-and-algebras)
3. [Mutation](#mutation)
4. [Equivalence, classes and invariants](#equivalence-classes-and-invariants)
5. [Shapes of relations on a line](#shapes-of-relations-on-a-line)
6. [Moves](#moves)
7. [Walks and verdicts](#walks-and-verdicts)
8. [The core census](#the-core-census)
9. [Sampling long lengths](#sampling-long-lengths)
10. [Searching and deduplication](#searching-and-deduplication)
11. [Runs, ledgers and records](#runs-ledgers-and-records)
12. [Equivalences that save work](#equivalences-that-save-work)

---

## Notation

**`n`, length.** The number of vertices of the line `1 -> 2 -> ... -> n`. "Length
12" and "`A_12`" mean the same line. It is *not* the number of arrows, which is
`n - 1`.

**Relation-length row**, `relLengths`, *row*. The standard way to name an LNA:
entry `i` (1-based) is the number of arrows in the relation starting at vertex
`i`, and `0` means no relation starts there. The row has `n - 2` entries,
because a relation needs at least two arrows. `0450000` at `n = 9` has a
four-arrow relation `2 -> ... -> 6` and a five-arrow relation `3 -> ... -> 8`.
Written as a tuple in code (`(0, 4, 5, 0, 0, 0, 0)`) and as a string in names.

**Class name**, e.g. `A7_22300`. An LNA's row prefixed with its length. Classes
in a classification table are named after the LNA that seeded them.

**Relation string**, e.g. `1;2;3|3;4;5;6`. Each relation as its vertices joined by
`;`, relations separated by `|`. The key of a classification table's rows.

**Mutation sequence**, e.g. `[4, 1]`, `[-5, -5]`. The vertices to mutate at, in
order. A positive entry is a **right** mutation, a negative one a **left**
mutation. See [Mutation](#mutation).

**Quipu notation**, `P^(m)_(k)`. The quipu with main string pieces `k_0 .. k_r+1`
and cord lengths `m_1 .. m_r`, as in arXiv:2305.06642. Not unique: one quipu can
be read along its main string in several ways, and
`quipuForms.quipuParameters` picks the lexicographically smallest reading.

**`C(p_1, ..., p_t)`**. The canonical algebra of weight type `(p_1, ..., p_t)`,
e.g. `C(2,4,4)`. In the hereditary-form column it means "has the Coxeter
polynomial of that canonical algebra", which is **not** a proof of anything -- see
*hereditary form*.

**Word**, *core word*. A row with its leading and trailing zeros removed, used
to name a configuration independently of where it sits: `45`, `504`, `3344`.
See [The core census](#the-core-census).

**Offset**. Where a word is placed in a row: the number of zeros in front of it.
`45` at offset 1 in a line of 11 is the row `045000000`.

**Slide**. A word's verdicts at every offset, read from the source to the sink,
as letters: `i` inside, `o` outside, `?` undecided, `.` not run. `45` at
`n = 13` is `iooooii`.

---

## Quivers, relations and algebras

**Quiver.** A directed graph, possibly with several arrows between the same two
vertices. Held as a `networkx.MultiDiGraph`.

**Path.** A sequence of composable arrows. In `pathAlgebra` a path is written
as its vertex list `[1, 2, 3]`; in `arrowPaths` it is a tuple of
`(tail, head, key)` arrows, which is the only form that can tell two parallel
arrows apart.

**Relation.** An element of the ideal a path algebra is divided by. A **zero
relation** (monomial relation) says one path is zero, `[[1,2,3]]`; a
**commutativity relation** says two paths from the same source to the same
target are equal up to sign, `[[1,2,4],[1,3,4]]`. In the storage format `rels`
a relation is a *set* of paths with no coefficients, which cannot express
`2p - 3q` and cannot tell `p - q` from `p + q` (F-005, F-006); `arrowRels` holds
the faithful combination.

**Length of a relation.** The number of **arrows** in it, not vertices. A
"relation of two arrows" is `i -> i+1 -> i+2 = 0`.

**Path algebra with relations**, `PathAlgebra`. A quiver and an ideal of
relations. Every algebra in the project is one of these.

**Admissible.** Two uses. (1) An *admissible ideal* is one whose relations all
have length two or more and form a minimal generating set; on a line this is
exactly "the relations' starts strictly increase and their ends strictly
increase" (`lnaMoves.isAdmissible`). (2) A mutation is *admissible* at a vertex
when the procedure is allowed to run there -- see *gate*.

**LNA**, linearly oriented Nakayama algebra. The path algebra of the line
`1 -> 2 -> ... -> n` with an admissible ideal of zero relations. There are
Catalan(n-1) of them: 429 at `n = 8`, 58786 at `n = 12`, 208012 at `n = 13`.
`nakayama.LinearNakayamaAlgebra`.

**Kupisch series.** The lengths of the indecomposable projective modules, one per
vertex. Another encoding of an LNA, equivalent to the row.

**Source, sink.** On the line, vertex `1` and vertex `n`. "The ends" means
these two. "The interior" means placements far from both.

**Quipu.** A tree of maximum degree three whose degree-three vertices all lie on
one path, the **main string**; the paths hanging off it are **cords**. A quipu
quiver orients the main string one way and each cord away from it.
`nakayama.QuipuAlgebra`. The trees of type `A`, `D` and `E` are quipus.

**Quipu with relations.** A quipu quiver that carries a nonzero ideal. The
candidate family for the classes the quipu theorem misses (F-034, H-014).

**Hereditary algebra.** A path algebra with no relations. For a tree, all
orientations are derived equivalent and the underlying tree is a complete
invariant.

**Parallel arrows**, *bundle*. Two or more arrows with the same tail and head.
The mutation procedure produces them (step 1 and step 3), so the model has to
state them; a *bundle* is the set of arrows between one pair of vertices. NOTES
"Parallel arrows", F-039, R-013.

**Kronecker quiver.** Two vertices, two parallel arrows. The smallest quiver
whose Coxeter polynomial the vertex model got wrong (F-039).

**Relation dual**, *opposite algebra*, *mirror*. Reverse every arrow. On a line
this sends vertex `v` to `n + 1 - v`, a relation `s -> s + l` to
`n + 1 - s - l -> n + 1 - s`, and exchanges right mutation for left. Takes the
`45` core at offset `o` to the `504` core at offset `n - 7 - o` (F-042).
`nakayama.LinearNakayamaAlgebra.relationDual`, `lnaMoves.dualRule`. It keeps an
LNA's derived class (F-026, F-032), and the move table is closed under it, so
it is a working symmetry of every walk as well (F-026, R-011).

---

## Mutation

**Tilting mutation.** Replacing one indecomposable summand of a tilting complex
to get another, which gives a derived equivalent algebra. The project computes
it combinatorially.

**The procedure**, *steps 1-7*. The combinatorial rule of arXiv:2112.08129 that
turns a quiver with relations into its mutation at a vertex. `procedure`.
Step 1 adds composite arrows, step 3 adds an arrow per relation out of the
vertex, steps 4 and 5 produce sums and differences, step 7 is a kernel
computation (R-007). The result is not reduced; see *reduction*.

**Right mutation, left mutation.** The two directions of mutation at a vertex.
Left mutation at `i` is computed as dual, right mutation, dual
(`leftQuiverMutationAtVertex`). The two are not inverse at the *same* vertex
label, because the procedure relabels.

**Reduction**, *cleanup*. `reducePathAlgebra`: substitute out every relation that
contains a single arrow, then drop relations lying in the ideal of the others.
F-008.

**Gate**, *admissibility gate*. `mutationIsPossibleAtVertex`: whether the
procedure may be run at a vertex. It rules mutation **out**, not in: the
theorem's condition is on the algebra, not the quiver, so passing the gate does
not prove the result is derived equivalent (F-016, F-038, R-012).

**Coxeter guard.** A search refusing any step whose Coxeter key differs from the
start's, since such a step cannot have been a derived equivalence (F-038).

**Double mutation.** `proposition:doubleMutation` of arXiv:2310.08346: two left
mutations at the end `t` of a relation `r: s -> t` shift the relations crossing
`r` by one. The mechanism most of the rule table was approximating (F-032).
`doubleMutation`.

**BGP reflection.** Reversing all arrows at a source or sink of a relation-free
quiver. Right mutation at a source of a relation-free tree is exactly this, so
all orientations of a tree are one **mutation** class (F-036). `reflections`.

**Sign gauge.** Rescaling an arrow by a nonzero scalar is an automorphism of the
path algebra, so `p = 0` and `-p = 0` present the same algebra. The procedure's
output depends on the choice, so a search must quotient by it (F-050).

---

## Equivalence, classes and invariants

**Derived equivalence.** Equivalence of bounded derived categories of modules.
What the classification is up to.

**Mutation equivalence**, *mutation class*. Connected by a sequence of
mutations. Implies derived equivalence; whether the converse holds for LNAs is
open (H-012). Anything proved with the free move is a derived statement only.

**Class.** An equivalence class of LNAs of one length. Which equivalence depends
on context: a classification table's classes are mutation classes as far as a
search found them, and merged by derived invariants.

**Quipu theorem**, `thm:QuipuToAn`. arXiv:2305.06642: an LNA whose relations are
*almost separate* is derived equivalent to the hereditary algebra of an explicit
quipu. F-003 inverts it, naming the class of such an LNA in constant time.

**Quipu class.** A derived class containing a hereditary quipu algebra;
equivalently, one containing an almost separate LNA.

**Piecewise hereditary.** Derived equivalent to a hereditary abelian category:
either a hereditary algebra or a canonical algebra. An LNA that is **not**
piecewise hereditary is in no quipu class, and `piecewiseHereditary` holds two
certificates of that (propositions A9 and A13 of arXiv:2310.08346). F-012 is
how a certificate propagates to longer LNAs.

**Canonical algebra**, *canonical type*, *weight type*. Ringel's canonical
algebras `C(p_1, ..., p_t)`. **Domestic** weight types are derived equivalent to
an extended Dynkin quiver (so to a quipu class); **tubular** and **wild** ones
are not. `C(2,4,4)` at `n = 9` is tubular.

**Coxeter polynomial.** The characteristic polynomial of the Coxeter matrix
`-C^T C^{-1}`. A derived invariant, not a complete one: it fails exactly at
cospectral quipus, first at order 9 (F-010). `invariants`.

**Coxeter key.** The Coxeter polynomial as an exact integer coefficient tuple,
hashable, computed by interpolation (`invariants.coxeterCoefficients`). What
tables are keyed by.

**Coxeter matrix up to Z-conjugacy.** A finer invariant than the polynomial;
separates the cospectral quipus the polynomial cannot (F-047).

**Cartan matrix.** Entry `(i, j)` counts the paths from `i` to `j` modulo the
relations. Two parallel arrows count twice (F-039).

**Euler form.** The bilinear form of the Cartan matrix. Indefinite plus periodic
Coxeter transformation certifies not piecewise hereditary (F-048).

**Cospectral quipus.** Non-isomorphic quipus with the same adjacency spectrum,
hence the same Coxeter polynomial, but not derived equivalent. One pair at
order 9, two at 10, four at 11, thirteen at 12.

**Hereditary form.** The relation-free quiver a search reaches from an LNA,
recorded as its underlying tree. A **complete** mutation invariant when it is a
tree (F-036). The classification column also carries two other values that are
*not* the same kind of thing: `not piecewise hereditary` (a certificate that
separates from every quipu class but may not merge) and `C(...)` (read off the
polynomial; may neither merge nor separate). NOTES "Naming a class", F-018.

**Brüstle's invariant.** A derived invariant of gentle and related algebras that
classifies the derived-tame LNAs outright (F-045). Not the
Avella-Alaminos-Geiss invariant, which does not apply directly (R-008).

---

## Shapes of relations on a line

**Overlap.** Two consecutive relations `(s, l)` and `(t, m)` share
`max(0, s + l - t)` arrows. `overlap.overlapProfile` lists it for each
consecutive pair; `maxOverlap` is the largest.

**Almost separate.** Every consecutive pair of relations shares **at most one**
arrow, i.e. `maxOverlap <= 1`. The hypothesis of the quipu theorem, so "what
the theorem names" and "almost separate" are the same set of LNAs (F-021).
`overlap.isAlmostSeparate`.

**Heavily overlapping**, *heavy overlap*. Two consecutive relations sharing two
or more arrows. Exactly what puts an LNA out of the quipu theorem's reach.

**Run**, *overlapping run*. A maximal sequence of relations linked by overlaps of
at least a threshold (`overlap.overlapRuns`).

**Heavy cluster.** A run at threshold two with at least two relations: a
connected block of heavy overlaps (F-040). "Single-cluster" and "two-cluster"
words are counted this way, from the row and not the word (`5600055` is one
cluster, because the six-arrow relation reaches across the gap).

**Free gap**, *free arrow between clusters*. Two heavy clusters with an arrow
that neither covers between them. Never seen outside a quipu class at `n <= 11`
(F-040).

**Barricade.** Two heavy clusters with a two-arrow relation walled in between
them. First fits at `n = 13` (F-040, E-037).

**Spectator.** A relation inside a rule's window that the rewrite leaves
unchanged. Most anchored rules need one (F-023). `spectatorMoves`.

**Frozen.** A heavily overlapping pair in the interior, which no bounded number
of mutations near it changes the overlap of (F-022).

**Reduced**, *stripped*. An LNA with no relation of two arrows. `stripLengthTwo`
deletes them all and gives the **reduced form**, the canonical representative
of the LNA's free-move class. `freeMoves.isReduced`, `freeMoves.reducedForms`.

---

## Moves

A **move** is a rewrite of a row into another row, known to preserve the class.
Every kind below comes with a statement of *which* class it preserves.

**Rule**, *move rule*, *the rule table*. A verified local rewrite of a window of
the row, carrying a mutation sequence: mutation equivalence. Found by discovery
and kept only if `verifyMove` confirms the predicted row, every mutation
admissible and the Coxeter polynomial fixed, at several lengths (R-005, R-009).
`lnaMoves.ALL_MOVES`.

**Window.** The stretch of arrows a rule reads and rewrites. `matchesAt` insists
no relation from outside reaches into it, which keeps rules honest and makes
some true families inexpressible (F-029, F-032).

**Floating rule.** A rule that holds at every position of the line.

**Anchored rule.** A rule that holds only against the source or the sink, where a
mutation can do what it cannot in the interior. Does most of the work (F-023,
F-024).

**Dual rule.** The relation dual of a rule. A rule's dual is always a rule
(F-026), so the table is closed under it.

**Edge moves.** The doubling and collapse at an end: a relation of `l` arrows at
the source with nothing starting at `2 .. l` gains a copy at vertex 2 under
`[-(l+1), -(l+1)]`, and dually at the sink. Inexpressible as a window rule
(F-029). `edgeMoves`.

**Pair slide**, *pair collapse*, *pair to triple*. Named families in the rule
table (F-013, F-022, F-030); most are single double mutations (F-032).

**Free move.** Adding or deleting a **relation of two arrows**, anywhere it fits.
`corollary:lengthtworelations` of arXiv:2310.08346: it never changes the derived
class. It has no mutation sequence, so it proves **derived** equivalence only
(F-028; whether it is also a mutation equivalence is H-012). It holds on a line
and on nothing else (F-035). A two-arrow relation never overlaps a neighbour by
more than one arrow, so the free move never changes whether an LNA is almost
separate: its whole value is **bridging** orbits, by removing or supplying the
spectator a rule needs.

**Plain walk**, `free = True`, `--walk plain` (the sampler's default). The free move as
every walk used it until 2026-09-22: delete every two-arrow relation, never add
one. One-way, so it can reach from `2404…` what it cannot reach from `0404…`
(F-052).

**Reduced walk**, `free = freeMoves.REDUCED`, `--walk reduced`. The walk on the
quotient by the free move: every state is a reduced row, and each step is a
move out of the reduced row or out of it with **one** two-arrow relation added,
stripped again (`freeMoves.reducedMovesFrom`). An LNA and its reduced form are
the same state, so they get the same verdict. Places strictly more than the
plain walk and costs 2x to 4x as much at `n = 11` and 12; its ledgers end in
`-reduced`. E-049, F-052.

**Shared walk**, `freeMoves.SharedWalk`, `--walk shared` (the census default).
Walks many rows of one length in one process and shares what each settles: a
union-find over classes (stripped row and mirror) carries "inside" to every row
of a class, and rows of closed orbits are never expanded twice. Plain first,
reduced only if that closes. Gives the reduced walk's verdicts at a thirtieth of
its cost at `n = 11` and 12, and wants **one worker per census**: split over
workers, each re-walks the shared orbits (E-050).

**Promotion.** A shared-walk unit recorded `outside` whose class a later unit
shows to be inside. The later unit lists it in `promotes`, and `--summary`
reads it as inside.

**Class-preserving operations** of `cor:EquivNakayamaAlgebras`: dropping short
relations, the exchanges at the first and last foot, and the relation dual
(`nakayama.LinearNakayamaAlgebra.classPreservingOrbit`). Checked against the
same symmetry on the tree side (F-014).

---

## Walks and verdicts

**Orbit**, *move orbit*. The set of rows the moves reach from a start.
`freeMoves.orbitOf` / `orbitReport`. **Forward only**: the moves are not
symmetric as rewrites, so an orbit need not contain the rows that reach it
(F-007). `derivedOrbits` is the symmetric version, a union-find over every row
of a length, affordable only up to about `n = 12`.

**Closed orbit.** The walk's frontier emptied: the returned set *is* the whole
forward orbit, and absence from it means something about the move set.

**Cap**, *orbit limit*. The row count at which a walk stops. A walk that stops
on it has measured the budget, not the moves (E-037).

**Join**, `movesJoin`. Walk out of two rows at once and stop where they meet.
Meeting proves equivalence; not meeting proves nothing.

**Meeting in the middle.** The same idea for mutation searches:
`search.meetingPoints` joins two algebras when their searches pass through the
same quiver, at twice the depth for the same cost (F-041).

**Certificate.** A row, path or proof that settles a question. An almost
separate row in an orbit certifies a quipu class; a failed criterion in
`piecewiseHereditary` certifies the opposite.

**Covered**, *placed*, *seeded*. An LNA is **seeded** when the quipu theorem
names it outright, and **covered** or **placed** when the moves carry it to a
seeded one. `freeMoves.coverage`.

**Inside, outside, undecided.** The three verdicts on one row (see
`batch.py` `CoresTask`):

| verdict | meaning |
|---|---|
| `inside` | the walk reached an almost separate row, or a join met one: a certificate of a quipu class |
| `outside` | the forward orbit **closed** without one: a statement about this move set, **not** about derived equivalence |
| `undecided` | the walk hit the cap and no join met: nothing was learned about the moves |

Never collapse `undecided` into `outside`. That is E-037.

**Leftover.** A row neither the theorem nor the moves place. In a sample it is
split into **closed** and **capped** leftovers, which are different facts.

**Outside band.** The run of offsets at which a core is outside, between its
head and tail (F-042, F-051).

---

## The core census

`python batch.py cores <n>`: put one small configuration into an otherwise empty
line at every offset, and give each placement a verdict.

**Core.** A short, heavily overlapping configuration of relations, named by its
word: `45`, `504`, `3344`. A catalogue word is not necessarily a core -- most
four-letter words are not LNAs at all (E-046).

**Catalogue.** Every word the digit ranges allow (`--max-word`, `--max-arrows`,
`--pair-word`, `--gaps`), sorted by length then value. Built the same way at
every length, so `--core-limit` cuts the same words everywhere. Under the
reduced walk it holds no word with a `2`, since that word is its reduced form at
another offset. Under either walk only the first of each **mirror pair** is
asked, unless `--no-mirror` (E-049).

**Placement.** A word at an offset in a line of length `n`: one unit of a
census, `45@1`.

**Alias.** A placement whose row holds a two-arrow relation. It is the same state
of the reduced walk as its stripped row, which is always another placement of
the same catalogue: `245@0` is `45@1`, `2045@k` and `2245@k` are `45@k+2`,
`2555@k` is `555@k+1`. Under the **plain** walk the two are not the same
computation and can get different verdicts, the alias being the better placed
(E-049), so the plain catalogue keeps them.

**Mirror pair.** A placement and its relation dual: `45@o` and `504@(n-7-o)`.
One verdict for both in every case measured (E-049). The census asks whichever
comes first in catalogue order, and `--summary` fills in the other.

**Single-cluster word, pair word.** A word with one heavy cluster, or two
single-cluster words joined by `--gaps` zeros. A gap counts **zeros**, not free
arrows: gaps 1 and 2 never separate the clusters (E-047).

**Head, tail.** How many offsets at the source end and at the sink end of a slide
are inside. H-020 says they belong to the core and not to the length.

**Interior verdict**, `P(c) v(c)^m S(c)`. H-020's amended form: a slide is a fixed
prefix word, the interior verdict repeated, and a fixed suffix word.

**Rescued.** A half of a pair word that is outside on its own but inside next to
its partner (E-047).

---

## Sampling long lengths

`python batch.py sample <n>`: draw LNAs uniformly and put each through the cheap
pipeline.

**Uniform draw.** Exactly uniform over the Catalan(n-1) LNAs, by a dynamic
program over the enumerator's recursion (`sampling.sampleRelations`). Each draw
is seeded by `"<seed>/<index>"`, so a run can stop and resume, and a draw can be
re-run alone.

**The cheap pipeline**, `sampling.probe`. Theorem, then moves, then leftover.

**Leftover rate.** The fraction of draws left over. An **upper** bound on the true
rate, because the walk is forward only and capped. H-019.

**Depth.** With `--depth` above 0, a deduplicated mutation search out of each
leftover. Much the most expensive part.

---

## Searching and deduplication

**Mutation search**, `mutationSearchDepthFirst`. Enumerate mutation sequences to
a depth from one algebra, collecting the LNAs and relation-free quivers reached.
Finds a **lower bound** on each class.

**Depth.** The number of mutations in a sequence.

**Node.** One algebra visited by a search. Many nodes are the same algebra
(F-049).

**Canonical key**, `fingerprint.canonicalKey`. An exact name for an algebra
reached in a walk: vertex labels do not move under mutation, so the only
ambiguities are parallel-arrow naming and the sign gauge, and both are
quotiented out (F-049, F-050).

**Digest.** A hash of the canonical key, for holding a large visited set. A
collision costs recall, never soundness.

**Visited set**, `fingerprint.Visited`. What makes the search walk each algebra
once.

**Relation-free sighting.** A quiver a search reaches with no relations left,
recorded with its graph and path (`search.relationFreeSightings`).

**Merge.** Deciding two classes of a table are one. Allowed only on a proof: a
mutation path, a shared hereditary tree, or the free move; never on an equal
Coxeter polynomial alone (F-018).

---

## Runs, ledgers and records

**Task.** A resumable long job in `batch.py`: `sample`, `cores`, and the older
ones `python batch.py --list` names.

**Unit.** One piece of a task's work: one draw, one placement.

**Ledger.** The append-only JSON-lines file under `logs/` a task writes, one line
per finished unit. Everything that changes what a unit's **answer** means is in
the file name (orbit limit, join limit, depth, seed, walk), so two different
instruments never share a ledger. Filters (`--cores`, `--core-limit`, `--count`)
are not, so two nights can split one ledger.

**Budget**, `--budget-hours`. Stops a run cleanly; exits 2, which is what
`overnight.py` restarts on.

**Core-hours.** CPU time summed across workers. Divide by `--jobs` for wall
clock.

**`--plan`.** Prints the units a run would do and how many are left. Always
first at a new length.

**Night.** One `overnight.py` invocation, from the menu in `OVERNIGHT.md`.

**Finding, hypothesis, retraction, experiment** (`F-`, `H-`, `R-`, `E-`). See
[`research/README.md`](research/README.md). Nothing is deleted; a wrong entry is
marked and corrected in a new one.

**Status** of a hypothesis: `OPEN`, `SUPPORTED`, `CONFIRMED → F-nnn`,
`REFUTED → R-nnn`, `PARKED`.

---

## Equivalences that save work

Every symmetry below makes two things the same thing. A run that treats them as
different does the work twice, and -- as E-049 found for the free move -- may
give the two different answers, which is worse. When you build an instrument,
go down this list and say for each one whether it is quotiented out.

| equivalence | what it identifies | kind | used by | not yet used by |
|---|---|---|---|---|
| **free move** (F-028) | an LNA and the same LNA with a two-arrow relation added or removed | derived | `derivedOrbits`, `coverage`, the quipu theorem's naming (F-003), the reduced and shared walks, the census (E-049, E-050) | the plain walk, still the sampler's default; the sampler's *draws* (on purpose: the rate is over all LNAs) |
| **relation dual** (F-026) | `45` at offset `o` and `504` at offset `n - 7 - o`; any row and its mirror | derived; the moves are closed under it | the rule table (`closeUnderDual`), edge moves, double mutation, the core census (one of each mirror pair, E-049) | `derivedOrbits` and the sampler; the orbit cache F-051 proposes |
| **move orbit** (F-051) | every placement in one closed forward orbit | move set | the shared walk's closed-orbit cache (E-050) | the plain and reduced walks, and the sampler unless `--walk shared` |
| **the class of a walk** (E-050) | every row a walk passes through, in both directions | derived | the shared walk's union-find | anything run with `--jobs` above 1, which splits it per worker |
| **vertex labels** (F-049) | nothing: labels do not move under mutation | exact | `fingerprint` | -- |
| **parallel-arrow naming** (F-049) | permutations within a bundle | exact | `fingerprint.canonicalKey` | -- |
| **sign gauge** (F-050) | `p` and `-p` in a relation | exact | `fingerprint.canonicalKey` | -- |
| **tree orientation** (F-036) | every orientation of a relation-free tree | mutation | `reflections`, the hereditary form | -- |
| **class-preserving operations** (F-014) | the exchanges at the first and last foot | derived | `nakayama`, `quipuForms` | -- |

**The free move and the census, worked through.** `245` at offset 0 and `45` at
offset 1 are one LNA up to a relation of two arrows, so one derived class.
Before E-049 the census ran both, and the plain walk -- which could delete the
relation but never add it -- gave some such pairs **different** verdicts: at
`n = 11`, `2404@0` reaches an almost separate row in nine steps while `404@1`
has no move at all. The derived truth is that both are inside. The reduced walk
gives both one state, so one verdict, and its catalogue asks it once. It is
not free: the reduced walk costs more per placement than the two placements
cost together under the plain walk, so the gain is in what is placed, not in
time. The saving in time came from the next step: the **shared walk** (E-050)
remembers every class it settles, so the alias, the mirror and every row the
walk passed through are answered for free, and it reaches the reduced walk's
verdicts at a ninth of the plain census's cost.
