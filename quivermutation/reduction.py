"""Cleaning up after a mutation.

The procedure's steps 1-7 stop short of the paper's own "Note": a relation
containing a path of length one says that arrow equals a combination of longer
paths, so the arrow and the relation both go and every other relation has that
arrow substituted out.  What is left after that has to be cut down to a minimal
set of generators.

`reducePathAlgebra` is both, and the work is `procedure.reduce`, which does them
on linear combinations of paths -- the substitution as arithmetic, the
minimality by deciding ideal membership exactly.  What used to be here was six
passes over the set-of-paths model that approximated the same thing by looking
for syntactic containment, and iterated to a fixed point because no single pass
was right; research F-008 and R-003 are what that cost.
"""

from . import procedure


def reducePathAlgebra(pathAlg):
    """The cleanup, to a fixed point.  Returns a new algebra."""
    return procedure.toPathAlgebra(*procedure.reduce(
        pathAlg.quiver, procedure.relationsFrom(pathAlg)))
