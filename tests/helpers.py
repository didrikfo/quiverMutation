"""Small helpers shared by the test modules."""

import contextlib
import io

import sympy

import pathAlgebraClass
import quiverMutation as qm

LAMBDA = sympy.Symbol("lambda")


def quiet(func, *args, **kwargs):
    """Run func, swallowing the progress chatter the library prints to stdout."""
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


def path_algebra(arrows, rels=()):
    """Build a PathAlgebra from a list of arrows and a list of relations."""
    pa = pathAlgebraClass.PathAlgebra()
    pa.add_arrows_from([list(a) for a in arrows])
    pa.add_rels_from([[list(p) for p in rel] for rel in rels])
    return pa


def line_algebra(length, rel_lengths):
    """An LNA on the linear quiver 1 -> 2 -> ... -> length.

    rel_lengths is the string/sequence notation used throughout the repo: entry
    i (0-based) is the number of arrows in the relation starting at vertex i+1,
    or 0 for no relation there.
    """
    if isinstance(rel_lengths, str):
        rel_lengths = [int(c) for c in rel_lengths]
    return qm.lineQuiverExample(length, list(rel_lengths))


def arrow_set(path_alg):
    """Arrows as a multiset-free set of (source, target) pairs."""
    return {(a[0], a[1]) for a in path_alg.quiver.edges}


def rel_set(path_alg):
    """Relations as a hashable, order-independent structure."""
    return {tuple(sorted(tuple(p) for p in rel)) for rel in path_alg.rels}


def coxeter_poly(path_alg):
    """Coxeter polynomial as a sympy expression in lambda."""
    return quiet(qm.coxeterPoly, path_alg).as_expr()


def dynkin_A_coxeter(n):
    """Coxeter polynomial of the path algebra of the Dynkin quiver A_n."""
    return sympy.expand(sum(LAMBDA ** i for i in range(n + 1)))


def dynkin_D_coxeter(n):
    """Coxeter polynomial of the path algebra of the Dynkin quiver D_n."""
    return sympy.expand((LAMBDA ** (n - 1) + 1) * (LAMBDA + 1))


def relation_string(rel_lengths):
    """Per-vertex relation lengths -> the 'a;b;c|d;e;f' string used in tables."""
    if isinstance(rel_lengths, str):
        rel_lengths = [int(c) for c in rel_lengths]
    paths = [
        ";".join(str(v) for v in range(start + 1, start + n + 2))
        for start, n in enumerate(rel_lengths)
        if n
    ]
    return "|".join(paths)
