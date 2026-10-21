"""Rerun of workshop/rounds/073/theorist_cycD.py at the current tree, where
procedure.mutateAtVertex refuses a cyclic D(vertex) (CyclicDescendantsError).
Wraps mutateAtVertex with allowCyclic=True (E-211 behaviour), then runs the 073 script unchanged."""
import functools, runpy, sys
sys.path.insert(0, ".")
from quivermutation import procedure as pr
_orig = pr.mutateAtVertex
pr.mutateAtVertex = functools.partial(_orig, allowCyclic=True)
runpy.run_path("workshop/rounds/073/theorist_cycD.py", run_name="__main__")
