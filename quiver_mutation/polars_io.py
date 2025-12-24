import json
from typing import Any, Dict, Iterable, List, Tuple

import polars as pl

from . import path_algebra_class


def path_algebra_to_dict(path_alg: path_algebra_class.PathAlgebra) -> Dict[str, Any]:
    return {
        "vertices": list(path_alg.quiver.nodes),
        "arrows": list(path_alg.quiver.edges),
        "relations": path_alg.rels,
    }


def dict_to_path_algebra(data: Dict[str, Any]) -> path_algebra_class.PathAlgebra:
    path_alg = path_algebra_class.PathAlgebra()
    vertices = data.get("vertices", [])
    arrows = data.get("arrows", [])
    relations = data.get("relations", [])
    path_alg.add_vertices_from(vertices)
    path_alg.add_arrows_from(arrows)
    path_alg.add_rels_from(relations)
    return path_alg


def _serialize(value: Any) -> str:
    return json.dumps(value, sort_keys=True)


def _deserialize(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, str):
        return json.loads(value)
    return value


def mutation_list_to_dataframe(
    mutation_list: Iterable[Tuple[path_algebra_class.PathAlgebra, List[int], Dict[int, int]]],
    serialize: bool = False,
) -> pl.DataFrame:
    rows = []
    for path_alg, mutation_vertices, vertex_relabeling in mutation_list:
        row = {
            "mutations": mutation_vertices,
            "numbering": vertex_relabeling,
            "vertices": list(path_alg.quiver.nodes),
            "arrows": list(path_alg.quiver.edges),
            "relations": path_alg.rels,
        }
        if serialize:
            row = {key: _serialize(value) for key, value in row.items()}
        rows.append(row)
    return pl.DataFrame(rows)


def _coerce_numbering_keys(numbering: Dict[Any, Any]) -> Dict[int, Any]:
    return {int(key): value for key, value in numbering.items()}


def dataframe_to_mutation_list(
    dataframe: pl.DataFrame,
    serialized: bool | None = None,
) -> List[Tuple[path_algebra_class.PathAlgebra, List[int], Dict[int, int]]]:
    if serialized is None:
        serialized = any(
            dtype == pl.Utf8
            for dtype in dataframe.schema.values()
        )
    mutation_list = []
    for row in dataframe.iter_rows(named=True):
        if serialized:
            row = {key: _deserialize(value) for key, value in row.items()}
        path_alg = path_algebra_class.PathAlgebra()
        path_alg.add_vertices_from(row.get("vertices", []))
        path_alg.add_arrows_from([tuple(arrow) for arrow in row.get("arrows", [])])
        path_alg.add_rels_from(row.get("relations", []))
        mutation_vertices = row.get("mutations", [])
        vertex_relabeling = row.get("numbering", {})
        if isinstance(vertex_relabeling, dict):
            vertex_relabeling = _coerce_numbering_keys(vertex_relabeling)
        mutation_list.append((path_alg, mutation_vertices, vertex_relabeling))
    return mutation_list


def write_mutation_dataframe_csv(dataframe: pl.DataFrame, file_name: str) -> None:
    dataframe.write_csv(file_name)


def read_mutation_dataframe_csv(file_name: str) -> pl.DataFrame:
    return pl.read_csv(file_name)


def write_mutation_list_csv(
    mutation_list: Iterable[Tuple[path_algebra_class.PathAlgebra, List[int], Dict[int, int]]],
    file_name: str,
) -> None:
    dataframe = mutation_list_to_dataframe(mutation_list, serialize=True)
    dataframe.write_csv(file_name)


def read_mutation_list_csv(
    file_name: str,
) -> List[Tuple[path_algebra_class.PathAlgebra, List[int], Dict[int, int]]]:
    dataframe = pl.read_csv(file_name)
    return dataframe_to_mutation_list(dataframe, serialized=True)
