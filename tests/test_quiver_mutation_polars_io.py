import tempfile
import unittest

from quiver_mutation import path_algebra_class
from quiver_mutation.quiver_mutation_polars_io import (
    dataframe_to_mutation_list,
    mutation_list_to_dataframe,
    read_mutation_list_csv,
    write_mutation_list_csv,
)


class TestQuiverMutationPolarsIO(unittest.TestCase):
    def setUp(self):
        path_alg = path_algebra_class.PathAlgebra()
        path_alg.add_vertices_from([1, 2, 3])
        path_alg.add_arrows_from([(1, 2), (2, 3)])
        path_alg.add_rels_from([[[1, 2, 3]]])
        self.mutation_list = [(path_alg, [2, -1], {1: 3, 2: 1, 3: 2})]

    def test_dataframe_roundtrip_serialized(self):
        dataframe = mutation_list_to_dataframe(self.mutation_list, serialize=True)
        rebuilt = dataframe_to_mutation_list(dataframe, serialized=True)

        self.assertEqual(len(rebuilt), 1)
        rebuilt_path_alg, rebuilt_mutations, rebuilt_numbering = rebuilt[0]
        original_path_alg, original_mutations, original_numbering = self.mutation_list[0]

        self.assertEqual(list(rebuilt_path_alg.quiver.nodes), list(original_path_alg.quiver.nodes))
        self.assertEqual(list(rebuilt_path_alg.quiver.edges), list(original_path_alg.quiver.edges))
        self.assertEqual(rebuilt_path_alg.rels, original_path_alg.rels)
        self.assertEqual(rebuilt_mutations, original_mutations)
        self.assertEqual(rebuilt_numbering, original_numbering)

    def test_csv_roundtrip_serialized(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = f"{tmpdir}/mutation_data.csv"
            write_mutation_list_csv(self.mutation_list, csv_path)
            rebuilt = read_mutation_list_csv(csv_path)

        self.assertEqual(len(rebuilt), 1)
        rebuilt_path_alg, rebuilt_mutations, rebuilt_numbering = rebuilt[0]
        original_path_alg, original_mutations, original_numbering = self.mutation_list[0]

        self.assertEqual(list(rebuilt_path_alg.quiver.nodes), list(original_path_alg.quiver.nodes))
        self.assertEqual(list(rebuilt_path_alg.quiver.edges), list(original_path_alg.quiver.edges))
        self.assertEqual(rebuilt_path_alg.rels, original_path_alg.rels)
        self.assertEqual(rebuilt_mutations, original_mutations)
        self.assertEqual(rebuilt_numbering, original_numbering)


if __name__ == "__main__":
    unittest.main()
