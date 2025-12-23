import copy
import unittest

import quiver_examples as qe
import quiver_mutation as qm


class TestQuiverMutation(unittest.TestCase):
    def test_all_rels_between_vertices_expands_relations(self):
        path_alg = copy.deepcopy(qe.A4_22)
        rels = qm.all_rels_between_vertices(path_alg, 1, 4)
        self.assertIn([[1, 2, 3, 4]], rels)

    def test_all_minimal_rels_between_vertices_scope(self):
        path_alg = copy.deepcopy(qe.A5_222)
        rels = qm.all_minimal_rels_between_vertices(path_alg, 1, 3)
        self.assertIn([[1, 2, 3]], rels)
        for rel in rels:
            for rel_path in rel:
                self.assertEqual(rel_path[0], 1)
                self.assertEqual(rel_path[-1], 3)

    def test_apply_rel_set_commutative_square(self):
        path_alg = copy.deepcopy(qe.oneCommutativeSqare)
        rel_set = path_alg.rels
        new_paths = qm.apply_rel_set_to_path([1, 2, 4], rel_set)
        self.assertIn([1, 3, 4], new_paths)
        self.assertTrue(all(path[0] == 1 and path[-1] == 4 for path in new_paths))

    def test_path_has_zero_relation(self):
        path_alg = copy.deepcopy(qe.twoSqaresOneZero)
        self.assertTrue(qm.path_has_zero_rel([1, 3, 4, 6, 7], path_alg.rels))
        self.assertFalse(qm.path_has_zero_rel([1, 2, 4, 5, 7], path_alg.rels))

    def test_remove_duplicate_rel_paths_and_rels(self):
        path_alg = copy.deepcopy(qe.oneCommutativeSqare)
        path_alg.rels[0].append(path_alg.rels[0][0][:])
        dedup_paths = qm.remove_duplicate_rel_paths(path_alg)
        self.assertEqual(len(dedup_paths.rels[0]), 2)

        path_alg = copy.deepcopy(qe.oneCommutativeSqare)
        path_alg.rels.append(copy.deepcopy(path_alg.rels[0]))
        dedup_rels = qm.remove_duplicate_rels(path_alg)
        self.assertEqual(len(dedup_rels.rels), 1)

    def test_number_of_paths_up_to_rels_line_quiver(self):
        path_alg = copy.deepcopy(qe.A5_000)
        self.assertEqual(qm.number_of_paths_up_to_rels(path_alg, 1, 5), 1)
        self.assertEqual(qm.number_of_paths_up_to_rels(path_alg, 1, 3), 1)
        self.assertEqual(qm.number_of_paths_up_to_rels(path_alg, 5, 1), 0)


if __name__ == "__main__":
    unittest.main()
