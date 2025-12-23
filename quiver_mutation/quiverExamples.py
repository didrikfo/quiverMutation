import networkx as nx
from . import pathAlgebraClass

# Automatically generate examples of line quivers with different sets of relations by using the function
# lineQuiverExample(lineLength , relationList )  where relationList is a list of length (lineLength - 2)
# whose entry i is the length of the relation starting at vertex i (0 if no relation starts at i).
# E.g.:
# A_6 = lineQuiverExample(lineLength = 6, relationList = [0000])
# A_6/Rad^3 = lineQuiverExample(lineLength = 6, relationList = [3333])
# etc.

# Alternatively, you can manually define the vertices, arrows and relations in your quiver.

# - Vertices are given as integers.
# - An arrows is given by a tuple containing 3 integers: (startVertex, endVertex, key).
#     - key is used to distinguish multiple arrows between the same pair of vertices. It is 0 if there is only one such
#       arrow, and it can usually be ignored.
# - A path is a list which contains all the vertices that path touches, in order.
# - A relation is given as a list of paths with the same start and end point, so a list of lists.
#   This is interpreted as the sum of these paths equalling zero. Coefficients are not supported.

# Below are a bunch of examples, both of lines and of more complicated quivers.
# Note that there are several different ways to construct the same quiver, use whichever you want.




A4_30 = pathAlgebraClass.PathAlgebra()
A4_30.add_vertices_from([1, 2, 3, 4])            # add_vertices_from() takes a list of integers, and adds them as vertices
A4_30.add_arrows_from([[1, 2], [2, 3], [3, 4]])  # add_arrows_from() takes a list of lists of 2 vertices, and adds an arrow for each list.
A4_30.add_rels_from([ [[1,2,3,4]] ])             # add_rels_from() takes a list of lists of lists, and adds a relation for each list of lists.

# It is also possible to add entire paths at once, using the function add_paths_from(). This will automatically add all
# vertices that appear in the path, and add an arrow for each pair of consecutive vertices in the path. It is possible
# to add more than one path to/from/through the same vertex, but adding two paths throguh the same arrow will add another
# copy of that arrow to the quiver.
A4_22 = pathAlgebraClass.PathAlgebra()
A4_22.add_paths_from([ [1,2,3,4] ]) # add_paths_from() takes a list of lists of vertices
A4_22.add_rels_from([ [[1,2,3]],[[2,3,4]] ])

A5_000 = pathAlgebraClass.PathAlgebra()
A5_000.add_vertices_from([1, 2, 3, 4, 5])
A5_000.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5]])

A5_222 = pathAlgebraClass.PathAlgebra()
A5_222.add_vertices_from([1, 2, 3, 4, 5])
A5_222.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5]])
A5_222.add_rels_from([ [[1,2,3]], [[2,3,4]], [[3,4,5]] ])

A5_230 = pathAlgebraClass.PathAlgebra()
A5_230.add_vertices_from([1, 2, 3, 4, 5])
A5_230.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5]])
A5_230.add_rels_from([ [[1,2,3]], [[2,3,4,5]] ])

A6_0000 = pathAlgebraClass.PathAlgebra()
A6_0000.add_vertices_from([1, 2, 3, 4, 5, 6])
A6_0000.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])

A6_2230 = pathAlgebraClass.PathAlgebra()
A6_2230.add_vertices_from([1, 2, 3, 4, 5, 6])
A6_2230.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
A6_2230.add_rels_from([[[1, 2, 3]], [[2, 3, 4]], [[3, 4, 5, 6]]])

A6_3030 = pathAlgebraClass.PathAlgebra()
A6_3030.add_vertices_from([1, 2, 3, 4, 5, 6])
A6_3030.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
A6_3030.add_rels_from([[[1, 2, 3, 4]], [[3, 4, 5, 6]]])

A7_00000 = pathAlgebraClass.PathAlgebra()
A7_00000.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_00000.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_00000.add_rels_from([])

A7_00400 = pathAlgebraClass.PathAlgebra()
A7_00400.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_00400.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_00400.add_rels_from([[[3, 4, 5, 6, 7]]])

A7_03030 = pathAlgebraClass.PathAlgebra()
A7_03030.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_03030.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_03030.add_rels_from([[[2, 3, 4, 5]], [[4, 5, 6, 7]]])

A7_04000 = pathAlgebraClass.PathAlgebra()
A7_04000.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_04000.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_04000.add_rels_from([[[2, 3, 4, 5, 6]]])

A7_05000 = pathAlgebraClass.PathAlgebra()
A7_05000.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_05000.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_05000.add_rels_from([[[2, 3, 4, 5, 6, 7]]])

A7_22230 = pathAlgebraClass.PathAlgebra()
A7_22230.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_22230.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_22230.add_rels_from([[[1, 2, 3]], [[2, 3, 4]], [[3, 4, 5]], [[4, 5, 6, 7]]])

A7_22300 = pathAlgebraClass.PathAlgebra()
A7_22300.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_22300.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_22300.add_rels_from([[[1, 2, 3]], [[2, 3, 4]], [[3, 4, 5, 6]]])

A7_22400 = pathAlgebraClass.PathAlgebra()
A7_22400.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_22400.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_22400.add_rels_from([[[1, 2, 3]], [[2, 3, 4]], [[3, 4, 5, 6, 7]]])

A7_23030 = pathAlgebraClass.PathAlgebra()
A7_23030.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_23030.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_23030.add_rels_from([[[1, 2, 3]], [[2, 3, 4, 5]], [[4, 5, 6, 7]]])

A7_23400 = pathAlgebraClass.PathAlgebra()
A7_23400.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_23400.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_23400.add_rels_from([[[1, 2, 3]], [[2, 3, 4, 5]], [[3, 4, 5, 6, 7]]])

A7_24000 = pathAlgebraClass.PathAlgebra()
A7_24000.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_24000.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_24000.add_rels_from([[[1, 2, 3]], [[2, 3, 4, 5, 6]]])

A7_30000 = pathAlgebraClass.PathAlgebra()
A7_30000.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_30000.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_30000.add_rels_from([[[1, 2, 3, 4]]])

A7_30300 = pathAlgebraClass.PathAlgebra()
A7_30300.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_30300.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_30300.add_rels_from([[[1, 2, 3, 4]], [[3, 4, 5, 6]]])

A7_30302 = pathAlgebraClass.PathAlgebra()
A7_30302.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_30302.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_30302.add_rels_from([[[1, 2, 3, 4]], [[3, 4, 5, 6]], [[5, 6, 7]]])

A7_30330 = pathAlgebraClass.PathAlgebra()
A7_30330.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_30330.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_30330.add_rels_from([[[1, 2, 3, 4]], [[3, 4, 5, 6]], [[4, 5, 6, 7]]])

A7_30400 = pathAlgebraClass.PathAlgebra()
A7_30400.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_30400.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_30400.add_rels_from([[[1, 2, 3, 4]], [[3, 4, 5, 6, 7]]])

A7_33030 = pathAlgebraClass.PathAlgebra()
A7_33030.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_33030.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_33030.add_rels_from([[[1, 2, 3, 4]],[[2, 3, 4, 5]], [[4, 5, 6, 7]]])

A7_33330 = pathAlgebraClass.PathAlgebra()
A7_33330.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_33330.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_33330.add_rels_from([[[1, 2, 3, 4]],[[2, 3, 4, 5]], [[3, 4, 5, 6]], [[4, 5, 6, 7]]])

A7_40030 = pathAlgebraClass.PathAlgebra()
A7_40030.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_40030.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_40030.add_rels_from([[[1, 2, 3, 4, 5]], [[4, 5, 6, 7]]])

A7_40400 = pathAlgebraClass.PathAlgebra()
A7_40400.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
A7_40400.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]])
A7_40400.add_rels_from([[[1, 2, 3, 4, 5]], [[3, 4, 5, 6, 7]]])

A8_222300 = pathAlgebraClass.PathAlgebra()
A8_222300.add_vertices_from([1, 2, 3, 4, 5, 6, 7, 8])
A8_222300.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]])
A8_222300.add_rels_from([[[1, 2, 3]], [[2, 3, 4]], [[3, 4, 5]], [[4, 5, 6, 7]]])

A8_333330 = pathAlgebraClass.PathAlgebra()
A8_333330.add_vertices_from([1, 2, 3, 4, 5, 6, 7, 8])
A8_333330.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]])
A8_333330.add_rels_from([[[1, 2, 3, 4]], [[2, 3, 4, 5]], [[3, 4, 5, 6]], [[4, 5, 6, 7]], [[5, 6, 7, 8]]])

A8_444400 = pathAlgebraClass.PathAlgebra()
A8_444400.add_vertices_from([1, 2, 3, 4, 5, 6, 7, 8])
A8_444400.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]])
A8_444400.add_rels_from([[[1, 2, 3, 4, 5]], [[2, 3, 4, 5, 6]], [[3, 4, 5, 6, 7]], [[4, 5, 6, 7, 8]]])

A8_700000 = pathAlgebraClass.PathAlgebra()
A8_700000.add_vertices_from([1, 2, 3, 4, 5, 6, 7, 8])
A8_700000.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]])
A8_700000.add_rels_from([[[1, 2, 3, 4, 5, 6, 7, 8]]])

A9_2222300 = pathAlgebraClass.PathAlgebra()
A9_2222300.add_vertices_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
A9_2222300.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8,9]])
A9_2222300.add_rels_from([[[1, 2, 3]], [[2, 3, 4]], [[3, 4, 5]], [[4, 5, 6]], [[5, 6, 7, 8]]])

A9_3345000 = pathAlgebraClass.PathAlgebra()
A9_3345000.add_vertices_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
A9_3345000.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8,9]])
A9_3345000.add_rels_from([[[1, 2, 3, 4]], [[2, 3, 4, 5]], [[3, 4, 5, 6, 7]], [[4, 5, 6, 7, 8,9]]])


A9_3030230 = pathAlgebraClass.PathAlgebra()
A9_3030230.add_vertices_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
A9_3030230.add_arrows_from([[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8,9]])
A9_3030230.add_rels_from([[[1, 2, 3, 4]], [[3, 4, 5, 6]], [[5, 6, 7]], [[6, 7, 8, 9]]])


oneCommutativeSqare = pathAlgebraClass.PathAlgebra()
oneCommutativeSqare.add_vertices_from([1, 2, 3, 4])
oneCommutativeSqare.add_arrows_from([[1,2], [1,3], [2,4], [3,4]])
oneCommutativeSqare.add_rels_from([[[1,2,4],[1,3,4]]])

twoSqares = pathAlgebraClass.PathAlgebra()
twoSqares.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
twoSqares.add_arrows_from([[1,2], [1,3], [2,4], [3,4], [4,5], [4,6], [5,7], [6,7]])
twoSqares.add_rels_from([[[1,2,4],[1,3,4]],[[4,5,7],[4,6,7]]])

twoSqaresOneZero = pathAlgebraClass.PathAlgebra()
twoSqaresOneZero.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
twoSqaresOneZero.add_arrows_from([[1,2], [1,3], [2,4], [3,4], [4,5], [4,6], [5,7], [6,7]])
twoSqaresOneZero.add_rels_from([[[1,2,4],[1,3,4]],[[4,5,7],[4,6,7]],[[3,4,6]]])

testingQuiver1 = pathAlgebraClass.PathAlgebra()
testingQuiver1.add_vertices_from([1,2,3,4,5])
testingQuiver1.add_arrows_from([[1,2],[1,5],[2,3],[3,4],[4,5]])
testingQuiver1.add_rels_from([[[2,3,4]]])

testingQuiver2 = pathAlgebraClass.PathAlgebra()
testingQuiver2.add_vertices_from([1, 2, 3, 4, 5])
testingQuiver2.add_arrows_from([(1, 5, 0), (2, 3, 0), (3, 1, 0), (3, 4, 0), (4, 5, 0)])
testingQuiver2.add_rels_from([[[3, 1, 5], [3, 4, 5]], [[2, 3, 1]]])

testingQuiver3 = pathAlgebraClass.PathAlgebra()
testingQuiver3.add_vertices_from([1, 2, 3, 4, 5])
testingQuiver3.add_arrows_from([(1, 2, 0), (2, 5, 0), (4, 3, 0), (5, 4, 0)])
testingQuiver3.add_rels_from([[[5, 4, 3]], [[2, 5, 4]]])

testingQuiver4 = pathAlgebraClass.PathAlgebra()
testingQuiver4.add_vertices_from([1, 2, 3, 4, 5, 6])
testingQuiver4.add_arrows_from([(1, 2, 0), (1, 3, 0), (1, 4, 0), (2, 5, 0), (3, 5, 0), (4, 5, 0), (5, 6, 0)])
testingQuiver4.add_rels_from([[[1, 2, 5], [1, 3, 5], [1, 4, 5]], [[2, 5, 6]], [[3, 5, 6]], [[4, 5, 6]]])

testingQuiver5 = pathAlgebraClass.PathAlgebra()
testingQuiver5.add_vertices_from([1, 2, 3, 4, 5, 6])
testingQuiver5.add_arrows_from([(1, 6, 0), (2, 1, 0), (3, 4, 0), (4, 2, 0), (4, 5, 0), (5, 6, 0)])
testingQuiver5.add_rels_from([[[4, 2, 1, 6], [4, 5, 6]], [[3, 4, 2]]])

testingQuiver6 = pathAlgebraClass.PathAlgebra()
testingQuiver6.add_vertices_from([1, 2, 3, 4, 5, 6])
testingQuiver6.add_arrows_from([(1, 2, 0), (2, 3, 0), (4, 6, 0), (5, 2, 0), (6, 5, 0), (6, 1, 0)])
testingQuiver6.add_rels_from([[[6, 1, 2], [6, 5, 2]], [[5, 2, 3]], [[4, 6, 5]], [[1, 2, 3]], [[4, 6, 1]]])

testingQuiver7 = pathAlgebraClass.PathAlgebra()
testingQuiver7.add_vertices_from([1, 2, 3, 4, 5, 6])
testingQuiver7.add_arrows_from([(1, 2, 0), (1, 5, 0), (2, 3, 0), (3, 4, 0), (4, 6, 0), (5, 4, 0)])
testingQuiver7.add_rels_from([[[1, 2, 3, 4], [1, 5, 4]], [[3, 4, 6]]])

testingQuiver8 = pathAlgebraClass.PathAlgebra()
testingQuiver8.add_vertices_from([1, 2, 3, 4, 5, 6])
testingQuiver8.add_arrows_from([(1, 2, 0), (1, 5, 0), (2, 3, 0), (3, 4, 0), (4, 6, 0), (5, 4, 0)])
testingQuiver8.add_rels_from([[[1, 2, 5], [1, 4, 5]], [[2, 3, 6], [2, 5, 6]], [[1, 2, 3, 6], [1, 4, 5, 6]]])

testingQuiver9 = pathAlgebraClass.PathAlgebra()
testingQuiver9.add_vertices_from([1, 2, 3, 4, 5, 6])
testingQuiver9.add_arrows_from([(1, 2, 0), (2, 3, 0), (2, 4, 0), (2, 5, 0), (3, 6, 0), (4, 6, 0), (5, 6, 0), (6, 7, 0)])
testingQuiver9.add_rels_from([[[1, 2, 5, 6, 7]], [[2, 3, 6], [2, 4, 6], [2, 5, 6]], [[2, 3, 6,], [2, 5, 6]], [[2, 4, 6], [2, 5, 6]]])

testingQuiver10 = pathAlgebraClass.PathAlgebra()
testingQuiver10.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
testingQuiver10.add_arrows_from([(1, 2, 0), (2, 6, 0), (2, 7, 0), (4, 3, 0), (5, 3, 0), (6, 5, 0), (7, 4, 0), (7, 5, 0)])
testingQuiver10.add_rels_from([[[2, 6, 5], [2, 7, 5]], [[2, 7, 4, 3]], [[7, 4, 3], [7, 5, 3]] ])

testingQuiver11 = pathAlgebraClass.PathAlgebra()
testingQuiver11.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
testingQuiver11.add_arrows_from([(1, 6, 0), (1, 7, 0), (2, 1, 0), (4, 3, 0), (5, 3, 0), (6, 5, 0), (7, 5, 0), (7, 4, 0)])
testingQuiver11.add_rels_from([[[1, 6, 5], [1, 7, 5]], [[2, 1, 6]], [[2, 1, 7]], [[7, 4, 3], [7, 5, 3]]])

testingQuiver12 = pathAlgebraClass.PathAlgebra()
testingQuiver12.add_vertices_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
testingQuiver12.add_arrows_from([(1, 9, 0), (2, 7, 0), (3, 5, 0), (4, 1, 0), (5, 4, 0), (5, 6, 0), (6, 9, 0), (8, 7, 0), (9, 8, 0), (9, 2, 0)])
testingQuiver12.add_rels_from([[[1, 9, 2]], [[3, 5, 4]], [[3, 5, 6, 9]], [[5, 4, 1, 9], [5, 6, 9]], [[6, 9, 2]], [[6, 9, 8]], [[9, 2, 7], [9, 8, 7]]])

testingQuiver13 = pathAlgebraClass.PathAlgebra()
testingQuiver13.add_vertices_from([1, 2, 3, 4, 5])
testingQuiver13.add_arrows_from([(1, 2, 0), (2, 3, 0), (2, 4, 0), (3, 5, 0), (4, 5, 0)])
testingQuiver13.add_rels_from([[[1, 2, 3]], [[2, 3, 5], [2, 4, 5]]])

testingQuiver14 = pathAlgebraClass.PathAlgebra()
testingQuiver14.add_vertices_from([1, 2, 3, 4, 5])
testingQuiver14.add_arrows_from([(1, 4, 0), (3, 2, 0), (4, 2, 0), (2, 5, 0)])
testingQuiver14.add_rels_from([[[1, 4, 2]]])

testingQuiver15 = pathAlgebraClass.PathAlgebra()
testingQuiver15.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
testingQuiver15.add_arrows_from([(1, 7, 0), (2, 1, 0), (3, 2, 0), (3, 5, 0), (5, 1, 0), (5, 6, 0), (6, 7, 0), (7, 4, 0)])
testingQuiver15.add_rels_from([[[3, 2, 1], [3, 5, 1]], [[5, 1, 7], [5, 6, 7]], [[2, 1, 7]], [[5, 6, 7, 4]]])

testingQuiver16 = pathAlgebraClass.PathAlgebra()
testingQuiver16.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
testingQuiver16.add_arrows_from([(1, 5, 0), (2, 1, 0), (3, 2, 0), (3, 6, 0), (5, 7, 0), (5, 4, 0), (6, 5, 0), (7, 4, 0)])
testingQuiver16.add_rels_from([[[3, 2, 1, 5], [3, 6, 5]], [[6, 5, 4], [6, 5, 7, 4]], [[1, 5, 4]], [[2, 1, 5, 7]]])

testingQuiver17 = pathAlgebraClass.PathAlgebra()
testingQuiver17.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
testingQuiver17.add_arrows_from([(1, 5, 0), (2, 1, 0), (3, 2, 0), (3, 6, 0), (5, 7, 0), (5, 4, 0), (6, 5, 0), (7, 4, 0)])
testingQuiver17.add_rels_from([[[3, 2, 1, 5], [3, 6, 5]], [[6, 5, 4], [6, 5, 7, 4]], [[1, 5, 4], [1, 5, 7, 4]], [[2, 1, 5, 7]]])

testingQuiver18 = pathAlgebraClass.PathAlgebra()
testingQuiver18.add_vertices_from([1, 2, 3, 4, 5, 6, 7])
testingQuiver18.add_arrows_from([(1, 5, 0), (2, 1, 0), (3, 2, 0), (3, 6, 0), (5, 7, 0), (6, 5, 0), (7, 4, 0)])
testingQuiver18.add_rels_from([[[3, 2, 1, 5], [3, 6, 5]], [[2, 1, 5, 7]]])
