import numpy as np
import networkx as nx
import sympy
from sympy.matrices import Matrix, eye
from quiver_mutation_relations import numberOfPathsUpToRels

def cartanMatrix(pathAlg):
    quiv = pathAlg.quiver
    vertices = quiv.nodes
    cartanMatrix = eye(len(vertices), len(vertices))
    for i in vertices:
        for j in vertices:
            if i == j:
                cartanMatrix[j-1,i-1] = numberOfPathsUpToRels(pathAlg, i, j) + 1
            else:
                cartanMatrix[j-1,i-1] = numberOfPathsUpToRels(pathAlg, i, j)
    return cartanMatrix

def coxeterPoly(pathAlg):
    cartanMat = cartanMatrix(pathAlg)
    print(np.matrix(cartanMat))
    cartanMatInvTrans = cartanMat.inv().transpose()
    coxeterMatrix = -cartanMatInvTrans*cartanMat
    coxeterPolynomial = coxeterMatrix.charpoly()
    return coxeterPolynomial

def generateAllCoxeterPolynomials(length):

    # generate Kupisch series
    kupisch = [[[1]]]
    for i in range(1, length):
        kupisch.append([])
        for s in kupisch[i - 1]:
            for x in range(2, s[0] + 2):
                kupisch[i].append([x] + s)
    print('Number of different possible sets of relations: ',len(kupisch[length - 1]))

    polynomials = []
    for s in kupisch[length - 1]:
        M = sympy.Matrix([[1 if (m >= n and m < n + s[n]) else 0 for m in range(length)] for n in range(length)])
        C = - M * M.inv().transpose()
        p = sympy.factor(C.charpoly(sympy.Symbol("x")).as_expr())
        if p not in polynomials:
            polynomials.append(p)
    print('Coxeter polynomials: ',polynomials)
    print('Number of different Coxeter polynomials: ' ,len(polynomials))
    return polynomials

def coxPolyOfTree(tree):
    adjMat = sympy.Matrix(nx.adjacency_matrix(tree).todense(), dtype=int)
    triuMat = np.triu(Matrix(adjMat))
    mat = sympy.Matrix(triuMat, dtype=int) + sympy.eye(len(tree))#sympy.Matrix(np.triu(Matrix(matrixAsList))) + sympy.eye(n)
    print(np.matrix(mat))
    matInvTrans = mat.inv().transpose()
    coxeterMatrix = -matInvTrans * mat
    coxeterPolynomial = coxeterMatrix.charpoly()
    #print(coxeterPolynomial.as_expr())
    return coxeterPolynomial

def cartanMatrixForCanonicalAlgebra(pathAlg):
    quiv = pathAlg.quiver
    vertices = quiv.nodes
    cartanMatrix = eye(len(vertices), len(vertices))
    for i in vertices:
        for j in vertices:
            if i == j:
                cartanMatrix[j-1,i-1] = 1
            else:
                allPaths = list(nx.all_simple_paths(pathAlg.quiver, i, j))
                cartanMatrix[j-1,i-1] = min(len(allPaths), 2)
    return cartanMatrix

def coxeterPolyForCanonicalAlgebra(pathAlg):
    cartanMat = cartanMatrixForCanonicalAlgebra(pathAlg)
    print(np.matrix(cartanMat))
    cartanMatInvTrans = cartanMat.inv().transpose()
    coxeterMatrix = -cartanMatInvTrans*cartanMat
    coxeterPolynomial = coxeterMatrix.charpoly()
    return coxeterPolynomial
