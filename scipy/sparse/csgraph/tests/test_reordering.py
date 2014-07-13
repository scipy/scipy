from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal
from scipy.sparse.csgraph import reverse_cuthill_mckee,\
        maximum_bipartite_matching
from scipy.sparse import diags, csr_matrix, csc_matrix,\
        coo_matrix

def test_graph_reverse_cuthill_mckee():
    A = np.array([[1, 0, 0, 0, 1, 0, 0, 0],
                [0, 1, 1, 0, 0, 1, 0, 1],
                [0, 1, 1, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 1, 0],
                [1, 0, 1, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 1, 0, 1],
                [0, 0, 0, 1, 0, 0, 1, 0],
                [0, 1, 0, 0, 0, 1, 0, 1]], dtype=int)
    
    graph = csr_matrix(A)
    perm = reverse_cuthill_mckee(graph)
    correct_perm = np.array([6, 3, 7, 5, 1, 2, 4, 0])
    assert_equal(perm, correct_perm)
    
    # Test int64 indices input
    graph.indices = graph.indices.astype('int64')
    graph.indptr = graph.indptr.astype('int64')
    perm = reverse_cuthill_mckee(graph, True)
    assert_equal(perm, correct_perm)


def test_graph_maximum_bipartite_matching():
    A = diags(np.ones(25), offsets=0, format='csc')
    rand_perm = np.random.permutation(25)
    rand_perm2 = np.random.permutation(25)

    Rrow = np.arange(25)
    Rcol = rand_perm
    Rdata = np.ones(25,dtype=int)
    Rmat = coo_matrix((Rdata,(Rrow,Rcol))).tocsc()

    Crow = rand_perm2
    Ccol = np.arange(25)
    Cdata = np.ones(25,dtype=int)
    Cmat = coo_matrix((Cdata,(Crow,Ccol))).tocsc()
    # Randomly permute identity matrix
    B = Rmat*A*Cmat
    
    # Row permute
    perm = maximum_bipartite_matching(B,perm_type='row')
    Rrow = np.arange(25)
    Rcol = perm
    Rdata = np.ones(25,dtype=int)
    Rmat = coo_matrix((Rdata,(Rrow,Rcol))).tocsc()
    C1 = Rmat*B
    
    # Column permute
    perm2 = maximum_bipartite_matching(B,perm_type='column')
    Crow = perm2
    Ccol = np.arange(25)
    Cdata = np.ones(25,dtype=int)
    Cmat = coo_matrix((Cdata,(Crow,Ccol))).tocsc()
    C2 = B*Cmat
    
    # Should get identity matrix back
    assert_equal(any(C1.diagonal() == 0), False)
    assert_equal(any(C2.diagonal() == 0), False)
    
    # Test int64 indices input
    B.indices = B.indices.astype('int64')
    B.indptr = B.indptr.astype('int64')
    perm = maximum_bipartite_matching(B,perm_type='row')
    Rrow = np.arange(25)
    Rcol = perm
    Rdata = np.ones(25,dtype=int)
    Rmat = coo_matrix((Rdata,(Rrow,Rcol))).tocsc()
    C3 = Rmat*B
    assert_equal(any(C3.diagonal() == 0), False)
    