import numpy as np
from scipy.sparse import csr_matrix, isspmatrix, isspmatrix_csc, isspmatrix_csr
from _tools import csgraph_to_dense, csgraph_from_dense

DTYPE = np.float64

def validate_graph(csgraph, directed, dtype=DTYPE,
                   csr_output=True, dense_output=True,
                   copy_if_dense=False, copy_if_sparse=False):
    """Routine for validation and conversion of csgraph inputs"""
    if not (csr_output or dense_output):
        raise ValueError("Internal: dense or csr output must be true")

    # if undirected and csc storage, then transposing is quicker than
    # converting to csr.
    if (not directed) and isspmatrix_csc(csgraph):
        csgraph = csgraph.T

    if isspmatrix(csgraph):
        if csr_output:
            if copy_if_sparse and isspmatrix_csr(csgraph):
                csgraph = csr_matrix(csgraph, dtype=DTYPE)
            else:
                csgraph = csgraph.tocsr().astype(dtype)
        else:
            csgraph = csgraph_to_dense(csgraph, null_value=np.inf)
    else:
        if dense_output:
            if copy_if_dense:
                csgraph = np.array(csgraph, dtype=dtype, order='C', copy=True)
            else:
                csgraph = np.asarray(csgraph, dtype=dtype, order='C')
        else:
            csgraph = csgraph_from_dense(csgraph)

    if csgraph.ndim != 2:
        raise ValueError("compressed-sparse graph must be two dimensional")

    if csgraph.shape[0] != csgraph.shape[1]:
        raise ValueError("compressed-sparse graph must be shape (N, N)")

    return csgraph
