__all__ =['inf_norm','diag_sparse']

import numpy,scipy,scipy.sparse,scipy.weave
from numpy import ravel,arange
from scipy.sparse import isspmatrix,isspmatrix_csr,isspmatrix_csc, \
                        csr_matrix,csc_matrix

def inf_norm(A):
    """
    Infinity norm of a sparse matrix (maximum absolute row sum).  This serves 
    as an upper bound on spectral radius.
    """
    
    if not isspmatrix_csr(A):
        return ValueError,'expected csr_matrix'
    
    abs_A = csr_matrix((abs(A.data),A.indices,A.indptr),dims=A.shape,check=False)
    return (abs_A * numpy.ones(A.shape[1],dtype=A.dtype)).max()

def diag_sparse(A):
    """
    If A is a sparse matrix (e.g. csr_matrix or csc_matrix)
       - return the diagonal of A as an array

    Otherwise
       - return a csr_matrix with A on the diagonal
    """
    
    if isspmatrix_csr(A) or isspmatrix_csc(A):
        n_row = len(A.indptr) - 1
        data,indices,indptr = A.data,A.indices,A.indptr

        diag = numpy.zeros(n_row,dtype=A.dtype)
        
        code = """
        #line 33 "sparse.py"

        for(int i = 0; i < n_row; i++){
          for(int jj = indptr(i); jj < indptr(i+1); jj++){
            if(indices(jj) == i){
              diag(i) = data(jj);
            }
          }
        }
        """

        err = scipy.weave.inline(code,
                                 ['data', 'indices', 'indptr', 'n_row', 'diag'],
                                 type_converters = scipy.weave.converters.blitz,
                                 compiler = 'gcc')
        return diag
    elif isspmatrix(A): 
        return ravel(array([float(A[i,i]) for i in range(min(A.shape))]))
    else:
        return csr_matrix((A,arange(len(A)),arange(len(A)+1)),(len(A),len(A)))


