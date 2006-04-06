from scipy.sparse import isspmatrix_csc, isspmatrix_csr, isspmatrix, spdiags
import _superlu

try:
    import umfpack
    isUmfpack = True
except:
    isUmfpack = False
useUmfpack = True

def use_solver( use ):
    """
    The default sparse solver is umfpack when available. This can be changed by
    passing "use = {'useUmfpack' : False}"
    which then causes the always present SuperLU based solver to be used.
    """
    for key, val in use.iteritems():
        globals()[key] = val

def _toCS_superLU( A ):
    if hasattr(A, 'tocsc') and not isspmatrix_csr( A ):
        mat = A.tocsc()
        csc = 1
    elif hasattr(A, 'tocsr'):
        mat = A.tocsr()
        csc = 0
    else:
        raise ValueError, "matrix cannot be converted to CSC/CSR"
    return mat, csc

def _toCS_umfpack( A ):
    if isspmatrix_csr( A ) or isspmatrix_csc( A ):
        mat = A
    else:
        if hasattr(A, 'tocsc'):
            mat = A.tocsc()
        elif hasattr(A, 'tocsr'):
            mat = A.tocsr()
        else:
            raise ValueError, "matrix cannot be converted to CSC/CSR"
    return mat

def spsolve(A, b, permc_spec=2):
    if not hasattr(A, 'tocsr') and not hasattr(A, 'tocsc'):
        raise ValueError, "sparse matrix must be able to return CSC format--"\
              "A.tocsc()--or CSR format--A.tocsr()"
    if not hasattr(A, 'shape'):
        raise ValueError, "sparse matrix must be able to return shape" \
                " (rows, cols) = A.shape"
    M, N = A.shape
    if (M != N):
        raise ValueError, "matrix must be square"

    if isUmfpack and useUmfpack:
        mat = _toCS_umfpack( A )

        if mat.dtype.char not in 'dD':
            raise ValueError, "convert matrix data to double, please, using"\
                  " .astype(), or set linsolve.useUmfpack = False"

        family = {'d' : 'di', 'D' : 'zi'}
        umf = umfpack.UmfpackContext( family[mat.dtype.char] )
        return umf.linsolve( umfpack.UMFPACK_A, mat, b, autoTranspose = True )

    else:
        mat, csc = _toCS_superLU( A )
        ftype, lastel, data, index0, index1 = \
               mat.ftype, mat.nnz, mat.data, mat.rowind, mat.indptr
        gssv = eval('_superlu.' + ftype + 'gssv')
        print "data-ftype: %s compared to data %s" % (ftype, data.dtype.char)
        print "Calling _superlu.%sgssv" % ftype
        return gssv(N, lastel, data, index0, index1, b, csc, permc_spec)[0]

def splu(A, permc_spec=2, diag_pivot_thresh=1.0,
         drop_tol=0.0, relax=1, panel_size=10):
    M, N = A.shape
    if (M != N):
        raise ValueError, "can only factor square matrices"

##     if isUmfpack:
##         print "UMFPACK is present - try umfpack.numeric and umfpack.solve instead!"

    csc = A.tocsc()
    gstrf = eval('_superlu.' + csc.ftype + 'gstrf')
    return gstrf(N, csc.nnz, csc.data, csc.rowind, csc.indptr, permc_spec,
                 diag_pivot_thresh, drop_tol, relax, panel_size)

def _testme():
    from scipy.sparse import csc_matrix, dok_matrix
    from numpy import transpose, array, arange
    
    print "Inverting a sparse linear system:"
    print "The sparse matrix (constructed from diagonals):"
    a = spdiags([[1, 2, 3, 4, 5], [6, 5, 8, 9, 10]], [0, 1], 5, 5)
    b = array([1, 2, 3, 4, 5])
    print "Solve: single precision complex:"
    globals()['useUmfpack'] = False
    a = a.astype('F')
    x = spsolve(a, b)
    print x
    print "Error: ", a*x-b

    print "Solve: double precision complex:"
    globals()['useUmfpack'] = True
    a = a.astype('D')
    x = spsolve(a, b)
    print x
    print "Error: ", a*x-b

    print "Solve: double precision:"
    a = a.astype('d')
    x = spsolve(a, b)
    print x
    print "Error: ", a*x-b

    print "Solve: single precision:"
    globals()['useUmfpack'] = False
    a = a.astype('f')
    x = spsolve(a, b.astype('f'))
    print x
    print "Error: ", a*x-b


if __name__ == "__main__":
    _testme()
