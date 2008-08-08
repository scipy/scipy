from warnings import warn

from numpy import asarray
from scipy.sparse import isspmatrix_csc, isspmatrix_csr, isspmatrix, \
        SparseEfficiencyWarning, csc_matrix

import _superlu

noScikit = False
try:
    import scikits.umfpack as umfpack
except ImportError:
    import umfpack
    noScikit = True

isUmfpack = hasattr( umfpack, 'UMFPACK_OK' )

if isUmfpack and noScikit:
    warn( 'scipy.sparse.linalg.dsolve.umfpack will be removed,'
          ' install scikits.umfpack instead', DeprecationWarning )

useUmfpack = True


__all__ = [ 'use_solver', 'spsolve', 'splu', 'factorized' ]

#convert numpy char to superLU char
superLU_transtabl = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}


def use_solver( **kwargs ):
    """
    Valid keyword arguments with defaults (other ignored):
      useUmfpack = True
      assumeSortedIndices = False

    The default sparse solver is umfpack when available. This can be changed by
    passing useUmfpack = False, which then causes the always present SuperLU
    based solver to be used.

    Umfpack requires a CSR/CSC matrix to have sorted column/row indices. If
    sure that the matrix fulfills this, pass assumeSortedIndices=True
    to gain some speed.
    """
    if 'useUmfpack' in kwargs:
        globals()['useUmfpack'] = kwargs['useUmfpack']

    if isUmfpack:
        umfpack.configure( **kwargs )


def spsolve(A, b, permc_spec=2):
    """Solve the sparse linear system Ax=b
    """
    if isspmatrix( b ):
        b = b.toarray()

    if b.ndim > 1:
        if max( b.shape ) == b.size:
            b = b.squeeze()
        else:
            raise ValueError, "rhs must be a vector (has shape %s)" % (b.shape,)

    if not (isspmatrix_csc(A) or isspmatrix_csr(A)):
        A = csc_matrix(A)
        warn('spsolve requires CSC or CSR matrix format', SparseEfficiencyWarning)

    A.sort_indices()
    A = A.asfptype()  #upcast to a floating point format

    M, N = A.shape
    if (M != N):
        raise ValueError, "matrix must be square (has shape %s)" % (M,N)
    if M != b.size:
        raise ValueError, "matrix - rhs size mismatch (%s - %s)"\
              % (A.shape, b.size)


    if isUmfpack and useUmfpack:
        if A.dtype.char not in 'dD':
            raise ValueError, "convert matrix data to double, please, using"\
                  " .astype(), or set linsolve.useUmfpack = False"

        b = asarray(b, dtype=A.dtype).reshape(-1)

        family = {'d' : 'di', 'D' : 'zi'}
        umf = umfpack.UmfpackContext( family[A.dtype.char] )
        return umf.linsolve( umfpack.UMFPACK_A, A, b,
                             autoTranspose = True )

    else:
        if isspmatrix_csc(A):
            flag = 1 # CSC format
        else:
            flag = 0 # CSR format

        ftype = superLU_transtabl[A.dtype.char]

        gssv = eval('_superlu.' + ftype + 'gssv')
        b = asarray(b, dtype=A.dtype)

        return gssv(N, A.nnz, A.data, A.indices, A.indptr, b, flag, permc_spec)[0]

def splu(A, permc_spec=2, diag_pivot_thresh=1.0,
         drop_tol=0.0, relax=1, panel_size=10):
    """
    A linear solver, for a sparse, square matrix A, using LU decomposition where
    L is a lower triangular matrix and U is an upper triagular matrix.

    Returns a factored_lu object. (scipy.sparse.linalg.dsolve._superlu.SciPyLUType)

    See scipy.sparse.linalg.dsolve._superlu.dgstrf for more info.
    """

    if not isspmatrix_csc(A):
        A = csc_matrix(A)
        warn('splu requires CSC matrix format', SparseEfficiencyWarning)

    A.sort_indices()
    A = A.asfptype()  #upcast to a floating point format

    M, N = A.shape
    if (M != N):
        raise ValueError, "can only factor square matrices" #is this true?

    ftype = superLU_transtabl[A.dtype.char]

    gstrf = eval('_superlu.' + ftype + 'gstrf')
    return gstrf(N, A.nnz, A.data, A.indices, A.indptr, permc_spec,
                 diag_pivot_thresh, drop_tol, relax, panel_size)

def factorized( A ):
    """
    Return a fuction for solving a sparse linear system, with A pre-factorized.

    Example:
      solve = factorized( A ) # Makes LU decomposition.
      x1 = solve( rhs1 ) # Uses the LU factors.
      x2 = solve( rhs2 ) # Uses again the LU factors.
    """
    if isUmfpack and useUmfpack:
        if not isspmatrix_csc(A):
            A = csc_matrix(A)
            warn('splu requires CSC matrix format', SparseEfficiencyWarning)

        A.sort_indices()
        A = A.asfptype()  #upcast to a floating point format

        if A.dtype.char not in 'dD':
            raise ValueError, "convert matrix data to double, please, using"\
                  " .astype(), or set linsolve.useUmfpack = False"

        family = {'d' : 'di', 'D' : 'zi'}
        umf = umfpack.UmfpackContext( family[A.dtype.char] )

        # Make LU decomposition.
        umf.numeric( A )

        def solve( b ):
            return umf.solve( umfpack.UMFPACK_A, A, b, autoTranspose = True )

        return solve
    else:
        return splu( A ).solve
