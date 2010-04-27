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

useUmfpack = True


__all__ = [ 'use_solver', 'spsolve', 'splu', 'spilu', 'factorized' ]

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


def spsolve(A, b, permc_spec=None, use_umfpack=True):
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

    use_umfpack = use_umfpack and useUmfpack

    if isUmfpack and use_umfpack:
        if noScikit:
            warn( 'scipy.sparse.linalg.dsolve.umfpack will be removed,'\
                    ' install scikits.umfpack instead', DeprecationWarning )
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
        elif isspmatrix_csr(A):
            flag = 0 # CSR format
        else:
            A = csc_matrix(A)
            flag = 1

        b = asarray(b, dtype=A.dtype)
        options = dict(ColPerm=permc_spec)
        return _superlu.gssv(N, A.nnz, A.data, A.indices, A.indptr, b, flag,
                             options=options)[0]

def splu(A, permc_spec=None, diag_pivot_thresh=None,
         drop_tol=None, relax=None, panel_size=None, options=dict()):
    """
    Compute the LU decomposition of a sparse, square matrix.

    Parameters
    ----------
    A
        Sparse matrix to factorize. Should be in CSR or CSC format.

    permc_spec : str, optional
        How to permute the columns of the matrix for sparsity preservation.
        (default: 'COLAMD')

        - ``NATURAL``: natural ordering.
        - ``MMD_ATA``: minimum degree ordering on the structure of A^T A.
        - ``MMD_AT_PLUS_A``: minimum degree ordering on the structure of A^T+A.
        - ``COLAMD``: approximate minimum degree column ordering

    diag_pivot_thresh : float, optional
        Threshold used for a diagonal entry to be an acceptable pivot.
        See SuperLU user's guide for details [SLU]_
    drop_tol : float, optional
        (deprecated) No effect.
    relax : int, optional
        Expert option for customizing the degree of relaxing supernodes.
        See SuperLU user's guide for details [SLU]_
    panel_size : int, optional
        Expert option for customizing the panel size.
        See SuperLU user's guide for details [SLU]_
    options : dict, optional
        Dictionary containing additional expert options to SuperLU.
        See SuperLU user guide [SLU]_ (section 2.4 on the 'Options' argument)
        for more details. For example, you can specify
        ``options=dict(Equil=False, IterRefine='SINGLE'))``
        to turn equilibration off and perform a single iterative refinement.

    Returns
    -------
    invA : scipy.sparse.linalg.dsolve._superlu.SciPyLUType
        Object, which has a ``solve`` method.

    See also
    --------
    spilu : incomplete LU decomposition

    Notes
    -----
    This function uses the SuperLU library.

    References
    ----------
    .. [SLU] SuperLU http://crd.lbl.gov/~xiaoye/SuperLU/

    """

    if not isspmatrix_csc(A):
        A = csc_matrix(A)
        warn('splu requires CSC matrix format', SparseEfficiencyWarning)

    A.sort_indices()
    A = A.asfptype()  #upcast to a floating point format

    M, N = A.shape
    if (M != N):
        raise ValueError, "can only factor square matrices" #is this true?

    _options = dict(DiagPivotThresh=diag_pivot_thresh, ColPerm=permc_spec,
                    PanelSize=panel_size, Relax=relax)
    if options is not None:
        _options.update(options)
    return _superlu.gstrf(N, A.nnz, A.data, A.indices, A.indptr,
                          ilu=False, options=_options)

def spilu(A, drop_tol=None, fill_factor=None, drop_rule=None, permc_spec=None,
          diag_pivot_thresh=None, relax=None, panel_size=None, options=None):
    """
    Compute an incomplete LU decomposition for a sparse, square matrix A.

    The resulting object is an approximation to the inverse of A.

    Parameters
    ----------
    A
        Sparse matrix to factorize

    drop_tol : float, optional
        Drop tolerance (0 <= tol <= 1) for an incomplete LU decomposition.
        (default: 1e-4)
    fill_factor : float, optional
        Specifies the fill ratio upper bound (>= 1.0) for ILU. (default: 10)
    drop_rule : str, optional
        Comma-separated string of drop rules to use.
        Available rules: ``basic``, ``prows``, ``column``, ``area``,
        ``secondary``, ``dynamic``, ``interp``. (Default: ``basic,area``)

        See SuperLU documentation for details.
    milu : str, optional
        Which version of modified ILU to use. (Choices: ``silu``,
        ``smilu_1``, ``smilu_2`` (default), ``smilu_3``.)

    Remaining other options
        Same as for `splu`

    Returns
    -------
    invA_approx : scipy.sparse.linalg.dsolve._superlu.SciPyLUType
        Object, which has a ``solve`` method.

    See also
    --------
    splu : complete LU decomposition

    Notes
    -----
    To improve the better approximation to the inverse, you may need to
    increase ``fill_factor`` AND decrease ``drop_tol``.

    This function uses the SuperLU library.

    References
    ----------
    .. [SLU] SuperLU http://crd.lbl.gov/~xiaoye/SuperLU/

    """

    if not isspmatrix_csc(A):
        A = csc_matrix(A)
        warn('splu requires CSC matrix format', SparseEfficiencyWarning)

    A.sort_indices()
    A = A.asfptype()  #upcast to a floating point format

    M, N = A.shape
    if (M != N):
        raise ValueError, "can only factor square matrices" #is this true?

    _options = dict(ILU_DropRule=drop_rule, ILU_DropTol=drop_tol,
                    ILU_FillFactor=fill_factor,
                    DiagPivotThresh=diag_pivot_thresh, ColPerm=permc_spec,
                    PanelSize=panel_size, Relax=relax)
    if options is not None:
        _options.update(options)
    return _superlu.gstrf(N, A.nnz, A.data, A.indices, A.indptr,
                          ilu=True, options=_options)

def factorized( A ):
    """
    Return a fuction for solving a sparse linear system, with A pre-factorized.

    Example:
      solve = factorized( A ) # Makes LU decomposition.
      x1 = solve( rhs1 ) # Uses the LU factors.
      x2 = solve( rhs2 ) # Uses again the LU factors.
    """
    if isUmfpack and useUmfpack:
        if noScikit:
            warn( 'scipy.sparse.linalg.dsolve.umfpack will be removed,'\
                    ' install scikits.umfpack instead', DeprecationWarning )

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
