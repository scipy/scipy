"""Matrix equation solver routines"""
# Author: Jeffrey Armstrong <jeff@approximatrix.com>
# February 24, 2012

# Modified: Chad Fulton <ChadFulton@gmail.com>
# June 19, 2014

# Modified: Ilhan Polat <ilhanpolat@gmail.com>
# September 13, 2016

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.linalg import inv, LinAlgError, norm, cond

from .basic import solve, solve_triangular
from .lapack import get_lapack_funcs
from .decomp_schur import schur
from .decomp_lu import lu
from .decomp_qr import qr
from ._decomp_qz import ordqz
from .decomp import _asarray_validated
from .special_matrices import kron

__all__ = ['solve_sylvester', 'solve_lyapunov', 'solve_discrete_lyapunov',
           'solve_continuous_are', 'solve_discrete_are']


def solve_sylvester(a, b, q):
    """
    Computes a solution (X) to the Sylvester equation :math:`AX + XB = Q`.

    Parameters
    ----------
    a : (M, M) array_like
        Leading matrix of the Sylvester equation
    b : (N, N) array_like
        Trailing matrix of the Sylvester equation
    q : (M, N) array_like
        Right-hand side

    Returns
    -------
    x : (M, N) ndarray
        The solution to the Sylvester equation.

    Raises
    ------
    LinAlgError
        If solution was not found

    Notes
    -----
    Computes a solution to the Sylvester matrix equation via the Bartels-
    Stewart algorithm.  The A and B matrices first undergo Schur
    decompositions.  The resulting matrices are used to construct an
    alternative Sylvester equation (``RY + YS^T = F``) where the R and S
    matrices are in quasi-triangular form (or, when R, S or F are complex,
    triangular form).  The simplified equation is then solved using
    ``*TRSYL`` from LAPACK directly.

    .. versionadded:: 0.11.0

    """

    # Compute the Schur decomp form of a
    r, u = schur(a, output='real')

    # Compute the Schur decomp of b
    s, v = schur(b.conj().transpose(), output='real')

    # Construct f = u'*q*v
    f = np.dot(np.dot(u.conj().transpose(), q), v)

    # Call the Sylvester equation solver
    trsyl, = get_lapack_funcs(('trsyl',), (r, s, f))
    if trsyl is None:
        raise RuntimeError('LAPACK implementation does not contain a proper '
                           'Sylvester equation solver (TRSYL)')
    y, scale, info = trsyl(r, s, f, tranb='C')

    y = scale*y

    if info < 0:
        raise LinAlgError("Illegal value encountered in "
                          "the %d term" % (-info,))

    return np.dot(np.dot(u, y), v.conj().transpose())


def solve_lyapunov(a, q):
    """
    Solves the continuous Lyapunov equation :math:`AX + XA^H = Q`.

    Uses the Bartels-Stewart algorithm to find :math:`X`.

    Parameters
    ----------
    a : array_like
        A square matrix

    q : array_like
        Right-hand side square matrix

    Returns
    -------
    x : array_like
        Solution to the continuous Lyapunov equation

    See Also
    --------
    solve_sylvester : computes the solution to the Sylvester equation

    Notes
    -----
    Because the continuous Lyapunov equation is just a special form of the
    Sylvester equation, this solver relies entirely on solve_sylvester for a
    solution.

    .. versionadded:: 0.11.0

    """

    return solve_sylvester(a, a.conj().transpose(), q)


def _solve_discrete_lyapunov_direct(a, q):
    """
    Solves the discrete Lyapunov equation directly.

    This function is called by the `solve_discrete_lyapunov` function with
    `method=direct`. It is not supposed to be called directly.
    """

    lhs = kron(a, a.conj())
    lhs = np.eye(lhs.shape[0]) - lhs
    x = solve(lhs, q.flatten())

    return np.reshape(x, q.shape)


def _solve_discrete_lyapunov_bilinear(a, q):
    """
    Solves the discrete Lyapunov equation using a bilinear transformation.

    This function is called by the `solve_discrete_lyapunov` function with
    `method=bilinear`. It is not supposed to be called directly.
    """
    eye = np.eye(a.shape[0])
    aH = a.conj().transpose()
    aHI_inv = inv(aH + eye)
    b = np.dot(aH - eye, aHI_inv)
    c = 2*np.dot(np.dot(inv(a + eye), q), aHI_inv)
    return solve_lyapunov(b.conj().transpose(), -c)


def solve_discrete_lyapunov(a, q, method=None):
    """
    Solves the discrete Lyapunov equation :math:`AXA^H - X + Q = 0`.

    Parameters
    ----------
    a, q : (M, M) array_like
        Square matrices corresponding to A and Q in the equation
        above respectively. Must have the same shape.

    method : {'direct', 'bilinear'}, optional
        Type of solver.

        If not given, chosen to be ``direct`` if ``M`` is less than 10 and
        ``bilinear`` otherwise.

    Returns
    -------
    x : ndarray
        Solution to the discrete Lyapunov equation

    See Also
    --------
    solve_lyapunov : computes the solution to the continuous Lyapunov equation

    Notes
    -----
    This section describes the available solvers that can be selected by the
    'method' parameter. The default method is *direct* if ``M`` is less than 10
    and ``bilinear`` otherwise.

    Method *direct* uses a direct analytical solution to the discrete Lyapunov
    equation. The algorithm is given in, for example, [1]_. However it requires
    the linear solution of a system with dimension :math:`M^2` so that
    performance degrades rapidly for even moderately sized matrices.

    Method *bilinear* uses a bilinear transformation to convert the discrete
    Lyapunov equation to a continuous Lyapunov equation :math:`(BX+XB'=-C)`
    where :math:`B=(A-I)(A+I)^{-1}` and
    :math:`C=2(A' + I)^{-1} Q (A + I)^{-1}`. The continuous equation can be
    efficiently solved since it is a special case of a Sylvester equation.
    The transformation algorithm is from Popov (1964) as described in [2]_.

    .. versionadded:: 0.11.0

    References
    ----------
    .. [1] Hamilton, James D. Time Series Analysis, Princeton: Princeton
       University Press, 1994.  265.  Print.
       http://www.scribd.com/doc/20577138/Hamilton-1994-Time-Series-Analysis
    .. [2] Gajic, Z., and M.T.J. Qureshi. 2008.
       Lyapunov Matrix Equation in System Stability and Control.
       Dover Books on Engineering Series. Dover Publications.

    """
    a = np.asarray(a)
    q = np.asarray(q)
    if method is None:
        # Select automatically based on size of matrices
        if a.shape[0] >= 10:
            method = 'bilinear'
        else:
            method = 'direct'

    meth = method.lower()

    if meth == 'direct':
        x = _solve_discrete_lyapunov_direct(a, q)
    elif meth == 'bilinear':
        x = _solve_discrete_lyapunov_bilinear(a, q)
    else:
        raise ValueError('Unknown solver %s' % method)

    return x


def solve_continuous_are(a, b, q, r):
    """
    Solves the continuous algebraic Riccati equation (CARE).

    The CARE is defined as

    .. math::
        A'X + XA - XB R^{-1} B'X + Q = 0

    It is solved directly using a Schur decomposition method.

    Parameters
    ----------
    a : (M, M) array_like
        Input
    b : (M, N) array_like
        Input
    q : (M, M) array_like
        Input
    r : (N, N) array_like
        Non-singular, square matrix

    Returns
    -------
    x : (M, M) ndarray
        Solution to the continuous algebraic Riccati equation

    See Also
    --------
    solve_discrete_are : Solves the discrete algebraic Riccati equation

    Notes
    -----
    Method is taken from [1]_.

    .. versionadded:: 0.11.0

    
    References
    ----------
    
    .. [1] Alan J Laub, "A Schur Method for Solving Algebraic Riccati 
       Equations.", Massachusetts Institute of Technology. Laboratory for 
       Information and Decision Systems. LIDS-R ; 859. Available online :
       http://hdl.handle.net/1721.1/1301
        
    """

    try:
        g = inv(r)
    except LinAlgError:
        raise ValueError('Matrix R in the algebraic Riccati equation solver '
                         'is ill-conditioned')

    g = np.dot(np.dot(b, g), b.conj().transpose())

    z11 = a
    z12 = -1.0*g
    z21 = -1.0*q
    z22 = -1.0*a.conj().transpose()

    z = np.vstack((np.hstack((z11, z12)), np.hstack((z21, z22))))

    # Note: we need to sort the upper left of s to have negative real parts,
    #       while the lower right is positive real components (Laub, p. 7)
    s, u, _ = schur(z, sort='lhp')

    (m, n) = u.shape

    u11 = u[0:m//2, 0:n//2]
    u21 = u[m//2:m, 0:n//2]
    u11i = inv(u11)

    return np.dot(u21, u11i)


def solve_discrete_are(a, b, q, r):
    """
    Solves the discrete algebraic Riccati equation (DARE).

    The DARE is defined as

    .. math::
        
        X = A'XA - (A'XB) (R+B'XB)^{-1} (B'XA) + Q

    The limitations for a solution to exist are :
        
        * All eigenvalues of :math:`A` outside the unit disc, should be 
          controllable. 

        * The associated symplectic pencil (See Notes), should have 
          eigenvalues sufficiently away from the unit circle. 
    
    Parameters
    ----------
    a : (M, M) array_like
        Square matrix
    b : (M, N) array_like
        Input
    q : (M, M) array_like
        Input
    r : (N, N) array_like
        Square matrix

    Returns
    -------
    x : ndarray
        Solution to the discrete algebraic Riccati equation.
        
    Raises
    ------
    LinAlgError
        For cases where the stable subspace of the pencil could not be 
        isolated. See Notes section and the references for details. 

    See Also
    --------
    solve_continuous_are : Solves the continuous algebraic Riccati equation

    Notes
    -----
    The equation is solved by forming the extended symplectic matrix pencil, 
    as described in [1]_, :math:`H - \lambda J` given by the block 
    matrices::
        
    
           [  A   0   B ]             [ I   0   B ] 
           [ -Q   I   0 ] - \lambda * [ 0  A^T  0 ]
           [  0   0   R ]             [ 0 -B^T  0 ]
        
    and using a QZ decomposition method.
    
    In this algorithm, the fail conditions are linked to the symmetrycity 
    of the product :math:`U_2 U_1^{-1}` and condition number of 
    :math:`U_1`. Here, :math:`U` is the 2m-by-m matrix that holds the 
    eigenvectors spanning the stable subspace with 2m rows and partitioned 
    into two m-row matrices. See [1]_ and [2]_ for more details.
    

    .. versionadded:: 0.11.0

    References
    ----------
    .. [1]  P. van Dooren , "A Generalized Eigenvalue Approach For Solving 
       Riccati Equations.", SIAM Journal on Scientific and Statistical 
       Computing, Vol.2(2), DOI: 10.1137/0902010
       
    .. [2] Alan J Laub, "A Schur Method for Solving Algebraic Riccati 
       Equations.", Massachusetts Institute of Technology. Laboratory for 
       Information and Decision Systems. LIDS-R ; 859. Available online :
       http://hdl.handle.net/1721.1/1301
       
    """
    
    a = np.atleast_2d(_asarray_validated(a))
    b = np.atleast_2d(_asarray_validated(b))
    q = np.atleast_2d(_asarray_validated(q))
    r = np.atleast_2d(_asarray_validated(r))
    
    # Get the correct data types otherwise Numpy complains about pushing 
    # complex numbers into real arrays.
    r_or_c = complex if np.iscomplexobj(b) else float

    for ind, mat in enumerate((a,q,r)):
        if np.iscomplexobj(mat):
            r_or_c = complex
            
        if not np.equal(*mat.shape):
            raise ValueError("Matrix {} should be square.".format("aqr"[ind]))
    
    # Shape consistency checks
    m, n = b.shape
    if m != a.shape[0]:
        raise ValueError("Matrix a and b should have the same number of rows.")
    if m != q.shape[0]:
        raise ValueError("Matrix a and q should have the same shape.")
    if n != r.shape[0]:
        raise ValueError("Matrix b and r should have the same number of cols.")

    # Check if the data matrices q, r are (sufficiently) hermitian
    for ind, mat in enumerate((q,r)):
        if norm(mat - mat.conj().T,1) > np.spacing(norm(mat,1))*100:
            raise ValueError("Matrix {} should be symmetric/hermitian."
                             "".format("qr"[ind]))

    # Form the matrix pencil        
    H = np.zeros((2*m + n, 2*m + n), dtype = r_or_c)
    H[:m, :m] = a
    H[:m, 2*m:] = b
    H[m:2*m, :m] = -q
    H[m:2*m, m:2*m] = np.eye(m)
    H[2*m:, 2*m:] = r

    J = np.zeros_like(H, dtype = r_or_c)
    J[:m, :m] = np.eye(m)
    J[m:2*m, m:2*m] = a.conj().T
    J[2*m:, m:2*m] = -b.conj().T

    # Deflate the pencil by the R column ala Ref. [1]
    q_of_qr, _ = qr(H[:,-n:])
    H = q_of_qr[:,n:].conj().T.dot(H[:,:2*m])
    J = q_of_qr[:,n:].conj().T.dot(J[:,:2*m])
    
    # Decide on which output type is needed for QZ 
    out_str = 'real' if r_or_c == float else 'complex'
    
    _, _, _, _, _, u = ordqz(H, J, sort = 'iuc', overwrite_a=True, 
                    overwrite_b=True, check_finite=False, output=out_str)
    
    # Get the relevant parts of the stable subspace basis
    u00 = u[:m, :m]
    u10 = u[m:, :m]
    
    # Solve via back-substituion after checking the condition of u00 
    up, ul, uu = lu(u00)

    if 1/cond(uu) < np.spacing(1.):
        raise LinAlgError('Failed to find a finite solution.')
    
    # Exploit the triangular structure
    x = solve_triangular(ul.conj().T,
                         solve_triangular(uu.conj().T, 
                                          u10.conj().T,
                                          lower=True),
                         unit_diagonal=True,
                         ).conj().T.dot(up.conj().T)

    # Check the deviation from symmetry for success
    u_sym = u00.conj().T.dot(u10)
    u_sym = u_sym - u_sym.conj().T
    sym_threshold = np.max(np.spacing([1000., norm(u_sym, 1)]))

    if np.max(np.abs(u_sym)) > sym_threshold:
        raise LinAlgError('The associated symplectic pencil has eigenvalues'
                          'too close to the unit circle')
    
    return (x + x.conj().T)/2
