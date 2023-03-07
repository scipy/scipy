import scipy
import numpy as np
import logging

def qr_solve(A,y,silent=True):
    """
    Solves over-determined and underdetemined linar systems :math:`y=Ax`.

    Parameters
    ----------
    A : (M, N) array_like
        Any matrix.
    
    y : ndarray
        A column vector to solve for with `M` rows.
        
    silent : {'True', 'False'}, optional
        To log if the solution is a true solution.

    Returns
    -------
    x : ndarray
        Solution to the linear system equation

    Notes
    -----
    This uses QR-Decomposition using the permutation matrix to find a soultion to the
    given linear system `y=Ax` through back substitution.

    References
    ----------
    .. [1] https://inst.eecs.berkeley.edu/~ee127/sp21/livebook/l_lineqs_solving.html

    Examples
    --------
    Given `A` is `m\times n`, `x \in \mathbb{R}^{n}` `y\in\mathbb{R}^m` in the equation `y=Ax`, solve for `x`:

    >>> import numpy as np
    >>> from scipy import linalg
    >>> A = np.random.rand( m, n )
    >>> x = np.random.rand( n )
    >>> y = A@x
    >>> x_0 = _qr_solve.qr_solve(A, y)
    >>> np.allclose(y, A@x_0)
    True

    """
    
    A = np.asarray(A)
    y = np.asarray(y)
    
    #Solving y=Ax for an x_0
    #A has m rows and n columns -> y has m rows and x has n rows
    m,n = A.shape
    
    #QR decomposition
    Q,R,P= scipy.linalg.qr(A,pivoting=True)
    #Q is an m by m orthogonal matrix
    #R is an m by n upper right triangle matrix
    
    #P is a permuntation matrix with n rows and n columns
    #P can order A by rank for *rank revealing*
    P = np.eye(len(P))[:,P]
    
    #Let r be the number of linearly independent columns & rows in A (the rank)
    rank = r = sum(np.around(np.diag(R),12)!=0)
    
    #Q is a m by m square orthogonal matrix
    #Q_1 has m rows and r columns
    Q_1 = Q[:,:rank]
    
    #R_1 is an r by r upper right triangular matrix
    R_1 = R[:rank,:rank]
    
    #R_2 is an r by m-r matrix
    R_2 = R[:,:-(rank-m)]
    
    #z_1 is a column vector with r elements
    z_1 = scipy.linalg.solve(R_1,Q_1.T@y)
    
    #We pad z_1 with r-m zeros at the bottom for a solution vector
    padding = np.zeros(n-r) #zero padding
    x_0 = P@np.hstack([z_1,padding]) #Add padding
    
    if not silent:
        #Log if there was a solution
        is_solution = np.allclose(y,A@x_0)
        logging.info("Solution Found!") if is_solution else loggin.info("No Solution!")
        
    return x_0