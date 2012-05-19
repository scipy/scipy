# Translated from Octave code at:
# http://www.ecs.shimane-u.ac.jp/~kyoshida/lpeng.htm
# and placed under MIT licence by Enzo Michelangeli with permission explicitly
# granted by the original author, Prof. Kazunobu Yoshida
#
# -----------------------------------------------------------------------------
# Copyright (c) 2010, Kazunobu Yoshida, Shimane University, and Enzo Michelangeli,
# IT Vision Limited
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# -----------------------------------------------------------------------------

from numpy import *
from optimize import Result

class LPResult(Result):
    """
    Solution to a linear programming problem

    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    fun : float
        The optimal value
    is_bounded : bool
        True if the solution is bounded; False if unbounded
    is_solvable : bool
        True if the problem is solvable; False if unsolvable
    basis : ndarray
        Indices of the basis of the solution.
    """
    _fields = ['x', 'min', 'is_bounded', 'solvable', 'basis']

def lp_solve(c, A, b, tol=1e-10):
    """
    Solves a linear programming problem using a two-phase method.

    The problem solved is::

        minimize    z = c' x
        subject to  A x = b, x >= 0

    Parameters
    ----------
    c : array_like
    A : array_like
    b : array_like
        Problem parameters
    tol : float, optional
        Tolerance

    Returns
    -------
    sol : LPResult
        A solution object, with the properties:
            x : ndarray
                An optimal solution.
            fun : float
                The optimal value.
            is_bounded : bool
                True if the solution is bounded; False if unbounded.
            is_solvable : bool
                True if the problem is solvable; False if unsolvable.
            basis : ndarray
                Indices of the basis of the solution.

    Notes
    -----
    The solution is found using the two phase method, where the
    simplex method is used at each stage.

    """
    c = asarray(c)
    A = asarray(A)
    b = asarray(b)

    m,n = A.shape   # m: number of constraints; n: number of variables
    for i in xrange(m):
        if b[i] < 0.0:
            A[i,:] = -A[i,:]
            b[i] = -b[i]
    d = -sum(A, axis=0)
    w0 = sum(b)
    H = vstack([     #  The initial simplex table of phase one
         hstack([A, array([b]).T]), # first m rows
         hstack([c, 0.]),   # last-but-one
         hstack([d, -w0])]) # last
    indx = range(n)
    basis = arange(n, n+m) # m elements from n to n+m-1
    is_bounded = _simplex(H, basis, indx, 1)
    if H[m+1,n] < -tol:   # last row, last column
        sol = False
        optx = None
        zmin = None
        is_bounded = None
    else:
        sol = True
        j = -1
        for i in xrange(n):
            j = j+1
            if H[m+1,j] > tol:
                H = delete(H, j, 1)     # H[:,j] = [] # delete column j from H
                del indx[j]
                j = j-1
        H = delete(H, m+1, 0)
        if indx > 0:
            # Phase two
            is_bounded = _simplex(H,basis,indx,2)
            if is_bounded:
                optx = zeros(n+m)
                for i in xrange(m):
                    optx[basis[i]] = H[i,-1]
                optx = optx[0:n]
                zmin = -H[-1,-1]    #  last row, last column
            else:
                optx = None
                zmin = -Inf
        else:
            optx = zeros(n)
            zmin = 0

    return LPResult(x=optx,
                    fun=zmin,
                    success=(sol and is_bounded),
                    is_bounded=is_bounded,
                    is_solvable=sol,
                    basis=basis)

def _simplex(H,basis,indx,s):
    '''
      [H1,basis,is_bounded] = _simplex(H,basis,indx,s)
      H: simplex table (MODIFIED).
      basis: the indices of basis (MODIFIED).
      indx: the indices of x.
      s: 1 for phase one; 2 for phase two.
      H1: new simplex table.
      is_bounded: True if the solution is bounded; False if unbounded.
    '''
    if s == 1:
        s0 = 2
    elif s == 2:
        s0 = 1
    n1 = H.shape[0]
    sol = False
    while not sol:
        q = H[-1, :-1] # last row, all columns but last
        jp = argmin(q)
        fm = q[jp]
        if fm >= 0:
            is_bounded = True    # bounded solution
            sol = True
        else:
            q = H[:-s0,jp]
            ip = argmax(q)
            hm = q[ip]
            if hm <= 0:
                is_bounded = False # unbounded solution
                sol = True
            else:
                h1 = zeros(n1-s0)
                for i in xrange(n1-s0):
                    if H[i,jp] > 0:
                        h1[i] = H[i,-1]/H[i,jp]
                    else:
                        h1[i] = Inf
                ip = argmin(h1)
                minh1 = h1[ip]
                basis[ip] = indx[jp]
                if not _pivot(H,ip,jp):
                    raise ValueError("the first parameter is a Singular matrix")
    return is_bounded

def _pivot(H,ip,jp):
    # H is MODIFIED
    n, m = H.shape
    piv = H[ip,jp]
    if piv == 0:
        #print('singular')
        return False
    else:
        H[ip,:] /= piv
        for i in xrange(n):
            if i != ip:
                H[i,:] -= H[i,jp]*H[ip,:]
    return True
