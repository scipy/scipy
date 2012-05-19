'''
Translated from Octave code at: http://www.ecs.shimane-u.ac.jp/~kyoshida/lpeng.htm
and placed under MIT licence by Enzo Michelangeli with permission explicitly
granted by the original author, Prof. Kazunobu Yoshida  

-----------------------------------------------------------------------------
Copyright (c) 2010, Kazunobu Yoshida, Shimane University, and Enzo Michelangeli, 
IT Vision Limited

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
-----------------------------------------------------------------------------

Usage:
 
 optx,zmin,is_bounded,sol,basis = lp(c,A,b)
 
  This program finds a solution of the standard linear programming problem:
    minimize    z = c'x
    subject to  Ax = b, x >= 0
  using the two phase method, where the simplex method is used at each stage.
  Returns the tuple:
    optx: an optimal solution.
    zmin: the optimal value. 
    is_bounded: True if the solution is bounded; False if unbounded.
    sol: True if the problem is solvable; False if unsolvable.
    basis: indices of the basis of the solution.
    
  All the non-scalar data types are numpy arrays.   
'''
from numpy import *

class NamedTuple(tuple):
    def __new__(cls, values, names):
        self = tuple.__new__(cls, values)
        for value, name in zip(values, names):
            setattr(self, name, value)
        return self

class Solution(NamedTuple):
    _fields = []
    def __new__(cls, *a, **kw):
        values = list(a)
        for name in cls._fields[len(values):]:
            values.append(kw.pop(name))
        if len(values) != len(cls._fields) or kw:
            raise ValueError("Invalid arguments")
        return NamedTuple.__new__(cls, values, cls._fields)

    def __repr__(self):
        return "%s%s" % (self.__class__.__name__, 
                         NamedTuple.__repr__(self))

class LPSolution(Solution):
    """
    Solution to a linear programming problem

    Attributes
    ----------
    x
        The optimal solution
    min
        The optimal value
    is_bounded
        True if the solution is bounded; False if unbounded
    solvable
        True if the problem is solvable; False if unsolvable
    basis
        Indices of the basis of the solution.
    """
    _fields = ['x', 'min', 'is_bounded', 'solvable', 'basis']

def lp(c, A, b, tol=1e-10):
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
    indx = arange(n)
    basis = arange(n, n+m) # m elements from n to n+m-1
    is_bounded = _simplex(H, basis, indx, 1)
    if H[m+1,n] < -tol:   # last row, last column
        sol = False
        #print('unsolvable')
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
                indx = delete(indx, j)  # indx[j] = [] #delete element j from indx
                j = j-1
        H = delete(H, m+1, 0)
        if size(indx) > 0:
        # Phase two
            is_bounded = _simplex(H,basis,indx,2);
            if is_bounded:
                optx = zeros(n+m);
                for i in xrange(m):
                    optx[basis[i]] = H[i,-1]
                optx = optx[0:n] 
                zmin = -H[-1,-1]    #  last row, last column
            else:
                optx = None
                zmin = -Inf
        else:
            optx = zeros(n)
            zmin = 0;
    return LPSolution(optx, zmin, is_bounded, sol, basis)  

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


######### Unit test section #########

from numpy.testing import *

def test_lp(prt=False):
    m1 = 20
    m2 = 50
    probs = [
        {
            'A': array([
                [2.,  5., 3., -1.,  0.,  0.],
                [3., 2.5, 8.,  0., -1.,  0.],
                [8.,10.,  4.,  0.,  0., -1.]]),
            'b': array([185., 155., 600.]),
            'c': array([4., 8., 3., 0., 0., 0.]),
            'result': [
                    array([ 66.25, 0., 17.5, 0., 183.75, 0.]),
                    317.5,
                    True,
                    True,
                    array([2, 0, 4])            
                ]
        },
        {        
            'A': array([
                [-1., -1., -1.,  0.,  0.,  0.],
                [ 0.,  0.,  0.,  1.,  1.,  1.],
                [ 1.,  0.,  0.,  1.,  0.,  0.],
                [ 0.,  1.,  0.,  0.,  1.,  0.],
                [ 0.,  0.,  1.,  0.,  0.,  1.]]),
            'b': array([-0.5, 0.4, 0.3, 0.3, 0.3]),
            'c': array([2.8, 6.3, 10.8, -2.8, -6.3, -10.8]),
            'result': [
                    array([0.3, 0.2, 0.0, 0.0, 0.1, 0.3]),
                    -1.77,
                    True,
                    True,
                    array([1, 7, 0, 4, 5])            
                ]
        },
        {   # with degeneracy
            'A': array([[cos(2*pi*i/(m1+1))-1., sin(2*pi*i/(m1+1))] for i in xrange(1,m1+1)]).T,
            'b': zeros(2).T,
            'c': -ones(m1).T,
            'result': [
                    zeros(m1),
                    0.,
                    True,
                    True,
                    array([0,19])
                ]
            
        },
        {   # with unboundedness (0 is a member of the convex hull of these vectors)
            'A': array([[cos(2*pi*i/(m2+1))-1., sin(2*pi*i/(m2+1))] for i in xrange(0,m2)]).T,
            'b': zeros(2).T,
            'c': -ones(m2).T,
            'result': [
                    None,   # unchecked when unbounded
                    -Inf,   # unchecked when unbounded
                    False,
                    True,
                    array([2, 49])
                ]
            
        }, 
        {   # Unsolvable
            'A': array([[cos(2*pi*i/(m2+1))-1., sin(2*pi*i/(m2+1))] for i in xrange(0,m2)]).T,
            'b': ones(2).T,
            'c': -ones(m2).T,
            'result': [
                    None,   # unchecked when unsolvable
                    None,   # unchecked when unsolvable
                    None,   # unchecked when unsolvable
                    False,
                    array([50, 1])
                ]
            
        }, # add other test cases here...
    ]


    for prob in probs:
        #optx, zmin, bounded, solvable, basis = lp(prob['c'],prob['A'],prob['b'])
        lpsol = lp(prob['c'],prob['A'],prob['b'])
        optx = lpsol.x
        zmin = lpsol.min
        bounded = lpsol.is_bounded
        solvable = lpsol.solvable
        basis = lpsol.basis
        if prt:
            print "A:\n",prob['A']
            print "b:",prob['b']
            print "c:",prob['c']
            print " ---->"
            print "optx:",optx
            print "zmin:",zmin
            print "bounded:",bounded
            print "solvable:",solvable
            print "basis:",basis
            print "-------------------------------------------"
        else:
            expected_res = prob['result']
            assert_equal(solvable, expected_res[3])
            assert_equal(basis, expected_res[4])
            if solvable:
                assert_equal(bounded, expected_res[2])
                if bounded:
                    assert_almost_equal(optx, expected_res[0])
                    assert_almost_equal(zmin, expected_res[1]) # when unbounded zmin == -Inf, but -Inf != -Inf so we won't check it...

if __name__ == "__main__":
    #test_lp(True)
    run_module_suite()
