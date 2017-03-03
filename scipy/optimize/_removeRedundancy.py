# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 22:08:23 2017

@author: Matt
"""


import numpy as np
import scipy as sp

def _rowCount(A):
    """
    Counts the number of nonzeros in each row of input array A.
    Nonzeros are defined as any element with absolute value greater than 
    tol = 1e-13. This value should probably be an input to the function.

    Parameters
    ----------
    A : 2-D array
        An array representing a matrix
        
    Returns
    -------
    rowcount : 1-D array
        Number of nonzeros in each row of A
        
    """
    tol = 1e-13
    return np.sum(np.abs(A) > tol, axis = 1)
    
def _getDensest(A, eligibleRows):
    """
    Returns the index of the densest row of A. Ignores rows that are not
    eligible for consideration.

    Parameters
    ----------
    A : 2-D array
        An array representing a matrix
    eligibleRows : 1-D logical array
        Values indicate whether the corresponding row of A is eligible 
        to be considered
        
    Returns
    -------
    i_densest : int
        Index of the densest row in A eligible for consideration
        
    """
    rowCounts = _rowCount(A)
    return np.argmax(rowCounts*eligibleRows)
    
def _removeZeroRows(A,b):
    """
    Eliminates trivial equations from system of equations defined by Ax = b
   and identifies trivial infeasibilities 

    Parameters
    ----------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations
        
    Returns
    -------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations
    status: int
        An integer indicating the status of the removal operation
        0: No infeasibility identified
        2: Trivially infeasible
        
    """
    status = 0
    i_zero = _rowCount(A) == 0
    A = A[np.logical_not(i_zero),:]
    if not(np.allclose(b[i_zero],0)):
        status = 2    
    b = b[np.logical_not(i_zero)]
    return A,b,status

def _removeRedundancy(A,b):
    """
    Eliminates redundant equations from system of equations defined by Ax = b
    and identifies infeasibilities.

    Parameters
    ----------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations
        
    Returns
    -------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations
    status: int
        An integer indicating the status of the system
        0: No infeasibility identified
        2: Trivially infeasible
        
    """

    A,b,status = _removeZeroRows(A,b)

    if status != 0:
        return A,b,status
        
    U,s,Vh = sp.linalg.svd(A)
    eps = np.finfo(float).eps
    tol = s.max() * max(A.shape) * eps
    
    m,n = A.shape
    s_min = s[-1] if m  <= n else 0 
    
    while abs(s_min) < tol:
        v = U[:,-1]
        eligibleRows = np.abs(v) > tol*10e6 # rows need to be represented in significant amount
        if not np.any(eligibleRows) or np.any(np.abs(v.dot(A)) > tol):
            status = 4
            break
#            print "Due to numerical issues, redundant equality constraints cannot be removed automatically."
        if np.any(np.abs(v.dot(b)) > tol):
            status = 2
            break
#            print "There is a linear combination of rows of A_eq that results in zero, suggesting a redundant constraint. However the same linear combination of b_eq is nonzero, suggesting that the constraints conflict and the problem is infeasible."
        i_remove = _getDensest(A,eligibleRows)
        A = np.delete(A,i_remove,axis = 0)
        b = np.delete(b,i_remove)
        U,s,Vh = sp.linalg.svd(A)
        m,n = A.shape
        s_min = s[-1] if m  <= n else 0 
        
    return A,b,status



#m = 5
##A_eq = np.eye(m)
#m,n = A_eq.shape
#A = np.hstack((np.eye(m),A_eq))
#indices = range(m)

#v = set(indices)
#b = set(indices)
#k = set(range(m+n))
#
#n = 0
#while n < len(indices):
#    i = indices[n]
#    print i
#    B = A[:,list(b)]
#    print B
#    e = np.zeros(m)
#    e[i] = 1
##    pi = np.linalg.solve(B.T,e)
#    ei = e[np.newaxis,:]
#    pi = e.dot(np.linalg.inv(B))
#    js = k-b.union(v)
#    dependent = True
#    for j in js:
#        if abs(pi.dot(A[:,j])) > 0:
#            print j
#            b.remove(i)
#            b.add(j)
#            indices.append(j)
#            dependent = False
#            break
#    if dependent:
#        print i, "is dependent"
#    n = n+1
#    

#v = set(indices)
#b = range(m)
#k = set(range(m+n))
#
#n = 0
#while n in v:
#    i = indices[n]
##    print i
#    B = A[:,b]
##    print B
#    e = np.zeros(m)
#    e[i] = 1
##    pi = np.linalg.solve(B.T,e)
#    ei = e[np.newaxis,:]
#    pi = e.dot(np.linalg.inv(B))
#    js = k-set(b).union(v)
#    dependent = True
#    for j in js:
#        if abs(pi.dot(A[:,j])) > 0:
##            print j
#            b[n] = j
#            indices.append(j)
#            dependent = False
#            break
#    if dependent:
#        print i, "is dependent"
#    n = n+1
    