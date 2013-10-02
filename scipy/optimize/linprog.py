"""
A basic linear programming function using a simplex method.

Functions
---------
.. autosummary::
   :toctree: generated/

    linprog

"""

from __future__ import division, print_function, absolute_import

import numpy as np


__all__ = ['linprog']


def linprog(c,A_eq=None,b_eq=None,A_lb=None,b_lb=None,A_ub=None,b_ub=None,lb=None,ub=None):
    """
    Solve a Linear Programming Problem in Standard Form.
    
    minimize:  dot(c,x)
    
    subject to:
               
               A_eq * x = b_eq
               A_lb * x >= b_lb
               A_ub * x <= b_ub
    """
    f = np.asarray(c)
    
    Aeq = np.asarray(A_eq) if not A_eq is None else None
    Alb = np.asarray(A_lb) if not A_lb is None else None
    Aub = np.asarray(A_ub) if not A_ub is None else None
    beq = np.asarray(b_eq) if not b_eq is None else None
    blb = np.asarray(b_lb) if not b_lb is None else None
    bub = np.asarray(b_ub) if not b_ub is None else None
    
    n = len(c)
    
    if not Aeq is None:
        m_eq, n_eq = Aeq.shape
    else:
        m_eq, n_eq = 0, n
        
    if not Alb is None:
        m_lb, n_lb = Alb.shape
    else:
        m_lb, n_lb = 0, n
        
    if not Aub is None:
        m_ub, n_ub = Aub.shape
    else:
        m_ub, n_ub = 0, n
    
    if not(n == n_eq == n_lb == n_ub ):
        # raise an exception
        print( """inconsistency in problem size.  
                 Vector c and matrices A should have the same number of columns""" )
        return None
        
    x = np.zeros([n])
    
    m = m_eq + m_lb + m_ub
       
    tableau = np.zeros([m+1,n+1])
    
    print(tableau)
    
    
def simplex_standard_maximization(c,A,b):
    """
    
    """
    AA = np.asarray(A) 
    bb = np.asarray(b) 
    cc = np.asarray(c) 
    
    m = len(b)
    n = len(c)
    
    A_rows, A_cols = AA.shape
    
    if not A_rows == m:
        # raise here
        print("Number of rows in A must be equal to the size of b")
    if not A_cols == n:
        # raise here
        print("Number of columsn in A must be equal to the size of c")
        
    # Add slack variables for our upper-bound constraints
    n_slack = m
    
    # Create the tableau
    T = np.zeros([m+1,n+n_slack+1])
    
    # Insert objective into tableau
    T[-1,:n] = -cc
    
    # Add A to the tableau
    T[:m,:n] = AA
    
    # Add the slack variables to the tableau
    T[:m,n:n+n_slack] = np.eye(n_slack)    
    
    # At b to the tableau
    T[:-1,-1] = bb
    
    # Tableau Finished
    print(T)
    
    count = 0
    
    error_flag = 0
    
    while np.any( T[-1,:-1] < 0 ):
        print("Iteration: ",count)
        
        count += 1
    
        # Find the most negative value in bottom row of tableau
        pivcol = np.argmin(T[-1,:-1])
    
        print(pivcol)
    
        # Find the minimizer of b[x_to_enter]/a[i,x_to_enter]
        test = bb/AA[:,pivcol]
        print(test)
        pivrow = np.argmin(test)
        
        print(pivrow)
        
        # Do pivot
        
        pivval = T[pivrow,pivcol]
        T[pivrow,:] = T[pivrow,:] / pivval
        for irow in range(T.shape[0]):
            if irow != pivrow:
                T[irow,:] = T[irow,:] - T[pivrow,:]*T[irow,pivcol]               
        
        print(T)
        
    if error_flag == 0:
        # if disp = True
        print("Optimal solution found:")
        x = np.linalg.solve( T[:-1,:n], T[:-1,-1] )
        print("x* = ", x)
        return x
        
    
    
if __name__ == "__main__":
    
    c = [1,1]
    A_ub = [[ 1, 2],
            [ 3, 1]]
    b_ub = [1,1]
    lb = [ 0, 0]

    #linprog(c=c,A_ub=A_ub,b_ub=b_ub,lb=lb)
    
    A = [[1,1],
         [2,3],
         [3,1]]
         
    b = [2,12,12]
    
    c = [3,1]
    
    # AAE550 Class 15 p 9
    
    c = [1,1]
    A = [[1,2],
         [3,1]]
    b = [1,1]
    
    simplex_standard_maximization(c,A,b)
    
    
    
    
    
    
    