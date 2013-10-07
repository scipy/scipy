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



def column_is_basic(j,tableau):
    """
    Determines if the jth column of the given simplex method tableau is canonical.
    
    A column is considered canonical if all rows within it contain value 0, and the
    remaining column contains a 1.
    
    If the column is canonical, return the row index whose value is 1.
    If the column is not canonical, return -1.
    """    
    if np.count_nonzero( tableau[:-1,j] ) == 1:
        if tableau[  np.flatnonzero( tableau[:-1,j] )[0], j ] == 1:
            return np.flatnonzero( tableau[:-1,j] )[0]
    return -1



def lpsimplex(tableau,n,n_surplus,n_slack,n_artificial,phase=2):
    """
    Solve a linear programming problem in "standard form" using the Simplex Method.
    
    Minimize :math:`f = c^T x`
    
    subject to
    
    .. math::
    
        Ax = b
        
        x_i >= 0
        
        b_j >= 0

    Status:
    
    0 : Optimal solution found - one feasible solution
    
    1 : Optimal solution found - iteration successful
    
    2 : 

    """
    print("Initial Tableau")
    print(tableau)
    T0 = np.empty_like(tableau)
    T0[:,:] = tableau
    
    # TODO FIGURE OUT WHEN TO JUST SOLVE AX=B
    t_rows, t_cols = tableau.shape
    
    if False: #t_cols == t_rows + 1:    
        # This is just Ax=b
        x = np.linalg.solve( tableau[:-1,:n], tableau[:-1,-1] )
        status = 0
    else:
        # Otherwise we need to solve it iteratively
    
        count = 0    
        status = 1
        cycle = 0
        
        while np.any( tableau[-1,:-1] < 0 ):
            if count > 2**n:
                if np.all(tableau == T0):
                    print("Cycling detected")
                    cycle += 1
            print("Iteration: ",count)
            
            count += 1
        
            # Find the most negative value in bottom row of tableau
            pivcol = np.argmin(tableau[-1,:-1])
        
            # Find the minimizer of b[x_to_enter]/a[i,x_to_enter]
            if phase == 1:
                row_quotient = list(zip( range( t_rows-2), tableau[:-2,-1] / tableau[:-2,pivcol] ))
            else:
                row_quotient = list(zip( range( t_rows-1), tableau[:-1,-1] / tableau[:-1,pivcol] ))
    
            print("Pivot col = ", pivcol)
            print("Row Quotients = ", row_quotient)
    
            # Sort row quotient, we want to use the minimizer of the quotient, as long
            # as that pivot element is positive (if negative, go to the next most negative quotient)
            row_quotient.sort(key= lambda tup: tup[1])
            
            for i in range(cycle, t_rows-1):
                if tableau[ row_quotient[i][0], pivcol ] > 0:
                    pivrow = row_quotient[i][0]
                    if cycle > 0:
                        cycle = 0
                    break
            else:
                print("No suitable pivot.  No solution.")
                status = 2
                break
                    
            # Do pivot
            print("Pivot is at",pivrow,pivcol)
            
            pivval = tableau[pivrow,pivcol]
            tableau[pivrow,:] = tableau[pivrow,:] / pivval
            for irow in range(tableau.shape[0]):
                if irow != pivrow:
                    tableau[irow,:] = tableau[irow,:] - tableau[pivrow,:]*tableau[irow,pivcol]               
            
            print(tableau)
    
    if status == 0:
        print("Optimal solution found - single feasible solution.")
        #print("x* = ", x)
        #print("f* = ", tableau[-1,-1])
        return x
      
    elif status == 1:
        # if disp = True
        print("Optimal solution found:")
        x = np.zeros(n)
        sl = np.zeros(n_slack)
        sp = np.zeros(n_surplus)
        for j in range( n ):
            row = column_is_basic( j, tableau )
            if row >= 0:
                x[j] = tableau[row,-1]
        for j in range( n_slack ):
            row = column_is_basic( n+j, tableau )
            if row >= 0:
                sl[j] = tableau[row,-1]
        print("x* = ", x)
        print("slack* = ", sl)
        print("surplus* = ", sp)
        print("f* = ", tableau[-1,-1])
        return x    
        
    elif status == 2:
        # if disp = True
        print("No solution")
        
    
  
    
def linprog(c,A_eq=None,b_eq=None,A_ub=None,b_ub=None,A_lb=None,b_lb=None,type='max'):
    """
    
    """
    Aeq = np.asarray(A_eq) if not A_eq is None else np.empty([0,0])
    Aub = np.asarray(A_ub) if not A_ub is None else np.empty([0,0])
    Alb = np.asarray(A_lb) if not A_lb is None else np.empty([0,0])
    
    beq = np.asarray(b_eq) if not b_eq is None else np.empty([0,0])
    bub = np.asarray(b_ub) if not b_ub is None else np.empty([0,0])
    blb = np.asarray(b_lb) if not b_lb is None else np.empty([0,0])
    
    cc = np.asarray(c) 
    
    mub = len(bub)
    mlb = len(blb)
    meq = len(beq)
    n = len(c)
    
    try:
        Aub_rows, Aub_cols = Aub.shape
    except ValueError:
        raise ValueError("A_ub must be two-dimensional")
        
    try:    
        Alb_rows, Alb_cols = Alb.shape
    except ValueError:
        raise ValueError("A_lb must be two-dimensional")
    
    if not Aub_rows + Alb_rows == mub + mlb + meq:
        # raise here
        print("Number of rows in A must be equal to the size of b")
    if Aub_cols and not Aub_cols == n:
        # raise here
        print("Number of columns in A_ub must be equal to the size of c")
    if Alb_cols and not Alb_cols == n:
        # raise here
        print(Alb_cols,n)
        print("Number of columns in A_lb must be equal to the size of c")
        
    # Add slack variables for our upper-bound constraints
    n_slack = mub
    
    # Add surplus variables for our lower-bound constraints
    n_surplus = mlb   
    
    # Keep track of the artificial variables added
    n_artificial = meq
    
    # Create the tableau
    T = np.zeros([meq+mlb+mub+1,n+n_slack+n_surplus])
    
    # Insert objective into tableau
    if type.lower()[:3] == "max":
        T[-1,:n] = -cc # maximize
    elif type.lower()[:3] == "min":
        T[-1,:n] = cc # minimize
    else:
        raise ValueError("Unknown type for linear programming.  Must be one of 'max' or 'min'.  Got " + type)
    
    # b will become the RHS of the tableau    
    b = np.zeros( mlb+mub+1 )
    
    # a_rows tracks those rows of the tableau with artificial variables
    a_rows = None
    
    if meq > 0:
    
        # Add Aeq to the tableau
        T[:meq,:n] = Aeq
        
        # Add beq to the tableau
        b[:meq] = beq
        
    
    if mlb > 0:
    
        # Add Alb to the tableau
        T[meq:meq+mlb,:n] = Alb
        
        # Add blb to the tableau
        ##T[:mlb,-1] = blb
        b[meq:meq+mlb] = blb
        
        # Add surplus variables to the tableau
        np.fill_diagonal( T[meq:meq+mlb,n:n+n_surplus], -1)
        
    if mub > 0:    
    
        # Add Aub to the tableau
        T[mlb:mlb+mub,:n] = Aub
    
        # At bub to the tableau
        ##T[mlb:mlb+mub,-1] = bub
        b[mlb:mlb+mub] = bub
    
        # Add the slack variables to the tableau
        np.fill_diagonal( T[mlb:mlb+mub,n+n_surplus:n+n_surplus+n_slack], 1) 
        
    # determine if artificial variables are needed
    # artificial variables are needed if
    # 1. the row is an equality constraint
    # 2. the row has a slack/surplus variable with a sign opposite of the RHS
    a_rows = []
    for i in range(n_slack + n_surplus):
        if b[meq+i] * T[meq+i,n+i] < 0:
            n_artificial += 1
            a_rows.append(i)
            print("Slack %i needs an artificial variable" % (i,))
            print(a_rows)
    # Concatenate the artificial variable columns to the tableau
    T = np.hstack( [ T, np.zeros( [ T.shape[0], n_artificial] ) ] )
    c = 0
    for i in range(mlb+mub):
        if i in a_rows:
            T[ i, n+n_slack+n_surplus+c ] = 1
            c += 1            
            
    # Add the RHS column (b)
    T = np.hstack( [ T, np.reshape(b,[len(b),1]) ] ) 
            
    # If we have artificial variables, solve the phase 1 problem first
    if n_artificial > 0:
        T = np.vstack( [ T, np.zeros([ 1, T.shape[1] ]) ])
        T[-1,n+n_slack+n_surplus:-1] = 1
        
        # Put into correct form by ensuring RHS is entirely positive
        for i in range(T.shape[0]):
            if T[i,-1] < 0:
                T[i,:] *= -1
                
        
            
        # Make the artificial variables basic feasible variables by subtracting each
        # row with an artificial variable from the Phase 1 objective
        for r in a_rows:
            T[-1,:] = T[-1,:] - T[r,:]
            
        lpsimplex(T,n,n_surplus,n_slack,n_artificial,phase=1)
        
        # if pseudo objective is zero, remove the last row from the tableau and
        # proceed to phase 2
        if T[-1,-1] == 0:
            # Remove the pseudo-objective row from the tableau
            T = T[:-1,:]
            # Remove the artificial variable columns from the tableau
            T= np.delete(T, np.s_[n+n_slack+n_surplus:n+n_slack+n_surplus+n_artificial], 1)
        else:
            print("Problem appears to be infeasible")
            return
    
    # Tableau Finished
    lpsimplex(T,n, n_surplus, n_slack, n_artificial, phase=2)
    

        
    
    
if __name__ == "__main__":
    
    
    #http://www.dam.brown.edu/people/huiwang/classes/am121/Archive/simplex_121_c.pdf
    
    c = [3,2]
    b_ub = [10,8,4]
    A_ub = [[2,1],
            [1,1],
            [1,0]]
            
    linprog(c,A_ub=A_ub,b_ub=b_ub,type='max')
    
    print("\n\n\n")
    
    c = [60, 30, 20 ]
    A_ub = [[ 8,   6,   1],
            [ 4,   2, 1.5],
            [ 2, 1.5, 0.5],
            [ 0,   1,   0]]
    b_ub = [48,20,8,5]        
    
    linprog(c,A_ub=A_ub,b_ub=b_ub,type='max')
    
    print("\n\n\n")
       
    
    # Cycling
    #http://math.la.asu.edu/~checkman/Beyond4_2.pdf
    
    c = [10,-57,-9,-24]
    A_ub = [[ 0.5, -11/2, -5/2, 9],
            [ 0.5, -3/2, -1/2,  1],
            [   1,    0,    0,  0]]
    b_ub = [0,0,1]
    
    linprog(c,A_ub=A_ub,b_ub=b_ub,type='max')
    
    print("\n\n\n")
    
    
    # Klee-Minty  http://www.math.ubc.ca/~israel/m340/kleemin3.pdf
    print("The Klee-Minty Example")
    c = [100,10,1]
    A_ub = [[  1, 0, 0],
            [ 20, 1, 0],
            [200,20, 1]]
            
    b_ub = [1,100,10000]
    
    linprog(c,A_ub=A_ub,b_ub=b_ub,type='max')
    
    print("\n\n\n")
    
    
    # AAE550 Class 15 p 9
    
#    c = [1,1]
#    A = [[1,2],
#         [3,1]]
#    b = [1,1]
#    linprog(c,A_ub=A,b_ub=b)
    
    print("AAE550 Class 16 p 9")
    
    c = [ 2, -4 ]
    A_lb = [[2,1]]
    b_lb = [2]
    A_ub = [[1,1]]
    b_ub = [3]
    linprog(c,A_ub=A_ub,b_ub=b_ub,A_lb=A_lb,b_lb=b_lb,type="min")
    
    print("\n\n\n")
    
    #http://www.statslab.cam.ac.uk/~ff271/teaching/opt/notes/notes8.pdf
    #http://staff.science.uva.nl/~walton/Notes/TwoPhaseSimplex.pdf
    print("Two Phase")
    
    c = [1,2]
    A_lb = [[ 1, 1],
            [ 1,-1],
            [-1, 2]]
            
    b_lb = [4,1,-1]
    
    
    linprog(c,A_lb=A_lb,b_lb=b_lb,type='min')
    
    
    
    
    