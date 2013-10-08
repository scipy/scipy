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
    
    A column is considered canonical if all but one row in it contain zeros.
    
    If the column is canonical, return the row index whose value is nonzero.
    If the column is not canonical, return -1.
    """    
    if np.count_nonzero( tableau[:-1,j] ) == 1:
        return np.nonzero(tableau[:-1,j])[0][0]
    return -1



def lpsimplex(tableau,n,n_surplus,n_slack,n_artificial,phase=2,disp=0):
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

    nit = 0
    status = 0
    message = ""

    # Make a copy of the tableau that will be used to detect cycling
    T0 = np.empty_like(tableau)
    T0[:,:] = tableau
    
    # TODO FIGURE OUT WHEN TO JUST SOLVE AX=B
    t_rows, t_cols = tableau.shape
    
    if phase == 2 and t_cols == t_rows:
        # This is just Ax=b
        x = np.linalg.solve( tableau[:-1,:-1], tableau[:-1,-1] )
        status = 0
        message = "Optimal solution found.  Problem has a single feasible solution."
        return x, status, message, 0
    else:    
        nit = 0
        status = 1
        cycle = 0
        
        while np.any( tableau[-1,:-1] < 0 ):

            if disp >= 2:
                print("Iteration: ",nit)
                print("Tableau")
                print(tableau)

            if nit > 2**n:
                if np.all(tableau == T0):
                    if disp >= 2:
                        print("Cycling detected")
                    cycle += 1



            # Find the most negative value in bottom row of tableau
            pivcol = np.argmin(tableau[-1,:-1])
        
            # Find the minimizer of b[x_to_enter]/a[i,x_to_enter]
            # temporarily disable divide by zero warnings
            err_save = np.geterr()['divide']
            np.seterr(divide='ignore')
            if phase == 1:
                row_quotient = np.empty([ t_rows-2, 2])
                row_quotient[:,0] = np.arange(t_rows-2)
                row_quotient[:,1] = tableau[:-2,-1] / tableau[:-2,pivcol]
            else:
                row_quotient = np.empty([ t_rows-1, 2])
                row_quotient[:,0] = np.arange(t_rows-1)
                row_quotient[:,1] = tableau[:-1,-1] / tableau[:-1,pivcol]
            np.seterr(divide=err_save)

            # Sort row quotient, we want to use the minimizer of the quotient, as long
            # as that pivot element is positive (if negative, go to the next most negative quotient)
            row_quotient = row_quotient[row_quotient[:,1].argsort()]
            
            for i in range(cycle, t_rows-1):
                if tableau[ row_quotient[i][0], pivcol ] > 0:
                    pivrow = row_quotient[i][0]
                    if cycle > 0:
                        cycle = 0
                    break
            else:
                message = "The problem appears to be unbounded."
                status = 3
                break
            
            pivval = tableau[pivrow,pivcol]
            tableau[pivrow,:] = tableau[pivrow,:] / pivval
            for irow in range(tableau.shape[0]):
                if irow != pivrow:
                    tableau[irow,:] = tableau[irow,:] - tableau[pivrow,:]*tableau[irow,pivcol]               


            nit += 1
        else:
            if disp >= 2:
                print("Iteration Complete.")
                print("Final Tableau:")
                print(tableau)
      
    if status == 1:
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
        
    
  
    
def linprog(c,A_eq=None,b_eq=None,A_lb=None,b_lb=None,A_ub=None,b_ub=None,type='max',disp=0):
    """
    Solve the following linear programming problem via a two-phase simplex algorithm.

    maximize:     c^T * x

    subject to:   A_eq * x == b_eq
                  A_lb * x >= b_lb
                  A_ub * x <= b_ub

    Parameters
    ----------
    c : array_like
        Coefficients of the linear objective function to be maximized.
    A_eq : array_like
        2-D array which, when matrix-multiplied by x, gives the values of the equality constraints at x.
    b_eq : array_like
        1-D array of values representing the RHS of each equality constraint (row) in A_eq.
    A_lb :
        2-D array which, when matrix-multiplied by x, gives the values of the lower-bound inequality constraints at x.
    b_lb : array_like
        1-D array of values representing the lower-bound of each inequality constraint (row) in A_lb.
    A_ub :
        2-D array which, when matrix-multiplied by x, gives the values of the upper-bound inequality constraints at x.
    b_ub : array_like
        1-D array of values representing the upper-bound of each inequality constraint (row) in A_ub.
    type : str
        The type of optimization performed.  Must be either 'max' (default) or 'min'
    disp : int
        The verbosity of the linprog routine::
        0 : Silent operation
        1 : Output a summary upon exit
        2 : Output the status of the optimization at each iteration

    Returns
    -------
    x : ndarray
        The independent variable vector which optimizes the linear programming problem.
    success : bool
        Returns True if the algorithm succeeded in finding an optimal solution.
    status : int
        An integer representing the exit status of the optimization::
       -1 : Invalid arguments
        0 : Optimization terminated successfully.
        1 : Optimization terminated successfully, single feasible solution.
        2 : Problem appears to be infeasible
        3 : Problem appears to be unbounded
    nit : int
        The number of iterations performed.
    message : str
        A string descriptor of the exit status of the optimization.
    bv : tuple
        The basic variables.
    nbv : tuple
        The nonbasic variables.
    """
    status = 0
    message = "Optimization Successful"

    Aeq = np.asarray(A_eq) if not A_eq is None else np.empty([0,0])
    Aub = np.asarray(A_ub) if not A_ub is None else np.empty([0,0])
    Alb = np.asarray(A_lb) if not A_lb is None else np.empty([0,0])
    
    beq = np.ravel(np.asarray(b_eq)) if not b_eq is None else np.empty([0])
    bub = np.ravel(np.asarray(b_ub)) if not b_ub is None else np.empty([0])
    blb = np.ravel(np.asarray(b_lb)) if not b_lb is None else np.empty([0])
    
    cc = np.asarray(c) 
    
    mub = len(bub)
    mlb = len(blb)
    meq = len(beq)
    n = len(c)
    
    try:
        Aub_rows, Aub_cols = Aub.shape
    except ValueError:
        status = -1
        message = "A_ub must be two-dimensional"
        
    try:    
        Alb_rows, Alb_cols = Alb.shape
    except ValueError:
        status = -1
        message = "A_lb must be two-dimensional"
        
    try:    
        Aeq_rows, Aeq_cols = Aeq.shape
    except ValueError:
        status = -1
        message = "A_eq must be two-dimensional"

    if Aeq_rows != meq:
        status = -1
        message = "The number of rows in A_eq must be equal to the number of values in b_eq"

    if Alb_rows != mlb:
        status = -1
        message = "The number of rows in A_lb must be equal to the number of values in b_lb"

    if Aub_rows != mub:
        status = -1
        message = "The number of rows in A_ub must be equal to the number of values in b_ub"

    if Aeq_cols != n:
        status = -1
        message = "Number of columns in A_eq must be equal to the size of c"

    if Alb_cols != n:
        status = -1
        message = "Number of columns in A_lb must be equal to the size of c"

    if Aub_cols != n:
        status = -1
        message = "Number of columns in A_ub must be equal to the size of c"

        
    # Add slack variables for our upper-bound constraints
    n_slack = mub
    
    # Add surplus variables for our lower-bound constraints
    n_surplus = mlb   
    
    # Add artificial variables for the equality constraints (more may be added later)
    n_artificial = meq
    
    # Create the tableau
    T = np.zeros([meq+mlb+mub+1,n+n_slack+n_surplus])
    
    # Insert objective into tableau
    if type.lower()[:3] == "max":
        T[-1,:n] = -cc # maximize
    elif type.lower()[:3] == "min":
        T[-1,:n] = cc # minimize
    else:
        status = -1
        message = "Unknown type for linear programming.  Must be one of 'max' or 'min'.  Got " + type


    if status == -1:
        x = np.empty([n])
        x.fill(np.nan)
        #return Result(status=status, success=False, message=message, x=x)


    # b will become the RHS of the tableau    
    b = np.zeros( meq+mlb+mub+1 )

    if meq > 0:
    
        # Add Aeq to the tableau
        T[:meq,:n] = Aeq
        
        # Add beq to the tableau
        b[:meq] = beq
        
    
    if mlb > 0:
    
        # Add Alb to the tableau
        T[meq:meq+mlb,:n] = Alb
        
        # Add blb to the tableau
        b[meq:meq+mlb] = blb
        
        # Add surplus variables to the tableau
        np.fill_diagonal( T[meq:meq+mlb,n:n+n_surplus], -1)
        
    if mub > 0:    
    
        # Add Aub to the tableau
        T[meq+mlb:meq+mlb+mub,:n] = Aub
    
        # At bub to the tableau
        ##T[mlb:mlb+mub,-1] = bub
        b[meq+mlb:meq+mlb+mub] = bub
    
        # Add the slack variables to the tableau
        np.fill_diagonal( T[meq+mlb:meq+mlb+mub,n+n_surplus:n+n_surplus+n_slack], 1) 
        
    # determine if artificial variables are needed
    # artificial variables are needed if
    # 1. the row is an equality constraint
    # 2. the row has a slack/surplus variable with a sign opposite of the RHS
    # a_rows tracks those rows of the tableau with artificial variables
    a_rows = range(meq)
    for i in range(n_slack + n_surplus):
        if b[meq+i] * T[meq+i,n+i] < 0:
            n_artificial += 1
            a_rows.append(i)

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
            
        lpsimplex(T,n,n_surplus,n_slack,n_artificial,phase=1,disp=disp)
        
        # if pseudo objective is zero, remove the last row from the tableau and
        # proceed to phase 2
        if T[-1,-1] == 0:
            # Remove the pseudo-objective row from the tableau
            T = T[:-1,:]
            # Remove the artificial variable columns from the tableau
            T= np.delete(T, np.s_[n+n_slack+n_surplus:n+n_slack+n_surplus+n_artificial], 1)
        else:
            print("Problem appears to be infeasible")
            status = 2
            message = "Unable to find a feasible starting point.  The problem appears to be infeasible"
            return
    
    # Tableau Finished
    lpsimplex(T,n, n_surplus, n_slack, n_artificial, phase=2,disp=disp)
    

        
    
    
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
    
    
    linprog(c,A_lb=A_lb,b_lb=b_lb,type='min',disp=2)
    
    print("\n\n\n")
    print("Equality constraint example")
    # http://www.sce.carleton.ca/faculty/chinneck/po/Chapter5.pdf
    
    c = [15,10]
    A_ub = [[1,0],
            [0,1]]
    b_ub = [2,3]
    A_eq = [[1,1]]
    b_eq = [4]
    
    linprog(c,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,type='max',disp=2)
    
    
    
    
    