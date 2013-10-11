"""
A basic linear programming function using a simplex method.

Functions
---------
.. autosummary::
   :toctree: generated/

    linprog
    lpsimplex
    verbose_callback
    terse_callback

"""

from __future__ import division, print_function, absolute_import

import numpy as np

from .optimize import Result, _check_unknown_options

__all__ = ['linprog','lpsimplex','verbose_callback','terse_callback']

__docformat__ = "restructuredtext en"


def verbose_callback(xk,**kwargs):
    """
    This is a sample callback for use with linprog, demonstrating the callback interface.
    This callback produces detailed output to sys.stdout before each iteration and after
    the final iteration of the simplex algorithm.

    Parameters
    ----------
    xk : array_like
        The current solution vector.
    **kwargs : dict
        A dictionary containing the following parameters:

        tableau : array_like
            The current tableau of the simplex algorithm.  Its structure is defined in lpsimplex.
        vars : tuple(str,...)
            Column headers for each column in tableau. "x[i]" for actual variables,
            "s[i]" for slack surplus variables, "a[i]" for artificial variables,
            and "RHS" for the constraint RHS vector
        phase : int
            The current Phase of the simplex algorithm (1 or 2)
        iter : int
            The current iteration number.
        pivot : tuple(int,int)
            The index of the tableau selected as the next pivot, or nan if no pivot exists
        basics : list[tuple(str,float)]
            A list of the current basic variables.  Each element contains the name of a basic variable and
            its value.
        complete : bool
            True if the simplex algorithm has completed (and this is the final call to callback), otherwise False.
    """
    tableau = kwargs["tableau"]
    columns = kwargs["vars"]
    iter = kwargs["iter"]
    pivrow,pivcol = kwargs["pivot"]
    phase = kwargs["phase"]
    basics = kwargs["basics"]
    complete = kwargs["complete"]

    t_rows,t_cols = tableau.shape

    saved_printoptions = np.get_printoptions()
    np.set_printoptions(linewidth=128,
                        formatter={'float':lambda x: "{: 12.4f}".format(x)})
    if complete:
        print("-------- Iteration Complete --------\n".format(iter))
    else:
        print("-------- Iteration {:d} --------\n".format(iter))
    print("Tableau:")
    print(" " + "".join(["{:>13}".format(columns[i]) for i in range(t_cols)]))
    print("" + str(tableau) + "\n")
    print("Pivot Element: T[{:.0f},{:.0f}]\n".format(pivrow,pivcol))
    print("Basic Variables:"),
    for basic_var in basics:
        print("{:<5s} = {: f}".format(basic_var[0],basic_var[1]))
    print()
    print("Current Solution:")
    print("x = ", xk)
    print()
    print("Current Objective Value:")
    print("f = ", tableau[-1,-1])
    print()
    np.set_printoptions(**saved_printoptions)


def terse_callback(xk, **kwargs):
    """
    This is a sample callback for use with linprog, demonstrating the callback interface.
    This callback produces brief output to sys.stdout before each iteration and after
    the final iteration of the simplex algorithm.

    Parameters
    ----------
    xk : array_like
        The current solution vector.
    **kwargs : dict
        A dictionary containing the following parameters:

        tableau : array_like
            The current tableau of the simplex algorithm.  Its structure is defined in lpsimplex.
        vars : tuple(str,...)
            Column headers for each column in tableau. "x[i]" for actual variables,
            "s[i]" for slack surplus variables, "a[i]" for artificial variables,
            and "RHS" for the constraint RHS vector
        phase : int
            The current Phase of the simplex algorithm (1 or 2)
        nit : int
            The current iteration number.
        pivot : tuple(int,int)
            The index of the tableau selected as the next pivot, or nan if no pivot exists
        basics : list[tuple(str,float)]
            A list of the current basic variables.  Each element contains the name of a basic variable and
            its value.
        complete : bool
            True if the simplex algorithm has completed (and this is the final call to callback), otherwise False.
    """
    iter = kwargs["iter"]

    if iter == 0:
        print("Iter:   X:")
    print("{: <5d}   ".format(iter),end="")
    print(xk)


def lpsimplex(tableau,n,n_slack,n_artificial,maxiter=1000,phase=2,callback=None,nit0=0):
    """
    Solve a linear programming problem in "standard maximization form" using the Simplex Method.

    Maximize :math:`f = c^T x`

    subject to

    .. math::

        Ax = b
        x_i >= 0
        b_j >= 0

    Parameters
    ----------
    tableau : array_like
        A 2-D array representing the simplex tableau corresponding to the maximization problem.
        It should have the form:

        [[A[0,0], A[0,1], ..., A[0,n_total], b[0]],
         [A[1,0], A[1,1], ..., A[1,n_total], b[1]],
         .
         .
         .
         [A[m,0], A[m,1], ..., A[m,n_total], b[m]],
         [c[0],   c[1], ...,   c[n_total],    0]]

        for a Phase 2 problem, or the form:

        [[A[0,0], A[0,1], ..., A[0,n_total], b[0]],
         [A[1,0], A[1,1], ..., A[1,n_total], b[1]],
         .
         .
         .
         [A[m,0], A[m,1], ..., A[m,n_total], b[m]],
         [c[0],   c[1], ...,   c[n_total],   0],
         [c'[0],  c'[1], ...,  c'[n_total],  0]]

         for a Phase 1 problem (a Problem in which a basic feasible solution is sought
         prior to maximizing the actual objective. n_total is the total of all variables:

    n : int
        The number of true variables in the problem.
    n_slack : int
        The number of slack/surplus variables in the problem.
    n_artificial : int
        The number of artificial variables in the problem.
    maxiter : int
        The maximum number of iterations to perform before aborting the optimization.
    phase : int
        The phase of the optimization being executed.  In phase 1 a basic feasible solution is sought and
        the tableau has an additional row representing an alternate objective function.
    callback : callable
        If a callback function is provided, it will be called within each iteration of the simplex algorithm.
        The callback must have the signature `callback(xk,**kwargs)` where xk is the current solution vector
        and kwargs is a dictionary containing the following::
        "tableau" : The current Simplex algorithm tableau
        "nit" : The current iteration.
        "pivot" : The pivot (row,column) used for the next iteration.
        "phase" : Whether the algorithm is in Phase 1 or Phase 2.
        "bv" : A structured array containing a string representation of each basic variable and its current value.
    nit0 : int
        The initial iteration number used to keep an accurate iteration total in a two-phase problem.

    Returns
    -------
    res : Result
        The optimization result represented as a ``Result`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. Possible
        values for the ``status`` attribute are:
        -1 : Invalid arguments
         0 : Optimization terminated successfully
         1 : Optimization terminated successfully, single feasible solution
         2 : Iteration limit reached
         3 : Problem appears to be infeasible
         4 : Problem appears to be unbounded

        See `Result` for a description of other attributes.
    """

    nit = nit0

    if phase not in (1,2):
        status = -1
        pass

    # Make a copy of the tableau that will be used to detect cycling
    tableau0 = np.empty_like(tableau)
    tableau0[:,:] = tableau

    t_rows, t_cols = tableau.shape

    pivrow = 0
    pivcol = 0

    if phase == 2 and t_cols == t_rows:
        # This is just Ax=b
        x = np.linalg.solve(tableau[:-1,:-1], tableau[:-1,-1])
        status = 1
        message = "Optimization terminated successfully.  (Single feasible solution)"
        return x, status, message, status
    else:
        status = -2
        message = ""
        cycle = 0
        x = np.zeros([n])

        if not callback is None:
            index_to_varname = {}
            for i in range(n):
                index_to_varname[i] = 'x[{:d}]'.format(i)
            for i in range(n_slack):
                index_to_varname[n+i] = 's[{:d}]'.format(i)
            for i in range(n_artificial):
                index_to_varname[n+n_slack+i] = 'a[{:d}]'.format(i)

            column_headers = [index_to_varname[i] for i in range(t_cols - 1)] + ["RHS"]

        while nit < maxiter and status == -2:

            if nit > 2**n:
                if np.all(tableau == tableau0):
                    cycle += 1

            # Find the most negative value in bottom row of tableau
            pivcol = np.argmin(tableau[-1,:-1])

            # sort the pivot column so that later we don't bother with rows where the value
            # in pivot column is negative
            if phase == 1:
                pivcol_sortorder = np.argsort(tableau[:-2, pivcol])[:,np.newaxis]
            else:
                pivcol_sortorder = np.argsort(tableau[:-1, pivcol])[:,np.newaxis]

            # Now row quotient has three columns: the original row index,
            # the value in the pivot column, and the RHS value
            row_quotient = np.hstack([pivcol_sortorder, tableau[pivcol_sortorder, pivcol], tableau[pivcol_sortorder, -1]])
            # Pare down row_quotient by removing those rows where the value in the pivot column is non-positive
            row_quotient = row_quotient[row_quotient[:,1] > 0]

            # Replace the 2nd column the ratio of the RHS value / pivot column value
            row_quotient[:,1] = row_quotient[:,2] / row_quotient[:,1]

            # Sort row quotient
            # With the negative column values removed, we now want to use the row with the minimum
            # quotient (in column 1) as the pivot row.
            try:
                row_quotient = row_quotient[row_quotient[:,1].argsort(),:1]
                pivrow = row_quotient[cycle][0]
            except IndexError:
                pivrow = np.nan
                if cycle > 0:
                    message = "Optimization failed. The problem appears to be unbounded. Unable to recover from cycling."
                status = 4
                message = "Optimization failed. The problem appears to be unbounded."
                break

            if cycle > 0:
                cycle = 0

            pivval = tableau[pivrow,pivcol]
            tableau[pivrow,:] = tableau[pivrow,:] / pivval
            for irow in range(tableau.shape[0]):
                if irow != pivrow:
                    tableau[irow,:] = tableau[irow,:] - tableau[pivrow,:]*tableau[irow,pivcol]

            if np.all(tableau[-1,:-1] >= 0):
                status = 0

            if not callback is None:
                bv_map = np.sum(tableau[:,:-1] != 0, 0)  # bv_map is the sum of the nonzero elements in each column
                basic_cols = np.where(bv_map == 1)[0]    # Columns with basic variables
                basic_vars = []
                x.fill(0.0)
                for i in basic_cols:
                    nonzero_row = np.nonzero(tableau[:-1,i])[0][0]
                    basic_vars.append((column_headers[i], tableau[nonzero_row,-1] / tableau[nonzero_row,i]))

                    if i < n:
                        x[i] = basic_vars[-1][1]

                callback(x, **{"tableau": tableau,
                                "iter":nit,
                                "vars":column_headers,
                                "pivot":(pivrow,pivcol),
                                "phase":phase,
                                "basics":basic_vars,
                                "complete":status != -2 and phase == 2})
            nit += 1

        else:
            if nit >= maxiter:
                message = "Iteration limit reached."
                status = 2
            else:
                message = "Optimization terminated successfully."
                status = 0

    bv_map = np.sum(tableau[:,:-1] != 0, 0)  # bv_map is the sum of the nonzero elements in each column
    basic_cols = np.where(bv_map == 1)[0]    # Columns with basic variables
    x.fill(0.0)
    for col in basic_cols:
        if col < n:
            nonzero_row = np.nonzero(tableau[:-1,col])[0][0]
            x[col] = tableau[nonzero_row,-1] / tableau[nonzero_row,col]

    return x, nit, status, message


def linprog(c,A_eq=None,b_eq=None,A_lb=None,b_lb=None,A_ub=None,b_ub=None,objtype='max',maxiter=1000,disp=False,callback=None):
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
    objtype : str
        The type of objective function represented by c.  Must be either 'max' (default) or 'min'
    maxiter : int
       The maximum number of iterations to perform.
    disp : bool
        If True, print exit status message to sys.stdout
    callback : callable
        If a callback function is provide, it will be called within each iteration of the simplex algorithm.
        The callback must have the signature `callback(xk,**kwargs)` where xk is the current solution vector
        and kwargs is a dictionary containing the following::
        "tableau" : The current Simplex algorithm tableau
        "nit" : The current iteration.
        "pivot" : The pivot (row,column) used for the next iteration.
        "phase" : Whether the algorithm is in Phase 1 or Phase 2.
        "bv" : A structured array containing a string representation of each basic variable and its current value.

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
    message = "Optimization terminated successfully"
    nit1 = 0  # Iterations used in Phase 1

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
        message = "Invalid input.  A_ub must be two-dimensional"

    try:
        Alb_rows, Alb_cols = Alb.shape
    except ValueError:
        status = -1
        message = "Invalid input.  A_lb must be two-dimensional"

    try:
        Aeq_rows, Aeq_cols = Aeq.shape
    except ValueError:
        status = -1
        message = "Invalid input.  A_eq must be two-dimensional"

    if Aeq_rows != meq:
        status = -1
        message = "Invalid input.  The number of rows in A_eq must be equal to the number of values in b_eq"

    if Alb_rows != mlb:
        status = -1
        message = "Invalid input.  The number of rows in A_lb must be equal to the number of values in b_lb"

    if Aub_rows != mub:
        status = -1
        message = "Invalid input.  The number of rows in A_ub must be equal to the number of values in b_ub"

    if Aeq_cols != n:
        status = -1
        message = "Invalid input.  Number of columns in A_eq must be equal to the size of c"

    if Alb_cols != n:
        status = -1
        message = "Invalid input.  Number of columns in A_lb must be equal to the size of c"

    if Aub_cols != n:
        status = -1
        message = "Invalid input.  Number of columns in A_ub must be equal to the size of c"

    # Add slack variables for our upper-bound constraints
    n_slack = mub

    # Add surplus variables for our lower-bound constraints
    n_surplus = mlb

    # Add artificial variables for the equality constraints (more may be added later)
    n_artificial = meq

    # Create the tableau
    T = np.zeros([meq+mlb+mub+1,n+n_slack+n_surplus])

    # Insert objective into tableau
    if objtype.lower()[:3] == "max":
        T[-1,:n] = -cc  # maximize
    elif objtype.lower()[:3] == "min":
        T[-1,:n] = cc  # minimize
    else:
        status = -1
        message = "Invalid input. Argument 'objtype' Must be one of 'max' or 'min'.  Got " + str(objtype)

    if status == -1:
        # b will become the RHS of the tableau
        b = np.zeros(meq+mlb+mub+1)

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
            np.fill_diagonal(T[meq:meq+mlb,n:n+n_surplus], -1)
        if mub > 0:
            # Add Aub to the tableau
            T[meq+mlb:meq+mlb+mub,:n] = Aub
            # At bub to the tableau
            ##T[mlb:mlb+mub,-1] = bub
            b[meq+mlb:meq+mlb+mub] = bub
            # Add the slack variables to the tableau
            np.fill_diagonal(T[meq+mlb:meq+mlb+mub,n+n_surplus:n+n_surplus+n_slack], 1)

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
        T = np.hstack([T, np.zeros([T.shape[0], n_artificial])])
        c = 0
        for i in range(mlb+mub):
            if i in a_rows:
                T[i, n+n_slack+n_surplus+c] = 1
                c += 1

        # Add the RHS column (b)
        T = np.hstack([T, np.reshape(b,[len(b),1])])

        # If we have artificial variables, solve the phase 1 problem first
        if n_artificial > 0:
            T = np.vstack([T, np.zeros([1, T.shape[1]])])
            T[-1,n+n_slack+n_surplus:-1] = 1

            # Put into correct form by ensuring RHS is entirely positive
            for i in range(T.shape[0]):
                if T[i,-1] < 0:
                    T[i,:] *= -1

            # Make the artificial variables basic feasible variables by subtracting each
            # row with an artificial variable from the Phase 1 objective
            for r in a_rows:
                T[-1,:] = T[-1,:] - T[r,:]

            x, nit1, status, message = lpsimplex(T,n,n_surplus+n_slack,n_artificial,phase=1,callback=callback,maxiter=maxiter)

            # if pseudo objective is zero, remove the last row from the tableau and
            # proceed to phase 2
            if T[-1,-1] == 0:
                # Remove the pseudo-objective row from the tableau
                T = T[:-1,:]
                # Remove the artificial variable columns from the tableau
                T = np.delete(T,np.s_[n+n_slack+n_surplus:n+n_slack+n_surplus+n_artificial],1)
            else:
                print("Optimization Failed.  The problem appears to be infeasible")
                status = 2
                message = "Unable to find a feasible starting point.  The problem appears to be infeasible"
                return

        # Tableau Finished
        x, nit2, status, message = lpsimplex(T,n,n_surplus+n_slack,n_artificial,maxiter=maxiter,
                                             phase=2,callback=callback,nit0=nit1)

    # Optimization complete at this point
    obj_mult = -1 if objtype == "min" else 1

    if status in (0,1):
        fstar = T[-1,-1] * obj_mult
        if disp:
            print(message)
            print("         Current function value: {: <12.6f}".format(fstar))
            print("         Iterations: {:d}".format(nit2))
    else:
        if disp:
            print(message)
            print("         Iterations: {:d}".format(nit2))

    return Result(x=x, fun=T[-1,-1]*obj_mult, nit=int(nit2), status=int(status),
                  message=message, success=(status in (0,1)))


if __name__ == "__main__":
    #http://www.dam.brown.edu/people/huiwang/classes/am121/Archive/simplex_121_c.pdf

    c = [3,2]
    b_ub = [10,8,4]
    A_ub = [[2,1],
            [1,1],
            [1,0]]

    print(linprog(c,A_ub=A_ub,b_ub=b_ub,objtype='max',disp=False))
    #
    #print("\n\n\n")
    #
    #c = [60, 30, 20]
    #A_ub = [[8,   6,   1],
    #        [4,   2, 1.5],
    #        [2, 1.5, 0.5],
    #        [0,   1,   0]]
    #b_ub = [48,20,8,5]
    #
    #linprog(c,A_ub=A_ub,b_ub=b_ub,type='max')
    #
    #print("\n\n\n")
    #
    # Cycling
    #http://math.la.asu.edu/~checkman/Beyond4_2.pdf
    #
    #c = [10,-57,-9,-24]
    #A_ub = [[0.5, -11/2, -5/2, 9],
    #        [0.5, -3/2, -1/2,  1],
    #        [1,    0,    0,  0]]
    #b_ub = [0,0,1]
    #
    #linprog(c,A_ub=A_ub,b_ub=b_ub,type='max',callback=verbose_callback)
    #
    #print("\n\n\n")
    #
    #import sys
    #
    # Two-Phase example
    c = [6,3]
    A_lb = [[1, 1],
            [2,-1]]
    b_lb = [1,1]
    A_ub = [[0,3]]
    b_ub = [2]

    #http://www.statslab.cam.ac.uk/~ff271/teaching/opt/notes/notes8.pdf
    linprog(c,A_ub=A_ub,b_ub=b_ub,A_lb=A_lb,b_lb=b_lb,objtype='min',callback=terse_callback,disp=True)

    import sys
    sys.exit(0)

    #unbounded example
    linprog(c,A_ub=A_ub,b_ub=b_ub,A_lb=A_lb,b_lb=b_lb,objtype='max',callback=verbose_callback,disp=True)

    # Klee-Minty  http://www.math.ubc.ca/~israel/m340/kleemin3.pdf
    print("The Klee-Minty Example")
    c = [100,10,1]
    A_ub = [[1, 0, 0],
            [20, 1, 0],
            [200,20, 1]]

    b_ub = [1,100,10000]

    linprog(c,A_ub=A_ub,b_ub=b_ub,objtype='max',disp=2)

    print("\n\n\n")

    # AAE550 Class 15 p 9

#    c = [1,1]
#    A = [[1,2],
#         [3,1]]
#    b = [1,1]
#    linprog(c,A_ub=A,b_ub=b)

    print("AAE550 Class 16 p 9")

    c = [2, -4]
    A_lb = [[2,1]]
    b_lb = [2]
    A_ub = [[1,1]]
    b_ub = [3]
    linprog(c,A_ub=A_ub,b_ub=b_ub,A_lb=A_lb,b_lb=b_lb,objtype="min",disp=2)

    print("\n\n\n")

    #http://www.statslab.cam.ac.uk/~ff271/teaching/opt/notes/notes8.pdf
    #http://staff.science.uva.nl/~walton/Notes/TwoPhaseSimplex.pdf
    print("Two Phase")

    c = [1,2]
    A_lb = [[1, 1],
            [1,-1],
            [-1,2]]

    b_lb = [4,1,-1]

    linprog(c,A_lb=A_lb,b_lb=b_lb,objtype='min',disp=2)

    print("\n\n\n")
    print("Equality constraint example")
    # http://www.sce.carleton.ca/faculty/chinneck/po/Chapter5.pdf

    c = [15,10]
    A_ub = [[1,0],
            [0,1]]
    b_ub = [2,3]
    A_eq = [[1,1]]
    b_eq = [4]

    linprog(c,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,objtype='max',disp=2)
