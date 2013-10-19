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
import numpy.ma as ma

from .optimize import Result, _check_unknown_options

__all__ = ['linprog','lpsimplex','liprog_verbose_callback','linprog_terse_callback']

__docformat__ = "restructuredtext en"


def linprog_verbose_callback(xk,**kwargs):
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
    np.set_printoptions(linewidth=256,
                        formatter={'float':lambda x: "{: 12.4f}".format(x)})
    if complete:
        print("--------- Iteration Complete - Phase {:d} -------\n".format(phase))
        print("Tableau:")
    elif iter == 0:
        print("--------- Initial Tableau - Phase {:d} ----------\n".format(phase))

    else:
        print("--------- Iteration {:d}  - Phase {:d} --------\n".format(iter,phase))
        print("Tableau:")

    if iter >= 0:
        print(" " + "".join(["{:>13}".format(columns[i]) for i in range(t_cols)]))
        print("" + str(tableau) + "\n")
        if not complete:
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


def linprog_terse_callback(xk, **kwargs):
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


def lpsimplex(tableau,n,n_slack,n_artificial,maxiter=1000,phase=2,callback=None,tol=1.0E-12,nit0=0):
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
    tol : float
        The tolerance which determines when a solution is "close enough" to zero in Phase 1 to be considered
        a basic feasible solution or close enough to positive to to serve as an optimal solution.
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
         1 : Iteration limit reached
         2 : Problem appears to be infeasible
         3 : Problem appears to be unbounded

        See `Result` for a description of other attributes.
    """

    nit = nit0

    if phase not in (1,2):
        status = -1
        message = "Invalid input to lpsimplex.  Phase must be 1 or 2."

    # Make a copy of the tableau that will be used to detect cycling
    tableau0 = np.empty_like(tableau)
    tableau0[:,:] = tableau

    t_rows, t_cols = tableau.shape

    pivrow = 0
    pivcol = 0

    status = None
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

    while nit < maxiter and status is None:

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
            status = 3
            message = "Optimization failed. The problem appears to be unbounded."
            break

        if cycle > 0:
            cycle = 0

        if np.all(tableau[-1,:-1] >= -tol):
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
                            "complete": (not status is None) and phase == 2})

        if status is None:
            pivval = tableau[pivrow,pivcol]
            tableau[pivrow,:] = tableau[pivrow,:] / pivval
            for irow in range(tableau.shape[0]):
                if irow != pivrow:
                    tableau[irow,:] = tableau[irow,:] - tableau[pivrow,:]*tableau[irow,pivcol]

            nit += 1

    else:
        if nit >= maxiter:
            message = "Iteration limit reached."
            status = 1
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


def linprog(c,A_eq=None,b_eq=None,A_lb=None,b_lb=None,A_ub=None,b_ub=None,
            bounds=None,objtype='max',maxiter=1000,disp=False,callback=None,tol=1.0E-12):
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
    bounds : array_like
        The bounds for each independent variable in the solution, which can take one of three forms::
        None : The default bounds, all variables are restricted to be non-negative.
        (lb,ub) : If a 2-element sequence is provided, the same lower bound (lb) and upper bound (ub) will be
                  applied to all variables.
        [(lb_0,ub_0),(lb_1,ub_1),...] : If an n x 2 sequence is provided, each variable x_i will be bounded by lb_i
                  and ub_i.
        Infinite bounds are specified using -np.inf (negative) or np.inf (positive).
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
    tol : float
        The tolerance which determines when a solution is "close enough" to zero in Phase 1 to be considered
        a basic feasible solution or close enough to positive to to serve as an optimal solution.

    Returns
    -------
    x : ndarray
        The independent variable vector which optimizes the linear programming problem.
    success : bool
        Returns True if the algorithm succeeded in finding an optimal solution.
    status : int
        An integer representing the exit status of the optimization::
        -1 : Invalid arguments
         0 : Optimization terminated successfully
         1 : Iteration limit reached
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
    message = ""
    nit1 = 0  # Iterations used in Phase 1
    nit2 = 0  # Iterations used after Phase 2
    have_floor_variable = False

    cc = np.asarray(c)

    # The initial value of the objective function element in the tableau
    f0 = 0

    # The number of variables as given by c
    n = len(c)

    Aeq_in = np.asarray(A_eq) if not A_eq is None else np.empty([0,len(cc)])
    Aub_in = np.asarray(A_ub) if not A_ub is None else np.empty([0,len(cc)])
    Alb_in = np.asarray(A_lb) if not A_lb is None else np.empty([0,len(cc)])

    beq_in = np.ravel(np.asarray(b_eq)) if not b_eq is None else np.empty([0])
    bub_in = np.ravel(np.asarray(b_ub)) if not b_ub is None else np.empty([0])
    blb_in = np.ravel(np.asarray(b_lb)) if not b_lb is None else np.empty([0])

    # Analyze the bounds and determine what modifications to me made to the constraints in order to accommodate them.
    # http://www.sce.carleton.ca/faculty/chinneck/po/Chapter5.pdf
    # 1. if lb == 0 and ub == inf:
    #        do nothing
    # 2. if abs(ub) != inf:
    # 3. if abs(lb) != inf:
    # 4. if lb is -inf:
    # 5. if lb == inf or ub == -inf:
    #        invalid inputs
    L = np.zeros(n,dtype=np.float64)
    U = np.ones(n,dtype=np.float64)*np.inf
    if bounds is None or len(bounds) == 0:
        pass
    elif len(bounds) == 1:
        # All bounds are the same
        L = np.asarray(n*bounds[0][0],dtype=np.float64)
        U = np.asarray(n*bounds[0][1],dtype=np.float64)
    else:
        if len(bounds) != n:
            status = -1
            message = "Invalid input.  Length of bounds is inconsistent with the length of c"
        else:
            try:
                L = np.asarray([val[0] for val in bounds],dtype=np.float64)
                U = np.asarray([val[1] for val in bounds],dtype=np.float64)
            except IndexError as err:
                status = -1
                message = "Invalid input.  bounds must be a n x 2 sequence/array where n = len(c)."

    for i in range(n):
        if L[i] is None:
            L[i] = -np.inf
        if U[i] is None:
            U[i] = np.inf

    if np.any(L == -np.inf):
        # If any lower-bound constraint is a free variable
        # add the first column variable as the "floor" variable which
        # accommodates the most negative variable in the problem.
        n = n + 1
        L = np.concatenate([np.array([0]),L])
        U = np.concatenate([np.array([np.inf]),U])
        cc = np.concatenate([np.array([0]),cc])
        Aeq_in = np.hstack([np.zeros([Aeq_in.shape[0],1]),Aeq_in])
        Alb_in = np.hstack([np.zeros([Alb_in.shape[0],1]),Alb_in])
        Aub_in = np.hstack([np.zeros([Aub_in.shape[0],1]),Aub_in])
        have_floor_variable = True

    # Now before we deal with any variables with lower bounds < 0,
    # deal with finite bounds which can be simply added as new constraints.
    # Also validate bounds inputs here.
    for i in range(n):
        if(L[i] > U[i]):
            status = -1
            message = "Invalid input.  Lower bound {:d} is greater than upper bound {:d}".format(i,i)

        if np.isinf(L[i]) and L[i] > 0:
            status = -1
            message = "Invalid input.  Lower bound may not be +infinity"

        if np.isinf(U[i]) and U[i] < 0:
            status = -1
            message = "Invalid input.  Upper bound may not be -infinity"

        if np.isfinite(L[i]) and L[i] > 0:
            # Add a new lower-bound constraint
            Alb_in = np.vstack([Alb_in, np.zeros(n)])
            Alb_in[-1,i] = 1
            blb_in = np.concatenate([blb_in,np.array([L[i]])])
            L[i] = 0

        if np.isfinite(U[i]):
            # Add a new upper-bound constraint
            Aub_in = np.vstack([Aub_in, np.zeros(n)])
            Aub_in[-1,i] = 1
            bub_in = np.concatenate([bub_in,np.array([U[i]])])
            U[i] = np.inf

    # Now find negative lower bounds (finite or infinite) which require a change of variables or free variables
    # and handle them appropriately
    for i in range(0,n):

        if L[i] < 0:
            if np.isfinite(L[i]) and L[i] < 0:
                # Add a change of variables for x[i]
                # For each row in the constraint matrices, we take the coefficient from column i in A,
                # and subtract the product of that and L[i] to the RHS b
                beq_in[:] = beq_in[:] - Aeq_in[:,i] * L[i]
                bub_in[:] = bub_in[:] - Aub_in[:,i] * L[i]
                blb_in[:] = blb_in[:] - Alb_in[:,i] * L[i]
                # We now have a nonzero initial value for the objective function as well.
                f0 = f0 - cc[i] * L[i]
            else:
                # This is an unrestricted variable, let x[i] = u[i] - v[0] where v is the first column in all matrices.
                Aeq_in[:,0] = Aeq_in[:,0] - Aeq_in[:,i]
                Alb_in[:,0] = Alb_in[:,0] - Alb_in[:,i]
                Aub_in[:,0] = Aub_in[:,0] - Aub_in[:,i]
                cc[0] = cc[0] - cc[i]

        if np.isinf(U[i]):
            if U[i] < 0:
                status = -1
                message = "Invalid input.  Upper bound may not be -inf."

    # Adjust the data if we have negative resource constraints in b_ub or b_lb.
    # Since the RHS of the tableau should be positive, we must multiply both
    # sides of these constraints by -1.  However, when doing so the inequality sign
    # flips and upper bounds become lower bounds, and vice-versa.
    pos_ub_A = Aub_in[np.where(bub_in >= 0),:][0] if not len(Aub_in) == 0 else np.empty([0,len(cc)])
    pos_ub_b = bub_in[np.where(bub_in >= 0),:][0] if not len(bub_in) == 0 else np.empty([0])
    neg_lb_A = Alb_in[np.where(blb_in < 0),:][0] if not len(Alb_in) == 0 else np.empty([0,len(cc)])
    neg_lb_b = blb_in[np.where(blb_in < 0),:][0] if not len(blb_in) == 0 else np.empty([0])
    Aub = np.vstack([pos_ub_A,-neg_lb_A])
    bub = np.concatenate([pos_ub_b,-neg_lb_b])

    pos_lb_A = Alb_in[np.where(blb_in >= 0),:][0] if not len(Alb_in) == 0 else np.empty([0,len(cc)])
    pos_lb_b = blb_in[np.where(blb_in >= 0),:][0] if not len(blb_in) == 0 else np.empty([0])
    neg_ub_A = Aub_in[np.where(bub_in < 0),:][0] if not len(Aub_in) == 0 else np.empty([0,len(cc)])
    neg_ub_b = bub_in[np.where(bub_in < 0),:][0] if not len(bub_in) == 0 else np.empty([0])
    Alb = np.vstack([pos_lb_A,-neg_ub_A])
    blb = np.concatenate([pos_lb_b,-neg_ub_b])

    # For equality constraints we can just multiply by -1
    pos_eq_A = Aeq_in[np.where(beq_in >= 0),:][0] if not len(Aeq_in) == 0 else np.empty([0,len(cc)])
    pos_eq_b = beq_in[np.where(beq_in >= 0),:][0] if not len(beq_in) == 0 else np.empty([0])
    neg_eq_A = Aeq_in[np.where(beq_in < 0),:][0] if not len(Aeq_in) == 0 else np.empty([0,len(cc)])
    neg_eq_b = beq_in[np.where(beq_in < 0),:][0] if not len(beq_in) == 0 else np.empty([0])
    Aeq = np.vstack([pos_eq_A,-neg_eq_A])
    beq = np.concatenate([pos_eq_b,-neg_eq_b])

    # The number of upper bound constraints (rows in A_ub and elements in b_ub)
    mub = len(bub)

    # The number of lower bound constraints (rows in A_lb and elements in b_lb)
    mlb = len(blb)

    # The number of equality constraints (rows in A_eq and elements in b_eq)
    meq = len(beq)

    # The number of slack variables (one for each of the upper-bound constraints)
    n_slack = mub

    # The number of surplus variables (one for each of the lower-bound constraints)
    n_surplus = mlb

    # The number of artificial variables (one for each equality constraint and lower-bound constraint)
    n_artificial = meq+mlb

    try:
        if not Aub is None:
            Aub_rows, Aub_cols = Aub.shape
        else:
            Aub_rows, Aub_cols = 0,0
    except ValueError:
        status = -1
        message = "Invalid input.  A_ub must be two-dimensional"

    try:
        if not Alb is None:
            Alb_rows, Alb_cols = Alb.shape
        else:
            Alb_rows, Alb_cols = 0,0
    except ValueError:
        status = -1
        message = "Invalid input.  A_lb must be two-dimensional"

    try:
        if not Aeq is None:
            Aeq_rows, Aeq_cols = Aeq.shape
        else:
            Aeq_rows, Aeq_cols = 0,0
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

    if Aeq_cols > 0 and Aeq_cols != n:
        status = -1
        message = "Invalid input.  Number of columns in A_eq must be equal to the size of c"

    if Alb_cols > 0 and Alb_cols != n:
        status = -1
        message = "Invalid input.  Number of columns in A_lb must be equal to the size of c"

    if Aub_cols > 0 and Aub_cols != n:
        status = -1
        message = "Invalid input.  Number of columns in A_ub must be equal to the size of c"

    # Create the tableau
    T = np.zeros([meq+mlb+mub+1,n+n_slack+n_surplus+n_artificial+1])

    # Insert objective into tableau
    if objtype.lower()[:3] == "max":
        T[-1,:n] = -cc  # maximize
        T[-1,-1] = -f0
    elif objtype.lower()[:3] == "min":
        T[-1,:n] = cc  # minimize
        T[-1,-1] = f0
    else:
        status = -1
        message = "Invalid input. Argument 'objtype' Must be one of 'max' or 'min'.  Got " + str(objtype)

    if status == 0:
        b = T[:-1,-1]

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
            b[meq+mlb:meq+mlb+mub] = bub
            # Add the slack variables to the tableau
            np.fill_diagonal(T[meq+mlb:meq+mlb+mub,n+n_surplus:n+n_surplus+n_slack], 1)

        # determine if artificial variables are needed
        # artificial variables are needed if
        # 1. the row is an equality constraint
        # 2. the row has a slack/surplus variable with a sign opposite of the RHS
        # a_rows tracks those rows of the tableau with artificial variables

        # We need artificial variables if:
        # 1. the row is an equality constraint
        # 2. the row has a slack/surplus variable with an opposite side of the RHS
        # Since the constraints have been reworked such that the RHS is entirely
        # nonnegative, the artificial variables apply to the first (meq+mlb)
        # rows.
        #
        T[:meq+mlb, n+n_slack+n_surplus:n+n_slack+n_surplus+n_artificial] = np.eye(meq+mlb)

        # If we have artificial variables, solve the phase 1 problem first
        if n_artificial > 0:
            T = np.vstack([T, np.zeros([1, T.shape[1]])])
            T[-1,n+n_slack+n_surplus:-1] = 1

            # Make the artificial variables basic feasible variables by subtracting each
            # row with an artificial variable from the Phase 1 objective
            for r in range(meq+mlb):
                T[-1,:] = T[-1,:] - T[r,:]

            x, nit1, status, message = lpsimplex(T,n,n_surplus+n_slack,n_artificial,
                                                 phase=1,callback=callback,maxiter=maxiter,tol=tol)

            # if pseudo objective is zero, remove the last row from the tableau and
            # proceed to phase 2
            if abs(T[-1,-1]) < tol:
                # Remove the pseudo-objective row from the tableau
                T = T[:-1,:]
                # Remove the artificial variable columns from the tableau
                T = np.delete(T,np.s_[n+n_slack+n_surplus:n+n_slack+n_surplus+n_artificial],1)
            else:
                status = 2
                message = "Optimization Failed.  Unable to find a feasible starting point."

        # Tableau Finished
        if status == 0:
            x, nit2, status, message = lpsimplex(T,n,n_surplus+n_slack,n_artificial,maxiter=maxiter,
                                                 phase=2,callback=callback,tol=tol,nit0=nit1)

        # For those variables with finite negative lower bounds, reverse the change of variables
        masked_L = ma.array(L,mask=np.isinf(L),fill_value=0.0)
        x = x + masked_L.data

        # For those variables with infinite negative lower bounds, take x[i] as the difference between x[i] and the floor variable.
        if have_floor_variable:
            for i in range(1,n):
                if np.isinf(L[i]):
                    x[i] -= x[0]

            x = x[1:]

        # Optimization complete at this point
        objmult = -1.0 if objtype == "min" else 1.0
        objval = T[-1,-1] * objmult

        if status in (0,1):
            if disp:
                print(message)
                print("         Current function value: {: <12.6f}".format(objval))
                print("         Iterations: {:d}".format(nit2))
        else:
            if disp:
                print(message)
                print("         Iterations: {:d}".format(nit2))

        return Result(x=x,fun=objval,nit=int(nit2),status=int(status),
                      message=message,success=(status == 0))
    else:
        # Invalid inputs provided
        print(message)
        return Result(x=np.zeros_like(cc),fun=0.0,nit=0,status=int(status),
                      message=message, success=False)
