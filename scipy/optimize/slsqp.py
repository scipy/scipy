from _slsqp import slsqp
from numpy import zeros, array, identity, linalg, rank, squeeze, append, \
                  asfarray,product
from sys import exit
from math import sqrt
from optimize import approx_fprime, wrap_function
import numpy

_epsilon = sqrt(numpy.finfo(float).eps)

def fmin_slsqp( func, x0 , eqcons=[], f_eqcons=None, ieqcons=[], f_ieqcons=None,
                bounds = [], fprime = None, fprime_cons=None,args = (),
                iter = 100, acc = 1.0E-6, iprint = 1, full_output = 0, 
                epsilon = _epsilon ):
    """
    Minimize a function using Sequential Least SQuares Programming
    
    Description:
        Python interface function for the SLSQP Optimization subroutine 
        originally implemented by Dieter Kraft.
    
    Inputs:
        func           - Objective function (in the form func(x, *args))
        x0             - Initial guess for the independent variable(s). 
        eqcons         - A list of functions of length n such that 
                         eqcons[j](x0,*args) == 0.0 in a successfully optimized 
                         problem.
        f_eqcons       - A function of the form f_eqcons(x, *args) that returns an
                         array in which each element must equal 0.0 in a 
                         successfully optimized problem.  If f_eqcons is
                         specified, eqcons is ignored.                    
        ieqcons        - A list of functions of length n such that 
                         ieqcons[j](x0,*args) >= 0.0 in a successfully optimized 
                         problem.
        f_ieqcons      - A function of the form f_ieqcons(x0, *args) that returns 
                         an array in which each element must be greater or equal 
                         to 0.0 in a successfully optimized problem.  If 
                         f_ieqcons is specified, ieqcons is ignored.                                               
        bounds         - A list of tuples specifying the lower and upper bound 
                         for each independent variable [(xl0, xu0),(xl1, xu1),...]
        fprime         - A function that evaluates the partial derivatives of func 
        fprime_cons    - A function of the form f(x, *args) that returns the
                         m by n array of constraint normals.  If not provided,
                         the normals will be approximated. Equality constraint
                         normals precede inequality constraint normals. The
                         array returned by fprime_cons should be sized as 
                         ( len(eqcons) + len(ieqcons), len(x0) ).  If 
                         fprime_cons is not supplied (normals are approximated)
                         then the constraints must be supplied via the eqcons
                         and ieqcons structures, not f_eqcons and f_ieqcons.
        args           - A sequence of additional arguments passed to func and fprime
        iter           - The maximum number of iterations (int)
        acc            - Requested accuracy (float)
        iprint         - The verbosity of fmin_slsqp.
                         iprint <= 0 : Silent operation
                         iprint == 1 : Print summary upon completion (default)
                         iprint >= 2 : Print status of each iterate and summary
        full_output    - If 0, return only the minimizer of func (default). 
                         Otherwise, output final objective function and summary 
                         information.
        epsilon        - The step size for finite-difference derivative estimates.
                     
    Outputs: ( x, { fx, gnorm, its, imode, smode })
        x            - The final minimizer of func.
        fx           - The final value of the objective function.
        its          - The number of iterations.
        imode        - The exit mode from the optimizer, as an integer.
        smode        - A string describing the exit mode from the optimizer. 
        
    Exit modes are defined as follows:
        -1 : Gradient evaluation required (g & a)
         0 : Optimization terminated successfully.
         1 : Function evaluation required (f & c)
         2 : More equality constraints than independent variables
         3 : More than 3*n iterations in LSQ subproblem
         4 : Inequality constraints incompatible
         5 : Singular matrix E in LSQ subproblem
         6 : Singular matrix C in LSQ subproblem
         7 : Rank-deficient equality constraint subproblem HFTI
         8 : Positive directional derivative for linesearch
         9 : Iteration limit exceeded

    """

    exit_modes = { -1 : "Gradient evaluation required (g & a)",
                    0 : "Optimization terminated successfully.", 
                    1 : "Function evaluation required (f & c)",
                    2 : "More equality constraints than independent variables",
                    3 : "More than 3*n iterations in LSQ subproblem",
                    4 : "Inequality constraints incompatible",
                    5 : "Singular matrix E in LSQ subproblem",
                    6 : "Singular matrix C in LSQ subproblem",
                    7 : "Rank-deficient equality constraint subproblem HFTI",
                    8 : "Positive directional derivative for linesearch",
                    9 : "Iteration limit exceeded" }
    
    # Wrap the functions
    # Wrap func
    feval, func = wrap_function(func, args)        
    if fprime:
        # Wrap fprime, if provided
        geval, fprime = wrap_function(fprime,args) 
    else:
        # Wrap approx_fprime, if fprime not provided
        geval, fprime = wrap_function(approx_fprime,(func,epsilon)) 
    if fprime_cons:
        approx_constraint_norms = False
        if f_eqcons:
            ceval, f_eqcons = wrap_function(f_eqcons,args)
        else:
            for i in range(len(eqcons)):
                if eqcons[i]:
                    ceval, eqcons[i] = wrap_function(eqcons[i],args)
        if f_ieqcons:
            ceval, f_ieqcons = wrap_function(f_ieqcons,args)
        else:
            for i in range(len(ieqcons)):
                if ieqcons[i]:
                    ceval, ieqcons[i] = wrap_function(ieqcons[i],args)
        geval, fprime_cons = wrap_function(fprime_cons,args)
    else:   
        approx_constraint_norms = True
        eqcons_prime = []
        for i in range(len(eqcons)):
            eqcons_prime.append(None)
            if eqcons[i]:
                ceval, eqcons[i] = wrap_function(eqcons[i],args)
                geval, eqcons_prime[i] = \
                wrap_function(approx_fprime, (eqcons[i],epsilon))
        ieqcons_prime = []
        for i in range(len(ieqcons)):
            ieqcons_prime.append(None)
            if ieqcons[i]:
                ceval, ieqcons[i] = wrap_function(ieqcons[i],args)
                geval, ieqcons_prime[i] = \
                wrap_function(approx_fprime, (ieqcons[i],epsilon))
    
    # Transform x0 into an array.  
    x = asfarray(x0).flatten()

    # Set the parameters that SLSQP will need
    meq = len(eqcons)        # meq = The number of equality constraints
    m = meq + len(ieqcons)   # m   = The total number of constraints
    la = array([1,m]).max()  # la  = 
    n = len(x)               # n   = The number of independent variables
    
    # Define the workspaces for SLSQP
    n1 = n+1
    mineq = m - meq + n1 + n1
    len_w = (3*n1+m)*(n1+1)+(n1-meq+1)*(mineq+2) + 2*mineq+(n1+mineq)*(n1-meq) \
            + 2*meq + n1 +(n+1)*n/2 + 2*m + 3*n + 3*n1 + 1
    len_jw = mineq
    w = zeros(len_w)
    jw = zeros(len_jw)     
    
    # Decompose bounds into xl and xu
    if len(bounds) == 0:
        bounds = [(-1.0E12, 1.0E12) for i in range(n)]
    if len(bounds) != n:
        raise IndexError, 'SLSQP Error:  If bounds is specified, len(bounds) == len(x0)'
    xl = array( [ b[0] for b in bounds ] )
    xu = array( [ b[1] for b in bounds ] )
    
    # Initialize the iteration counter and the mode value        
    mode = array(0,int)
    acc = array(acc,float)    
    majiter = array(iter,int)
    majiter_prev = 0

    # Print the header if iprint >= 2
    if iprint >= 2:
        print "%5s %5s %16s %16s" % ("NIT","FC","OBJFUN","GNORM")

    while 1:    
        if mode == 0 or mode == 1: # objective and constraint evaluation requird    
            # Compute objective function
            fx = func(x)            
            # Compute the constraints
            if f_eqcons:
               ceq = f_eqcons(x)
            else:
               ceq = [ eqcons[i](x) for i in range(meq) ]
            if f_ieqcons:
               cieq = f_ieqcons(x)
            else:
               cieq = [ ieqcons[i](x) for i in range(len(ieqcons)) ] 
            c = numpy.concatenate( (ceq,cieq), 1)
            #c = array ( [ eqcons[i](x) for i in range(meq) ] + 
            #            [ ieqcons[i](x) for i in range(len(ieqcons)) ] )                    
        if mode == 0 or mode == -1: # gradient evaluation required
            # Compute the derivatives of the objective function
            # For some reason SLSQP wants g dimensioned to n+1
            g = append(fprime(x),0.0)         
        
            # Compute the normals of the constraints
            if approx_constraint_norms:
                a = zeros([la,n1])
                for i in range(meq):
                    a[i] = append(eqcons_prime[i](x),0.0)
                for i in range(meq+1,m):
                    a[i] = append(ieqcons_prime[i](x),0.0)
            else:
                a = numpy.concatenate( (fprime_cons(x),zeros([la,1])),1)                
        
        # Call SLSQP
        slsqp(m, meq, x, xl, xu, fx, c, g, a, acc, majiter, mode, w, jw)
        
        # Print the status of the current iterate if iprint > 2 and the 
        # major iteration has incremented
        if iprint >= 2 and majiter > majiter_prev:
            print "%5i %5i % 16.6E % 16.6E" % (majiter,feval[0],
                                               fx,linalg.norm(g))    
        
        # If exit mode is not -1 or 1, slsqp has completed
        if abs(mode) != 1:
            break
            
        majiter_prev = int(majiter)
        
    # Optimization loop complete.  Print status if requested
    if iprint >= 1:
        print exit_modes[int(mode)] + "    (Exit mode " + str(mode) + ')'  
        print "            Current function value:", fx 
        print "            Iterations:", majiter
        print "            Function evaluations:", feval[0]
        print "            Gradient evaluations:", geval[0]
  
    if full_output == 0:
        return x
    else: 
        return [list(x), 
                float(fx), 
                int(majiter), 
                int(mode), 
                exit_modes[int(mode)] ]
        
