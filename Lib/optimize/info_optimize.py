"""
Optimization Tools
==================

 A collection of general-purpose optimization routines.

   fmin        --  Nelder-Mead Simplex algorithm
                     (uses only function calls)
   fmin_powell --  Powell's (modified) level set method (uses only 
                     function calls)
   fmin_cg     --  Non-linear (Polak-Rubiere) conjugate gradient algorithm
                     (can use function and gradient).
                     
   fmin_bfgs   --  Quasi-Newton method (can use function and gradient)
   fmin_ncg    --  Line-search Newton Conjugate Gradient (can use
                     function, gradient and hessian).
   leastsq     --  Minimize the sum of squares of M equations in
                     N unknowns given a starting estimate.

  Constrained Optimizers (multivariate)

   fmin_l_bfgs_b -- Zhu, Byrd, and Nocedal's L-BFGS-B constrained optimizer
                      (if you use this please quote their papers -- see help)

   fmin_tnc      -- Truncated Newton Code originally written by Stephen Nash and
                      adapted to C by Jean-Sebastien Roy.

   fmin_cobyla   -- Contrained Optimization BY Linear Approximation

  Global Optimizers

   anneal      --  Simulated Annealing
   brute       --  Brute Force searching Optimizer

  Scalar function minimizers
 
   fminbound   --  Bounded minimization of a scalar function.
   brent       --  1-D function minimization using Brent method.
   golden      --  1-D function minimization using Golden Section method
   bracket     --  Bracket a minimum (given two starting points)
  
 Also a collection of general_purpose root-finding routines.

   fsolve      --  Non-linear multi-variable equation solver.

  Scalar function solvers

   brentq      --  quadratic interpolation Brent method
   brenth      --  Brent method (modified by Harris with
                     hyperbolic extrapolation)
   ridder      --  Ridder's method
   bisect      --  Bisection method
   newton      --  Secant method or Newton's method

   fixed_point -- Single-variable fixed-point solver.

 Utility Functions

   line_search -- Return a step that satisfies the strong Wolfe conditions.
   check_grad  -- Check the supplied derivative using finite difference
                   techniques.
"""

postpone_import = 1
