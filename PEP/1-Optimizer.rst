Scientific PEP -- Introduction of Optimizer and Function classes
================================================================

.. outline::

   * Abstract
   * Introduction
       * Here's what minimization does...
           * It minimizes a function
           * These are -or should be- fairly independent -- functions and optimizers are not tied together.
       * Point to users of...
           * Minimization in general
           * scipy.optimize.minimize (many users, do a github search)
   * Motivation
       * No standard interface for optimizers or functions.
           * Have to explain why minimize isn't a standard interface.
               * Currently there is a hotch potch of warn_flag numbers that indicate problems when a minimizer stops. Using
                 an Optimizer class could standardise these. See #7819 for discussion on this. The Optimizer class could
                 return an 
           * Minimization is a common problem and implemented in many places
           * Providing a standard interface for this could help unify users and libraries
       * Current API needs improvement
          * minimize is trying to be a class
                 * method: should be subclasses
                 * show_options: show method-specific args
                 * some options specific to method (jac, hess, hessp, contraints, options, bounds)
                 * OptimizeResult: trying to expose what should be properties of class
                 * callback: not adequate (only sends one arg, not any internal state)
                      * only sends `x`, not the potentially expensive `f(x), g(x), h(x)`.
                          **the opposing argument here is that we could just add extra solver state information to the**
                          **callback. ironically the easiest way to achieve this by using Optimizer objects, where**
                          **once you've implemented a change to the base class all Optimizers access the benefits.**
                      * What if some internal state is wanted?
          * function arg is trying to be a class
              * jac, hess, hessp
              * args (kwargs?)
          * there is no separation of concerns between function and minimizer
              * meaning the minimizer is carrying out numerical gradient calculations.
              * The correct place for grad computation belongs with the function, not the minimizer.
              * Mixing of function arguments with optimization arguments (plus, there are too many arguments)
              * no kwargs for func, only args
          * scipy.optimize.minimize is a black box (have to explain why)
              * hides all details. Some are literal black boxes and implemented in Fortran/C.
                  * e.g., what if want to change step size? Choosing an initial step size is difficult. There's theoritical
                    bounds, but these are not known in practice.
                  * if the user doesn't provide a gradient function the minimizers currently use the same absolute step size
                      for numerical differentiation for the duration of the minimization. However, the fd-step size should
                      be relative to parameter value as it changes. Not easy to fix this in current implementation without 
                      placing the onus on the user to write their own grad function, this is the job of the library.
                      The new Function object will offer more options for numerical differentiation (absolute step, relative
                      step, 2-point/3-point/complex step, bounds). Of course, the user can still provide their own gradient
                      implementation if preferred.
                  * would like ability to proceed stepwise through iteration
                      * What if running some web server, and don't have time to wait for minimization to finish?
                      * There's no easy way of halting minimization and still returning a solution. With the Optimizer
                        approach one can simply stop on the current iteration, if you're doing the stepping, and you
                        retain access to the current best solution. You can then restart at a later point. Moreover
                        if you are using the Optimizer.solve method that runs to convergence you can simply halt at anytime
                        by raising a StopIteration exception, either in the 'callback', or in your Function evaluation.
                        This could be done for current Optimizers, but only by amending all minimizers.
                      * user can use their own convergence criteria, don't need to depend on minimizer to halt.
                  * would like to access solver state
                      * e.g., current value of f(x)
                      * e.g., for coding gradients
                  * can't access solver state or hyper parameters, and change on fly
                     * e.g. gradient coding as example
                     * e.g. change convergence tolerances as we're going
                     * e.g. change mutation constant during differential evolution.
          * addition of new features to minimizers leads to lengthy functions and lots of duplicate code.
              * Classes => inheritance. Base class improves => all improve. For example, placing numerical differentiation in 
              the Function class allows either absolute or relative delta change to be made easily, and in one place. To do 
              that for all minimizers would require modifications and extra keywords to all minimizer functions with the 
              attendant risk of introducing bugs in lots of places. Testing those changes is a lot harder.
              * With Optimizer objects testing can be made a lot easier. If the base class is tested thoroughly then 
              subclasses with inherited methods are by definition covered. This is not the case for a multiplicity of 
              minimizer functions.
              * Unix philisophy, small sharp tools for one job and one job only. Not many dull tools for the same job.
   * The following open issues/PRs would be significantly easier to be addressed (or tackled by the user themselves) with
     subclassing of an Optimizer base class. That there are many signifies the level of difficulty implementing a coherent 
     solution across scipy.optimize.
      # 5832 grad.T should be returned but not documented
      # 7819 WIP: Basin hopping improvements. **discusses behaviour of how a minimizer should signify success/failure, e.g.**
        **if a constraint is violated**
      # 7425 ENH: optimize: more complete callback signature. **easily achieved, Optimizer base class calls the callback**
      # 6907 differential_evolution: improve callback **easily achieved, Optimizer base class calls the callback**
      # 4384 ENH: optimize, returning True from callback function halts minimization **callback could return a StopIteration**
        **which would simply stop at the current iteration in Optimizer.solve(), the optimization could then be restarted if**
        **if desired**.
      # 8375 optimize - check that maxiter is not exceeded **correct implementation is inherited by all Optimizers.**
        **testing is simple for all Optimizers**
      # 8419 (comment): "some optimize.minimize methods modify the parameter vector in-place", **is inherited by all** 
        **Optimizers**
      # 8031 Scipy optimize.minimize maxfun has confusing behavior **maxfun behaviour is implemented by Optimizer base**
        **class. Documentation in one place should make things clear**
      # 8373 "scipy.optimize has broken my trust." mismatch between callback x and displayed output from L-BFGS-B
      # 6019 "minimize_scalar doesn't honor disp option". **Optimizer base class can standardise iteration by iteration**
        **displaying, and end of solve displaying. Inheriting Optimizers can override if absolutely necessary**
      # 7854: "BUG: L-BFGS-B does two more iterations than specified in maxiter" **More easily tested with Optimizer class**
      # 6673, "return value of scipy.optimize.minimize not consistent for 1D", **This can be standardised more easily**
      # 7306 "any way of stopping optimization?". **Easily implemented by Optimizer. Either by raising StopIteration,**
        **or by controlling the iteration yourself on a stepwise basis** One comment in this issue: "Beyond a pre-specified
        iteration limit, I always wanted some way of gracefully terminating an optimization routine during execution. I was
        working on problems that took a very long time to solve and sometimes I wanted to see what was going on when the
        algorithm seemed close to a solution but never seemed to achieve the termination conditions.
      # 6878 differential_evolution: make callback receive fun(xk) **User has full access to Optimizer, this is available**
        **during stepwise iteration. Otherwise it should be straightforward to introduce an expanded callback**
        **in a standardised fashion**
      # 6026 Replace approx_grad with _numdiff.approx_derivative in scipy.optimize **all numerical differentiation done in**
        **Function class, fix is only needed in one place. Optimizers don't need to know.**.
      # 6019 minimize_scalar doesn't seem to honor "disp" option
      # 5481 "1D root-finding interface and documentation could be improved" **Asking for a standardised approach to root**
        **finding. May be possible to inherit Optimizer class for root finding to standardise behaviour.**
      # 5161 Optimizers reporting success when the minimum is NaN. **this would be standardised to make success False**
      # 4921 scipy.optimize maxiter option not working as expected **Optimizer.solve standardises for all subclasses**
      # 3816 wrap_function seems not to be working when wrapper_args is a one element list **fix in Optimizer, fix in all**
        *subclasses**
      
 
   * Existing work
       * Class defs: PyTorch, skopt
       * Functional class wrapper around minimize: statsmodels, astropy, scikits.fitting
       * Functional defs: sklearn, daskml, skimage
   * Proposed solution
       * Classes (idea: `Function` and `Optimizer` class)
           * `Optimizer` - takes care of minimization and stepping
           * `Function` - takes care of evaluating function, gradient, and hessian. Takes care of numerical differentiation
             for grad and hess if required. Can be overridden if the user wishes to define their own grad/hess
             implementations. This pattern is intrinsic, and is sort of **already in use** in scipy at
             scipy/benchmarks/benchmarks/test_functions.py.
       * Goal:
           * provide minimal class interface
           * preserve backwards compatibility
           * targetted at minimization of scalar functions to start with, although the Optimizer class and its methods should
             be a suitable base class for implementing for class based root and least-squares solvers. For example, both of
             those examples need to iterate, they both finish up with an OptimizeResult, they both have convergence criteria,
             etc.
       * Give an example
   * Enhancements
       * Provide standard interface
           * for enhancements to sklearn, dask-ml, etc. Possibly PyTorch. **Would those projects be prepared to state that?**
           * it would provide a standard way to operate the object, but all the classes would still have different names
           * give example of how sklearn could revamp (ask the developers how they'd use it)
       * Provide class features
           * object interaction. Useful for experts, intermediates.
           * expose alg hyperparameters (grid search, etc)
           * keyboard interrupts
           * introduction of context manager enables easy setup of cleanup actions
              * would make it easier have wholesale introduction of things like multiprocessing.
              * We should think about multiprocessing or multithreaded algorithms like Hogwild!. How will these be used?
      * Clean up minimize API (it's complicated right now)
         * Require fewer arguments to minimize, and separate them
   * Implementation
       * List functions, attributes in more depth
       * Existing code
           * How would it work with C/Fortran optimizers?
           * What interface are we proposing? See proposed code below
       * Speed
         * will be benchmarked to check that performance is not damaged. Class based system is easy to convert to cython.
           **Using asv it's about a 25% extra time penalty for bfgs, lbfgsb, fmin (e.g. 252us to 310us). However,**
           **those benchmarks use really quick functions. If one of the benchmarks was on much slower function**
           **the overhead will be relatively minor compared to that going to an Optimizer class**
       * Backwards compatibility
         * backwards compatibility is a focus
         * the functionality will remain but rely on the solver objects. Should be able to remove `_minimize_lbfgsb`, etc.
         * new solver objects can be used by themselves.

*Abstract*

Introduction
============

Optimization is extremely common and often critical in many
applications. Imaging, machine learning and regression problems 
all depend on optimization. As such, it has been adopted by
libraries including SciPy and many related libraries (e.g.,
scikit-learn). Optimization has received significant attention from
industry as well -- Google, Facebook, Amazon and Microsoft have
developed Tensorflow, PyTorch, MXNet and CNTK respectively, all of which
use optimization, have Python bindings and are open source.

Optimization is the minimization or maximization (though typically
minimization) of a certain function. Minimization tries to find which
argument yields the smallest function value, or in pseudo-code,

.. code:: python

    import numpy as np
    from scipy.optimize import minimize
    
    def f(x):
        return (x - 1) ** 2

    result = minimize(f, x0=np.random.randn())
    assert np.allclose(result.x, 1) and np.allclose(result.fun, 0)

The SciPy ``minimize`` function has been widely used. Over 17,000
results for "``from scipy.optimize import minimize``" appear from a
GitHub search, and ``minimize`` is included in many libraries including
scikit-learn, scikit-image, statsmodels and astropy. Preserving
backwards compatibility to keep this code functional is a priority.
However, we believe that we can improve upon SciPy's minimization API.
We believe implementation of this will allow easier use, enable more
widespread use and unify various interfaces.

Motivation
==========
No standard interface
---------------------

Current API needs improvement
-----------------------------

`minimize` has many class features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`minimize`'s `func` argument has many class features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`minimize` is a black box
^^^^^^^^^^^^^^^^^^^^^^^^^
Separation of function and minimizer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Repeated code
^^^^^^^^^^^^^
Open bugs
^^^^^^^^^

Existing work
=============

Proposed solution
=================
Goals
-----
Embodiment
----------
Example
-------

Enhancements
============
Standard interface
------------------
Class features
--------------
API cleaning
------------



Implementation
==============
Definition
----------
Existing code
-------------
Backward compatibility
----------------------

Proposed code
-------------


.. code-block:: python

    def func(x, *args):
        return x**2 + args[0]
    def grad(x, *args):
        return 2 * x

    def callback(x): print(x)

    # existing call has lots of parameters, mixing optimizer args with func args
    # it might be nice to have **kwds as well, but not possible with current approach
    result = minimize(func, x0, args=(2,), jac=grad, method='BFGS', maxiter=10, callback=callback)

    # proposed

    function = Function(func=func, args=(2,), kwargs=kwargs, grad=grad)
    opt = BFGS(function, x0)
    result = opt.solve(maxiter=10, callback=callback)

    # could also have
    result = BFGS(function, x0).solve(maxiter=10, callback=callback)

    # alternatively control how iteration occurs
    d = opt.hyper_parameters
    for i, v in enumerate(opt):
      x, f = v print(i, f, x)
      d['my_hyper_parameter'] = np.inf

    # use function classes encapsulates the whole function and offers the potential for more sophisticated calculation.

    class Quad(Function):
        def __init__(self, bkg):
            super(Quad, self).__init__(self)
            self.bkg = bkg

        def func(self, x):
            return (x**2 + args[0])

        def grad(self, x):
            return 2*x

        def hess(self, x):
            return 2

    opt = BFGS(function, x0).solve(maxiter=10)

    # context managers offer the chance for cleanup actions, for example multiprocessing.

    with DifferentialEvolutionSolver(function, bounds, workers=2) as opt:
        # the __entry__ and __exit__ in the solver can create and close
        # multiprocessing pools.
        res = opt.solve()
