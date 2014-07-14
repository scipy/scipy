# Original Author: Travis Oliphant 2002
# Bug-fixes in 2006 by Tim Leslie

from __future__ import division, print_function, absolute_import

import numpy
from numpy import (asarray, tan, exp, ones, squeeze, sign,
        all, log, sqrt, pi, shape, array, minimum, where, random, deprecate)
from .optimize import OptimizeResult, _check_unknown_options
from scipy.lib.six import xrange

__all__ = ['anneal']

_double_min = numpy.finfo(float).min
_double_max = numpy.finfo(float).max


class base_schedule(object):
    def __init__(self):
        self.dwell = 20
        self.learn_rate = 0.5
        self.lower = -10
        self.upper = 10
        self.Ninit = 50
        self.accepted = 0
        self.tests = 0
        self.feval = 0
        self.k = 0
        self.T = None

    def init(self, **options):
        self.__dict__.update(options)
        self.lower = asarray(self.lower)
        self.lower = where(self.lower == numpy.NINF, -_double_max, self.lower)
        self.upper = asarray(self.upper)
        self.upper = where(self.upper == numpy.PINF, _double_max, self.upper)
        self.k = 0
        self.accepted = 0
        self.feval = 0
        self.tests = 0

    def getstart_temp(self, best_state):
        """ Find a matching starting temperature and starting parameters vector
        i.e. find x0 such that func(x0) = T0.

        Parameters
        ----------
        best_state : _state
            A _state object to store the function value and x0 found.

        Returns
        -------
        x0 : array
            The starting parameters vector.
        """

        assert(self.dims is not None)
        lrange = self.lower
        urange = self.upper
        fmax = _double_min
        fmin = _double_max
        for _ in range(self.Ninit):
            x0 = random.uniform(size=self.dims)*(urange-lrange) + lrange
            fval = self.func(x0, *self.args)
            self.feval += 1
            if fval > fmax:
                fmax = fval
            if fval < fmin:
                fmin = fval
                best_state.cost = fval
                best_state.x = array(x0)

        self.T0 = (fmax-fmin)*1.5
        return best_state.x

    def accept_test(self, dE):
        T = self.T
        self.tests += 1
        if dE < 0:
            self.accepted += 1
            return 1
        p = exp(-dE*1.0/self.boltzmann/T)
        if (p > random.uniform(0.0, 1.0)):
            self.accepted += 1
            return 1
        return 0

    def update_guess(self, x0):
        pass

    def update_temp(self, x0):
        pass


#  A schedule due to Lester Ingber
class fast_sa(base_schedule):
    def init(self, **options):
        self.__dict__.update(options)
        if self.m is None:
            self.m = 1.0
        if self.n is None:
            self.n = 1.0
        self.c = self.m * exp(-self.n * self.quench)

    def update_guess(self, x0):
        x0 = asarray(x0)
        u = squeeze(random.uniform(0.0, 1.0, size=self.dims))
        T = self.T
        y = sign(u-0.5)*T*((1+1.0/T)**abs(2*u-1)-1.0)
        xc = y*(self.upper - self.lower)
        xnew = x0 + xc
        return xnew

    def update_temp(self):
        self.T = self.T0*exp(-self.c * self.k**(self.quench))
        self.k += 1
        return


class cauchy_sa(base_schedule):
    def update_guess(self, x0):
        x0 = asarray(x0)
        numbers = squeeze(random.uniform(-pi/2, pi/2, size=self.dims))
        xc = self.learn_rate * self.T * tan(numbers)
        xnew = x0 + xc
        return xnew

    def update_temp(self):
        self.T = self.T0/(1+self.k)
        self.k += 1
        return


class boltzmann_sa(base_schedule):
    def update_guess(self, x0):
        std = minimum(sqrt(self.T) * ones(self.dims),
                      (self.upper - self.lower) / 3.0 / self.learn_rate)
        x0 = asarray(x0)
        xc = squeeze(random.normal(0, 1.0, size=self.dims))

        xnew = x0 + xc*std*self.learn_rate
        return xnew

    def update_temp(self):
        self.k += 1
        self.T = self.T0 / log(self.k+1.0)
        return


class _state(object):
    def __init__(self):
        self.x = None
        self.cost = None


# TODO:
#     allow for general annealing temperature profile
#     in that case use update given by alpha and omega and
#     variation of all previous updates and temperature?

# Simulated annealing

@deprecate(message='Deprecated in scipy 0.14.0, use basinhopping instead')
def anneal(func, x0, args=(), schedule='fast', full_output=0,
           T0=None, Tf=1e-12, maxeval=None, maxaccept=None, maxiter=400,
           boltzmann=1.0, learn_rate=0.5, feps=1e-6, quench=1.0, m=1.0, n=1.0,
           lower=-100, upper=100, dwell=50, disp=True):
    """
    Minimize a function using simulated annealing.

    Uses simulated annealing, a random algorithm that uses no derivative
    information from the function being optimized. Other names for this
    family of approaches include: "Monte Carlo", "Metropolis",
    "Metropolis-Hastings", `etc`. They all involve (a) evaluating the
    objective function on a random set of points, (b) keeping those that
    pass their randomized evaluation critera, (c) cooling (`i.e.`,
    tightening) the evaluation critera, and (d) repeating until their
    termination critera are met.  In practice they have been used mainly in
    discrete rather than in continuous optimization.

    Available annealing schedules are 'fast', 'cauchy' and 'boltzmann'.

    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        `f(x, *args)`, where `x` is the argument in the form of a 1-D array
        and `args` is a  tuple of any additional fixed parameters needed to
        completely specify the function.
    x0: 1-D array
        An initial guess at the optimizing argument of `func`.
    args : tuple, optional
        Any additional fixed parameters needed to completely
        specify the objective function.
    schedule : str, optional
        The annealing schedule to use.  Must be one of 'fast', 'cauchy' or
        'boltzmann'.  See `Notes`.
    full_output : bool, optional
        If `full_output`, then return all values listed in the Returns
        section. Otherwise, return just the `xmin` and `status` values.
    T0 : float, optional
        The initial "temperature".  If None, then estimate it as 1.2 times
        the largest cost-function deviation over random points in the
        box-shaped region specified by the `lower, upper` input parameters.
    Tf : float, optional
        Final goal temperature.  Cease iterations if the temperature
        falls below `Tf`.
    maxeval : int, optional
        Cease iterations if the number of function evaluations exceeds
        `maxeval`.
    maxaccept : int, optional
        Cease iterations if the number of points accepted exceeds `maxaccept`.
        See `Notes` for the probabilistic acceptance criteria used.
    maxiter : int, optional
        Cease iterations if the number of cooling iterations exceeds `maxiter`.
    learn_rate : float, optional
        Scale constant for tuning the probabilistc acceptance criteria.
    boltzmann : float, optional
        Boltzmann constant in the probabilistic acceptance criteria
        (increase for less stringent criteria at each temperature).
    feps : float, optional
        Cease iterations if the relative errors in the function value over the
        last four coolings is below `feps`.
    quench, m, n : floats, optional
        Parameters to alter the `fast` simulated annealing schedule.
        See `Notes`.
    lower, upper : floats or 1-D arrays, optional
        Lower and upper bounds on the argument `x`.  If floats are provided,
        they apply to all components of `x`.
    dwell : int, optional
        The number of times to execute the inner loop at each value of the
        temperature.  See `Notes`.
    disp : bool, optional
        Print a descriptive convergence message if True.

    Returns
    -------
    xmin : ndarray
        The point where the lowest function value was found.
    Jmin : float
        The objective function value at `xmin`.
    T : float
        The temperature at termination of the iterations.
    feval : int
        Number of function evaluations used.
    iters : int
        Number of cooling iterations used.
    accept : int
        Number of tests accepted.
    status : int
        A code indicating the reason for termination:

        - 0 : Points no longer changing.
        - 1 : Cooled to final temperature.
        - 2 : Maximum function evaluations reached.
        - 3 : Maximum cooling iterations reached.
        - 4 : Maximum accepted query locations reached.
        - 5 : Final point not the minimum amongst encountered points.

    See Also
    --------
    basinhopping : another (more performant) global optimizer
    brute : brute-force global optimizer

    Notes
    -----
    Simulated annealing is a random algorithm which uses no derivative
    information from the function being optimized. In practice it has
    been more useful in discrete optimization than continuous
    optimization, as there are usually better algorithms for continuous
    optimization problems.

    Some experimentation by trying the different temperature
    schedules and altering their parameters is likely required to
    obtain good performance.

    The randomness in the algorithm comes from random sampling in numpy.
    To obtain the same results you can call `numpy.random.seed` with the
    same seed immediately before calling `anneal`.

    We give a brief description of how the three temperature schedules
    generate new points and vary their temperature.  Temperatures are
    only updated with iterations in the outer loop.  The inner loop is
    over loop over ``xrange(dwell)``, and new points are generated for
    every iteration in the inner loop.  Whether the proposed new points
    are accepted is probabilistic.

    For readability, let ``d`` denote the dimension of the inputs to func.
    Also, let ``x_old`` denote the previous state, and ``k`` denote the
    iteration number of the outer loop.  All other variables not
    defined below are input variables to `anneal` itself.

    In the 'fast' schedule the updates are::

        u ~ Uniform(0, 1, size = d)
        y = sgn(u - 0.5) * T * ((1 + 1/T)**abs(2*u - 1) - 1.0)

        xc = y * (upper - lower)
        x_new = x_old + xc

        c = n * exp(-n * quench)
        T_new = T0 * exp(-c * k**quench)

    In the 'cauchy' schedule the updates are::

        u ~ Uniform(-pi/2, pi/2, size=d)
        xc = learn_rate * T * tan(u)
        x_new = x_old + xc

        T_new = T0 / (1 + k)

    In the 'boltzmann' schedule the updates are::

        std = minimum(sqrt(T) * ones(d), (upper - lower) / (3*learn_rate))
        y ~ Normal(0, std, size = d)
        x_new = x_old + learn_rate * y

        T_new = T0 / log(1 + k)

    References
    ----------
    [1] P. J. M. van Laarhoven and E. H. L. Aarts, "Simulated Annealing: Theory
        and Applications", Kluwer Academic Publishers, 1987.

    [2] W.H. Press et al., "Numerical Recipies: The Art of Scientific Computing",
        Cambridge U. Press, 1987.

    Examples
    --------
    *Example 1.* We illustrate the use of `anneal` to seek the global minimum
    of a function of two variables that is equal to the sum of a positive-
    definite quadratic and two deep "Gaussian-shaped" craters.  Specifically,
    define the objective function `f` as the sum of three other functions,
    ``f = f1 + f2 + f3``.  We suppose each of these has a signature
    ``(z, *params)``, where ``z = (x, y)``, ``params``, and the functions are
    as defined below.

    >>> params = (2, 3, 7, 8, 9, 10, 44, -1, 2, 26, 1, -2, 0.5)
    >>> def f1(z, *params):
    ...     x, y = z
    ...     a, b, c, d, e, f, g, h, i, j, k, l, scale = params
    ...     return (a * x**2 + b * x * y + c * y**2 + d*x + e*y + f)

    >>> def f2(z, *params):
    ...     x, y = z
    ...     a, b, c, d, e, f, g, h, i, j, k, l, scale = params
    ...     return (-g*np.exp(-((x-h)**2 + (y-i)**2) / scale))

    >>> def f3(z, *params):
    ...     x, y = z
    ...     a, b, c, d, e, f, g, h, i, j, k, l, scale = params
    ...     return (-j*np.exp(-((x-k)**2 + (y-l)**2) / scale))

    >>> def f(z, *params):
    ...     x, y = z
    ...     a, b, c, d, e, f, g, h, i, j, k, l, scale = params
    ...     return f1(z, *params) + f2(z, *params) + f3(z, *params)

    >>> x0 = np.array([2., 2.])     # Initial guess.
    >>> from scipy import optimize
    >>> np.random.seed(555)   # Seeded to allow replication.
    >>> res = optimize.anneal(f, x0, args=params, schedule='boltzmann',
                              full_output=True, maxiter=500, lower=-10,
                              upper=10, dwell=250, disp=True)
    Warning: Maximum number of iterations exceeded.
    >>> res[0]  # obtained minimum
    array([-1.03914194,  1.81330654])
    >>> res[1]  # function value at minimum
    -3.3817...

    So this run settled on the point [-1.039, 1.813] with a minimum function
    value of about -3.382.  The final temperature was about 212. The run used
    125301 function evaluations, 501 iterations (including the initial guess as
    a iteration), and accepted 61162 points. The status flag of 3 also
    indicates that `maxiter` was reached.

    This problem's true global minimum lies near the point [-1.057, 1.808]
    and has a value of about -3.409.  So these `anneal` results are pretty
    good and could be used as the starting guess in a local optimizer to
    seek a more exact local minimum.

    *Example 2.* To minimize the same objective function using
    the `minimize` approach, we need to (a) convert the options to an
    "options dictionary" using the keys prescribed for this method,
    (b) call the `minimize` function with the name of the method (which
    in this case is 'Anneal'), and (c) take account of the fact that
    the returned value will be a `OptimizeResult` object (`i.e.`, a dictionary,
    as defined in `optimize.py`).

    All of the allowable options for 'Anneal' when using the `minimize`
    approach are listed in the ``myopts`` dictionary given below, although
    in practice only the non-default values would be needed.  Some of their
    names differ from those used in the `anneal` approach.  We can proceed
    as follows:

    >>> myopts = {
            'schedule'     : 'boltzmann',   # Non-default value.
            'maxfev'       : None,  # Default, formerly `maxeval`.
            'maxiter'      : 500,   # Non-default value.
            'maxaccept'    : None,  # Default value.
            'ftol'         : 1e-6,  # Default, formerly `feps`.
            'T0'           : None,  # Default value.
            'Tf'           : 1e-12, # Default value.
            'boltzmann'    : 1.0,   # Default value.
            'learn_rate'   : 0.5,   # Default value.
            'quench'       : 1.0,   # Default value.
            'm'            : 1.0,   # Default value.
            'n'            : 1.0,   # Default value.
            'lower'        : -10,   # Non-default value.
            'upper'        : +10,   # Non-default value.
            'dwell'        : 250,   # Non-default value.
            'disp'         : True   # Default value.
            }
    >>> from scipy import optimize
    >>> np.random.seed(777)  # Seeded to allow replication.
    >>> res2 = optimize.minimize(f, x0, args=params, method='Anneal',
                                 options=myopts)
    Warning: Maximum number of iterations exceeded.
    >>> res2
      status: 3
     success: False
      accept: 61742
        nfev: 125301
           T: 214.20624873839623
         fun: -3.4084065576676053
           x: array([-1.05757366,  1.8071427 ])
     message: 'Maximum cooling iterations reached'
     nit: 501

    """

    opts = {'schedule': schedule,
            'T0': T0,
            'Tf': Tf,
            'maxfev': maxeval,
            'maxaccept': maxaccept,
            'maxiter': maxiter,
            'boltzmann': boltzmann,
            'learn_rate': learn_rate,
            'ftol': feps,
            'quench': quench,
            'm': m,
            'n': n,
            'lower': lower,
            'upper': upper,
            'dwell': dwell,
            'disp': disp}

    res = _minimize_anneal(func, x0, args, **opts)

    if full_output:
        return res['x'], res['fun'], res['T'], res['nfev'], res['nit'], \
            res['accept'], res['status']
    else:
        return res['x'], res['status']


def _minimize_anneal(func, x0, args=(),
                     schedule='fast', T0=None, Tf=1e-12, maxfev=None,
                     maxaccept=None, maxiter=400, boltzmann=1.0,
                     learn_rate=0.5, ftol=1e-6, quench=1.0, m=1.0, n=1.0,
                     lower=-100, upper=100, dwell=50, disp=False,
                     **unknown_options):
    """

    Minimization of scalar function of one or more variables using the
    simulated annealing algorithm.

    Options for the simulated annealing algorithm are:
        disp : bool
            Set to True to print convergence messages.
        schedule : str
            Annealing schedule to use. One of: 'fast', 'cauchy' or
            'boltzmann'.
        T0 : float
            Initial Temperature (estimated as 1.2 times the largest
            cost-function deviation over random points in the range).
        Tf : float
            Final goal temperature.
        maxfev : int
            Maximum number of function evaluations to make.
        maxaccept : int
            Maximum changes to accept.
        maxiter : int
            Maximum number of iterations to perform.
        boltzmann : float
            Boltzmann constant in acceptance test (increase for less
            stringent test at each temperature).
        learn_rate : float
            Scale constant for adjusting guesses.
        ftol : float
            Relative error in ``fun(x)`` acceptable for convergence.
        quench, m, n : float
            Parameters to alter fast_sa schedule.
        lower, upper : float or ndarray
            Lower and upper bounds on `x`.
        dwell : int
            The number of times to search the space at each temperature.

    This function is called by the `minimize` function with
    `method=anneal`. It is not supposed to be called directly.
    """
    _check_unknown_options(unknown_options)
    maxeval = maxfev
    feps = ftol

    x0 = asarray(x0)
    lower = asarray(lower)
    upper = asarray(upper)

    schedule = eval(schedule+'_sa()')
    #   initialize the schedule
    schedule.init(dims=shape(x0), func=func, args=args, boltzmann=boltzmann,
                  T0=T0, learn_rate=learn_rate, lower=lower, upper=upper,
                  m=m, n=n, quench=quench, dwell=dwell)

    current_state, last_state, best_state = _state(), _state(), _state()
    if T0 is None:
        x0 = schedule.getstart_temp(best_state)
    else:
        best_state.x = None
        best_state.cost = numpy.Inf

    last_state.x = asarray(x0).copy()
    fval = func(x0, *args)
    schedule.feval += 1
    last_state.cost = fval
    if last_state.cost < best_state.cost:
        best_state.cost = fval
        best_state.x = asarray(x0).copy()
    schedule.T = schedule.T0
    fqueue = [100, 300, 500, 700]
    iters = 0
    while 1:
        for n in xrange(dwell):
            current_state.x = schedule.update_guess(last_state.x)
            current_state.cost = func(current_state.x, *args)
            schedule.feval += 1

            dE = current_state.cost - last_state.cost
            if schedule.accept_test(dE):
                last_state.x = current_state.x.copy()
                last_state.cost = current_state.cost
                if last_state.cost < best_state.cost:
                    best_state.x = last_state.x.copy()
                    best_state.cost = last_state.cost
        schedule.update_temp()
        iters += 1
        # Stopping conditions
        # 0) last saved values of f from each cooling step
        #     are all very similar (effectively cooled)
        # 1) Tf is set and we are below it
        # 2) maxeval is set and we are past it
        # 3) maxiter is set and we are past it
        # 4) maxaccept is set and we are past it

        fqueue.append(squeeze(last_state.cost))
        fqueue.pop(0)
        af = asarray(fqueue)*1.0
        if all(abs((af-af[0])/af[0]) < feps):
            retval = 0
            if abs(af[-1]-best_state.cost) > feps*10:
                retval = 5
                if disp:
                    print("Warning: Cooled to %f at %s but this is not"
                          % (squeeze(last_state.cost),
                             str(squeeze(last_state.x)))
                          + " the smallest point found.")
            break
        if (Tf is not None) and (schedule.T < Tf):
            retval = 1
            break
        if (maxeval is not None) and (schedule.feval > maxeval):
            retval = 2
            break
        if (iters > maxiter):
            if disp:
                print("Warning: Maximum number of iterations exceeded.")
            retval = 3
            break
        if (maxaccept is not None) and (schedule.accepted > maxaccept):
            retval = 4
            break

    result = OptimizeResult(x=best_state.x, fun=best_state.cost,
                            T=schedule.T, nfev=schedule.feval, nit=iters,
                            accept=schedule.accepted, status=retval,
                            success=(retval <= 1),
                            message={0: 'Points no longer changing',
                                     1: 'Cooled to final temperature',
                                     2: 'Maximum function evaluations',
                                     3: 'Maximum cooling iterations reached',
                                     4: 'Maximum accepted query locations reached',
                                     5: 'Final point not the minimum amongst '
                                        'encountered points'}[retval])
    return result


if __name__ == "__main__":
    from numpy import cos
    # minimum expected at ~-0.195
    func = lambda x: cos(14.5 * x - 0.3) + (x + 0.2) * x
    print(anneal(func, 1.0, full_output=1, upper=3.0, lower=-3.0,
                 feps=1e-4, maxiter=2000, schedule='cauchy'))
    print(anneal(func, 1.0, full_output=1, upper=3.0, lower=-3.0,
                 feps=1e-4, maxiter=2000, schedule='fast'))
    print(anneal(func, 1.0, full_output=1, upper=3.0, lower=-3.0,
                 feps=1e-4, maxiter=2000, schedule='boltzmann'))

    # minimum expected at ~[-0.195, -0.1]
    func = lambda x: (cos(14.5 * x[0] - 0.3) + (x[1] + 0.2) * x[1] +
                      (x[0] + 0.2) * x[0])
    print(anneal(func, [1.0, 1.0], full_output=1,
                 upper=[3.0, 3.0], lower=[-3.0, -3.0],
                 feps=1e-4, maxiter=2000, schedule='cauchy'))
    print(anneal(func, [1.0, 1.0], full_output=1,
                 upper=[3.0, 3.0], lower=[-3.0, -3.0],
                 feps=1e-4, maxiter=2000, schedule='fast'))
    print(anneal(func, [1.0, 1.0], full_output=1,
                 upper=[3.0, 3.0], lower=[-3.0, -3.0],
                 feps=1e-4, maxiter=2000, schedule='boltzmann'))
