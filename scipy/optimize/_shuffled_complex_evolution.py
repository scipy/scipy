#!/usr/bin/env python
"""
The Shuffled Complex Evolution (SCE) global optimizer

This code is based on a Fortran program of Qingyun Duan (2004), ported to
Python by Stijn Van Hoey (2011). It was taken up, debugged, enhanced and is
maintained by Matthias Cuntz while at Department of Computational Hydrosystems
(CHS), Helmholtz Centre for Environmental Research - UFZ, Leipzig, Germany, and
continued while at Institut National de Recherche pour l'Agriculture,
l'Alimentation et l'Environnement (INRAE), Nancy, France.

:copyright: Copyright 2004-2023 Qingyun Duan [1]_, Stijn Van Hoey, Matthias Cuntz, see AUTHORS.rst for details.
:license: MIT License, see LICENSE for details.

References
----------
.. [1] Duan, Sorooshian and Gupta (1992) Effective and efficient global
       optimization for conceptual rainfall-runoff models, Water Resour Res 28,
       1015-1031, https://doi.org/10.1029/91WR02985

.. moduleauthor:: Matthias Cuntz

The following functions are provided:

.. autosummary::
   shuffled_complex_evolution

History
    * Written in Fortran by Q Duan, Sep 2004
    * Ported to Python by Stijn Van Hoey, 2011
      https://github.com/stijnvanhoey/Optimization_SCE
    * Synchronised with enhanced Fortran version of CHS,
      Oct 2013, Matthias Cuntz
    * Added functionality to call external executable, Nov 2016, Matthias Cuntz
    * Treat NaN and Inf in function output, Nov 2016, Matthias Cuntz
    * Possibility to exclude (mask) parameters from optimisation,
      Nov 2016, Matthias Cuntz
    * Added restart possibility, Nov 2016, Matthias Cuntz
    * Return also function value of best parameter set if maxit==True,
      Nov 2016, Matthias Cuntz
    * Removed functionality to call external executable,
      Dec 2017, Matthias Cuntz
    * Print out number of function evaluations with printit=1,
      Mar 2018, Matthias Cuntz
    * Mask parameters with degenerated ranges, e.g. upper<lower bound,
      Mar 2018, Matthias Cuntz
    * Use only masked parameters in calculation of geometric range,
      Mar 2018, Matthias Cuntz
    * Removed bug that calculated the size of the complexes using all
      parameters and not only the masked parameters, Mar 2018, Matthias Cuntz
    * Fixed bug where masked parameters were always out of bounds,
      Mar 2018, Matthias Cuntz
    * Allow scalar bounds, which will be taken for all parameters,
      Mar 2018, Matthias Cuntz
    * Get number of parameters from lower boundary in SampleInputMatrix,
      Mar 2018, Matthias Cuntz
    * Define number of parameters by inquiring x0 in case of restart,
      May 2018, Matthias Cuntz
    * Removed multiplication with one hundred in criter_change,
      regarded a bug compared to Fortran code, May 2018, Matthias Cuntz
    * Removed exec command to make restart work with Python 3,
      May 2020, Matthias Cuntz
    * Use underscore before private routines, May 2020, Matthias Cuntz
    * Code refactoring, Sep 2021, Matthias Cuntz
    * Added keywords args and kwargs to pass to function,
      Apr 2022, Matthias Cuntz
    * Pass RandomState to _SampleInputMatrix, Jul 2022, Matthias Cuntz
    * Different sampling of input parameters, Jul 2022, Matthias Cuntz
    * Use helper class to pass arguments to objective function,
      Dec 2022, Matthias Cuntz
    * Copy strtobool from distutils.util because distutils is deprecated,
      Dec 2022, Matthias Cuntz
    * Include number of complexes ngs in restartfile, Dec 2022, Matthias Cuntz
    * No restart files written by default, Dec 2022, Matthias Cuntz
    * Rewrite into class SCESolver, Dec 2022, Matthias Cuntz
    * Output OptimizeResult class, Dec 2022, Matthias Cuntz
    * Polish results with L-BFGS-B, Dec 2022, Matthias Cuntz
    * Warn only if lb > ub, simply set mask if lb == ub,
      May 2023, Matthias Cuntz
    * Rename sce to shuffled_complex_evolution and SCESolver to
      ShuffledComplexEvolutionSolver, May 2023, Matthias Cuntz
    * Exit if initial population failed twice, May 2023, Matthias Cuntz
    * random_sample(1)[0] to assure scalar, Jul 2023, Matthias Cuntz
    * call_func method to assure scalar output and max or min,
      Jul 2023, Matthias Cuntz
    * Require keyword names after mask, Jul 2023, Matthias Cuntz

"""
import warnings
import numpy as np
from scipy.optimize import OptimizeResult, minimize
from scipy._lib._util import check_random_state
from scipy.optimize._constraints import Bounds

# ToDo:
# - write tmp/population files (of Fortran code)

__all__ = ['shuffled_complex_evolution']


# extended scipy/_util.py version
class _FunctionWrapper:
    """
    Wrap user function with arguments and keywords, allowing picklability

    Parameters
    ----------
    func : callable
        Function in the form ``func(x, *args, **kwargs)``, where ``x`` are
        the parameters in the form of an iterable.
        ``args`` and ``kwargs`` are passed to the function via the usual
        unpacking operators.

    """
    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs
        self.nargs = len(args)
        self.nkwargs = len(kwargs)
        if self.nkwargs == 0:
            if ( (len(self.args) == 2) and
                 isinstance(self.args[-1], dict) and
                 (len(self.args[-1]) == 0) ):
                # if kwargs={} then **kwargs={} and hence counted as args
                self.args = self.args[0]
                self.nargs = len(args)

    def __call__(self, x):
        if (self.nargs > 0) and (self.nkwargs > 0):
            return self.func(x, *self.args, **self.kwargs)
        elif (self.nargs > 0) and (self.nkwargs == 0):
            return self.func(x, *self.args)
        elif (self.nargs == 0) and (self.nkwargs > 0):
            return self.func(x, **self.kwargs)
        else:
            return self.func(x)


# from distutils.util (Python 3.11.1)
def _strtobool(val):
    """
    Convert a string representation of truth to true (1) or false (0).
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.

    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value {!r}".format(val))


def shuffled_complex_evolution(
        func, x0, lb, ub=None,
        mask=None, *,
        args=(), kwargs={},
        sampling='half-open',
        maxn=1000, kstop=10, pcento=0.0001, peps=0.001,
        ngs=2, npg=0, nps=0, nspl=0, mings=0,
        seed=None, iniflg=True,
        alpha=0.8, beta=0.45, maxit=False, printit=2,
        polish=True,
        restart=False, restartfile1='',
        restartfile2=''):
    """
    Shuffled Complex Evolution algorithm for finding the minimum of a
    multivariate function

    The SCE or SCE-UA method is a general purpose global optimization, which
    can be used with high-dimensional problems. It was used successfully,
    for example, to calibrate hydrologic models with more than 50 parameters
    [5]_.
    The algorithm has been described in detail in Duan et al. [1]_ and [2]_.
    Another paper of Duan et al. [3]_ discusses how to use the method
    effectively. The implementation here also includes the recommendations of
    Behrangi et al. [4]_.

    Parameters
    ----------
    func : callable
        Function in the form ``func(x, *args, **kwargs)``, where ``x`` are
        the parameters in the form of an iterable.
        ``args`` and ``kwargs`` are passed to the function via the usual
        unpacking operators.
    x0 : array_like
        Parameter values with `mask==1` used in initial complex if
        `iniflg==1`.
    lb : array_like, sequence or `Bounds`
        Lower bounds of parameters if ``ub`` given.
        If `ub` is not given, then `lb` are the bounds of the parameters,
        either 1) as an instance of the `Bounds` class, or 2) as ``(min, max)``
        pairs for each element in ``x``, defining the finite lower and upper
        bounds for the parameters of `func`.
    ub : array_like, optional
        Upper bounds of parameters.
    mask : array_like, optional
        Include (1, True) or exclude (0, False) parameters in minimization
        (default: include all parameters). The number of parameters ``nopt`` is
        ``sum(mask)``.
    args : tuple, optional
        Extra arguments passed to the function *func*. Note that ``args`` must
        be iterable. `args=int`and `args=(int)` are not valid (with int being
        any scalar variable) but should be `args=(int,)`.
    kwargs : dict, optional
        Extra keyword arguments passed to the function `func`.
    sampling : string or array_like of strings, optional
        Options for sampling random numbers. Options can be on of:

            - 'half-open': same as 'right-half-open' [lb, ub)
            - 'left-half-open': sample random floats in half-open
              interval (lb, ub]
            - 'right-half-open': sample random floats in half-open
              interval [lb, ub)
            - 'open': sample random floats in open interval (lb, ub)
            - 'log': sample half-open interval [log(lb), log(ub)), which
              samples better intervals spanning orders or magnitude such as
              `lb=1e-9` and `ub=1e-4`

        The default is 'half-open'.
    maxn : int, optional
        Maximum number of function evaluations allowed during minimization
        (without polishing) (default: 1000).
    kstop : int, optional
        Maximum number of evolution loops before convergence (default: 10).
    pcento : float, optional
        Percentage change allowed in kstop loops before convergence
        (default: 0.0001).
    peps : float, optional
        Value of normalised geometric range needed for convergence
        (default: 0.001).
    ngs : int, optional
        Number of complexes (default: 2).
    npg : int, optional
        Number of points in each complex (default: `2*nopt+1`).
    nps : int, optional
        Number of points in each sub-complex (default: `nopt+1`).
    mings : int, optional
        Minimum number of complexes required if the number of complexes is
        allowed to reduce as the optimization proceeds (default: `ngs`).
    nspl : int, optional
        Number of evolution steps allowed for each complex before complex
        shuffling (default: `2*nopt+1`).
    seed : {None, int, `numpy.random.Generator`, `numpy.random.RandomState`}, optional
        If `seed` is None (or `numpy.random`), the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new `RandomState` instance is used,
        seeded with `seed`.
        If `seed` is already a `Generator` or `RandomState` instance then
        that instance is used.
        Specify `seed` for repeatable results.
    iniflg : bool, optional
        If True: include initial parameters ``x0`` in initial population
        (default: True).
    alpha : float, optional
        Parameter for reflection of points in complex (default: 0.8).
    beta : float, optional
        Parameter for contraction of points in complex (default: 0.45).
    maxit : bool, optional
        If True: maximize instead of minimize `func` (default: False).
    printit : int, optional
        Controlling print-out (default: 2):

            - 0: print information for the best point of the population
            - 1: print information for each function evaluation
            - 2: no printing.

        The default is 2.
    polish : bool, optional
        If True (default), then `scipy.optimize.minimize` is used with the
        `L-BFGS-B` method to polish the result at the end, which
        can improve the minimization slightly. For large problems, polishing
        can take a long time due to the computation of the Jacobian.
    restart : bool, optional
        If True, continue from saved state in `restartfile1` and
        `restartfile2` (default: False).
    restartfile1 : str, optional
        Filename for saving state of array variables of optimizer
        (default: '').
        If `restart==True` and `restartfile1==''` then
        `restartfile1='sce.restart.npz'` will be taken.
    restartfile2 : int, optional
        Filename for saving state of non-array variables of optimizer
        (default: `restartfile1 + '.txt'`).

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a `OptimizeResult` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully (0) and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes. If `polish`
        was employed, and a lower minimum was obtained by the polishing, then
        OptimizeResult also contains the ``jac`` attribute.

    References
    ----------
    .. [1] Duan QY, Sorooshian S, and Gupta VK,
           Effective and efficient global optimization for conceptual
           rainfall-runoff models, Water Resources Research 28(4), 1015-1031,
           1992, https://doi.org/10.1029/91WR02985
    .. [2] Duan QY, Gupta VK, and Sorooshian S,
           Shuffled Complex Evolution approach for effective and efficient
           global minimization, Journal of Optimization Theory and its
           Applications 76(3), 501-521, 1993,
           https://doi.org/10.1007/BF00939380
    .. [3] Duan QY, Gupta VK, and Sorooshian S,
           Optimal use of the SCE-UA global optimization method for calibrating
           watershed models, Journal of Hydrology 158, 265-284, 1994,
           https://doi.org/10.1016/0022-1694(94)90057-4
    .. [4] Behrangi A, Khakbaz B, Vrugt JA, Duan Q, and Sorooshian S,
           Comment on "Dynamically dimensioned search algorithm for
           computationally efficient watershed model calibration" by
           Bryan A. Tolson and Christine A. Shoemaker, Water Resources
           Research 44, W12603, 2008, http://doi.org/10.1029/2007WR006429
    .. [5] Cuntz M, Mai J, Zink M, Thober S, Kumar R, Schäfer D, Schrön M,
           Craven J, Rakovec O, Spieler D, Prykhodko V, Dalmasso G, Musuuza J,
           Langenberg B, Attinger S, and Samaniego L,
           Computationally inexpensive identification of noninformative model
           parameters by sequential screening. Water Resources Research 51,
           6417–6441, 2015, https://doi.org/10.1002/2015WR016907

    Examples
    --------
    Search the minimum of the Rosenbrock function, implemented in `rosen`
    in `scipy.optimize`.

    >>> import numpy as np
    >>> from scipy.optimize import rosen

    One has to provide the function `rosen`, an initial guess of the
    parameters, and the lower and upper limits of the parameters.
    The 2D version is:

    >>> lb = np.array([-5., -2.])
    >>> ub = np.array([5., 8.])
    >>> x0 = np.array([-2., 7.])
    >>> res = shuffled_complex_evolution(rosen, x0, lb, ub, seed=1, maxn=1000,
    ...     printit=2)
    >>> print(res.nfev)
    298
    >>> print('{:.3f}'.format(res.fun))
    0.001

    A 10-dimensional version using `(min, max)` pairs for parameter bounds
    as well as setting a number of keyword parameters for the SCE algorithm is:

    >>> nopt = 10
    >>> lb = np.full(10, -5.)
    >>> ub = np.full(10, 5.)
    >>> x0 = np.full(10, 0.5)
    >>> res = shuffled_complex_evolution(rosen, x0, zip(lb, ub), maxn=30000,
    ...     kstop=10, pcento=0.0001,
    ...     seed=12358, ngs=5, npg=5*nopt+1, nps=nopt+1, nspl=5*nopt+1,
    ...     mings=2, iniflg=True, printit=2, alpha=0.8, beta=0.45)
    >>> print(res.nfev)
    30228
    >>> print('{:.3g}'.format(res.fun))
    3.38e-12

    """
    # using a context manager means that any created Pool objects are
    # cleared up.
    ret = None
    with ShuffledComplexEvolutionSolver(
            func, x0, lb, ub=ub,
            mask=mask, args=args, kwargs=kwargs, sampling=sampling,
            maxn=maxn, kstop=kstop, pcento=pcento,
            ngs=ngs, npg=npg, nps=nps, nspl=nspl, mings=mings,
            peps=peps, seed=seed, iniflg=iniflg,
            alpha=alpha, beta=beta, maxit=maxit, printit=printit,
            polish=polish,
            restart=restart, restartfile1=restartfile1,
            restartfile2=restartfile2) as solver:
        ret = solver.solve()

    return ret


class ShuffledComplexEvolutionSolver:
    """
    This class implements the Shuffled-Complex-Evolution algorithm

    Parameters
    ----------
    func : callable
        Function in the form ``func(x, *args, **kwargs)``, where ``x`` are
        the parameters in the form of an iterable.
        ``args`` and ``kwargs`` are passed to the function via the usual
        unpacking operators.
    x0 : array_like
        Parameter values with `mask==1` used in initial complex if
        `iniflg==1`.
    lb : array_like, sequence or `Bounds`
        Lower bounds of parameters if ``ub`` given.
        If `ub` is not given, then `lb` are the bounds of the parameters,
        either 1) as an instance of the `Bounds` class, or 2) as ``(min, max)``
        pairs for each element in ``x``, defining the finite lower and upper
        bounds for the parameters of `func`.
    ub : array_like
        Upper bounds of parameters.
    mask : array_like, optional
        Include (1, True) or exclude (0, False) parameters in minimization
        (default: include all parameters). The number of parameters ``nopt`` is
        ``sum(mask)``.
    args : tuple, optional
        Extra arguments passed to the function *func*.
    kwargs : dict, optional
        Extra keyword arguments passed to the function `func`.
    sampling : string or array_like of strings, optional
        Options for sampling random numbers. Options can be on of:

            - 'half-open': same as 'right-half-open' [lb, ub)
            - 'left-half-open': sample random floats in half-open
              interval (lb, ub]
            - 'right-half-open': sample random floats in half-open
              interval [lb, ub)
            - 'open': sample random floats in open interval (lb, ub)
            - 'log': sample half-open interval [log(lb), log(ub)), which
              samples better intervals spanning orders or magnitude such as
              `lb=1e-9` and `ub=1e-4`

        The default is 'half-open'.
    maxn : int, optional
        Maximum number of function evaluations allowed during minimization
        (without polishing) (default: 1000).
    kstop : int, optional
        Maximum number of evolution loops before convergence (default: 10).
    pcento : float, optional
        Percentage change allowed in kstop loops before convergence
        (default: 0.0001).
    peps : float, optional
        Value of normalised geometric range needed for convergence
        (default: 0.001).
    ngs : int, optional
        Number of complexes (default: 2).
    npg : int, optional
        Number of points in each complex (default: `2*nopt+1`).
    nps : int, optional
        Number of points in each sub-complex (default: `nopt+1`).
    mings : int, optional
        Minimum number of complexes required if the number of complexes is
        allowed to reduce as the optimization proceeds (default: `ngs`).
    nspl : int, optional
        Number of evolution steps allowed for each complex before complex
        shuffling (default: `2*nopt+1`).
    seed : {None, int, `numpy.random.Generator`, `numpy.random.RandomState`}, optional
        If `seed` is None (or `numpy.random`), the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new `RandomState` instance is used,
        seeded with `seed`.
        If `seed` is already a `Generator` or `RandomState` instance then
        that instance is used.
        Specify `seed` for repeatable results.
    iniflg : bool, optional
        If True: include initial parameters ``x0`` in initial population
        (default: True).
    alpha : float, optional
        Parameter for reflection of points in complex (default: 0.8).
    beta : float, optional
        Parameter for contraction of points in complex (default: 0.45).
    maxit : bool, optional
        If True: maximize instead of minimize `func` (default: False).
    printit : int, optional
        Controlling print-out (default: 2):

            - 0: print information for the best point of the population
            - 1: print information for each function evaluation
            - 2: no printing.

        The default is 2.
    polish : bool, optional
        If True (default), then `scipy.optimize.minimize` is used with the
        `L-BFGS-B` method to polish the result at the end, which
        can improve the minimization slightly. For large problems, polishing
        can take a long time due to the computation of the Jacobian.
    restart : bool, optional
        If True, continue from saved state in `restartfile1` and
        `restartfile2` (default: False).
    restartfile1 : str, optional
        Filename for saving state of array variables of optimizer
        (default: '').
        If `restart==True` and `restartfile1==''` then
        `restartfile1='sce.restart.npz'` will be taken.
    restartfile2 : int, optional
        Filename for saving state of non-array variables of optimizer
        (default: `restartfile1 + '.txt'`).

    """

    def __init__(self, func, x0, lb, ub=None,
                 mask=None, args=(), kwargs={},
                 sampling='half-open',
                 maxn=1000, kstop=10, pcento=0.0001,
                 ngs=2, npg=0, nps=0, nspl=0, mings=0,
                 peps=0.001, seed=None, iniflg=True,
                 alpha=0.8, beta=0.45, maxit=False, printit=2,
                 polish=True,
                 restart=False, restartfile1='',
                 restartfile2=''):

        # function to minimize
        self.func = _FunctionWrapper(func, *args, **kwargs)
        # seed random number generator - will be overwritten by restart
        # self.rnd = np.random.RandomState(seed=seed)
        self.rnd = check_random_state(seed)
        # parameters for initial run and for restart
        self.sampling = sampling
        self.maxn = maxn
        self.kstop = kstop
        self.pcento = pcento
        self.peps = peps
        self.alpha = alpha
        self.beta = beta
        self.maxit = maxit
        self.printit = printit
        self.polish = polish
        self.restartfile1 = restartfile1
        self.restartfile2 = restartfile2
        if restart and (not self.restartfile1):
            self.restartfile1 = 'sce.restart.npz'
        if self.restartfile1 and (not self.restartfile2):
            self.restartfile2 = self.restartfile1 + '.txt'
        if not restart:
            # initialize SCE parameters
            self.nn    = len(x0)
            self.mask  = np.ones(self.nn, dtype=bool) if mask is None else mask
            self.nopt  = np.sum(self.mask)
            self.ngs   = ngs if ngs > 0 else 2
            self.npg   = npg if npg > 0 else 2 * self.nopt + 1
            self.nps   = nps if nps > 0 else self.nopt + 1
            self.nspl  = nspl if nspl > 0 else 2 * self.nopt + 1
            self.mings = mings if mings > 0 else self.ngs
            self.npt   = self.npg * self.ngs

            # assure lb and ub are numpy arrays
            if ub is None:
                if isinstance(lb, Bounds):
                    self.lb = lb.lb
                    self.ub = lb.ub
                else:
                    self.lb, self.ub = zip(*lb)
            else:
                self.lb = lb
                self.ub = ub
            self.lb = np.array(self.lb)
            self.ub = np.array(self.ub)

            # same bounds for all parameters
            try:
                if len(self.lb) == 1:
                    self.lb = np.full(self.nn, self.lb[0])
            except TypeError:  # could be size 0 array if lb was a scalar
                self.lb = np.full(self.nn, self.lb)
            try:
                if len(self.ub) == 1:
                    self.ub = np.full(self.nn, self.ub[0])
            except TypeError:
                self.ub = np.full(self.nn, self.ub)
            bound = self.ub - self.lb
            # degenerated bounds
            if np.any(bound == 0.):
                ii = np.where(bound == 0.)[0]
                self.mask[ii] = False
            if np.any(bound < 0.):
                ii = np.where(bound < 0.)[0]
                warnings.warn(
                    f'shuffled_complex_evolution: found lower bound lb >'
                    f' upper bound ub for parameter(s) {ii} with'
                    f' lb={self.lb[ii]} and ub={self.ub[ii]} => masking'
                    f' the parameter(s).',
                    UserWarning, stacklevel=2)
                self.mask[ii] = False
                self.ub[ii] = self.lb[ii]

            # set 'large' for function runs that return NaN
            self.large = 0.5 * np.finfo(float).max

            # create an initial population to fill array x(npt, nparams)
            self.x = self.sample_input_matrix(self.npt)
            for i in range(self.npt):
                self.x[i, :] = np.where(self.mask, self.x[i, :], x0)
            if iniflg == 1:
                self.x[0, :] = x0

            self.icall = 0
            self.xf = np.zeros(self.npt)
            for i in range(self.npt):
                self.xf[i] = self.call_func(self.x[i, :])
                self.icall += 1
                if self.printit == 1:
                    print('  i, f, X: ', self.icall, self.xf[i], self.x[i, :])

            # redo an initial population if all runs failed
            if not np.any(np.isfinite(self.xf)):
                if self.printit < 2:
                    print('Redo initial population because all failed')
                self.x = self.sample_input_matrix(self.npt)
                for i in range(self.npt):
                    self.x[i, :] = np.where(self.mask, self.x[i, :], x0)
                for i in range(self.npt):
                    self.xf[i] = self.call_func(self.x[i, :])
                    self.icall += 1
                    if self.printit == 1:
                        print('  i, f, X: ', self.icall, self.xf[i],
                              self.x[i, :])
            if not np.any(np.isfinite(self.xf)):
                raise ValueError(
                    'Did not succeed to produce initial population:'
                    ' all function evaluations failed. Give an initial'
                    ' value x0 that works and set iniflg=True (default).')

            # remember large for treating of NaNs
            self.large = self.xf[np.isfinite(self.xf)].max()
            self.large = (1.1 * self.large if self.large > 0. else
                          0.9 * self.large)

            # sort the population in order of increasing function values
            # and report best point
            self.xf = np.where(np.isfinite(self.xf), self.xf, self.large)
            idx = np.argsort(self.xf)
            self.xf = self.xf[idx]
            self.x  = self.x[idx, :]
            self.bestx  = self.x[0, :]
            self.bestf  = self.xf[0]

            # compute the normalized geometric range of the parameters
            self.gnrng = self.calc_gnrng()

            # initialise evolution loops
            self.nloop = 0
            self.criter = []
            self.criter_change = 1e+5

            # save restart
            self.write_restartfiles()

        else:  # if not restart
            self.nn = len(x0)
            self.read_restartfiles()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return self

    # def __iter__(self):
    #     return self

    def __next__(self):
        # loop on complexes (sub-populations)
        for igs in range(self.ngs):
            k1 = np.array(range(self.npg))
            k2 = k1 * self.ngs + igs

            # partition the population into complexes (sub-populations)
            cx = np.zeros((self.npg, self.nn))
            cf = np.zeros((self.npg))
            k1 = np.array(range(self.npg))
            k2 = k1 * self.ngs + igs
            cx[k1, :] = self.x[k2, :]
            cf[k1]    = self.xf[k2]

            # evolve sub-population igs for nspl steps:
            for loop in range(self.nspl):
                # select simplex by sampling the complex according to a
                # linear probability distribution
                lcs = np.zeros(self.nps, dtype=int)
                lcs[0] = 1
                for k3 in range(1, self.nps):
                    for i in range(1000):
                        lpos = int(np.floor(
                            self.npg + 0.5 -
                            np.sqrt((self.npg + 0.5)**2 -
                                    self.npg * (self.npg + 1) *
                                    self.rnd.random_sample(1)[0]) ))
                        # check if element was already chosen
                        idx = (lcs[0:k3] == lpos).nonzero()
                        if idx[0].size == 0:
                            break
                    lcs[k3] = lpos
                lcs.sort()

                # construct the simplex
                s  = np.zeros((self.nps, self.nn))
                s  = cx[lcs, :]
                sf = cf[lcs]

                # remember large for treating of NaNs
                self.large = cf[np.isfinite(cf)].max()
                self.large = (1.1 * self.large if self.large > 0. else
                              0.9 * self.large)

                snew, fnew, icall = self.cce(s, sf)
                # replace the worst point in simplex with the new point
                s[-1, :] = snew
                sf[-1]   = fnew
                # self.icall += icall
                self.icall += icall

                # reinsert the simplex into the complex
                cx[lcs, :] = s
                cf[lcs]    = sf

                # sort the complex
                cf  = np.where(np.isfinite(cf), cf, self.large)
                idx = np.argsort(cf)
                cf  = cf[idx]
                cx  = cx[idx, :]
                # end of inner loop for competitive evolution of simplexes:
                #     for loop in range(self.nspl):
                # i.e. end of evolve sub-population igs for nspl steps

            self.x[k2, :] = cx[k1, :]
            self.xf[k2]   = cf[k1]

        return

    def calc_criter_change(self):
        """
        Computes the percentage change in function output

        """
        if self.nloop >= self.kstop:
            criter_change = np.abs(self.criter[self.nloop - 1] -
                                   self.criter[self.nloop - self.kstop])
            criter_change = (
                criter_change /
                np.maximum( 1e-15, np.mean( np.abs(
                    self.criter[self.nloop - self.kstop:self.nloop]) ) ) )
        else:
            criter_change = self.criter_change

        return criter_change

    def calc_gnrng(self):
        """
        Computes the normalized geometric range of the parameters

        """
        bound = self.ub - self.lb
        rrange = (np.ma.array(self.x.max(axis=0) - self.x.min(axis=0),
                              mask=~self.mask) /
                  np.ma.array(bound, mask=~self.mask))
        gnrng = np.ma.exp(np.ma.mean(np.ma.log(rrange)))

        return gnrng

    def call_func(self, x):
        """
        Call function `func` asserting scalar output and maximum or minimum

        """
        fuc = self.func(x)
        if isinstance(fuc, np.ndarray):
            if fuc.size > 1:
                raise RuntimeError(
                    'func(x, *args, **kwargs) must return a'
                    ' scalar value.')
            fuc = fuc[0]
        fuc = -fuc if self.maxit else fuc
        return fuc

    def cce(self, s, sf):
        """
        Generate a new point in a simplex

        Parameters
        ----------
        s : 2D-array
            The sorted simplex in order of increasing function values
        sf : 1D-array
            Function values in increasing order

        Results
        -------
        new parameter set, function value of new set,
        number of function evaluations performed

        """
        # best point
        sb = s[0, :]
        # fb = sf[0]
        # worst point and function value
        sw = s[-1, :]
        fw = sf[-1]

        # centroid of the simplex excluding worst point
        ce = np.mean(s[:-1, :], axis=0)

        # attempt a reflection point
        snew = ce + self.alpha * (ce - sw)
        # sb should have initial params at mask==False
        snew = np.where(self.mask, snew, sb)

        # check if new point is outside bounds
        ibound = 0
        if np.ma.any(np.ma.array(snew - self.lb, mask=~self.mask) < 0.):
            ibound = 1
        if np.ma.any(np.ma.array(self.ub - snew, mask=~self.mask) < 0.):
            ibound = 2
        if ibound >= 1:
            snew = self.sample_input_matrix(1)[0, :]
            snew = np.where(self.mask, snew, sb)

        icall = 0
        # calc function for reflection point
        fnew = self.call_func(snew)
        icall += 1
        if self.printit == 1:
            print('  i, f, X: ', self.icall + icall, fnew, snew)

        # reflection failed: attempt a contraction point
        if fnew > fw:
            snew = sw + self.beta * (ce - sw)
            snew = np.where(self.mask, snew, sb)
            fnew = self.call_func(snew)
            icall += 1
            if self.printit == 1:
                print('  i, f, X: ', self.icall + icall, fnew, snew)

        # both reflection and contraction have failed, attempt a random point
        if fnew > fw:
            snew = self.sample_input_matrix(1)[0, :]
            snew = np.where(self.mask, snew, sb)
            fnew = self.call_func(snew)
            icall += 1
            if self.printit == 1:
                print('  i, f, X: ', self.icall + icall, fnew, snew)

        # end of cce
        return snew, fnew, icall

    def check_number_calls(self):
        """
        Check for maximum number of function calls reached

        """
        if self.icall < self.maxn:
            return True
        else:
            if self.printit < 2:
                print(f'Optimisation terminated because trial number'
                      f' {self.maxn} reached maximum number of trials'
                      f' {self.icall}.')
            return False

    def check_criter_change(self):
        """
        Check if percentage change in function output is below pcento

        """
        if self.criter_change < self.pcento:
            if self.printit < 2:
                print(f'The best point has improved by less then {self.pcento}'
                      f' in the last {self.kstop} loops.')
            return True
        else:
            return False

    def check_geometric_range(self):
        """
        Check if normalized geometric range of the parameters is below peps

        """
        self.gnrng = self.calc_gnrng()
        if self.gnrng < self.peps:
            if self.printit < 2:
                print(f'The population has converged to a small parameter'
                      f' space {self.gnrng} (<{self.peps}).')
            return True
        else:
            return False

    def read_restartfiles(self):
        """
        Read from restart files

        """
        with open(self.restartfile1, 'rb') as p1:
            pp = np.load(p1)
            self.lb       = pp['lb']
            self.ub       = pp['ub']
            self.mask     = pp['mask']
            self.criter   = pp['criter']
            self.x        = pp['x']
            self.xf       = pp['xf']
            rs2           = pp['rs2']

        with open(self.restartfile2, 'r') as p2:
            (self.ngs, self.nopt, self.npg, self.nps, self.nspl,
             self.mings, self.npt, self.nloop, self.icall, rs3, rs4) = [
                 int(inn) for inn in p2.readline().rstrip().split(',') ]
            (self.gnrng, self.criter_change, self.bestf, rs5) = [
                float(inn) for inn in p2.readline().rstrip().split(',') ]
            self.maxit = bool(_strtobool(p2.readline().rstrip()))
            rs1 = p2.readline().rstrip()

        self.rnd.set_state((rs1, rs2, rs3, rs4, rs5))

    def write_restartfiles(self):
        """
        Write restart files if demanded

        """
        if self.restartfile1:
            # compute the normalized geometric range of the parameters
            self.gnrng = self.calc_gnrng()
            # compute percentage change in function output
            self.criter_change = self.calc_criter_change()
            # get the current state of the random number generator
            rs1, rs2, rs3, rs4, rs5 = self.rnd.get_state()

            # only arrays into savez_compressed - restartfile1
            np.savez_compressed(self.restartfile1, lb=self.lb, ub=self.ub,
                                mask=self.mask,
                                criter=self.criter, x=self.x, xf=self.xf,
                                bestx=self.bestx, rs2=rs2)

            # scalars into text file - restartfile2
            with open(self.restartfile2, 'w') as p:
                iout = [self.ngs, self.nopt, self.npg, self.nps, self.nspl,
                        self.mings, self.npt, self.nloop, self.icall,
                        rs3, rs4]
                iform = '{:d},' * (len(iout) - 1) + '{:d}'
                print(iform.format(*iout), file=p)
                fout = [self.gnrng, self.criter_change, self.bestf, rs5]
                fform = '{:.15g},' * (len(fout) - 1) + '{:.15g}'
                print(fform.format(*fout), file=p)
                print(self.maxit, file=p)
                print(rs1, file=p)

        return

    def sample_input_matrix(self, npt=0):
        """
        Create input parameter matrix (npt, npars) for
        npt simulations and npars parameters with bounds lb and ub

        Returns
        -------
        parameter matrix (npt, npars) with new parameter samples

        """
        if npt < 1:
            npt = self.npt
        npars = len(self.lb)
        if isinstance(self.sampling, str):
            isampling = [self.sampling] * npars
        else:
            isampling = self.sampling
        assert len(isampling) == npars, (
            f'sampling must be string or list of strings'
            f' with {npars} entries')
        x = np.zeros((npt, npars))
        bound = self.ub - self.lb
        for i in range(npt):
            irnd = self.rnd.random_sample(npars)
            for j in range(npars):
                opt = isampling[j].lower()
                if (opt == 'half-open') or (opt == 'right-half-open'):
                    x[i, j] = self.lb[j] + irnd[j] * bound[j]
                elif opt == 'left-half-open':
                    irnd[j] = 1. - irnd[j]
                    x[i, j] = self.lb[j] + irnd[j] * bound[j]
                elif opt == 'open':
                    iirnd = irnd[j]
                    while not (iirnd > 0.):
                        iirnd = self.rnd.random_sample(1)[0]
                    x[i, j] = self.lb[j] + iirnd * bound[j]
                elif opt == 'log':
                    # x must be > 0. for ln(x)
                    xshift = 0.
                    if (self.lb[j] * self.ub[j]) < 0.:
                        # lb < 0 and ub > 0 -> shift both > 0
                        xshift = 2. * np.maximum(np.abs(self.lb[j]),
                                                 np.abs(self.ub[j]))
                    elif (self.lb[j] * self.ub[j]) == 0.:
                        if self.lb[j] == 0.:
                            # lb == 0 and ub > 0 -> shift to [ub, 2*ub)
                            xshift = self.ub[j]
                        if self.ub[j] == 0.:
                            # lb < 0 and ub == 0 -> shift to [-lb, -2*lb) > 0.
                            xshift = -2. * self.lb[j]
                    elif (self.lb[j] < 0.) and (self.ub[j] < 0.):
                        # lb < 0 and ub < 0 -> shift both > 0
                        xshift = 2. * np.maximum(np.abs(self.lb[j]),
                                                 np.abs(self.ub[j]))
                    lnlb = np.log(self.lb[j] + xshift)
                    lnub = np.log(self.ub[j] + xshift)
                    x[i, j] = np.exp(lnlb + irnd[j] * (lnub - lnlb)) - xshift
                else:
                    raise ValueError(f'unknown sampling option'
                                     f' {isampling[j]}.\n'
                                     f'Known samplings are: half-open,'
                                     f' left-half-open, right-half-open, open,'
                                     f' log')

        return x

    def solve(self):
        while (self.check_number_calls() and
               (not self.check_geometric_range()) and
               (not self.check_criter_change())):
            self.nloop += 1

            # loop on complexes (sub-populations)
            next(self)

            # shuffle the complexes
            idx = np.argsort(self.xf)
            self.xf = self.xf[idx]
            self.x  = self.x[idx, :]

            # record best and worst points
            self.bestx = self.x[0, :]
            self.bestf = self.xf[0]

            if self.printit < 2:
                print(f'\nEvolution loop {self.nloop}, trials {self.icall}.'
                      f' Best f: {self.bestf:.2f}, worst f:'
                      f' {self.xf[-1]:.2f}\n'
                      f'  best x:\n', self.bestx, '\n')

            self.criter = np.append(self.criter, self.bestf)
            self.criter_change = self.calc_criter_change()

            # save restart
            self.write_restartfiles()

        # finish up
        if self.printit < 2:
            print('Search stopped at trial number {0:d} with normalized'
                  ' geometric range {1:f}. '.format(self.icall, self.gnrng))
            print('The best point has improved by {:f} in the last'
                  ' {:d} loops.'.format(self.criter_change, self.kstop))

        if not self.check_number_calls():
            success = False
            status  = 2
            message = f'Reached maximum number of trials {self.maxn}'
        if self.check_criter_change():
            success = True
            status  = 1
            message = (f'Best point improved less than {self.pcento}'
                       f' in last {self.kstop} loops')
        if self.check_geometric_range():
            success = True
            status  = 0
            message = (f'Normalized geometric range of parameters'
                       f' {self.gnrng} < {self.peps}')

        if self.maxit:
            self.bestf *= -1.

        sce_result = OptimizeResult(
            x=self.bestx,
            success=success,
            status=status,
            message=message,
            fun=self.bestf,
            nfev=self.icall,
            nit=self.nloop)

        # polish results, only works if no mask
        if self.polish and (len(self.lb) == self.nopt):
            polish_method = 'L-BFGS-B'
            result = minimize(self.func,
                              np.copy(sce_result.x),
                              method=polish_method,
                              bounds=zip(self.lb, self.ub))

            self.icall += result.nfev
            sce_result.nfev = self.icall

            # Polishing solution is only accepted if there is an improvement in
            # the cost function, the polishing was successful and the solution
            # lies within the bounds.
            if ( (result.fun < sce_result.fun) and
                 result.success and
                 np.all(result.x <= self.ub) and
                 np.all(self.lb <= result.x) ):
                sce_result.fun = result.fun
                sce_result.x = result.x
                sce_result.jac = result.jac

        return sce_result


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
