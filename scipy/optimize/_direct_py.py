# -*- coding: utf-8 -*-
r"""
scipydirect - A python wrapper to the DIRECT algorithm.
=======================================================

DIRECT is a method to solve global bound constraint optimization problems and
was originally developed by D. R. Jones, C. D. Perttunen and B. E. Stuckmann.
It is designed to find **global** solutions of mathematical optimization problems of the from

.. math::

       \min_ {x \in R^n} f(x)

subject to

.. math::

       x_L \leq  x  \leq x_U

Where :math:`x` are the optimization variables (with upper and lower
bounds), :math:`f(x)` is the objective function.

The DIRECT package uses the Fortran implementation of DIRECT written by
Joerg.M.Gablonsky, DIRECT Version 2.0.4. More information on the DIRECT
algorithm can be found in Gablonsky's `thesis <http://repository.lib.ncsu.edu/ir/bitstream/1840.16/3920/1/etd.pdf>`_.

.. codeauthor:: Andreas Mayer <andimscience@gmail.com>, Amit Aides <amitibo@tx.technion.ac.il>
"""

from __future__ import print_function
import numpy as np
from ._direct import direct
from .optimize import OptimizeResult

ERROR_MESSAGES = (
    "u[i] < l[i] for some i",
    "maxfun is too large",
    "Initialization failed",
    "There was an error in the creation of the sample points",
    "An error occured while the function was sampled",
    "Maximum number of levels has been reached.",
)

SUCCESS_MESSAGES = (
    "Number of function evaluations done is larger then maxfun",
    "Number of iterations is equal to maxiter",
    "The best function value found is within fglper of the (known) global optimum",
    "The volume of the hyperrectangle with best function value found is below volper",
    "The volume of the hyperrectangle with best function value found is smaller then volper",
)


def _minimize_direct(
    func,
    bounds=None,
    nvar=None,
    args=(),
    disp=False,
    eps=1e-4,
    maxfun=20000,
    maxiter=6000,
    method=0,
    fglobal=-1e100,
    fglper=0.01,
    volper=-1.0,
    sigmaper=-1.0,
):
    r"""

    Solve an optimization problem using the DIRECT (Dividing Rectangles) algorithm.
    It can be used to solve general nonlinear programming problems of the form:

    .. math::

           \min_ {x \in R^n} f(x)

    subject to

    .. math::

           x_L \leq  x  \leq x_U

    Where :math:`x` are the optimization variables (with upper and lower
    bounds), :math:`f(x)` is the objective function.

    Parameters
    ----------
    func : objective function
        called as `func(x, *args)`; does not need to be defined everywhere,
        raise an Exception where function is not defined
    bounds : array_like, optional
            ``(min, max)`` pairs for each element in ``x``, defining
            the bounds on that parameter.
    nvar : int, optional
        Dimensionality of x (only needed if `bounds` is not defined)
    eps : float, optional
        Ensures sufficient decrease in function value when a new potentially
        optimal interval is chosen.
    maxfun : int, optional
        Approximate upper bound on objective function evaluations.
        Maximal allowed value is 90000 see documentation of Fortran library.
    maxiter : int, optional
        Maximum number of iterations.
        Maximal allowed value is 6000 see documentation of Fortran library.
    method : integer, optional
        Whether to use the original or modified DIRECT algorithm. Possible values:

        * ``method=0`` - use the original DIRECT algorithm
        * ``method=1`` - use the modified DIRECT-l algorithm
    fglobal : float, optional
        Function value of the global optimum. If this value is not known set this
        to a very large negative value.
    fglper : float, optional
        Terminate the optimization when the percent error satisfies:

        .. math::

            100*(f_{min} - f_{global})/\max(1, |f_{global}|) \leq f_{glper}
    volper : float, optional
        Terminate the optimization once the volume of a hyperrectangle is less
        than volper percent of the original hyperrectangel.
    sigmaper : float, optional
        Terminate the optimization once the measure of the hyperrectangle is less
        than sigmaper.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.
    """

    if bounds is None:
        l = np.zeros(nvar, dtype=np.float64)
        u = np.ones(nvar, dtype=np.float64)
    else:
        bounds = np.asarray(bounds)
        l = bounds[:, 0]
        u = bounds[:, 1]

    def _objective_wrap(x, iidata, ddata, cdata, n, iisize, idsize, icsize):
        r"""
        Wrap the python objective to comply with the signature required by the
        Fortran library.

        Returns the function value and a flag indicating whether function is defined.
        If function is not defined return np.nan
        """
        try:
            return func(x, *args), 0
        except:
            return np.nan, 1

    #
    # Dummy values so that the python wrapper will comply with the required
    # signature of the fortran library.
    #
    iidata = np.ones(0, dtype=np.int32)
    ddata = np.ones(0, dtype=np.float64)
    cdata = np.ones([0, 40], dtype=np.uint8)

    #
    # Call the DIRECT algorithm
    #
    x, fun, ierror = direct(
        _objective_wrap,
        eps,
        maxfun,
        maxiter,
        l,
        u,
        method,
        "dummylogfile",
        fglobal,
        fglper,
        volper,
        sigmaper,
        iidata,
        ddata,
        cdata,
        disp,
    )

    if ierror > 0:
        message = SUCCESS_MESSAGES[ierror - 1]
    else:
        message = ERROR_MESSAGES[abs(ierror) - 1]

    return OptimizeResult(
        x=x, fun=fun, status=ierror, success=ierror > 0, message=message
    )
