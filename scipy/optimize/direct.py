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
try:
    from ._direct import direct
except ImportError:
    print('Fortran code not compiled, module not functional')
    direct = None


__version_info__ = ('1', '1')
__version__ = '.'.join(__version_info__)

ERROR_MESSAGES = (
    'u[i] < l[i] for some i',
    'maxf is too large',
    'Initialization failed',
    'There was an error in the creation of the sample points',
    'An error occured while the function was sampled',
    'Maximum number of levels has been reached.',
)

SUCCESS_MESSAGES = (
    'Number of function evaluations done is larger then maxf',
    'Number of iterations is equal to maxT',
    'The best function value found is within fglper of the (known) global optimum',
    'The volume of the hyperrectangle with best function value found is below volper',
    'The volume of the hyperrectangle with best function value found is smaller then volper'
)

# Class for returning the result of an optimization algorithm (copied from
# scipy.optimize)
class OptimizeResult(dict):
    r""" Represents the optimization result.

    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    status : int
        Termination status of the optimizer. Its value depends on the
        underlying solver. Refer to `message` for details.
    message : str
        Description of the cause of the termination.
    fun, jac, hess, hess_inv : ndarray
        Values of objective function, Jacobian, Hessian or its inverse (if
        available). The Hessians may be approximations, see the documentation
        of the function in question.
    nfev, njev, nhev : int
        Number of evaluations of the objective functions and of its
        Jacobian and Hessian.
    nit : int
        Number of iterations performed by the optimizer.
    maxcv : float
        The maximum constraint violation.

    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in self.items()])
        else:
            return self.__class__.__name__ + "()"

def _minimize_direct(func, bounds=None, nvar=None, args=(), disp=False,
             eps=1e-4,
             maxf=20000,
             maxT=6000,
             algmethod=0,
             fglobal=-1e100,
             fglper=0.01,
             volper=-1.0,
             sigmaper=-1.0,
             **kwargs
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
    
    bounds : array-like
            ``(min, max)`` pairs for each element in ``x``, defining
            the bounds on that parameter.

    nvar: integer
        Dimensionality of x (only needed if `bounds` is not defined)
        
    eps : float
        Ensures sufficient decrease in function value when a new potentially
        optimal interval is chosen.

    maxf : integer
        Approximate upper bound on objective function evaluations.
        
        .. note::
        
            Maximal allowed value is 90000 see documentation of Fortran library.
    
    maxT : integer
        Maximum number of iterations.
        
        .. note::
        
            Maximal allowed value is 6000 see documentation of Fortran library.
        
    algmethod : integer
        Whether to use the original or modified DIRECT algorithm. Possible values:
        
        * ``algmethod=0`` - use the original DIRECT algorithm
        * ``algmethod=1`` - use the modified DIRECT-l algorithm
    
    fglobal : float
        Function value of the global optimum. If this value is not known set this
        to a very large negative value.
        
    fglper : float
        Terminate the optimization when the percent error satisfies:
        
        .. math::

            100*(f_{min} - f_{global})/\max(1, |f_{global}|) \leq f_{glper}
        
    volper : float
        Terminate the optimization once the volume of a hyperrectangle is less
        than volper percent of the original hyperrectangel.
        
    sigmaper : float
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
                        maxf,
                        maxT,
                        l,
                        u,
                        algmethod,
                        'dummylogfile', 
                        fglobal,
                        fglper,
                        volper,
                        sigmaper,
                        iidata,
                        ddata,
                        cdata,
                        disp
                        )

    return OptimizeResult(x=x,fun=fun, status=ierror, success=ierror>0,
                          message=SUCCESS_MESSAGES[ierror-1] if ierror>0 else ERROR_MESSAGES[abs(ierror)-1])

