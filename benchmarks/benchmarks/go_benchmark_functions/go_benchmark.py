import numpy as np
from numpy import abs, asarray

from ..common import safe_import  # noqa:F401


class Benchmark:

    """
    Defines a global optimization benchmark problem.

    This abstract class defines the basic structure of a global
    optimization problem. Subclasses should implement the ``fun`` method
    for a particular optimization problem.

    Attributes
    ----------
    N : int
        The dimensionality of the problem.
    bounds : sequence
        The lower/upper bounds to be used for minimizing the problem.
        This a list of (lower, upper) tuples that contain the lower and upper
        bounds for the problem.  The problem should not be asked for evaluation
        outside these bounds. ``len(bounds) == N``.
    xmin : sequence
        The lower bounds for the problem
    xmax : sequence
        The upper bounds for the problem
    fglob : float
        The global minimum of the evaluated function.
    global_optimum : sequence
        A list of vectors that provide the locations of the global minimum.
        Note that some problems have multiple global minima, not all of which
        may be listed.
    nfev : int
        the number of function evaluations that the object has been asked to
        calculate.
    change_dimensionality : bool
        Whether we can change the benchmark function `x` variable length (i.e.,
        the dimensionality of the problem)
    custom_bounds : sequence
        a list of tuples that contain lower/upper bounds for use in plotting.
    """
    change_dimensionality = False

    def __init__(self, dimensions):
        """
        Initialises the problem

        Parameters
        ----------

        dimensions : int
            The dimensionality of the problem
        """

        self._dimensions = dimensions
        self.nfev = 0
        self.fglob = np.nan
        self.global_optimum = None
        self.custom_bounds = None

    def __str__(self):
        return f'{self.__class__.__name__} ({self.N} dimensions)'

    def __repr__(self):
        return self.__class__.__name__

    def initial_vector(self):
        """
        Random initialisation for the benchmark problem.

        Returns
        -------
        x : sequence
            a vector of length ``N`` that contains random floating point
            numbers that lie between the lower and upper bounds for a given
            parameter.
        """

        return asarray([np.random.uniform(l, u) for l, u in self.bounds])

    def success(self, x, tol=1.e-5):
        """
        Tests if a candidate solution at the global minimum.
        The default test is

        Parameters
        ----------
        x : sequence
            The candidate vector for testing if the global minimum has been
            reached. Must have ``len(x) == self.N``
        tol : float
            The evaluated function and known global minimum must differ by less
            than this amount to be at a global minimum.

        Returns
        -------
        bool : is the candidate vector at the global minimum?
        """
        val = self.fun(asarray(x))
        if abs(val - self.fglob) < tol:
            return True

        # the solution should still be in bounds, otherwise immediate fail.
        bounds = np.asarray(self.bounds, dtype=np.float64)
        if np.any(x > bounds[:, 1]):
            return False
        if np.any(x < bounds[:, 0]):
            return False

        # you found a lower global minimum.  This shouldn't happen.
        if val < self.fglob:
            raise ValueError("Found a lower global minimum",
                             x,
                             val,
                             self.fglob)

        return False

    def fun(self, x):
        """
        Evaluation of the benchmark function.

        Parameters
        ----------
        x : sequence
            The candidate vector for evaluating the benchmark problem. Must
            have ``len(x) == self.N``.

        Returns
        -------
        val : float
              the evaluated benchmark function
        """

        raise NotImplementedError

    def change_dimensions(self, ndim):
        """
        Changes the dimensionality of the benchmark problem

        The dimensionality will only be changed if the problem is suitable

        Parameters
        ----------
        ndim : int
               The new dimensionality for the problem.
        """

        if self.change_dimensionality:
            self._dimensions = ndim
        else:
            raise ValueError('dimensionality cannot be changed for this'
                             'problem')

    @property
    def bounds(self):
        """
        The lower/upper bounds to be used for minimizing the problem.
        This a list of (lower, upper) tuples that contain the lower and upper
        bounds for the problem.  The problem should not be asked for evaluation
        outside these bounds. ``len(bounds) == N``.
        """
        if self.change_dimensionality:
            return [self._bounds[0]] * self.N
        else:
            return self._bounds

    @property
    def N(self):
        """
        The dimensionality of the problem.

        Returns
        -------
        N : int
            The dimensionality of the problem
        """
        return self._dimensions

    @property
    def xmin(self):
        """
        The lower bounds for the problem

        Returns
        -------
        xmin : sequence
            The lower bounds for the problem
        """
        return asarray([b[0] for b in self.bounds])

    @property
    def xmax(self):
        """
        The upper bounds for the problem

        Returns
        -------
        xmax : sequence
            The upper bounds for the problem
        """
        return asarray([b[1] for b in self.bounds])
