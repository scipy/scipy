"""rbf - Radial basis functions for interpolation/smoothing scattered Nd data.

Written by John Travers <jtravs@gmail.com>, February 2007
Based closely on Matlab code by Alex Chirokov
Additional, large, improvements by Robert Hetland
Some additional alterations by Travis Oliphant

Permission to use, modify, and distribute this software is given under the
terms of the SciPy (BSD style) license.  See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.

Copyright (c) 2006-2007, Robert Hetland <hetland@tamu.edu>
Copyright (c) 2007, John Travers <jtravs@gmail.com>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

    * Neither the name of Robert Hetland nor the names of any
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from __future__ import division, print_function, absolute_import

from itertools import chain, combinations_with_replacement as combinations_w_r
import sys
import numpy as np

from scipy import linalg, sparse
from scipy._lib.six import callable, get_method_function, get_function_code
from scipy.misc import comb
from scipy.special import xlogy, jn

__all__ = ['Rbf']


def combinations(n_features, degree):
    # Same functionality as sklearn.preprocessing.PolynomialFeatures
    return chain.from_iterable(combinations_w_r(range(n_features), i)
                               for i in range(degree + 1))


class Rbf(object):
    """
    Rbf(*args)

    A class for radial basis function approximation/interpolation of
    n-dimensional scattered data.

    For better interpolation accuracy, we generate an augmented matrix with
    polynomial coefficients as described in Section 2.1 of:
    http://www.sciencedirect.com/science/article/pii/S0021999116301632

    Parameters
    ----------
    *args : arrays
        x, y, z, ..., d, where x, y, z, ... are the coordinates of the nodes
        and d is the array of values at the nodes
    function : str or callable, optional
        The radial basis function, based on the radius, r, given by the norm
        (default is Euclidean distance); the default is 'multiquadric'::

            'multiquadric': sqrt((r/self.epsilon)**2 + 1)
            'inverse': 1.0/sqrt((r/self.epsilon)**2 + 1)
            'gaussian': exp(-(r/self.epsilon)**2)
            'linear': r
            'cubic': r**3
            'quintic': r**5
            'thin_plate': r**2 * log(r)
            'bessel': jn(self.p/2-1, 1.0/self.epsilon*r)/(arg**(self.p/2-1))
            'polyharmonic': r**(2*self.m + 1)

        If callable, then it must take 2 arguments (self, r).  The epsilon
        parameter will be available as self.epsilon.  Other keyword
        arguments passed in will be available as well.
    degree : integer, optional
        Parameter denoting the maximum polynomial degree to use for
        generating additional terms from the data to use for interpolation.
        The higher the degree, the more accurate the interpolation
        coefficients.
        - defaults to None, which performs the default RBF interpolation
        without any additional generated terms.
    epsilon : float, optional
        Parameter used for interpolation. Only used if function is 'gaussian',
        'bessel', 'multiquadric', or 'inverse'.
        - defaults to approximate average distance between nodes (which is
        a good start).
    m : float, optional
        Parameter used for interpolation. Only used if function is 'polyharmonic'.
        Denotes the power of the polyharmonic function.
    p : float, optional
        Parameter used for interpolation. Only used if function is 'bessel'.
        The order of the Bessel function to be used is p/2-1.
    smooth : float, optional
        Values greater than zero increase the smoothness of the
        approximation.  0 is for interpolation (default), the function will
        always go through the nodal points in this case.
    norm : callable, optional
        A function that returns the 'distance' between two points, with
        inputs as arrays of positions (x, y, z, ...), and an output as an
        array of distance.  E.g, the default::

            def euclidean_norm(x1, x2):
                return sqrt( ((x1 - x2)**2).sum(axis=0) )

        which is called with ``x1 = x1[ndims, newaxis, :]`` and
        ``x2 = x2[ndims, : ,newaxis]`` such that the result is a matrix of the
        distances from each point in ``x1`` to each point in ``x2``.

    Examples
    --------
    >>> from scipy.interpolate import Rbf
    >>> x, y, z, d = np.random.rand(4, 50)
    >>> rbfi = Rbf(x, y, z, d)  # radial basis function interpolator instance
    >>> xi = yi = zi = np.linspace(0, 1, 20)
    >>> di = rbfi(xi, yi, zi)   # interpolated values
    >>> di.shape
    (20,)

    """

    def _euclidean_norm(self, x1, x2):
        return np.sqrt(((x1 - x2)**2).sum(axis=0))

    def _h_multiquadric(self, r):
        return np.sqrt((1.0/self.epsilon*r)**2 + 1)

    def _h_inverse_multiquadric(self, r):
        return 1.0/np.sqrt((1.0/self.epsilon*r)**2 + 1)

    def _h_gaussian(self, r):
        return np.exp(-(1.0/self.epsilon*r)**2)

    def _h_bessel(self, r):
        arg = 1.0/self.epsilon*r
        return jn(self.p/2-1, arg)/(arg**(self.p/2-1))

    def _h_inverse_quadratic(self, r):
        return 1.0/(1.0 + (1.0/self.epsilon*r)**2)

    def _h_polyharmonic(self, r):
        return r**(2*self.m+1)

    def _h_log_polyharmonic(self, r):
        return xlogy(r**self.m, r)

    def _h_linear(self, r):
        return r

    def _h_cubic(self, r):
        return r**3

    def _h_quintic(self, r):
        return r**5

    def _h_thin_plate(self, r):
        return xlogy(r**2, r)

    # Setup self._function and do smoke test on initial r
    def _init_function(self, r):
        if isinstance(self.function, str):
            self.function = self.function.lower()
            _mapped = {'inverse': 'inverse_multiquadric',
                       'inverse multiquadric': 'inverse_multiquadric',
                       'thin-plate': 'thin_plate'}
            if self.function in _mapped:
                self.function = _mapped[self.function]

            func_name = "_h_" + self.function
            if hasattr(self, func_name):
                self._function = getattr(self, func_name)
            else:
                functionlist = [x[3:] for x in dir(self) if x.startswith('_h_')]
                raise ValueError("function must be a callable or one of " +
                                 ", ".join(functionlist))
            self._function = getattr(self, "_h_"+self.function)
        elif callable(self.function):
            allow_one = False
            if hasattr(self.function, 'func_code') or \
                    hasattr(self.function, '__code__'):
                val = self.function
                allow_one = True
            elif hasattr(self.function, "im_func"):
                val = get_method_function(self.function)
            elif hasattr(self.function, "__call__"):
                val = get_method_function(self.function.__call__)
            else:
                raise ValueError("Cannot determine number of arguments "
                                 "to function")

            argcount = get_function_code(val).co_argcount
            if allow_one and argcount == 1:
                self._function = self.function
            elif argcount == 2:
                if sys.version_info[0] >= 3:
                    self._function = self.function.__get__(self, Rbf)
                else:
                    import new
                    self._function = new.instancemethod(self.function, self,
                                                        Rbf)
            else:
                raise ValueError("Function argument must take 1 or 2 "
                                 "arguments.")

        a0 = self._function(r)
        if a0.shape != r.shape:
            raise ValueError("Callable must take array and return array "
                             "of the same shape")
        return a0

    def __init__(self, *args, **kwargs):
        self.xi = np.asarray([np.asarray(a, dtype=np.float64).flatten()
                             for a in args[:-1]])
        self.N = self.xi.shape[-1]
        self.di = np.asarray(args[-1]).flatten()

        if not all([x.size == self.di.size for x in self.xi]):
            raise ValueError("All arrays must be equal length.")

        self.norm = kwargs.pop('norm', self._euclidean_norm)
        self.degree = kwargs.pop('degree', None)
        self.epsilon = kwargs.pop('epsilon', None)
        self.m = kwargs.pop('m', None)
        self.p = kwargs.pop('p', None)
        if self.epsilon is None:
            # default epsilon is the "the average distance between nodes" based
            # on a bounding hypercube
            ximax = np.amax(self.xi, axis=1)
            ximin = np.amin(self.xi, axis=1)
            edges = ximax-ximin
            edges = edges[np.nonzero(edges)]
            self.epsilon = np.power(np.prod(edges)/self.N, 1.0/edges.size)
        self.smooth = kwargs.pop('smooth', 0.0)

        self.function = kwargs.pop('function', 'multiquadric')

        # attach anything left in kwargs to self
        #  for use by any user-callable function or
        #  to save on the object returned.
        for item, value in kwargs.items():
            setattr(self, item, value)

        self._compute_nodes()

    @property
    def A(self):
        # this only exists for backwards compatibility: self.A was available
        # and, at least technically, public.
        r = self._call_norm(self.xi, self.xi)
        return self._init_function(r) - np.eye(self.N)*self.smooth

    def _call_norm(self, x1, x2):
        if len(x1.shape) == 1:
            x1 = x1[np.newaxis, :]
        if len(x2.shape) == 1:
            x2 = x2[np.newaxis, :]
        x1 = x1[..., :, np.newaxis]
        x2 = x2[..., np.newaxis, :]
        return self.norm(x1, x2)

    def __call__(self, *args):
        """
        Instances of this class are callable with points passed in.
        We generate the upper part of the augmented matrix with the
        arguments and existing data points, and solve this the same way
        as in self._compute_nodes().
        """
        args = [np.asarray(x) for x in args]
        if not all([x.shape == y.shape for x in args for y in args]):
            raise ValueError("Array lengths must be equal")
        shp = args[0].shape
        xa = np.asarray([a.flatten() for a in args], dtype=np.float64)
        r = self._call_norm(xa, self.xi)
        if self.degree:
            n_points, n_dims = shp[0] if len(shp) else 1, len(args)
            POLY = np.empty((n_points, self.n_coeffs), dtype=np.float64)
            combs = combinations(n_dims, self.degree)
            for i, c in enumerate(combs):
                POLY[:, i] = xa[c, :].prod(axis=0)
            aug = np.hstack((self._function(r), POLY))
            return np.dot(aug, self.nodes).reshape(shp)
        else:
            return np.dot(self._function(r), self.nodes).reshape(shp)

    def _compute_nodes(self):
        """
        First, we generate all possible combinations of different powers
        from each input dimension such that the sum of the powers of each
        term is equal to self.degree. Denote this matrix as C.
        C is stacked with the A matrix to create the augmented matrix as
        shown below.

             ---------
             | A | C |
        M =  ---------
             |C.T| 0 |
             ---------

        where C is a matrix of size (n, n_coeffs) generated from the
        A matrix.

        We then solve the equation Mx = self.di to find x.
        """
        n_dims, n_points = self.xi.shape
        if self.degree:
            # closed-form formula for number of polynomial coefficients
            self.n_coeffs = int(comb(self.degree + n_dims, n_dims))

            # allocate output data
            POLY = np.empty((n_points, self.n_coeffs), dtype=np.float64)

            combs = combinations(n_dims, self.degree)
            for i, c in enumerate(combs):
                POLY[:, i] = self.xi[c, :].prod(axis=0)
            A_aug = sparse.bmat([[self.A, POLY], [POLY.T, None]], format='csr')
            di_aug = np.hstack((self.di, np.zeros(self.n_coeffs)))
            self.nodes = sparse.linalg.spsolve(A_aug, di_aug)
        else:
            self.nodes = linalg.solve(self.A, self.di)
