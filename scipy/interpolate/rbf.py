# coding: utf-8
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
met::

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

import sys
import numpy as np

from scipy import linalg
from scipy._lib.six import callable, get_method_function, get_function_code
from scipy.special import xlogy

__all__ = ['Rbf']


class Rbf(object):
    r"""
    Rbf(*args)

    A class for radial basis function approximation/interpolation of
    `n`-dimensional scattered data.

    Parameters
    ----------
    *args : arrays
        `x`, `y`, `z`, ..., `f`, where `x`, `y`, `z`, ... are the vectors
        of the coordinates of the nodes and `f` is the array of values at
        the nodes
    function : {'multiquadric', 'inverse_multiquadric', 'gaussian', 'thin_plate'} or callable
        The radial basis function (see below for details).
        If callable, then it must take 2 arguments `(self, r)`.  The `epsilon`
        parameter will be available as `self.epsilon`.  Other keyword
        arguments passed in will be available as well.
    epsilon : float, optional
        Adjustable constant for gaussian or multiquadrics functions ---
        defaults to approximate average distance between adjacent nodes
        (which is a reasonable starting value).
    smooth : float, optional
        Values greater than zero increase the smoothness of the
        approximation.  The default value of zero gives interpolation, i.e.,
        the function will always go through the nodal points.
    poly_order : 1, 0 or None, optional
        The order of a (low-order) polynomial to be fitted and added to the
        input data; the default of 1 corresponds to a linear term, 0 to a
        constant, `None` for no polynomial part.
    norm : callable, optional
        A function that returns the 'distance' between two points, with an
        input of arrays of positions (`x`, `y`, `z`, ...), and an output of
        an array of distances between those points.

    Returns
    -------
    s : callable
        A function which performs the interpolation.

    Raises
    ------
    ValueError
        In the case that the vector arguments `x`, `y`, `z`,..., `f` are
        not of the same length,

    Notes
    -----
    **Interpolation**

    Let :math:`f` be a function which has been sampled at distinct
    points :math:`\mathbf{x}_1, \mathbf{x}_2, \ldots`,

    .. math::

        f_i = f(\mathbf{x}_1, \mathbf{x}_2, \ldots)

    and we wish to create an interpolant (or approximant) of :math:`f`
    from the data :math:`f_i`, :math:`\mathbf{x}_i`.  If :math:`\phi`
    is a non-negative function defined on :math:`[0, \infty)` then

    .. math::

        \phi_i(\mathbf{x}) = \
        \phi\left(\left\| \mathbf{x} - \mathbf{x}_i \right\| \right)

    is a radially-symmetric function centred on :math:`\mathbf{x}_i`.  A
    *radial basis function interpolant* :math:`s` for :math:`f` is a linear
    combination of such :math:`\phi_i`, i.e.,

    .. math::

        s(\mathbf{x}) = \sum_{i=1}^{n} \alpha_i \phi_i(\mathbf{x})

    where the coefficients :math:`\alpha_i` are determined by the requirement
    that :math:`s` be an interpolant:

    .. math::

        f_j = s(\mathbf{x}_j) = \sum_{i=1}^{n} \alpha_i \phi_i(\mathbf{x}_j)

    which is a linear system which can be written in matrix form

    .. math::
        :label: interp

        \mathbf{f} = \Phi \mathbf{a}

    where :math:`\mathbf{f}` is the vector of :math:`f_i`,
    :math:`\mathbf{a}` is the vector of :math:`\alpha_i`, and :math:`\Phi` is
    the matrix with entries
    :math:`\phi\left(\left\| \mathbf{x}_j - \mathbf{x}_i \right\| \right)`:
    one solves the linear system :eq:`interp` to find the :math:`\alpha_i`.

    This interpolation method was introduced by Hardy [2]_ for applications
    in geophysics, and has since become one of the most popular methods for
    multidimensional scattered data interpolation.  The texts by Buhmann [1]_
    or Wendland [4]_ should be consulted for a detailed treatment of the
    theory.

    **Approximation**

    In the case that the :math:`f_i` are contaminated by noise, one
    typically seeks an *approximation* rather than an interpolation, and
    this can be achieved by penalising the size of the function
    :math:`s` by finding

    .. math::

        \min\left\{ \
        \sum_{j=1}^{N}\left| s(\mathbf{x}_j) - f_j \right| + \
        \lambda\left\| s \right\|
        \right\}

    where :math:`\lambda` is a *smoothing parameter* which balances the
    smoothness and accuracy of the approximation.  It can be shown that
    the minimiser of this problem is the solution to the linear
    system

    .. math::
        :label: approx

        \mathbf{f} = (\Phi + \lambda I)\mathbf{a}

    where :math:`I` is the identity matrix of the same order as
    :math:`\Phi`. Clearly the system :eq:`approx` reduces to
    the interpolation system :eq:`interp` as
    :math:`\lambda\rightarrow 0`.

    See Wendland's paper [3]_ for a fuller discussion and references
    on the practicalites of radial basis function approximation.

    **The polynomial part**

    The quality of the interpolation (or approximation) is often improved
    by combining a `low-order polynomial` interpolant with that of the
    radial basis functions; in other words one seeks an interpolant
    :math:`s` of the form

    .. math::

        s(\mathbf{x}) = \
        \sum_{i=1}^{n} \alpha_i \phi_i(\mathbf{x}) + p(\mathbf{x}).

    This module implements only linear (order 1) and constant (order 0)
    polynomials.

    **The radial basis functions implemented**

    The `multiquadric` (conditionally negative definite):

    .. math::

        \phi(r) = \sqrt{(r/\epsilon)^2 + 1},

    `inverse_multiquadric`, alias `inverse` (positive definte):

    .. math::

        \phi(r) = \frac{1}{\sqrt{(r/\epsilon)^2 + 1}},

    `gaussian` (positive definite):

    .. math::

        \phi(r) = \exp\left(-(r/\epsilon)^2\right),

    and the `thin_plate` (conditionally positive definite)

    .. math::

        \phi(r) = r^2 \log(r).

    (The deprecated functions `linear`, `cubic` and `quintic` are, fairly
    obviously, :math:`r`, :math:`r^3` and :math:`r^5` respecively.)

    It is common to find variations in these definitions in the literature
    (:math:`\epsilon` replaced by :math:`1/\epsilon`, for example).

    References
    ----------
    .. [1] Martin D. Buhmann,
       "Radial Basis Functions: Theory and Implementations",
       Cambridge University Press, 2003.

    .. [2] R. L. Hardy,
       "Multiquadric equations of topography and other irregular surfaces",
       Journal of Geophysical Research, 76(8):1905--1915, 1971.

    .. [3] Holger Wendland,
       "Computational Aspects of Radial Basis Function Approximation",
       in K. Jetter *et al.* (eds.)
       Topics in Multivariate Approximation and Interpolation,
       Elsevier B.V., 2005, 231--256.

    .. [4] Holger Wendland,
       "Scattered Data Approximation"
       Cambridge University Press, 2004.

    Examples
    --------
    For a set of `n` points in 3-dimensional space with coordinates
    in the `n`-vectors `x`, `y` and `z`; and with `f` being a
    `n`-vector of the data from which to interpolate, the interpolant
    `s` is created with

    >>> from scipy.interpolate import Rbf
    >>> x, y, z, f = np.random.rand(4, 50)
    >>> s = Rbf(x, y, z, f)

    One can then evaluate the interpolant at arbitrary points
    `xi`, `yi`, `zi` with

    >>> xi = yi = zi = np.linspace(0, 1, 20)
    >>> si = s(xi, yi, zi)
    >>> si.shape
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

    def _h_thin_plate(self, r):
        return xlogy(r**2, r)

    @np.deprecate(message="not useful as an RBF")
    def _h_linear(self, r):
        return r

    @np.deprecate(message="not useful as an RBF")
    def _h_cubic(self, r):
        return r**3

    @np.deprecate(message="not useful as an RBF")
    def _h_quintic(self, r):
        return r**5

    def _init_function(self):
        if isinstance(self.function, str):
            self.function = self.function.lower()
            _mapped = {
                'inverse': 'inverse_multiquadric',
                'inverse multiquadric': 'inverse_multiquadric',
                'thin-plate': 'thin_plate'
            }
            if self.function in _mapped:
                self.function = _mapped[self.function]

            func_name = "_h_" + self.function
            if hasattr(self, func_name):
                self._function = getattr(self, func_name)
            else:
                functionlist = [x[3:] for x in dir(self) if x.startswith('_h_')]
                raise ValueError("function must be a callable or one of " +
                                 ", ".join(functionlist))
            self._function = getattr(self, "_h_" + self.function)
        elif callable(self.function):
            allow_one = False
            if hasattr(self.function, 'func_code') or hasattr(self.function, '__code__'):
                val = self.function
                allow_one = True
            elif hasattr(self.function, "im_func"):
                val = get_method_function(self.function)
            elif hasattr(self.function, "__call__"):
                val = get_method_function(self.function.__call__)
            else:
                raise ValueError("Cannot determine number of arguments to function")

            argcount = get_function_code(val).co_argcount
            if allow_one and argcount == 1:
                self._function = self.function
            elif argcount == 2:
                if sys.version_info[0] >= 3:
                    self._function = self.function.__get__(self, Rbf)
                else:
                    import new
                    self._function = new.instancemethod(self.function, self, Rbf)
            else:
                raise ValueError("Function argument must take 1 or 2 arguments")

    def _init_default_epsilon(self):
        # default epsilon is the "the average distance between adjacent nodes"
        ximax = np.amax(self.xi, axis=1)
        ximin = np.amin(self.xi, axis=1)
        self.epsilon = np.power(np.prod(ximax - ximin)/self.N, 1.0/self.dim)

    def _interpolation_matrix(self):
        R = self._call_norm(self.xi, self.xi)
        A = self._function(R)
        if R.shape != A.shape:
            raise ValueError("Callable must take array and return array of the same shape")
        if self.smooth != 0.0:
            eigidx = 0 if self.function == 'multiquadric' else self.N - 1
            eig = linalg.eigh(A, eigvals_only=True, eigvals=(eigidx, eigidx))[0]
            A += np.eye(self.N) * self.smooth * eig
        return A

    def __init__(self, *args, **kwargs):
        # the data points and the value of the function to be interpolated;
        self.xi = np.array([np.asarray(a, dtype=np.float_).flatten() for a in args[:-1]])
        self.fi = np.asarray(args[-1]).flatten()

        if not all([x.size == self.fi.size for x in self.xi]):
            raise ValueError("All arrays must be equal length")

        # the dimension of the interpolation space, and the number of
        # interpolation samples
        self.dim = self.xi.shape[0]
        self.N = self.xi.shape[-1]

        # the norm
        self.norm = kwargs.pop('norm', self._euclidean_norm)

        # the shape parameter for the RBF
        self.epsilon = kwargs.pop('epsilon', None)
        if self.epsilon is None:
            self._init_default_epsilon()

        # the smoothing parameter, if zero then the RBF will interpolate,
        # if non-zero then it will approximate.  For the standard (non-
        # deprecated) RBFs, the value should be positve.
        self.smooth = kwargs.pop('smooth', 0.0)

        # order of polynomial term ('None' for no polynomial term)
        self.poly_order = kwargs.pop('poly_order', 1)
        if self.poly_order not in set([None, 0, 1]):
            raise ValueError("poly_order should be None, 0 or 1")

        # the radial basis function, and then attach anything left in kwargs
        # to self for use by any user-callable function
        self.function = kwargs.pop('function', 'multiquadric')
        for item, value in kwargs.items():
            setattr(self, item, value)

        self._init_function()

        # the interpolation (without the polynomial term)
        A = self._interpolation_matrix()

        def solve_system(A, f):
            kwargs = {
                'overwrite_a': True,
                'overwrite_b': False,
                'check_finite': False
            }
            coef, _, _, _ = linalg.lstsq(A, f, **kwargs)
            return coef

        if self.poly_order is None:

            self.nodes = solve_system(A, self.fi)
            self.poly = None

        else:

            # construct and augment the polynomial coefficients
            P = np.ones((1, self.N))
            if self.poly_order >= 1:
                P = np.vstack((P, self.xi))
            k = P.shape[0]
            Z = np.zeros((k, k))
            A = np.vstack((np.hstack((A, P.T)),
                           np.hstack((P, Z))))
            f = np.hstack((self.fi, np.zeros((k,))))

            # solve
            coefs = solve_system(A, f)

            # save node coefficients
            self.nodes = coefs[:self.N]

            # save polynomial coefficients
            self.poly = (coefs[self.N],)
            if self.poly_order >= 1:
                self.poly += (coefs[self.N+1:],)

    def _call_norm(self, x1, x2):
        if len(x1.shape) == 1:
            x1 = x1[np.newaxis, :]
        if len(x2.shape) == 1:
            x2 = x2[np.newaxis, :]
        x1 = x1[..., :, np.newaxis]
        x2 = x2[..., np.newaxis, :]
        return self.norm(x1, x2)

    def __call__(self, *args):
        args = [np.asarray(x) for x in args]
        if len(args) == 0:
            raise ValueError("Need at least one argument")
        shp = args[0].shape
        if not all([arg.shape == shp for arg in args]):
            raise ValueError("Array lengths must be equal")

        # points at which to interpolate
        x = np.asarray([a.flatten() for a in args], dtype=np.float_)

        # the polynomial part
        poly = 0
        if self.poly_order >= 0:
            poly += self.poly[0]
            if self.poly_order >= 1:
                poly += np.dot(self.poly[1], x)

        # the RBF part
        r = self._call_norm(x, self.xi)

        return (np.dot(self._function(r), self.nodes) + poly).reshape(shp)
