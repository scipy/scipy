"""rbf - Radial basis functions for interpolation/smoothing scattered Nd data.

Written by John Travers <jtravs@gmail.com>, February 2007
Based closely on Matlab code by Alex Chirokov
Additional, large, improvements by Robert Hetland

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

from numpy import sqrt, log, asarray, newaxis, all, dot, float64, exp, eye
from scipy import linalg

class Rbf(object):
    """ A class for radial basis function approximation/interpolation of
        n-dimensional scattered data.
    """

    def _euclidean_norm(self, x1, x2):
        return sqrt( ((x1 - x2)**2).sum(axis=0) )

    def _function(self, r):
        if self.function.lower() == 'multiquadric':
            return sqrt((1.0/self.epsilon*r)**2 + 1)
        elif self.function.lower() == 'inverse multiquadric':
            return 1.0/sqrt((1.0/self.epsilon*r)**2 + 1)
        elif self.function.lower() == 'gaussian':
            return exp(-(self.epsilon*r)**2)
        elif self.function.lower() == 'cubic':
            return r**3
        elif self.function.lower() == 'quintic':
            return r**5
        elif self.function.lower() == 'thin-plate':
            return r**2 * log(r)
        else:
            raise ValueError, 'Invalid basis function name'

    def __init__(self, *args, **kwargs):
        """ Constructor for Rbf class.

        Parameters
        ----------
        *args : arrays
            x, y, z, ..., d, where x, y, z, ... are the coordinates of the nodes
            and d is the array of values at the nodes
        function : str, optional
            The radial basis function, based on the radius, r, given by the norm
            (defult is Euclidean distance); the default is 'multiquadratic'.

            ::
                'multiquadric': sqrt((self.epsilon*r)**2 + 1)
                'inverse multiquadric': 1.0/sqrt((self.epsilon*r)**2 + 1)
                'gaussian': exp(-(self.epsilon*r)**2)
                'cubic': r**3
                'quintic': r**5
                'thin-plate': r**2 * log(r)
        epsilon : float, optional
            Adjustable constant for gaussian or multiquadrics functions
            - defaults to approximate average distance between nodes (which is
            a good start).
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

            which is called with x1=x1[ndims,newaxis,:] and
            x2=x2[ndims,:,newaxis] such that the result is a symetric, square
            matrix of the distances between each point to each other point.

        Returns
        -------
        rbf : Rbf
            Interpolator object that returns interpolated values at new positions.

        Examples
        --------
        >>> rbfi = Rbf(x, y, z, d)      # radial basis function interpolator instance
        >>> di = rbfi(xi, yi, zi)       # interpolated values
        """
        self.xi = asarray([asarray(a, dtype=float64).flatten() for a in args[:-1]])
        self.N = self.xi.shape[-1]
        self.di = asarray(args[-1], dtype=float64).flatten()

        assert [x.size==self.di.size for x in self.xi], \
               'All arrays must be equal length'

        self.norm = kwargs.pop('norm', self._euclidean_norm)
        r = self._call_norm(self.xi, self.xi)
        self.epsilon = kwargs.pop('epsilon', r.mean())
        self.function = kwargs.pop('function', 'multiquadric')
        self.smooth = kwargs.pop('smooth', 0.0)

        self.A = self._function(r) - eye(self.N)*self.smooth
        self.nodes = linalg.solve(self.A, self.di)

    def _call_norm(self, x1, x2):
        if len(x1.shape) == 1:
            x1 = x1[newaxis, :]
        if len(x2.shape) == 1:
            x2 = x2[newaxis, :]
        x1 = x1[..., :, newaxis]
        x2 = x2[..., newaxis, :]
        return self.norm(x1, x2)

    def __call__(self, *args):
        assert all([x.shape == y.shape \
                    for x in args \
                    for y in args]), 'Array lengths must be equal'
        shp = args[0].shape
        self.xa = asarray([a.flatten() for a in args], dtype=float64)
        r = self._call_norm(self.xa, self.xi)
        return dot(self._function(r), self.nodes).reshape(shp)
