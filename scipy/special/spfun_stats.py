#! /usr/bin/env python
# Last Change: Sat Mar 21 02:00 PM 2009 J

# Copyright (c) 2001, 2002 Enthought, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   a. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#   b. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#   c. Neither the name of the Enthought nor the names of its contributors
#      may be used to endorse or promote products derived from this software
#      without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

"""Some more special functions which may be useful for multivariate statistical
analysis."""

import numpy as np
from scipy.special import gammaln as loggam

__all__ = ['multigammln']


def multigammaln(a, d):
    """returns the log of multivariate gamma, also sometimes called the
    generalized gamma.

    Parameters
    ----------
    a : ndarray
        the multivariate gamma is computed for each item of a
    d : int
        the dimension of the space of integration.

    Returns
    -------
    res : ndarray
        the values of the log multivariate gamma at the given points a.

    Note
    ----
    The formal definition of the multivariate gamma of dimension d for a real a
    is :

    \Gamma_d(a) = \int_{A>0}{e^{-tr(A)\cdot{|A|}^{a - (m+1)/2}dA}}

    with the condition a > (d-1)/2, and A>0 being the set of all the positive
    definite matrices of dimension s. Note that a is a scalar: the integrand
    only is multivariate, the argument is not (the function is defined over a
    subset of the real set).

    This can be proven to be equal to the much friendler equation:

    \Gamma_d(a) = \pi^{d(d-1)/4}\prod_{i=1}^{d}{\Gamma(a - (i-1)/2)}.

    Notes
    -----
    Reference:

    R. J. Muirhead, Aspects of multivariate statistical theory (Wiley Series in
    probability and mathematical statistics). """
    a = np.asarray(a)
    if not np.isscalar(d) or (np.floor(d) != d):
        raise ValueError("d should be a positive integer (dimension)")
    if np.any(a <= 0.5 * (d - 1)):
        raise ValueError("condition a (%f) > 0.5 * (d-1) (%f) not met" \
                         % (a, 0.5 * (d-1)))

    res = (d * (d-1) * 0.25) * np.log(np.pi)
    if a.size == 1:
        axis = -1
    else:
        axis = 0
    res += np.sum(loggam([(a - (j - 1.)/2) for j in range(1, d+1)]), axis)
    return res
