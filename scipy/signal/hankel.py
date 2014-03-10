"""
hankel -- A numerical Hankel transform using the digital filter method.
"""

# Hankel transform for scipy
# by Tom Grydeland (tom.grydeland@norut.no)
# No copyright notices on the original code, nor on this

# The Python routines in this directory implement Hankel transforms of order
# 0 and 1 by digital filtering.  These are adaptations of a MATLAB code by
# Prof. Brian Borchers at New Mexico Tech, which again are based on a
# FORTRAN program by Walt Anderson

# The program is independent of the actual evaluation points and filter
# coefficients used, and the filters are provided separately.

import numpy as np


class HankelTransform(object):
    """An object implementing Hankel transformations

    Parameters
    ----------
    args : arguments
        The 'HankelTransform' class can be instantiated with a string argument
        (defaulting to 'Anderson801', or with three vectors.  A string must
        indicate which filter coefficients to use, while the three vectors
        gives the user the possibility of providing filter coefficients from
        other sources, e.g. Fan-Nian Kong's 61- 121- or 241-point filters[1].

        When three vectors are given, these are `ybase`, `wt0` and `wt1`, respectively.

        [1] Kong, F. N. (2007). Hankel transform filters for dipole antenna
        radiation in a conductive medium. Geophysical Prospecting, 55(1),
        83–89. doi:10.1111/j.1365-2478.2006.00585.x
    """
    def __init__(self, ybase='Anderson801', wt0=None, wt1=None):

        # Asking for a specific implementation by name
        if isinstance(ybase, str):
            if ybase.lower() == 'anderson801':
                import anderson801 as A
                self.__init__(A.YBASE, A.WT0, A.WT1)
                return

            raise ValueError('Hankel Transform of type ' + ybase + \
                    ' not known')

        self.ybase = ybase
        self.wt0   = wt0
        self.wt1   = wt1

    def hankel(self, func, B, order=[0]):
        """
        zeroth- and first-order Hankel transforms.

        Parameters
        ----------
        func : function to transform.  Must deal with vector arguments
        B    : point(s) at which to transform
        order : list of orders to evaluate (i.e. 0, 1, or [0, 1])

        Returns a list of Hankel transforms of `func`, evaluated at points `B`,
        conforming to the order(s) given as inputs.
        """

        try:
            iter(B)
        except TypeError:
            B = [B]

        try:
            iter(order)
        except TypeError:
            order = [order]

        retval = []
        if 0 in order:
            HB0 = np.zeros(B.shape)
            retval.append(HB0)
        if 1 in order:
            HB1 = np.zeros(B.shape)
            retval.append(HB1)

        for i in range(len(B)):
            bval = B[i]
            fy = func(self.ybase/bval)
            if 0 in order:
                HB0[i] = np.dot(self.wt0, fy) / bval
            if 1 in order:
                HB1[i] = np.dot(self.wt1, fy) / bval

        return retval

    def hankel0(self, func, B):
        """
        zeroth-order Hankel transform.

        Parameters
        ----------
        func : function to transform.  Must deal with vector arguments
        B    : point(s) at which to transform

        Returns
        -------
        z0 : zeroth-order Hankel transforms of func,
        evaluated at points B.
        """
        return self.hankel(func, B, [0])[0]

    def hankel1(self, func, B):
        """
            first-order Hankel transform.

        Parameters
        ----------
        func : function to transform.  Must deal with vector arguments
        B    : point(s) at which to transform

        Returns
        -------
        z1 : first-order Hankel transforms of func,
        evaluated at points B.
        """
        return self.hankel(func, B, [1])[0]

    def hankel01(self, func, B):
        """
        zeroth- and first-order Hankel transforms.

        Parameters
        ----------
        func : function to transform.  Must deal with vector arguments
        B    : point(s) at which to transform
        order : list of orders to evaluate (i.e. 0, 1, or [0, 1])

        Returns
        -------
        [z0, z1] : zeroth- and first-order Hankel transforms of func,
        evaluated at points B.
        """
        return self.hankel(func, B, [0, 1])

