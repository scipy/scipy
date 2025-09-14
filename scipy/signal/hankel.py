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


def hankel(func, B, order=[0, 1], 
           method='Anderson801', ybase=None, wt0=None, wt1=None):
    """Zeroth- and first-order Hankel transforms.

    This routine implements zeroth- and first-order Hankel transforms through
    digital filtering [1]_.  No attempts are made at being efficient here, so
    by default each avaluation of the transformed function (each element of
    `B`) requires 801 evaluations of the input function.

    Parameters
    ----------
    func : callable
        The function to transform.  Must deal with vector arguments.
    B    : array_like
        point(s) at which to transform
    order : {0, 1, [0, 1]}, optional
        order or list of orders to evaluate
    method : string, optional
        A name which identifies the filter coefficients used.  Currently, the
        default 'Anderson801' is the only name recognized.  To specify other
        filters, use the `ybase`, `wt0` and `wt1` keywords.
    ybase, wt0, wt1 : array_like, optional
        The scaling factor and weights (for zeroth- and first-order
        transforms), respectively.  For use when transforming with other
        filters than the ones provided with SciPy.  If `ybase` is provided, 
        `wt0` and/or `wt1` must also be.

    Returns
    -------
    list
        a list of Hankel transforms of `func`, evaluated at points `B`,
        conforming to the order(s) given as inputs. The list will have one or
        two elements, each of which is the same size as the `B` parameter.

    Notes
    -----

    The default transformation filters are the 801-point filters from [1]_.
    Other filters have been described in the literature, e.g. the 61- 121- or
    241-point filters due to Fan-Nian Kong [2]_.

    .. [1] Anderson, W. L., "Computer Program Numerical Integration of Related
       Hankel Transforms of Orders 0 and 1 by Adaptive Digital Filtering,"
       Geophysics, vol. 44(7), pp. 1287-1305, 1979.

    .. [2] Kong, F. N., "Hankel transform filters for dipole antenna radiation
       in a conductive medium," Geophysical Prospecting, vol. 55(1), pp. 83-89,
       doi:10.1111/j.1365-2478.2006.00585.x, 2007.

    Examples
    --------

    The zeroth-order Hankel transform of `g*np.exp(-g**2)` is `np.exp(-b**2/4)/2`

    >>> B = np.r_[0.0001, 0.001, 0.01:0.05:0.01,  0.05:5.01:0.05]
    >>> f = lambda g: g*np.exp(-g**2)
    >>> hf = hankel0(f, B)
    >>> ht = lambda b: np.exp(-b**2/4)/2
    >>> if np.allclose(hf, ht(B), atol=1e-7): print "success!"

    """

    if not ybase:
        # find specific filters by name
        if method.lower() == 'anderson801':
            import _anderson801 as A
            ybase, wt0, wt1 = A.YBASE, A.WT0, A.WT1
        ## add more named filters here
        # elif method.lower() == 'kong61':
        #     import _kong61 as K
        #     ybase, wt0, wt1 = K.YBASE, K.WT0, K.WT1
        else:
            raise ValueError('Hankel Transform of type ' + method + \
                    ' not known')

    B = np.asarray(B)

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
        fy = func(ybase/bval)
        if 0 in order:
            HB0[i] = np.dot(wt0, fy) / bval
        if 1 in order:
            HB1[i] = np.dot(wt1, fy) / bval

    return retval

def hankel0(func, B, method='Anderson801', ybase=None, wt0=None, wt1=None):
    """Zeroth-order Hankel transform.

    This is a convenience wrapper around `hankel`, which computes and returns
    only the zeroth-order transform.

    Returns
    -------
    array_like
        The transform of `func`, evaluated at points `B`.

    """
    return hankel(func, B, [0], method, ybase, wt0, wt1)[0]

def hankel1(self, func, B):
    """First-order Hankel transform.

    This is a convenience wrapper around `hankel`, which computes and returns
    only the first-order transform.

    Returns
    -------
    array_like
        The transform of `func`, evaluated at points `B`.
    """
    return self.hankel(func, B, [1], method, ybase, wt0, wt1)[0]

