
# Adding documentation to the logit and expit ufuncs

from _logit import logit, expit

try:
    from numpy.lib import add_newdoc_ufunc

    add_newdoc_ufunc(logit,
        """
        Logit ufunc for ndarrays.

        The logit function is defined as logit(p) = log(p/(1-p)).
        Note that logit(0) = -inf, logit(1) = inf, and logit(p)
        for p<0 or p>1 yields nan.

        Parameters
        ----------
        x : ndarray
            The ndarray to apply logit to element-wise.

        Returns
        -------
        out : ndarray
            An ndarray of the same shape as x. Its entries
            are logit of the corresponding entry of x.

        Notes
        -----
        As a ufunc logit takes a number of optional
        keywork arguments. For more information
        see `ufuncs <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_
        """)

    add_newdoc_ufunc(expit,
        """
        Expit ufunc for ndarrays.

        The expit function is defined as expit(x) = 1/(1+exp(-x)).
        Note that expit is the inverse logit function.

        Parameters
        ----------
        x : ndarray
            The ndarray to apply expit to element-wise.

        Returns
        -------
        out : ndarray
            An ndarray of the same shape as x. Its entries
            are expit of the corresponding entry of x.

        Notes
        -----
        As a ufunc logit takes a number of optional
        keywork arguments. For more information
        see `ufuncs <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_
        """)

except ImportError:
    pass






