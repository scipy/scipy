import operator
import numpy as np

from scipy._lib._util import prod

from . import _bspl

__all__ = ["NdBSpline"]


def _get_dtype(dtype):
    """Return np.complex128 for complex dtypes, np.float64 otherwise."""
    if np.issubdtype(dtype, np.complexfloating):
        return np.complex_
    else:
        return np.float_


class NdBSpline:
    """Tensor product spline object.

    c[i1, i2, ..., id] * B(x1, i1) * B(x2, i2) * ... * B(xd, id)

    Parameters
    ----------
    c : ndarray, shape (n1, n2, ..., nd, ...)
        b-spline coefficients
    t : tuple of 1D ndarrays
        knot vectors in directions 1, 2, ... d
        ``len(t[i]) == n[i] + k + 1``
    k : int or length-d tuple of integers
        spline degrees.
    extrapolate : bool, optional
        Whether to exrapolate based on first and last intervals in each
        dimension, or return `nan`. Default is to extrapolate.
    """
    def __init__(self, t, c, k, extrapolate=None):
        ndim = len(t)
        assert ndim <= len(c.shape)

        try:
            len(k)
        except TypeError:
            # make k a tuple
            k = (k,)*ndim

        if extrapolate is None:
            extrapolate = True
        self.extrapolate = bool(extrapolate)

        self.k = tuple(operator.index(ki) for ki in k)
        self.t = tuple(np.ascontiguousarray(ti, dtype=float) for ti in t)
        self.c = np.asarray(c)

        if len(k) != ndim:
            raise ValueError(f"len(t) = {ndim} != {len(k)} = len(k)")

        for d in range(ndim):
            td = self.t[d]
            kd = self.k[d]
            n = td.shape[0] - kd - 1
            if kd < 0:
                raise ValueError(f"Spline degree in dimension {d} cannot be"
                                 f" negative.")
            if td.ndim != 1:
                raise ValueError(f"Knot vector in dimension {d} must be"
                                 f" one-dimensional.")
            if n < kd + 1:
                raise ValueError(f"Need at least {2*kd + 2} knots for degree"
                                 f" {kd} in dimension {d}.")
            if (np.diff(td) < 0).any():
                raise ValueError(f"Knots in dimension {d} must be in a"
                                 f" non-decreasing order.")
            if len(np.unique(td[kd:n + 1])) < 2:
                raise ValueError(f"Need at least two internal knots in"
                                 f" dimension {d}.")
            if not np.isfinite(td).all():
                raise ValueError(f"Knots in dimension {d} should not have"
                                 f" nans or infs.")
            if self.c.ndim < ndim:
                raise ValueError(f"Coefficients must be at least"
                                 f" {d}-dimensional.")
            if self.c.shape[d] != n:
                raise ValueError(f"Knots, coefficients and degree in dimension"
                                 f" {d} are inconsistent:"
                                 f" got {self.c.shape[d]} coefficients for"
                                 f" {len(td)} knots, need at least {n} for"
                                 f" k={k}.")

        dt = _get_dtype(self.c.dtype)
        self.c = np.ascontiguousarray(self.c, dtype=dt)

        # tabulate the flat indices for iterating over the (k+1)**ndim subarray
        shape = tuple(kd + 1 for kd in self.k)
        indices = np.unravel_index(np.arange(prod(shape)), shape)
        self._indices_k1d = np.asarray(indices).T

        # replacement for np.ravel_multi_index for indexing of `c1`:
        c1 = self.c.reshape(self.c.shape[:ndim] + (-1,))
        strides_c1 = [1]*(ndim + 1)
        for d in range(ndim-1, -1, -1):
            strides_c1[d] = strides_c1[d+1] * c1.shape[d+1]

        assert strides_c1 == [_//c1.dtype.itemsize for _ in c1.strides]
        assert strides_c1[-1] == 1
        self._strides_c1 = np.asarray(strides_c1)

    def __call__(self, xi, nu=None, extrapolate=None):
        """Evaluate the tensor product b-spline at `xi`.

        Parameters
        ----------
        xi : array_like, shape(..., ndim)
            The coordinates to evaluate the interpolator at.
            This can be a list or tuple of ndim-dimensional points
            or an array with the shape (num_points, ndim).
        nu : array_like, optional, shape (ndim,)
            Orders of derivatives to evaluate. Each must be non-negative.
        extrapolate : book, optional
            Whether to exrapolate based on first and last intervals in each
            dimension, or return `nan`. Default is to ``self.extrapolate`.

        Returns
        -------
        values : ndarray, shape xi.shape[:-1] + self.c.shape[ndim:]
            Interpolated values at xi
        """
        ndim = len(self.t)

        if extrapolate is None:
            extrapolate = self.extrapolate
        extrapolate = bool(extrapolate)

        if nu is None:
            nu = np.zeros((ndim,), dtype=np.intc)
        else:
            nu = np.asarray(nu, dtype=np.intc)
            if nu.ndim != 1 or nu.shape[0] != ndim:
                raise ValueError("invalid number of derivative orders nu")

        # prepare xi : shape (..., m1, ..., md) -> (1, m1, ..., md)
        xi = np.asarray(xi, dtype=float)
        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])
        xi = np.ascontiguousarray(xi)

        if xi_shape[-1] != ndim:
            raise ValueError(f"Shapes: xi.shape={xi_shape} and ndim={ndim}")
        assert xi_shape[-1] == xi.shape[-1]

        # prepare the coefficients: flatten the trailing dimensions
        c1 = self.c.reshape(self.c.shape[:ndim] + (-1,))
        c1r = c1.ravel()
        assert c1r.flags.c_contiguous

        num_c_tr = c1.shape[-1]  # # of trailing coefficients
        out = np.empty(xi.shape[:-1] + (num_c_tr,), dtype=c1.dtype)

        _bspl.evaluate_ndbspline(xi,
                                 self.t,
                                 self.k,
                                 nu,
                                 extrapolate,
                                 c1r,
                                 num_c_tr,
                                 self._strides_c1,
                                 self._indices_k1d,
                                 out,)

        return out.reshape(xi_shape[:-1] + self.c.shape[ndim:])
