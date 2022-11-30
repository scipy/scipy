import operator
import itertools
import numpy as np

from scipy._lib._util import prod

__all__ = ["NdBspline"]

def _get_dtype(dtype):
    """Return np.complex128 for complex dtypes, np.float64 otherwise."""
    if np.issubdtype(dtype, np.complexfloating):
        return np.complex_
    else:
        return np.float_

# XXX: remove
def B(x, k, i, t):
    if k == 0:
        return 1.0 if t[i] <= x < t[i+1] else 0.0
    if t[i+k] == t[i]:
        c1 = 0.0
    else:
        c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, k-1, i, t)
    if t[i+k+1] == t[i+1]:
        c2 = 0.0
    else:
        c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, k-1, i+1, t)
    return c1 + c2
    

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
    """
    def __init__(self, t, c, k=3):
        ndim = len(t)
        assert ndim <= len(c.shape)
        
        try:
            len(k)
        except TypeError:
            # make k a tuple
            k = (k,)*ndim

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
                raise ValueError(f"Spline degree in dimension {d} cannot be negative.")
            if td.ndim != 1:
                raise ValueError(f"Knot vector in dimension {d} must be one-dimensional.")           
            if n < kd + 1:
                raise ValueError(f"Need at least {2*kd + 2} knots for degree {kd}"
                                 f" in dimension {d}.")
            if (np.diff(td) < 0).any():
                raise ValueError(f"Knots in dimension {d} must be in a non-decreasing order.")
            if len(np.unique(td[kd:n + 1])) < 2:
                raise ValueError(f"Need at least two internal knots in dimension {d}.")
            if not np.isfinite(td).all():
                raise ValueError(f"Knots in dimension {d} should not have nans or infs.")
            if self.c.ndim < ndim:
                raise ValueError(f"Coefficients must be at least {d}-dimensional.")
            if self.c.shape[d] < n:
                raise ValueError(f"Knots, coefficients and degree in dimension {d} are inconsistent:"
                                 f" got {self.c.shape[d]} coefficients for {len(td)} knots,"
                                 f" need at least {n} for k={k}.")            

        dt = _get_dtype(self.c.dtype)
        self.c = np.ascontiguousarray(self.c, dtype=dt)
                                 
        # tabulate the flat indices for iterating over the (k+1)**ndim subarray
        shape = tuple(kd + 1 for kd in self.k)
        indices = np.unravel_index( np.arange(prod(shape)), shape)
        self._indices_k1d = np.asarray(indices).T
        
    def __call__(self, xi):
        """Evaluate the tensor product b-spline at coordinates.
        
        Parameters
        ----------
        xi : array_like, shape(..., ndim)
            The coordinates to evaluate the interpolator at.
            This can be a list or tuple of ndim-dimensional points
            or an array with the shape (num_points, ndim).
            
        Returns
        -------
        values : ndarray, shape xi.shape[:-1] + self.c.shape[ndim:]
            Interpolated values at xi
        """
        ndim = len(self.t)

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
        out = np.empty(xi.shape[:-1] + (num_c_tr,), dtype=float)    

        # replacement for np.ravel_multi_index for indexing of `c1`
        strides_c1 = [1]*(ndim + 1)
        for d in range(ndim-1, -1, -1):
            strides_c1[d] = strides_c1[d+1] * c1.shape[d+1]
        assert strides_c1 == [_//8 for _ in c1.strides]
        assert strides_c1[-1] == 1

        # 'intervals': indices for a point in xi into the knot arrays t
        i = [-101,]*ndim
        # container for non-zero b-splines at each point in xi
        b = np.empty((ndim, max(self.k)+1), dtype=float) * np.nan

        ### Finally, iterate over the data points
        for j in range(xi.shape[0]):           
            x = xi[j]

            # get the indices in an ndim-dimensional vector
            for d in range(ndim):
                td, xd = self.t[d], x[d]
                k = self.k[d]

                # find the index for x[d]
                if xd == td[k]:
                    i[d] = k
                else:
                    i[d] = np.searchsorted(td, xd) - 1
                assert td[i[d]] <= xd <= td[i[d]+1]
                assert i[d] >= k and i[d] < len(td) - k

                # (k+1) b-splines which are non-zero at x[d] 
                b[d, :k+1] = [B(xd, k, j, td) for j in range(i[d]-k, i[d]+1)]   

            # iterate over the dimensions, form linear combinations of
            # products B(x_1) * B(x_2) * ... B(x_N) of (k+1)**N b-splines
            # which are non-zero at `i = (i_1, i_2, ..., i_N)`.
            result = np.zeros(num_c_tr, dtype=float)
            iters = [range(i[d] - self.k[d], i[d] + 1) for d in range(ndim)]
            for idx in itertools.product(*iters):
                factor = prod(b[d, idx[d] - i[d] + self.k[d]] for d in range(ndim))
                # loop over the trailing values of self.c explicitly
                for i_c in range(num_c_tr):
                    result[i_c] += c1[idx + (i_c,)] * factor
            
            
            ### The above is in principle enough for a prototype.  ###
            ### Replicate it with flat indexing below.             ###
            result_flat = np.zeros(num_c_tr, dtype=float)
            
            volume = prod([k+1 for k in self.k])  # (k+1)**ndim
            for iflat in range(volume):
                # idx_b = np.unravel_index(iflat, (k+1,)*ndim)   # equiv below
                idx_b = self._indices_k1d[iflat, :]
                assert all(idx_b == np.unravel_index(iflat, [kd+1 for kd in self.k])) # (k+1,)*ndim))
                
                # 1. Shift the subblock indices into indices into c1.ravel()
                # 2. Collect the product of non-zero b-splines at this value of the $x$ vector
                # 3. Compute the base index for iterating over the c1 array
                idx_c = [-101]*ndim
                idx_cflat_base = 0
                factor = 1.0
                for d in range(ndim):
                    factor *= b[d, idx_b[d]]
                    idx_c[d] = idx_b[d] + i[d] - self.k[d]             # XXX: remove later
                    idx_cflat_base += idx_c[d] * strides_c1[d]

                ##### double-check the above loop -- XXX: remove later
                # Product of non-zero bsplines at this value of the $x$ vector.
                factor2 = prod(b[d, idx_b[d]] for d in range(ndim))
        
                # shift the indices into c1.ravel()
                idx_c2 = list(idx_b)
                for d in range(ndim):
                    idx_c2[d] += i[d] - self.k[d] 

                idx_cflat_base2 = sum(idx_c[d] * strides_c1[d] for d in range(ndim))
                    
                assert idx_c2 == idx_c
                assert factor == factor2
                assert idx_cflat_base == idx_cflat_base2
                ####### end double-check
                    
                ### collect linear combinations of coef * factor
                for i_c in range(num_c_tr):
                    # this is equivalent to 
                    # idx_cflat = np.ravel_multi_index(tuple(idx_c) + (i_c,), c1.shape)
                    # we pre-computed the first ndim strides of `c1r` array and use the
                    # fact that the `c1r` array is C-ordered by construction
                    result_flat[i_c] += c1r[idx_cflat_base + i_c] * factor
                    out[j, i_c] = result_flat[i_c]

            # XXX: remove later
            from numpy.testing import assert_allclose
            assert_allclose(result, out[j, :], atol=1e-15)
            
            # copy the result over: this is in fact a loop over num_c_tr
            #out[j, ...] = result_flat
            """            
            from gobbledegook2 import _evaluate_ndbspline
            out2 = np.empty_like(out)
            _evaluate_ndbspline(xi,
                                self.t,
                                self.k,
                                c1r,
                                num_c_tr,
                                out2,
                                c1,
                                np.asarray(self._indices_k1d),
                                np.asarray(strides_c1)
                               )
            '''
            def _evaluate_ndbspline(const double[:, ::1] xi,
                        tuple t,
                        tuple ktuple,
                        const double[::1] c1r,
                        npy_intp num_c_tr,
                        double[:, ::1] out,
                        c1,        # XXX: remove
                        const npy_intp[:, ::] indices_k1d,
                        const npy_intp[::1] strides_c1,
            '''
            assert_allclose(out, out2, atol=1e-14)
            """
        return out.reshape(xi_shape[:-1] + self.c.shape[ndim:])
