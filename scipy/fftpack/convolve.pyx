from scipy.fft._pocketfft.pypocketfft import r2r_fftpack
import numpy as np
cimport numpy as np
cimport cython

np.import_array()

__all__ = ['destroy_convolve_cache', 'convolve', 'convolve_z',
           'init_convolution_kernel']


def destroy_convolve_cache():
    pass  # We don't cache anything, needed for compatibility


@cython.boundscheck(False)
@cython.wraparound(False)
def convolve(inout, omega, swap_real_imag=False, overwrite_x=False):
    """y = convolve(x,omega,[swap_real_imag,overwrite_x])

    Wrapper for ``convolve``.

    Parameters
    ----------
    x : input rank-1 array('d') with bounds (n)
    omega : input rank-1 array('d') with bounds (n)

    Other Parameters
    ----------------
    overwrite_x : input int, optional
        Default: 0
    swap_real_imag : input int, optional
         Default: 0

    Returns
    -------
    y : rank-1 array('d') with bounds (n) and x storage
    """
    cdef:
        np.ndarray[np.float64_t, ndim=1] X_arr, w_arr
        double [:] w, X
        double c
        size_t n, i

    X = X_arr = np.array(inout, np.float64, copy=not overwrite_x)
    w = w_arr = np.asarray(omega, np.float64)
    n = X_arr.shape[0]
    if X_arr.ndim != 1 or w.ndim != 1 or w.shape[0] != n:
        raise ValueError(
            "inout and omega must be 1-dimensional arrays of the same length")

    r2r_fftpack(X_arr, None, True, True, out=X_arr)

    if swap_real_imag:
        # Swap packed real and imag components
        X, w = X_arr, w_arr

        X[0] *= w[0];

        for i in range(1, n - 1, 2):
            c = X[i] * w[i]
            X[i] = X[i + 1] * w[i + 1]
            X[i + 1] = c

        if (n % 2) == 0:
            X[n - 1] *= w[n - 1]
    else:
        X_arr *= w_arr

    r2r_fftpack(X_arr, None, False, False, out=X_arr)
    return X_arr

@cython.boundscheck(False)
@cython.wraparound(False)
def convolve_z(inout, omega_real, omega_imag, overwrite_x=False):
    """y = convolve_z(x,omega_real,omega_imag,[overwrite_x])

    Wrapper for ``convolve_z``.

    Parameters
    ----------
    x : input rank-1 array('d') with bounds (n)
    omega_real : input rank-1 array('d') with bounds (n)
    omega_imag : input rank-1 array('d') with bounds (n)

    Other Parameters
    ----------------
    overwrite_x : input int, optional
        Default: 0

    Returns
    -------
    y : rank-1 array('d') with bounds (n) and x storage
    """
    cdef:
        np.ndarray[np.float64_t, ndim=1] X_arr
        double [:] wr, wi, X
        size_t n, i
        double c

    X = X_arr = np.array(inout, np.float64, copy=not overwrite_x)
    wr = np.asarray(omega_real, np.float64)
    wi = np.asarray(omega_imag, np.float64)

    n = X_arr.shape[0]
    if (X_arr.ndim != 1 or wr.ndim != 1 or wr.shape[0] != n
        or wi.ndim != 1 or wi.shape[0] != n):
        raise ValueError(
            "inout and omega must be 1-dimensional arrays of the same length")

    r2r_fftpack(X_arr, None, True, True, out=X_arr)

    X[0] *= wr[0] + wi[0]
    if (n % 2) == 0:
        X[n - 1] *= wr[n - 1] + wi[n - 1]

    for i in range(1, n - 1, 2):
        c = X[i + 1] * wr[i + 1] + X[i] * wi[i]
        X[i] = X[i] * wr[i] + X[i + 1] * wi[i + 1]
        X[i + 1] = c

    r2r_fftpack(X_arr, None, False, False, out=X_arr)
    return X_arr


@cython.boundscheck(False)
@cython.wraparound(False)
def init_convolution_kernel(size_t n, object kernel_func,
                            ssize_t d=0, zero_nyquist=None,
                            tuple kernel_func_extra_args=()):
    """omega = init_convolution_kernel(n,kernel_func,[d,zero_nyquist,kernel_func_extra_args])

    Wrapper for ``init_convolution_kernel``.

    Parameters
    ----------
    n : input int
    kernel_func : call-back function

    Other Parameters
    ----------------
    d : input int, optional
        Default: 0
    kernel_func_extra_args : input tuple, optional
        Default: ()
    zero_nyquist : input int, optional
        Default: d%2

    Returns
    -------
    omega : rank-1 array('d') with bounds (n)

    Notes
    -----
    Call-back functions::

      def kernel_func(k): return kernel_func
      Required arguments:
        k : input int
      Return objects:
        kernel_func : float
    """
    cdef:
        np.ndarray[np.float64_t, ndim=1] omega_arr
        double [::1] omega
        size_t j, k, l
        double scale_real, scale_imag, x

    if zero_nyquist is None:
        zero_nyquist = (d % 2 != 0)

    omega = omega_arr = np.empty(n, np.float64)
    l = n if n % 2 != 0 else n - 1

    # omega[k] = pow(sqrt(-1),d) * kernel_func(k)
    # omega[0] = kernel_func(0)
    # conjugate(omega[-k]) == omega[k]

    x = kernel_func(0, *kernel_func_extra_args)
    omega[0] = x / n

    d %= 4
    scale_real = 1./n if d == 0 or d == 1 else -1./n
    scale_imag = 1./n if d == 0 or d == 3 else -1./n

    k = 1
    for j in range(1, l, 2):
        x = kernel_func(k, *kernel_func_extra_args)
        omega[j] = scale_real * x
        omega[j + 1] = scale_imag * x
        k += 1

    if (n % 2) == 0:
        if zero_nyquist:
            omega[n - 1] = 0.
        else:
            x = kernel_func(k, *kernel_func_extra_args)
            omega[n - 1] = scale_real * x

    return omega_arr
