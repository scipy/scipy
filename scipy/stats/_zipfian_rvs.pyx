#cython: language_level=3
"""
This module defines the function zipfian_rvs(bitgen_t gen, a, n, size=None)
that generates variates from the Zipfian distribution (scipy.stats.zipfian)
using a rejection method.
"""

# C standard library...
from libc.math cimport trunc, pow
from libc.stdint cimport int64_t

# Cython...
cimport cython

# CPython...
from cpython.pycapsule cimport PyCapsule_IsValid, PyCapsule_GetPointer
from cpython cimport PyFloat_AsDouble

# NumPy...
import numpy as np
cimport numpy as cnp
from numpy.random cimport bitgen_t

# SciPy Cython API from scipy.special...
from scipy.special.cython_special cimport boxcox, inv_boxcox

cnp.import_array()

@cython.boundscheck(False)
@cython.wraparound(False)
cdef verify_equal_shapes(cnp.broadcast it, cnp.ndarray out):
    # Check that the shape of the iterator `it` is the same as that of
    # the ndarray `out`.

    cdef bint valid_shapes = True
    cdef cnp.npy_intp i
    cdef cnp.npy_intp it_ndim = cnp.PyArray_MultiIter_NDIM(it)
    cdef cnp.npy_intp out_ndim = cnp.PyArray_NDIM(out)

    if it_ndim != out_ndim:
        valid_shapes = False
    else:
        for i in range(it_ndim):
            if cnp.PyArray_MultiIter_DIMS(it)[i] != cnp.PyArray_DIMS(out)[i]:
                valid_shapes = False
                break

    if not valid_shapes:
        # Create Python tuples to include the shapes in the error message.
        it_shape = tuple(cnp.PyArray_MultiIter_DIMS(it)[i] for i in range(it_ndim))
        out_shape = tuple(cnp.PyArray_DIMS(out)[i] for i in range(out_ndim))
        bad_shapes_msg = (f"The output size {out_shape} must exactly match the "
                          "result of broadcasting the input shapes with that "
                          f"size, but the broadcast shape is {it_shape}.")
        raise ValueError(bad_shapes_msg)


#
# Functions used in the implementation of the rejection method for
# the Zipfian distribution.
#

cdef inline double g(double x, double a, int64_t n) noexcept nogil:
    # g is the *nonnormalized* PDF of the dominating distribution.
    # It is assumed that 1 <= x < n + 1.
    if x < 2:
        return 1
    return pow(x - 1, -a)


cdef inline double g_rv(bitgen_t *bit_generator,
                        double a, int64_t n, double maxG) noexcept nogil:
    # Generate a random variate from the dominating distribution using
    # inversion.
    # G is proportional to the CDF of the dominating distribution, but G is
    # not normalized, so to implement the inversion method, the random
    # uniform input is drawn from the interval [0, maxG], where
    # maxG = max(G) = G(n + 1, a, n) = boxcox(n, 1 - a) + 1.
    # We save a few CPU cycles by computing maxG once outside of this function
    # instead of every call here.
    cdef double y = maxG * bit_generator.next_double(bit_generator.state)

    # The following computes Ginv(y, a, n), the inverse with respect to the
    # first parameter of G(x, a, n).
    if 0 <= y <= 1:
        return y + 1
    return 1 + inv_boxcox(y - 1, 1 - a)


cdef int zipfian_mrvs(bitgen_t *bit_generator,
                      double a, int64_t n,
                      size_t m, int64_t *out) noexcept nogil:
    # Generate m variates from the Zipfian distribution using the rejection method.
    #
    # The Cython-level lock of the wrapper of this bitgen_t instance must be locked
    # when this function is called.
    #
    # The target distribution for continuous x has the nonnormalized PDF floor(x)**-a,
    # which we implement here as pow(trunc(x), -a).
    #
    # The normal return value is 0.  If -1 is returned, it means the maximum
    # number of rejections was reached at least once.  The value of the random
    # variate in the `out` array will be -1 in that case.  This shouldn't happen,
    # and in extensive testing it has not been observed.  The check is implemented
    # out of an overabundance of caution.
    #
    # Compute maxG here instead of computing it in every call of g_rv() to save a
    # few CPU cycles.
    cdef double maxG = boxcox(n, 1 - a) + 1
    cdef double x
    cdef size_t i
    cdef int num_rejections
    cdef int max_rejections = 100
    cdef int status = 0

    for i in range(m):
        num_rejections = 0
        while num_rejections <= max_rejections:
            x = g_rv(bit_generator, a, n, maxG)
            # The dominating function g and the target function f coincide on the
            # interval 1 <= x < 2, so a candidate variate in that interval is never
            # rejected.
            if x < 2 or (bit_generator.next_double(bit_generator.state) * g(x, a, n)
                         <= pow(trunc(x), -a)):
                out[i] = <int64_t>(trunc(x))
                break
            num_rejections += 1
        else:
            # Too many rejections (this should not happen).
            out[i] = <int64_t>(-1)
            status = -1
    return status


def _zipfian_value_error(varname, badvalue):
    # varname is either 'a' or 'n'.
    if varname == 'a':
        constraint = "nonnegative"
    else:
        constraint = "greater than 0"
    return ValueError(f'zipfian: {varname} must be {constraint}, got {badvalue}')


cdef inline bitgen_t* get_bit_generator(bit_generator):
    # Get the C bitgen_t from the Python bit_generator instance.
    # The input bit_generator must be a Python capsule holding the pointer
    # to an underlying bitgen_t.
    cdef const char *capsule_name = "BitGenerator"
    capsule = bit_generator.capsule
    if not PyCapsule_IsValid(capsule, capsule_name):
        raise ValueError("bit_generator has an invalid capsule")
    return <bitgen_t *> PyCapsule_GetPointer(capsule, capsule_name)


@cython.boundscheck(False)
@cython.wraparound(False)
def zipfian_rvs(bit_generator, *, a, n, size=None):
    """
    Generate random variates from the "Zipfian" distribution.

    The distribution is also known as the generalized Zipf distribution.
    It is a discrete distribution with finite support {1, 2, ..., n}.  The
    probability of integer k is proportional to k**-a.

    Note: This implementation assumes that this function is being called
    from scipy.stats.zipfian._rvs() because a user called the .rvs() method.
    This means `a` and `n` are assumed to be numpy arrays, with `a` having dtype
    np.float64 and `n` having dtype `int64`.

    Parameters
    ----------
    bit_generator: NumPy BitGenerator instance
    a: np.float64
        The probability of integer k in the support is proportional to k**-a.
        `a` must be nonnegative.  When `a` is 0, the distribution is uniform.
    n: np.int64
        Determines the support {1, 2, ..., n} of the distribution.
        Must be at least 1.
    size: int or tuple of int
        Number of variates to generate.

    Returns
    -------
    variates: int64 or array of int64
        Zipfian random variates.
        The implementation uses a rejection method to generate the random variates.
        Theoretically, the rejection rate for this implementation should be very
        low; the average number of rejections per variate should be less than 1.
        As a precaution against extremely unlikely events (and against bugs in the
        code), the algorithm will return -1 if the number of rejections reaches 100.

    Examples
    --------
    >>> import numpy as np
    >>> from _zipfian_rvs import zipfian_rvs

    >>> bitgen = np.random.PCG64(121263137472525314065)
    >>> zipfian_rvs(bitgen, a=1.05, n=400, size=13)
    array([  1,   1, 374,   8,  57,   1,   2,  82,  45,   4, 122, 236,   2])

    The parameters broadcast:

    >>> a = np.array([0.5, 1.25])
    >>> n = np.array([[20], [100], [500]])

    `a` has shape (2,) and `n` has shape (3, 1).  The broadcast shape is (3, 2):

    >>> zipfian_rvs(bitgen, a=a, n=n)
    array([[ 13,   1],
           [  4,   2],
           [146,   2]])
    """
    cdef Py_ssize_t i
    cdef bitgen_t *bitgen
    cdef bint is_scalar
    # `variates` is the output array of random variates.
    cdef cnp.ndarray variates
    cdef cnp.int64_t *variates_data
    # `variate` is the return value when `a` and `n` are scalars and `size` is None.
    cdef cnp.int64_t variate
    cdef cnp.broadcast it
    cdef double a1
    cdef int64_t n1
    cdef int64_t *out
    cdef cnp.ndarray a_arr, n_arr

    # Get the C bitgen_t from the Python bit_generator instance.
    bitgen = get_bit_generator(bit_generator)

    # Assume that `a` is a numpy array with dtype float64, and
    # `n` is a numpy array with dtype int64.
    a_arr = <cnp.ndarray> a
    n_arr = <cnp.ndarray> n

    is_scalar = cnp.PyArray_NDIM(a_arr) == 0 and cnp.PyArray_NDIM(n_arr) == 0

    if not is_scalar:
        # At least one of 'a' and 'n' has ndim >= 1.

        # These checks are commented out; we assume these conditions were
        # checked by the caller. (Delete this eventually.)
        # if np.any(a_arr < 0):
        #     raise _zipfian_value_error("a", a_arr[a_arr < 0].item(0))
        # if np.any(n_arr < 1):
        #     raise _zipfian_value_error("n", n_arr[n_arr < 1].item(0))

        if size is not None:
            variates = <cnp.ndarray> np.empty(size, np.int64)
        else:
            # `size` was not given, so the output shape is determined by
            # broadcasting a_arr and n_arr.
            variates = <cnp.ndarray> np.empty(np.broadcast(a_arr, n_arr).shape,
                                              np.int64)

        # This call will catch most shape mismatches.
        it = cnp.PyArray_MultiIterNew3(a_arr, n_arr, variates)

        if size is not None:
            # If `size` was given, the output shape (variates.shape) must be
            # the same as the shape of broadcasting a_arr, n_arr and variates.
            # This is an additional requirement of random variate generators
            # that is not caught when `it` is created above.
            verify_equal_shapes(it, variates)

        with bit_generator.lock, nogil:
            # This is the main loop over the input and output arrays.
            for i in range(cnp.PyArray_SIZE(variates)):
                a1 = (<double*> cnp.PyArray_MultiIter_DATA(it, 0))[0]
                n1 = (<int64_t*> cnp.PyArray_MultiIter_DATA(it, 1))[0]
                out = <int64_t*> cnp.PyArray_MultiIter_DATA(it, 2)
                zipfian_mrvs(bitgen, a1, n1, 1, out)
                cnp.PyArray_MultiIter_NEXT(it)

        return variates

    # a and n are scalars...

    a1 = PyFloat_AsDouble(a)
    if a1 < 0:
        raise _zipfian_value_error("a", a)
    n1 = <int64_t> n
    if n1 < 1:
        raise _zipfian_value_error("n", n)

    if size is None:
        with bit_generator.lock:
            zipfian_mrvs(bitgen, a1, n1, 1, &variate)
        return variate

    # a and n are scalars, size is not None...

    variates = <cnp.ndarray> np.empty(size, np.int64)
    nvars = cnp.PyArray_SIZE(variates)
    variates_data = <cnp.int64_t *> cnp.PyArray_DATA(variates)

    with bit_generator.lock, nogil:
        # variates was just created above, so we know it will have contiguous
        # data, no need to worry about strides.
        zipfian_mrvs(bitgen, a1, n1, nvars, variates_data)

    return variates
