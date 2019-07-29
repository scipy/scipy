import numpy as np


# TODO: Build with OpenMp and add configuration support
_default_workers = 1


def _init_nd_shape_and_axes(x, shape, axes):
    """Handles shape and axes arguments for nd transforms"""
    noshape = shape is None
    noaxes = axes is None

    if noaxes:
        axes = range(x.ndim)
    else:
        axes = np.atleast_1d(axes)

        if axes.size == 0:
            axes = axes.astype(np.intc)

        if not axes.ndim == 1:
            raise ValueError("when given, axes values must be a scalar or vector")
        if not np.issubdtype(axes.dtype, np.integer):
            raise ValueError("when given, axes values must be integers")

        axes = [a + x.ndim if a < 0 else a for a in axes]

        if any(a >= x.ndim or a < 0 for a in axes):
            raise ValueError("axes exceeds dimensionality of input")
        if len(set(axes)) != len(axes):
            raise ValueError("all axes must be unique")

    if not noshape:
        shape = np.atleast_1d(shape)

        if shape.size == 0:
            shape = shape.astype(np.intc)

        if shape.ndim != 1:
            raise ValueError("when given, shape values must be a scalar or vector")
        if not np.issubdtype(shape.dtype, np.integer):
            raise ValueError("when given, shape values must be integers")
        if len(axes) != len(shape):
            raise ValueError("when given, axes and shape arguments"
                             " have to be of the same length")

        shape = [x.shape[a] if s == -1 else s for s, a in zip(shape, axes)]
    elif noaxes:
        shape = list(x.shape)
    else:
        shape = [x.shape[a] for a in axes]

    if any(s < 1 for s in shape):
        raise ValueError(
            "invalid number of data points ({0}) specified".format(shape))

    return shape, axes


def _asfarray(x):
    """
    Convert to array with floating or complex dtype.

    float16 values are also promoted to float32.
    """
    if not hasattr(x, "dtype"):
        x = np.asarray(x)

    if x.dtype == np.float16:
        return np.asarray(x, np.float32)
    elif x.dtype.kind not in 'fc':
        return np.asarray(x, np.float64)

    # Always align input
    return np.array(x, copy=not x.flags['ALIGNED'])

def _datacopied(arr, original):
    """
    Strict check for `arr` not sharing any data with `original`,
    under the assumption that arr = asarray(original)
    """
    if arr is original:
        return False
    if not isinstance(original, np.ndarray) and hasattr(original, '__array__'):
        return False
    return arr.base is None


def _fix_shape(x, shape, axes):
    """Internal auxiliary function for _raw_fft, _raw_fftnd."""
    must_copy = False

    # Build an nd slice with the dimensions to be read from x
    index = [slice(None)]*x.ndim
    for n, ax in zip(shape, axes):
        if x.shape[ax] >= n:
            index[ax] = slice(0, n)
        else:
            index[ax] = slice(0, x.shape[ax])
            must_copy = True

    index = tuple(index)

    if not must_copy:
        return x[index], False

    s = list(x.shape)
    for n, axis in zip(shape, axes):
        s[axis] = n

    z = np.zeros(s, x.dtype)
    z[index] = x[index]
    return z, True


def _fix_shape_1d(x, n, axis):
    if n < 1:
        raise ValueError(
            "invalid number of data points ({0}) specified".format(n))

    return _fix_shape(x, (n,), (axis,))


def _normalization(norm, forward):
    """Returns the pypocketfft normalization mode from the norm argument"""

    if norm is None:
        return 0 if forward else 2

    if norm == 'ortho':
        return 1

    raise ValueError(
        "Invalid norm value {}, should be None or \"ortho\".".format(norm))
