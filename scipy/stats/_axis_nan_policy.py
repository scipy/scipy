# Many scipy.stats functions support `axis` and `nan_policy` parameters.
# When the two are combined, it can be tricky to get all the behavior just
# right. This file contains utility functions useful for scipy.stats functions
# that support `axis` and `nan_policy`, including a decorator that
# automatically adds `axis` and `nan_policy` arguments to a function.

import numpy as np
import scipy.stats
import scipy.stats._stats_py
from functools import wraps
from scipy._lib._docscrape import FunctionDoc, Parameter
import inspect


def _broadcast_array_shapes_remove_axis(arrays, axis=None):
    """
    Broadcast shapes of arrays, dropping specified axes

    Given a sequence of arrays `arrays` and an integer or tuple `axis`, find
    the shape of the broadcast result after consuming/dropping `axis`.
    In other words, return output shape of a typical hypothesis test on
    `arrays` vectorized along `axis`.

    Examples
    --------
    >>> a = np.zeros((5, 2, 1))
    >>> b = np.zeros((9, 3))
    >>> _broadcast_array_shapes((a, b), 1)
    (5, 3)
    """
    # Note that here, `axis=None` means do not consume/drop any axes - _not_
    # ravel arrays before broadcasting.
    shapes = [arr.shape for arr in arrays]
    return _broadcast_shapes_remove_axis(shapes, axis)


def _broadcast_shapes_remove_axis(shapes, axis=None):
    """
    Broadcast shapes, dropping specified axes

    Same as _broadcast_array_shapes, but given a sequence
    of array shapes `shapes` instead of the arrays themselves.
    """
    n_dims = max([len(shape) for shape in shapes])
    new_shapes = np.ones((len(shapes), n_dims), dtype=int)
    for row, shape in zip(new_shapes, shapes):
        row[len(row)-len(shape):] = shape  # can't use negative indices (-0:)
    if axis is not None:
        new_shapes = np.delete(new_shapes, axis, axis=1)
    new_shape = np.max(new_shapes, axis=0)
    new_shape *= new_shapes.all(axis=0)
    if np.any(~((new_shapes == 1) | (new_shapes == new_shape))):
        raise ValueError("Array shapes are incompatible for broadcasting.")
    return tuple(new_shape)


def _broadcast_concatenate(xs, axis):
    """Concatenate arrays along an axis with broadcasting."""
    # prepend 1s to array shapes as needed
    ndim = max([x.ndim for x in xs])
    xs = [x.reshape([1]*(ndim-x.ndim) + list(x.shape)) for x in xs]
    # move the axis we're concatenating along to the end
    xs = [np.swapaxes(x, axis, -1) for x in xs]
    # determine final shape of all but the last axis
    shape = _broadcast_array_shapes_remove_axis(xs, axis=-1)
    # broadcast along all but the last axis
    xs = [np.broadcast_to(x, shape + (x.shape[-1],)) for x in xs]
    # concatenate along last axis
    res = np.concatenate(xs, axis=-1)
    # move the last axis back to where it was
    res = np.swapaxes(res, axis, -1)
    return res


# TODO: add support for `axis` tuples
def _remove_nans(samples, paired):
    "Remove nans from paired or unpaired samples"
    # potential optimization: don't copy arrays that don't contain nans
    if not paired:
        return [sample[~np.isnan(sample)] for sample in samples]

    # for paired samples, we need to remove the whole pair when any part
    # has a nan
    nans = np.isnan(samples[0])
    for sample in samples[1:]:
        nans = nans | np.isnan(sample)
    not_nans = ~nans
    return [sample[not_nans] for sample in samples]


def _check_empty_inputs(samples, axis):
    """
    Check for empty sample; return appropriate output for a vectorized hypotest
    """
    # if none of the samples are empty, we need to perform the test
    if not any((sample.size == 0 for sample in samples)):
        return None
    # otherwise, the statistic and p-value will be either empty arrays or
    # arrays with NaNs. Produce the appropriate array and return it.
    output_shape = _broadcast_array_shapes_remove_axis(samples, axis)
    output = np.ones(output_shape) * np.nan
    return output


# Standard docstring / signature entries for `axis` and `nan_policy`
_name = 'axis'
_type = "int or None, default: 0"
_desc = (
    """If an int, the axis of the input along which to compute the statistic.
The statistic of each axis-slice (e.g. row) of the input will appear in a
corresponding element of the output.
If ``None``, the input will be raveled before computing the statistic."""
    .split('\n'))
_axis_parameter_doc = Parameter(_name, _type, _desc)
_axis_parameter = inspect.Parameter(_name,
                                    inspect.Parameter.KEYWORD_ONLY,
                                    default=0)

_name = 'nan_policy'
_type = "{'propagate', 'omit', 'raise'}"
_desc = (
    """Defines how to handle input NaNs.

- ``propagate``: if a NaN is present in the axis slice (e.g. row) along
  which the  statistic is computed, the corresponding entry of the output
  will be NaN.
- ``omit``: NaNs will be omitted when performing the calculation.
  If insufficient data remains in the axis slice along which the
  statistic is computed, the corresponding entry of the output will be
  NaN.
- ``raise``: if a NaN is present, a ``ValueError`` will be raised."""
    .split('\n'))
_nan_policy_parameter_doc = Parameter(_name, _type, _desc)
_nan_policy_parameter = inspect.Parameter(_name,
                                          inspect.Parameter.KEYWORD_ONLY,
                                          default='propagate')


def _axis_nan_policy_factory(result_object, default_axis=0,
                             n_samples=1, paired=False,
                             result_unpacker=None, too_small=0):
    """Factory for a wrapper that adds axis/nan_policy params to a function.

    Parameters
    ----------
    result_object : callable
        Callable that returns an object of the type returned by the function
        being wrapped (e.g. the namedtuple or dataclass returned by a
        statistical test) provided the separate components (e.g. statistic,
        pvalue).
    default_axis : int, default: 0
        The default value of the axis argument. Standard is 0 except when
        backwards compatibility demands otherwise (e.g. `None`).
    n_samples : int or callable, default: 1
        The number of data samples accepted by the function
        (e.g. `mannwhitneyu`), a callable that accepts a dictionary of
        parameters passed into the function and returns the number of data
        samples (e.g. `wilcoxon`), or `None` to indicate an arbitrary number
        of samples (e.g. `kruskal`).
    paired : {False, True}
        Whether the function being wrapped treats the samples as paired (i.e.
        corresponding elements of each sample should be considered as different
        components of the same sample.)
    result_unpacker : callable, optional
        Function that unpacks the results of the function being wrapped into
        a tuple. This is essentially the inverse of `result_object`. Default
        is `None`, which is appropriate for statistical tests that return a
        statistic, pvalue tuple (rather than, e.g., a non-iterable datalass).
    too_small : int, default: 0
        The largest unnacceptably small sample for the function being wrapped.
        For example, some functions require samples of size two or more or they
        raise an error. This argument prevents the error from being raised when
        input is not 1D and instead places a NaN in the corresponding element
        of the result.
    """

    if result_unpacker is None:
        def result_unpacker(res):
            return res[..., 0], res[..., 1]

    def is_too_small(samples):
        for sample in samples:
            if len(sample) <= too_small:
                return True
        return False

    def axis_nan_policy_decorator(hypotest_fun_in):
        @wraps(hypotest_fun_in)
        def axis_nan_policy_wrapper(*args, _no_deco=False, **kwds):

            if _no_deco:  # for testing, decorator does nothing
                return hypotest_fun_in(*args, **kwds)

            # We need to be flexible about whether position or keyword
            # arguments are used, but we need to make sure users don't pass
            # both for the same parameter. To complicate matters, some
            # functions accept samples with *args, and some functions already
            # accept `axis` and `nan_policy` as positional arguments.
            # The strategy is to make sure that there is no duplication
            # between `args` and `kwds`, combine the two into `kwds`, then
            # the samples, `nan_policy`, and `axis` from `kwds`, as they are
            # dealt with separately.

            # Check for intersection between positional and keyword args
            params = list(inspect.signature(hypotest_fun_in).parameters)
            if n_samples is None:
                # Give unique names to each positional sample argument
                # Note that *args can't be provided as a keyword argument
                params = [f"arg{i}" for i in range(len(args))] + params[1:]

            d_args = dict(zip(params, args))
            intersection = set(d_args) & set(kwds)
            if intersection:
                message = (f"{hypotest_fun_in.__name__}() got multiple values "
                           f"for argument '{list(intersection)[0]}'")
                raise TypeError(message)

            # Consolidate other positional and keyword args into `kwds`
            kwds.update(d_args)

            # rename avoids UnboundLocalError
            if callable(n_samples):
                n_samp = n_samples(kwds)
            else:
                n_samp = n_samples or len(args)

            # Extract the things we need here
            samples = [np.atleast_1d(kwds.pop(param))
                       for param in params[:n_samp]]
            vectorized = True if 'axis' in params else False
            axis = kwds.pop('axis', default_axis)
            nan_policy = kwds.pop('nan_policy', 'propagate')
            del args  # avoid the possibility of passing both `args` and `kwds`

            if axis is None:
                samples = [sample.ravel() for sample in samples]
                axis = 0
            elif axis != int(axis):
                raise ValueError('`axis` must be an integer')
            axis = int(axis)

            # if axis is not needed, just handle nan_policy and return
            ndims = np.array([sample.ndim for sample in samples])
            if np.all(ndims <= 1):
                # Addresses nan_policy == "raise"
                contains_nans = []
                for sample in samples:
                    contains_nan, _ = (
                        scipy.stats._stats_py._contains_nan(sample, nan_policy))
                    contains_nans.append(contains_nan)

                # Addresses nan_policy == "propagate"
                # Consider adding option to let function propagate nans, but
                # currently the hypothesis tests this is applied to do not
                # propagate nans in a sensible way
                if any(contains_nans) and nan_policy == 'propagate':
                    return result_object(np.nan, np.nan)

                # Addresses nan_policy == "omit"
                if any(contains_nans) and nan_policy == 'omit':
                    # consider passing in contains_nans
                    samples = _remove_nans(samples, paired)

                # ideally, this is what the behavior would be, but some
                # existing functions raise exceptions, so overriding it
                # would break backward compatibility.
                # if is_too_small(samples):
                #     return result_object(np.nan, np.nan)

                return hypotest_fun_in(*samples, **kwds)

            # check for empty input
            # ideally, move this to the top, but some existing functions raise
            # exceptions for empty input, so overriding it would break
            # backward compatibility.
            empty_output = _check_empty_inputs(samples, axis)
            if empty_output is not None:
                statistic = empty_output
                pvalue = empty_output.copy()
                return result_object(statistic, pvalue)

            # otherwise, concatenate all samples along axis, remembering where
            # each separate sample begins
            lengths = np.array([sample.shape[axis] for sample in samples])
            split_indices = np.cumsum(lengths)
            x = _broadcast_concatenate(samples, axis)

            # Addresses nan_policy == "raise"
            contains_nan, _ = (
                scipy.stats._stats_py._contains_nan(x, nan_policy))

            if vectorized and not contains_nan:
                return hypotest_fun_in(*samples, axis=axis, **kwds)

            # Addresses nan_policy == "omit"
            if contains_nan and nan_policy == 'omit':
                def hypotest_fun(x):
                    samples = np.split(x, split_indices)[:n_samp]
                    samples = _remove_nans(samples, paired)
                    if is_too_small(samples):
                        return result_object(np.nan, np.nan)
                    return hypotest_fun_in(*samples, **kwds)

            # Addresses nan_policy == "propagate"
            elif contains_nan and nan_policy == 'propagate':
                def hypotest_fun(x):
                    if np.isnan(x).any():
                        return result_object(np.nan, np.nan)
                    samples = np.split(x, split_indices)[:n_samp]
                    return hypotest_fun_in(*samples, **kwds)

            else:
                def hypotest_fun(x):
                    samples = np.split(x, split_indices)[:n_samp]
                    return hypotest_fun_in(*samples, **kwds)

            x = np.moveaxis(x, axis, -1)
            res = np.apply_along_axis(hypotest_fun, axis=-1, arr=x)
            return result_object(*result_unpacker(res))

        doc = FunctionDoc(axis_nan_policy_wrapper)
        parameter_names = [param.name for param in doc['Parameters']]
        if 'axis' in parameter_names:
            doc['Parameters'][parameter_names.index('axis')] = (
                _axis_parameter_doc)
        else:
            doc['Parameters'].append(_axis_parameter_doc)
        if 'nan_policy' in parameter_names:
            doc['Parameters'][parameter_names.index('nan_policy')] = (
                _nan_policy_parameter_doc)
        else:
            doc['Parameters'].append(_nan_policy_parameter_doc)
        doc = str(doc).split("\n", 1)[1]  # remove signature
        axis_nan_policy_wrapper.__doc__ = str(doc)

        sig = inspect.signature(axis_nan_policy_wrapper)
        parameters = sig.parameters
        parameter_list = list(parameters.values())
        if 'axis' not in parameters:
            parameter_list.append(_axis_parameter)
        if 'nan_policy' not in parameters:
            parameter_list.append(_nan_policy_parameter)
        sig = sig.replace(parameters=parameter_list)
        axis_nan_policy_wrapper.__signature__ = sig

        return axis_nan_policy_wrapper
    return axis_nan_policy_decorator
