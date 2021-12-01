# Many scipy.stats functions support `axis` and `nan_policy` parameters.
# When the two are combined, it can be tricky to get all the behavior just
# right. This file contains a suite of common tests for scipy.stats functions
# that support `axis` and `nan_policy` and additional tests for some associated
# functions in stats._util.

from itertools import product, combinations_with_replacement
import re
import pickle
import pytest

import numpy as np
from numpy.lib import NumpyVersion
from numpy.testing import assert_allclose, assert_equal
from scipy import stats


axis_nan_policy_cases = [
    # function, args, kwds, number of samples, paired, unpacker function
    # args, kwds typically aren't needed; just showing that they work
    (stats.kruskal, tuple(), dict(), 3, False, None),  # 4 samples is slow
    (stats.ranksums, ('less',), dict(), 2, False, None),
    (stats.mannwhitneyu, tuple(), {'method': 'asymptotic'}, 2, False, None),
    (stats.wilcoxon, ('pratt',), {'mode': 'auto'}, 2, True, None),
    (stats.wilcoxon, tuple(), dict(), 1, True, None),
    ]

# If the message is one of those expected, put nans in
# appropriate places of `statistics` and `pvalues`
too_small_messages = {"The input contains nan",  # for nan_policy="raise"
                      "Degrees of freedom <= 0 for slice",
                      "x and y should have at least 5 elements",
                      "Data must be at least length 3",
                      "The sample must contain at least two",
                      "x and y must contain at least two",
                      "division by zero",
                      "Mean of empty slice",
                      "Data passed to ks_2samp must not be empty",
                      "Not enough test observations",
                      "Not enough other observations",
                      "At least one observation is required",
                      "zero-size array to reduction operation maximum",
                      "`x` and `y` must be of nonzero size.",
                      "The exact distribution of the Wilcoxon test"}


def _mixed_data_generator(n_samples, n_repetitions, axis, rng,
                          paired=False):
    # generate random samples to check the response of hypothesis tests to
    # samples with different (but broadcastable) shapes and various
    # nan patterns (e.g. all nans, some nans, no nans) along axis-slices

    data = []
    for i in range(n_samples):
        n_patterns = 6  # number of distinct nan patterns
        n_obs = 20 if paired else 20 + i  # observations per axis-slice
        x = np.ones((n_repetitions, n_patterns, n_obs)) * np.nan

        for j in range(n_repetitions):
            samples = x[j, :, :]

            # case 0: axis-slice with all nans (0 reals)
            # cases 1-3: axis-slice with 1-3 reals (the rest nans)
            # case 4: axis-slice with mostly (all but two) reals
            # case 5: axis slice with all reals
            for k, n_reals in enumerate([0, 1, 2, 3, n_obs-2, n_obs]):
                # for cases 1-3, need paired nansw  to be in the same place
                indices = rng.permutation(n_obs)[:n_reals]
                samples[k, indices] = rng.random(size=n_reals)

            # permute the axis-slices just to show that order doesn't matter
            samples[:] = rng.permutation(samples, axis=0)

        # For multi-sample tests, we want to test broadcasting and check
        # that nan policy works correctly for each nan pattern for each input.
        # This takes care of both simultaneosly.
        new_shape = [n_repetitions] + [1]*n_samples + [n_obs]
        new_shape[1 + i] = 6
        x = x.reshape(new_shape)

        x = np.moveaxis(x, -1, axis)
        data.append(x)
    return data


def _homogeneous_data_generator(n_samples, n_repetitions, axis, rng,
                                paired=False, all_nans=True):
    # generate random samples to check the response of hypothesis tests to
    # samples with different (but broadcastable) shapes and homogeneous
    # data (all nans or all finite)
    data = []
    for i in range(n_samples):
        n_obs = 20 if paired else 20 + i  # observations per axis-slice
        shape = [n_repetitions] + [1]*n_samples + [n_obs]
        shape[1 + i] = 2
        x = np.ones(shape) * np.nan if all_nans else rng.random(shape)
        x = np.moveaxis(x, -1, axis)
        data.append(x)
    return data


def nan_policy_1d(hypotest, data1d, unpacker, *args,
                  nan_policy='raise', paired=False, _no_deco=True, **kwds):
    # Reference implementation for how `nan_policy` should work for 1d samples

    if nan_policy == 'raise':
        for sample in data1d:
            if np.any(np.isnan(sample)):
                raise ValueError("The input contains nan values")

    elif nan_policy == 'propagate':
        # For all hypothesis tests tested, returning nans is the right thing.
        # But many hypothesis tests don't propagate correctly (e.g. they treat
        # np.nan the same as np.inf, which doesn't make sense when ranks are
        # involved) so override that behavior here.
        for sample in data1d:
            if np.any(np.isnan(sample)):
                return np.nan, np.nan

    elif nan_policy == 'omit':
        # manually omit nans (or pairs in which at least one element is nan)
        if not paired:
            data1d = [sample[~np.isnan(sample)] for sample in data1d]
        else:
            nan_mask = np.isnan(data1d[0])
            for sample in data1d[1:]:
                nan_mask = np.logical_or(nan_mask, np.isnan(sample))
            data1d = [sample[~nan_mask] for sample in data1d]

    return unpacker(hypotest(*data1d, *args, _no_deco=_no_deco, **kwds))


@pytest.mark.parametrize(("hypotest", "args", "kwds", "n_samples", "paired",
                          "unpacker"), axis_nan_policy_cases)
@pytest.mark.parametrize(("nan_policy"), ("propagate", "omit", "raise"))
@pytest.mark.parametrize(("axis"), (1,))
@pytest.mark.parametrize(("data_generator"), ("mixed",))
def test_axis_nan_policy_fast(hypotest, args, kwds, n_samples, paired,
                              unpacker, nan_policy, axis,
                              data_generator):
    _axis_nan_policy_test(hypotest, args, kwds, n_samples, paired,
                          unpacker, nan_policy, axis, data_generator)


@pytest.mark.slow
@pytest.mark.parametrize(("hypotest", "args", "kwds", "n_samples", "paired",
                          "unpacker"), axis_nan_policy_cases)
@pytest.mark.parametrize(("nan_policy"), ("propagate", "omit", "raise"))
@pytest.mark.parametrize(("axis"), range(-3, 3))
@pytest.mark.parametrize(("data_generator"),
                         ("all_nans", "all_finite", "mixed"))
def test_axis_nan_policy_full(hypotest, args, kwds, n_samples, paired,
                              unpacker, nan_policy, axis,
                              data_generator):
    _axis_nan_policy_test(hypotest, args, kwds, n_samples, paired,
                          unpacker, nan_policy, axis, data_generator)


def _axis_nan_policy_test(hypotest, args, kwds, n_samples, paired,
                          unpacker, nan_policy, axis, data_generator):
    # Tests the 1D and vectorized behavior of hypothesis tests against a
    # reference implementation (nan_policy_1d with np.ndenumerate)

    # Some hypothesis tests return a non-iterable that needs an `unpacker` to
    # extract the statistic and p-value. For those that don't:
    if not unpacker:
        def unpacker(res):
            return res

    if NumpyVersion(np.__version__) < '1.18.0':
        pytest.xfail("Generator `permutation` method doesn't support `axis`")
    rng = np.random.default_rng(0)

    # Generate multi-dimensional test data with all important combinations
    # of patterns of nans along `axis`
    n_repetitions = 3  # number of repetitions of each pattern
    data_gen_kwds = {'n_samples': n_samples, 'n_repetitions': n_repetitions,
                     'axis': axis, 'rng': rng, 'paired': paired}
    if data_generator == 'mixed':
        inherent_size = 6  # number of distinct types of patterns
        data = _mixed_data_generator(**data_gen_kwds)
    elif data_generator == 'all_nans':
        inherent_size = 2  # hard-coded in _homogeneous_data_generator
        data_gen_kwds['all_nans'] = True
        data = _homogeneous_data_generator(**data_gen_kwds)
    elif data_generator == 'all_finite':
        inherent_size = 2  # hard-coded in _homogeneous_data_generator
        data_gen_kwds['all_nans'] = False
        data = _homogeneous_data_generator(**data_gen_kwds)

    output_shape = [n_repetitions] + [inherent_size]*n_samples

    # To generate reference behavior to compare against, loop over the axis-
    # slices in data. Make indexing easier by moving `axis` to the end and
    # broadcasting all samples to the same shape.
    data_b = [np.moveaxis(sample, axis, -1) for sample in data]
    data_b = [np.broadcast_to(sample, output_shape + [sample.shape[-1]])
              for sample in data_b]
    statistics = np.zeros(output_shape)
    pvalues = np.zeros(output_shape)

    for i, _ in np.ndenumerate(statistics):
        data1d = [sample[i] for sample in data_b]
        with np.errstate(divide='ignore', invalid='ignore'):
            try:
                res1d = nan_policy_1d(hypotest, data1d, unpacker, *args,
                                      nan_policy=nan_policy, paired=paired,
                                      _no_deco=True, **kwds)

                # Eventually we'll check the results of a single, vectorized
                # call of `hypotest` against the arrays `statistics` and
                # `pvalues` populated using the reference `nan_policy_1d`.
                # But while we're at it, check the results of a 1D call to
                # `hypotest` against the reference `nan_policy_1d`.
                res1db = unpacker(hypotest(*data1d, *args,
                                           nan_policy=nan_policy, **kwds))
                assert_equal(res1db[0], res1d[0])
                if len(res1db) == 2:
                    assert_equal(res1db[1], res1d[1])

            # When there is not enough data in 1D samples, many existing
            # hypothesis tests raise errors instead of returning nans .
            # For vectorized calls, we put nans in the corresponding elements
            # of the output.
            except (RuntimeWarning, ValueError, ZeroDivisionError) as e:

                # whatever it is, make sure same error is raised by both
                # `nan_policy_1d` and `hypotest`
                with pytest.raises(type(e), match=re.escape(str(e))):
                    nan_policy_1d(hypotest, data1d, unpacker, *args,
                                  nan_policy=nan_policy, paired=paired,
                                  _no_deco=True, **kwds)
                with pytest.raises(type(e), match=re.escape(str(e))):
                    hypotest(*data1d, *args, nan_policy=nan_policy, **kwds)

                if any([str(e).startswith(message)
                        for message in too_small_messages]):
                    res1d = np.nan, np.nan
                else:
                    raise e
        statistics[i] = res1d[0]
        if len(res1d) == 2:
            pvalues[i] = res1d[1]

    # Perform a vectorized call to the hypothesis test.
    # If `nan_policy == 'raise'`, check that it raises the appropriate error.
    # If not, compare against the output against `statistics` and `pvalues`
    if nan_policy == 'raise' and not data_generator == "all_finite":
        message = 'The input contains nan values'
        with pytest.raises(ValueError, match=message):
            hypotest(*data, axis=axis, nan_policy=nan_policy, *args, **kwds)

    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            res = unpacker(hypotest(*data, axis=axis, nan_policy=nan_policy,
                                    *args, **kwds))

        assert_equal(res[0], statistics)
        assert_equal(res[0].dtype, statistics.dtype)
        if len(res) == 2:
            assert_equal(res[1], pvalues)
            assert_equal(res[1].dtype, pvalues.dtype)


@pytest.mark.parametrize(("hypotest", "args", "kwds", "n_samples", "paired",
                          "unpacker"), axis_nan_policy_cases)
@pytest.mark.parametrize(("nan_policy"), ("propagate", "omit", "raise"))
@pytest.mark.parametrize(("data_generator"),
                         ("all_nans", "all_finite", "mixed", "empty"))
def test_axis_nan_policy_axis_is_None(hypotest, args, kwds, n_samples, paired,
                                      unpacker, nan_policy, data_generator):
    # check for correct behavior when `axis=None`

    if not unpacker:
        def unpacker(res):
            return res

    if NumpyVersion(np.__version__) < '1.18.0':
        pytest.xfail("Generator `permutation` method doesn't support `axis`")
    rng = np.random.default_rng(0)

    if data_generator == "empty":
        data = [rng.random((2, 0)) for i in range(n_samples)]
    else:
        data = [rng.random((2, 20)) for i in range(n_samples)]

    if data_generator == "mixed":
        masks = [rng.random((2, 20)) > 0.9 for i in range(n_samples)]
        for sample, mask in zip(data, masks):
            sample[mask] = np.nan
    elif data_generator == "all_nans":
        data = [sample * np.nan for sample in data]

    data_raveled = [sample.ravel() for sample in data]

    if nan_policy == 'raise' and data_generator not in {"all_finite", "empty"}:
        message = 'The input contains nan values'

        # check for correct behavior whether or not data is 1d to begin with
        with pytest.raises(ValueError, match=message):
            hypotest(*data, axis=None, nan_policy=nan_policy,
                     *args, **kwds)
        with pytest.raises(ValueError, match=message):
            hypotest(*data_raveled, axis=None, nan_policy=nan_policy,
                     *args, **kwds)

    else:
        # behavior of reference implementation with 1d input, hypotest with 1d
        # input, and hypotest with Nd input should match, whether that means
        # that outputs are equal or they raise the same exception

        ea_str, eb_str, ec_str = None, None, None
        with np.errstate(divide='ignore', invalid='ignore'):
            try:
                res1da = nan_policy_1d(hypotest, data_raveled, unpacker, *args,
                                       nan_policy=nan_policy, paired=paired,
                                       _no_deco=True, **kwds)
            except (RuntimeWarning, ValueError, ZeroDivisionError) as ea:
                ea_str = str(ea)

            try:
                res1db = unpacker(hypotest(*data_raveled, *args,
                                           nan_policy=nan_policy, **kwds))
            except (RuntimeWarning, ValueError, ZeroDivisionError) as eb:
                eb_str = str(eb)

            try:
                res1dc = unpacker(hypotest(*data, *args, axis=None,
                                           nan_policy=nan_policy, **kwds))
            except (RuntimeWarning, ValueError, ZeroDivisionError) as ec:
                ec_str = str(ec)

            if ea_str or eb_str or ec_str:
                assert any([str(ea_str).startswith(message)
                            for message in too_small_messages])
                assert ea_str == eb_str == ec_str
            else:
                assert_equal(res1db, res1da)
                assert_equal(res1dc, res1da)


@pytest.mark.parametrize(("axis"), (0, 1, 2))
def test_axis_nan_policy_decorated_positional_axis(axis):
    # Test for correct behavior of function decorated with
    # _axis_nan_policy_decorator whether `axis` is provided as positional or
    # keyword argument
    if NumpyVersion(np.__version__) < '1.18.0':
        pytest.xfail("Avoid test failures due to old version of NumPy")

    shape = (8, 9, 10)
    rng = np.random.default_rng(0)
    x = rng.random(shape)
    y = rng.random(shape)
    res1 = stats.mannwhitneyu(x, y, True, 'two-sided', axis)
    res2 = stats.mannwhitneyu(x, y, True, 'two-sided', axis=axis)
    assert_equal(res1, res2)

    message = "mannwhitneyu() got multiple values for argument 'axis'"
    with pytest.raises(TypeError, match=re.escape(message)):
        stats.mannwhitneyu(x, y, True, 'two-sided', axis, axis=axis)


def test_axis_nan_policy_decorated_positional_args():
    # Test for correct behavior of function decorated with
    # _axis_nan_policy_decorator when function accepts *args
    if NumpyVersion(np.__version__) < '1.18.0':
        pytest.xfail("Avoid test failures due to old version of NumPy")

    shape = (3, 8, 9, 10)
    rng = np.random.default_rng(0)
    x = rng.random(shape)
    x[0, 0, 0, 0] = np.nan
    stats.kruskal(*x)

    message = "kruskal() got an unexpected keyword argument 'args'"
    with pytest.raises(TypeError, match=re.escape(message)):
        stats.kruskal(args=x)

    with pytest.raises(TypeError, match=re.escape(message)):
        stats.kruskal(*x, args=x)


def test_axis_nan_policy_decorated_keyword_samples():
    # Test for correct behavior of function decorated with
    # _axis_nan_policy_decorator whether samples are provided as positional or
    # keyword arguments
    if NumpyVersion(np.__version__) < '1.18.0':
        pytest.xfail("Avoid test failures due to old version of NumPy")

    shape = (2, 8, 9, 10)
    rng = np.random.default_rng(0)
    x = rng.random(shape)
    x[0, 0, 0, 0] = np.nan
    res1 = stats.mannwhitneyu(*x)
    res2 = stats.mannwhitneyu(x=x[0], y=x[1])
    assert_equal(res1, res2)

    message = "mannwhitneyu() got multiple values for argument"
    with pytest.raises(TypeError, match=re.escape(message)):
        stats.mannwhitneyu(*x, x=x[0], y=x[1])


@pytest.mark.parametrize(("hypotest", "args", "kwds", "n_samples", "paired",
                          "unpacker"), axis_nan_policy_cases)
def test_axis_nan_policy_decorated_pickled(hypotest, args, kwds, n_samples,
                                           paired, unpacker):
    if NumpyVersion(np.__version__) < '1.18.0':
        rng = np.random.RandomState(0)
    else:
        rng = np.random.default_rng(0)

    # Some hypothesis tests return a non-iterable that needs an `unpacker` to
    # extract the statistic and p-value. For those that don't:
    if not unpacker:
        def unpacker(res):
            return res

    data = rng.uniform(size=(n_samples, 2, 30))
    pickled_hypotest = pickle.dumps(hypotest)
    unpickled_hypotest = pickle.loads(pickled_hypotest)
    res1 = unpacker(hypotest(*data, *args, axis=-1, **kwds))
    res2 = unpacker(unpickled_hypotest(*data, *args, axis=-1, **kwds))
    assert_allclose(res1, res2, rtol=1e-12)


def test_check_empty_inputs():
    # Test that _check_empty_inputs is doing its job, at least for single-
    # sample inputs. (Multi-sample functionality is tested below.)
    # If the input sample is not empty, it should return None.
    # If the input sample is empty, it should return an array of NaNs or an
    # empty array of appropriate shape. np.mean is used as a reference for the
    # output because, like the statistics calculated by these functions,
    # it works along and "consumes" `axis` but preserves the other axes.
    for i in range(5):
        for combo in combinations_with_replacement([0, 1, 2], i):
            for axis in range(len(combo)):
                samples = (np.zeros(combo),)
                output = stats._axis_nan_policy._check_empty_inputs(samples,
                                                                    axis)
                if output is not None:
                    with np.testing.suppress_warnings() as sup:
                        sup.filter(RuntimeWarning, "Mean of empty slice.")
                        sup.filter(RuntimeWarning, "invalid value encountered")
                        reference = samples[0].mean(axis=axis)
                    np.testing.assert_equal(output, reference)


def _check_arrays_broadcastable(arrays, axis):
    # https://numpy.org/doc/stable/user/basics.broadcasting.html
    # "When operating on two arrays, NumPy compares their shapes element-wise.
    # It starts with the trailing (i.e. rightmost) dimensions and works its
    # way left.
    # Two dimensions are compatible when
    # 1. they are equal, or
    # 2. one of them is 1
    # ...
    # Arrays do not need to have the same number of dimensions."
    # (Clarification: if the arrays are compatible according to the criteria
    #  above and an array runs out of dimensions, it is still compatible.)
    # Below, we follow the rules above except ignoring `axis`

    n_dims = max([arr.ndim for arr in arrays])
    if axis is not None:
        # convert to negative axis
        axis = (-n_dims + axis) if axis >= 0 else axis

    for dim in range(1, n_dims+1):  # we'll index from -1 to -n_dims, inclusive
        if -dim == axis:
            continue  # ignore lengths along `axis`

        dim_lengths = set()
        for arr in arrays:
            if dim <= arr.ndim and arr.shape[-dim] != 1:
                dim_lengths.add(arr.shape[-dim])

        if len(dim_lengths) > 1:
            return False
    return True


@pytest.mark.slow
@pytest.mark.parametrize(("hypotest", "args", "kwds", "n_samples", "paired",
                          "unpacker"), axis_nan_policy_cases)
def test_empty(hypotest, args, kwds, n_samples, paired, unpacker):
    # test for correct output shape when at least one input is empty

    def small_data_generator(n_samples, n_dims):

        def small_sample_generator(n_dims):
            # return all possible "small" arrays in up to n_dim dimensions
            for i in n_dims:
                # "small" means with size along dimension either 0 or 1
                for combo in combinations_with_replacement([0, 1, 2], i):
                    yield np.zeros(combo)

        # yield all possible combinations of small samples
        gens = [small_sample_generator(n_dims) for i in range(n_samples)]
        for i in product(*gens):
            yield i

    n_dims = [2, 3]
    for samples in small_data_generator(n_samples, n_dims):

        # this test is only for arrays of zero size
        if not any((sample.size == 0 for sample in samples)):
            continue

        max_axis = max((sample.ndim for sample in samples))

        # need to test for all valid values of `axis` parameter, too
        for axis in range(-max_axis, max_axis):

            try:
                # After broadcasting, all arrays are the same shape, so
                # the shape of the output should be the same as a single-
                # sample statistic. Use np.mean as a reference.
                concat = stats._stats_py._broadcast_concatenate(samples, axis)
                with np.testing.suppress_warnings() as sup:
                    sup.filter(RuntimeWarning, "Mean of empty slice.")
                    sup.filter(RuntimeWarning, "invalid value encountered")
                    expected = np.mean(concat, axis=axis) * np.nan

                res = hypotest(*samples, *args, axis=axis, **kwds)

                if hasattr(res, 'statistic'):
                    assert_equal(res.statistic, expected)
                    assert_equal(res.pvalue, expected)
                else:
                    assert_equal(res, expected)

            except ValueError:
                # confirm that the arrays truly are not broadcastable
                assert not _check_arrays_broadcastable(samples, axis)

                # confirm that _both_ `_broadcast_concatenate` and `hypotest`
                # produce this information.
                message = "Array shapes are incompatible for broadcasting."
                with pytest.raises(ValueError, match=message):
                    stats._stats_py._broadcast_concatenate(samples, axis)
                with pytest.raises(ValueError, match=message):
                    hypotest(*samples, *args, axis=axis, **kwds)
