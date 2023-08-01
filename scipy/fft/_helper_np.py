from functools import update_wrapper, lru_cache

from ._pocketfft import helper as _helper


def next_fast_len(target, real=False):
    pass


# Directly wrap the c-function good_size but take the docstring etc., from the
# next_fast_len function above
next_fast_len = update_wrapper(lru_cache(_helper.good_size), next_fast_len)
next_fast_len.__wrapped__ = _helper.good_size


def _init_nd_shape_and_axes(x, shape, axes):
    return _helper._init_nd_shape_and_axes(x, shape, axes)
