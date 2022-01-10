import functools
import numpy as np
from scipy._lib.uarray import Dispatchable, all_of_type, create_multimethod
from scipy.linalg import _api
from scipy.linalg._backend import scalar_tuple_callable_array

__all__ = [
    # sketches
    'clarkson_woodruff_transform'
]


_create_linalg = functools.partial(
                    create_multimethod,
                    domain="numpy.scipy.linalg"
                 )


_mark_scalar_tuple_callable_array = functools.partial(
                            Dispatchable,
                            dispatch_type=scalar_tuple_callable_array,
                            coercible=True
                        )


def _get_docs(func):
    func.__doc__ = getattr(_api, func.__name__).__doc__
    return func


def _input_sketch_seed_replacer(args, kwargs, dispatchables):
    def self_method(input_matrix, sketch_size, seed=None, *args, **kwargs):
        return (dispatchables[0], dispatchables[1],
                dispatchables[2]) + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_input_sketch_seed_replacer)
@all_of_type(np.ndarray)
@_get_docs
def clarkson_woodruff_transform(input_matrix, sketch_size, seed=None):
    return (input_matrix, Dispatchable(sketch_size, int),
            _mark_scalar_tuple_callable_array(seed))
