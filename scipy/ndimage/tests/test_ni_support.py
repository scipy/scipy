import pytest

import numpy as np
from .._ni_support import _get_output


@pytest.mark.parametrize(
    ('output', 'in_spec', 'shape', 'cplx_out', 'expect_spec', 'error'),
    [
        # Different ways to spell float32
        ('f', ((2, 3), 'float32'), None, False,
         ((2, 3), 'float32'), None),
        ('f4', ((2, 3), 'float32'), None, False,
         ((2, 3), 'float32'), None),
        ('float32', ((2, 3), 'float32'), None, False,
         ((2, 3), 'float32'), None),
        (np.float32, ((2, 3), 'float32'), None, False,
         ((2, 3), 'float32'), None),
        (np.dtype('float32'), ((2, 3), 'float32'), None, False,
         ((2, 3), 'float32'), None),
        # Float64
        ('float', ((2, 3), 'float32'), None, False,
         ((2, 3), 'float64'), None),
        ('float64', ((2, 3), 'float32'), None, False,
         ((2, 3), 'float64'), None),
        (np.float64, ((2, 3), 'float32'), None, False,
         ((2, 3), 'float64'), None),
        (float, ((2, 3), 'float32'), None, False,
         ((2, 3), 'float64'), None),
        # Complex output
        ('complex64', ((2, 3), 'float32'), None, True,
         ((2, 3), 'complex64'), None),
        ('F', ((2, 3), 'float32'), None, True,
         ((2, 3), 'complex64'), None),
        (np.complex64, ((2, 3), 'float32'), None, True,
         ((2, 3), 'complex64'), None),
        (complex, ((2, 3), 'float32'), None, True,
         ((2, 3), 'complex128'), None),
        # Pre-allocated arrays - shape must match input shape or explicit shape
        (np.zeros((2, 3), 'float32'), ((2, 3), 'float32'), None, False,
         ((2, 3), 'float32'), None),
        (np.zeros((2, 3), 'float64'), ((2, 3), 'float32'), None, False,
         ((2, 3), 'float64'), None),
        (np.zeros((3, 2), 'float32'), ((2, 3), 'float32'), (3, 2), False,
         ((3, 2), 'float32'), None),
        # None output follows input, shape and complex_output
        (None, ((2, 3), 'float32'), None, False,
         ((2, 3), 'float32'), None),
        (None, ((2, 3), 'float32'), (3, 2), False,
         ((3, 2), 'float32'), None),
        (None, ((2, 3), 'float32'), None, True,
         ((2, 3), 'complex64'), None),
        (None, ((2, 3), 'float64'), None, True,
         ((2, 3), 'complex128'), None),
        # Error cases
        ('float32', ((2, 3), 'float32'), None, True,
         None, "output must have complex dtype"),
        ('void', ((2, 3), 'float32'), None, False,
         None, "output must have numeric dtype"),
        (np.zeros((3, 2)), ((2, 3), 'float32'), None, False,
         None, "shape not correct"),
        (np.zeros((2, 3)), ((2, 3), 'float32'), None, True,
         None, "output must have complex dtype"),
    ],
)
def test_get_output_defaults(output, in_spec, shape, cplx_out, expect_spec, error):
    input = np.zeros(*in_spec)
    if expect_spec is None:
        with pytest.raises(RuntimeError, match=error):
            _get_output(output, input, shape=shape, complex_output=cplx_out)
    else:
        result = _get_output(output, input, shape=shape, complex_output=cplx_out)
        assert result.shape == expect_spec[0]
        assert result.dtype == np.dtype(expect_spec[1])
