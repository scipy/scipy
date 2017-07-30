from __future__ import division, print_function, absolute_import

import gc
import itertools

import numpy as np
from numpy.testing import assert_, assert_array_equal, assert_raises
from scipy._lib._numpy_compat import assert_raises_regex
from scipy.ndimage._nd_image import _test_converters


def make_array(shape, dtype=np.intp, byteswapped=False, misaligned=False,
               readonly=False):
    """Creates an array with the required characteristics."""
    dtype = np.dtype(dtype)
    if ((byteswapped and not dtype.isnative) or
            (not byteswapped and dtype.isnative)):
        dtype = dtype.newbyteorder()
    size = np.product(shape)
    if misaligned and dtype.itemsize > 1:
        base_array = np.arange(size + 1).astype(dtype).view(np.uint8)
        base_array = base_array[dtype.itemsize - 1:-1].view(dtype)
    else:
        base_array = np.arange(size).astype(dtype)
    array = base_array.reshape(shape)
    if readonly:
        array.flags.writeable = False
    return array


def check_converted_input_array(array, converted_array):
    # Converted array must have same contents as original.
    assert_array_equal(array, converted_array)
    # Converted array must be aligned and not swapped.
    assert_(converted_array.flags.aligned)
    assert_(converted_array.dtype.isnative)


def check_converted_output_array(array, converted_array):
    # Converted array must be writeable, aligned and not swapped.
    assert_(converted_array.flags.writeable)
    assert_(converted_array.flags.aligned)
    assert_(converted_array.dtype.isnative)
    # Either the array is passed unchanged...
    if converted_array is not array:
        # ...or it should use the updateifcopy mechanism.
        assert_(converted_array.flags.updateifcopy)
    # Changes to the converted array should be reflected in the
    # original one after deletion of the former.
    converted_array[:] = 1


def check_converted_input_output_array(array, converted_array):
    # Converted array must have same contents as original.
    assert_array_equal(array, converted_array)
    # Everything else is the same as a converted output array.
    check_converted_output_array(array, converted_array)


class TestNIConverters:

    def setUp(self):
        self.shape = (3, 5, 7)
        all_types = '?' + np.typecodes['AllInteger'] + 'fd'
        f_t = (False, True)
        self.all_args = list(itertools.product(all_types, f_t, f_t, f_t))

    def test_input_array(self):
        other_array1 = make_array(self.shape)
        other_array2 = make_array(self.shape)
        origin = (1,) * len(self.shape)
        for args in self.all_args:
            in_array = make_array(self.shape, *args)
            check_converted_input_array(
                in_array,
                _test_converters(in_array, None, other_array1, other_array2,
                                 None, origin)[0])

    def test_optional_input_array(self):
        other_array1 = make_array(self.shape)
        other_array2 = make_array(self.shape)
        other_array3 = make_array(self.shape)
        origin = (1,) * len(self.shape)
        # When passing None the return should also be None.
        assert_(_test_converters(other_array1, None, other_array2,
                                 other_array3, None, origin)[1]
                is None)
        for args in self.all_args:
            opt_in_array = make_array(self.shape, *args)
            check_converted_input_array(
                opt_in_array,
                _test_converters(other_array1, opt_in_array, other_array2,
                                 other_array3, None, origin)[1])

    def test_input_output_array(self):
        other_array1 = make_array(self.shape)
        other_array2 = make_array(self.shape)
        origin = (1,) * len(self.shape)
        # Passing a non-array should raise an error.
        assert_raises_regex(
            TypeError, 'input-output array must be an array',
            _test_converters, other_array1, None, None, other_array2, None,
            origin)
        for args in self.all_args:
            dt, bs, ma, ro = args
            in_out_array = make_array(self.shape, *args)
            if ro:
                assert_raises_regex(
                    ValueError, 'input-output array is read-only',
                    _test_converters, other_array1, None, in_out_array,
                    other_array2, None, origin)
            else:
                check_converted_input_output_array(
                    in_out_array,
                    _test_converters(other_array1, None, in_out_array,
                                     other_array2, None, origin)[2])
                # Changes to the converted array should be reflected in the
                # original one after deletion of the former.
                assert_(np.all(in_out_array == 1))
                # In the end, the original array must again be writeable.
                assert_(in_out_array.flags.writeable)

    def test_output_array(self):
        other_array1 = make_array(self.shape)
        other_array2 = make_array(self.shape)
        origin = (1,) * len(self.shape)
        # Passing a non-array should raise an error.
        assert_raises_regex(
            TypeError, 'output array must be an array',
            _test_converters, other_array1, None, other_array2, None, None,
            origin)
        for args in self.all_args:
            dt, bs, ma, ro = args
            out_array = make_array(self.shape, *args)
            if ro:
                assert_raises_regex(
                    ValueError, 'output array is read-only',
                    _test_converters, other_array1, None, other_array2,
                    out_array, None, origin)
            else:
                check_converted_output_array(
                    out_array,
                    _test_converters(other_array1, None, other_array2,
                                     out_array, None, origin)[3])
                # Changes to the converted array should be reflected in the
                # original one after deletion of the former.
                assert_(np.all(out_array == 1))
                # In the end, the original array must again be writeable.
                assert_(out_array.flags.writeable)

    def test_optional_output_array(self):
        other_array1 = make_array(self.shape)
        other_array2 = make_array(self.shape)
        other_array3 = make_array(self.shape)
        origin = (1,) * len(self.shape)
        # When passing None the return should also be None.
        assert_(_test_converters(other_array1, None, other_array2,
                                 other_array3, None, origin)[4]
                is None)
        # Passing a non-array should raise an error.
        assert_raises_regex(
            TypeError, 'output array must be an array',
            _test_converters, other_array1, None, other_array2, None,
            'not an array', origin)
        for args in self.all_args:
            dt, bs, ma, ro = args
            opt_out_array = make_array(self.shape, *args)
            if ro:
                assert_raises_regex(
                    ValueError, 'output array is read-only',
                    _test_converters, other_array1, None, other_array2,
                    other_array3, opt_out_array, origin)
            else:
                check_converted_output_array(
                    opt_out_array,
                    _test_converters(other_array1, None, other_array2,
                                     other_array3, opt_out_array, origin)[4])
                # Changes to the converted array should be reflected in the
                # original one after deletion of the former.
                assert_(np.all(opt_out_array == 1))
                # In the end, the original array must again be writeable.
                assert_(opt_out_array.flags.writeable)

    def test_origin_and_new_array(self):
        array1 = make_array(self.shape)
        array2 = make_array(self.shape)
        array3 = make_array(self.shape)
        assert_raises(
            TypeError,
            _test_converters, array1, None, array2, array3, None, 'letters')
        assert_raises_regex(
            ValueError, "invalid 2 element 'origin' sequence",
            _test_converters, array1, None, array2, array3, None, (11, 13))
        origin = (17, 19, 23)
        assert_array_equal(
            origin,
            _test_converters(array1, None, array2, array3, None, origin)[5])
        assert_array_equal(
            (0, 0, 0),
            _test_converters(array1, None, array2, array3, None, origin)[6])
