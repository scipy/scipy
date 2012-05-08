from os import path
import warnings

DATA_PATH = path.join(path.dirname(__file__), 'data')

import numpy as np
from numpy.compat import asbytes_nested, asbytes
from numpy.testing import assert_equal, assert_array_equal, run_module_suite
from numpy.testing.utils import WarningManager
from nose.tools import assert_true

from scipy.io.idl import readsav


def object_array(*args):
    '''Constructs a numpy array of objects'''
    array = np.empty(len(args), dtype=np.object)
    for i in range(len(args)):
        array[i] = args[i]
    return array


def assert_identical(a, b):
    '''Assert whether value AND type are the same'''
    assert_equal(a, b)
    if type(b) is np.str:
        assert_equal(type(a), type(b))
    else:
        assert_equal(np.asarray(a).dtype.type, np.asarray(b).dtype.type)


def assert_array_identical(a, b):
    '''Assert whether values AND type are the same'''
    assert_array_equal(a, b)
    assert_equal(a.dtype.type, b.dtype.type)


# Define vectorized ID function for pointer arrays
vect_id = np.vectorize(id)


class TestIdict:
    '''Test the idict= argument to read'''

    def test_idict(self):
        custom_dict = {'a': np.int16(999)}
        original_id = id(custom_dict)
        s = readsav(path.join(DATA_PATH, 'scalar_byte.sav'), idict=custom_dict, verbose=False)
        assert_equal(original_id, id(s))
        assert_true('a' in s)
        assert_identical(s['a'], np.int16(999))
        assert_identical(s['i8u'], np.uint8(234))


class TestScalars:
    '''Test that scalar values are read in with the correct value and type'''

    def test_byte(self):
        s = readsav(path.join(DATA_PATH, 'scalar_byte.sav'), verbose=False)
        assert_identical(s.i8u, np.uint8(234))

    def test_int16(self):
        s = readsav(path.join(DATA_PATH, 'scalar_int16.sav'), verbose=False)
        assert_identical(s.i16s, np.int16(-23456))

    def test_int32(self):
        s = readsav(path.join(DATA_PATH, 'scalar_int32.sav'), verbose=False)
        assert_identical(s.i32s, np.int32(-1234567890))

    def test_float32(self):
        s = readsav(path.join(DATA_PATH, 'scalar_float32.sav'), verbose=False)
        assert_identical(s.f32, np.float32(-3.1234567e+37))

    def test_float64(self):
        s = readsav(path.join(DATA_PATH, 'scalar_float64.sav'), verbose=False)
        assert_identical(s.f64, np.float64(-1.1976931348623157e+307))

    def test_complex32(self):
        s = readsav(path.join(DATA_PATH, 'scalar_complex32.sav'), verbose=False)
        assert_identical(s.c32, np.complex64(3.124442e13-2.312442e31j))

    def test_bytes(self):
        s = readsav(path.join(DATA_PATH, 'scalar_string.sav'), verbose=False)
        assert_identical(s.s, np.bytes_("The quick brown fox jumps over the lazy python"))

    def test_structure(self):
        pass

    def test_complex64(self):
        s = readsav(path.join(DATA_PATH, 'scalar_complex64.sav'), verbose=False)
        assert_identical(s.c64, np.complex128(1.1987253647623157e+112-5.1987258887729157e+307j))

    def test_heap_pointer(self):
        pass

    def test_object_reference(self):
        pass

    def test_uint16(self):
        s = readsav(path.join(DATA_PATH, 'scalar_uint16.sav'), verbose=False)
        assert_identical(s.i16u, np.uint16(65511))

    def test_uint32(self):
        s = readsav(path.join(DATA_PATH, 'scalar_uint32.sav'), verbose=False)
        assert_identical(s.i32u, np.uint32(4294967233))

    def test_int64(self):
        s = readsav(path.join(DATA_PATH, 'scalar_int64.sav'), verbose=False)
        assert_identical(s.i64s, np.int64(-9223372036854774567))

    def test_uint64(self):
        s = readsav(path.join(DATA_PATH, 'scalar_uint64.sav'), verbose=False)
        assert_identical(s.i64u, np.uint64(18446744073709529285))


class TestCompressed(TestScalars):
    '''Test that compressed .sav files can be read in'''

    def test_compressed(self):
        warn_ctx = WarningManager()
        warn_ctx.__enter__()
        try:
            warnings.filterwarnings('ignore', message="warning: empty strings")
            s = readsav(path.join(DATA_PATH, 'various_compressed.sav'), verbose=False)
        finally:
            warn_ctx.__exit__()

        assert_identical(s.i8u, np.uint8(234))
        assert_identical(s.f32, np.float32(-3.1234567e+37))
        assert_identical(s.c64, np.complex128(1.1987253647623157e+112-5.1987258887729157e+307j))
        assert_equal(s.array5d.shape, (4, 3, 4, 6, 5))
        assert_identical(s.arrays.a[0], np.array([1, 2, 3], dtype=np.int16))
        assert_identical(s.arrays.b[0], np.array([4., 5., 6., 7.], dtype=np.float32))
        assert_identical(s.arrays.c[0], np.array([np.complex64(1+2j), np.complex64(7+8j)]))
        assert_identical(s.arrays.d[0], np.array(asbytes_nested(["cheese", "bacon", "spam"]), dtype=np.object))


class TestArrayDimensions:
    '''Test that multi-dimensional arrays are read in with the correct dimensions'''

    def test_1d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_1d.sav'), verbose=False)
        assert_equal(s.array1d.shape, (123, ))

    def test_2d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_2d.sav'), verbose=False)
        assert_equal(s.array2d.shape, (22, 12))

    def test_3d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_3d.sav'), verbose=False)
        assert_equal(s.array3d.shape, (11, 22, 12))

    def test_4d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_4d.sav'), verbose=False)
        assert_equal(s.array4d.shape, (4, 5, 8, 7))

    def test_5d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_5d.sav'), verbose=False)
        assert_equal(s.array5d.shape, (4, 3, 4, 6, 5))

    def test_6d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_6d.sav'), verbose=False)
        assert_equal(s.array6d.shape, (3, 6, 4, 5, 3, 4))

    def test_7d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_7d.sav'), verbose=False)
        assert_equal(s.array7d.shape, (2, 1, 2, 3, 4, 3, 2))

    def test_8d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_8d.sav'), verbose=False)
        assert_equal(s.array8d.shape, (4, 3, 2, 1, 2, 3, 5, 4))


class TestStructures:
    '''Test that structures are correctly read in'''

    def test_scalars(self):
        s = readsav(path.join(DATA_PATH, 'struct_scalars.sav'), verbose=False)
        assert_identical(s.scalars.a, np.array(np.int16(1)))
        assert_identical(s.scalars.b, np.array(np.int32(2)))
        assert_identical(s.scalars.c, np.array(np.float32(3.)))
        assert_identical(s.scalars.d, np.array(np.float64(4.)))
        assert_identical(s.scalars.e, np.array(asbytes_nested(["spam"]), dtype=np.object))
        assert_identical(s.scalars.f, np.array(np.complex64(-1.+3j)))

    def test_scalars_replicated(self):
        s = readsav(path.join(DATA_PATH, 'struct_scalars_replicated.sav'), verbose=False)
        assert_identical(s.scalars_rep.a, np.repeat(np.int16(1), 5))
        assert_identical(s.scalars_rep.b, np.repeat(np.int32(2), 5))
        assert_identical(s.scalars_rep.c, np.repeat(np.float32(3.), 5))
        assert_identical(s.scalars_rep.d, np.repeat(np.float64(4.), 5))
        assert_identical(s.scalars_rep.e, np.repeat(asbytes("spam"), 5).astype(np.object))
        assert_identical(s.scalars_rep.f, np.repeat(np.complex64(-1.+3j), 5))

    def test_scalars_replicated_3d(self):
        s = readsav(path.join(DATA_PATH, 'struct_scalars_replicated_3d.sav'), verbose=False)
        assert_identical(s.scalars_rep.a, np.repeat(np.int16(1), 24).reshape(4, 3, 2))
        assert_identical(s.scalars_rep.b, np.repeat(np.int32(2), 24).reshape(4, 3, 2))
        assert_identical(s.scalars_rep.c, np.repeat(np.float32(3.), 24).reshape(4, 3, 2))
        assert_identical(s.scalars_rep.d, np.repeat(np.float64(4.), 24).reshape(4, 3, 2))
        assert_identical(s.scalars_rep.e, np.repeat(asbytes("spam"), 24).reshape(4, 3, 2).astype(np.object))
        assert_identical(s.scalars_rep.f, np.repeat(np.complex64(-1.+3j), 24).reshape(4, 3, 2))

    def test_arrays(self):
        s = readsav(path.join(DATA_PATH, 'struct_arrays.sav'), verbose=False)
        assert_array_identical(s.arrays.a[0], np.array([1, 2, 3], dtype=np.int16))
        assert_array_identical(s.arrays.b[0], np.array([4., 5., 6., 7.], dtype=np.float32))
        assert_array_identical(s.arrays.c[0], np.array([np.complex64(1+2j), np.complex64(7+8j)]))
        assert_array_identical(s.arrays.d[0], np.array(asbytes_nested(["cheese", "bacon", "spam"]), dtype=np.object))

    def test_arrays_replicated(self):

        s = readsav(path.join(DATA_PATH, 'struct_arrays_replicated.sav'), verbose=False)

        # Check column types
        assert_true(s.arrays_rep.a.dtype.type is np.object_)
        assert_true(s.arrays_rep.b.dtype.type is np.object_)
        assert_true(s.arrays_rep.c.dtype.type is np.object_)
        assert_true(s.arrays_rep.d.dtype.type is np.object_)

        # Check column shapes
        assert_equal(s.arrays_rep.a.shape, (5, ))
        assert_equal(s.arrays_rep.b.shape, (5, ))
        assert_equal(s.arrays_rep.c.shape, (5, ))
        assert_equal(s.arrays_rep.d.shape, (5, ))

        # Check values
        for i in range(5):
            assert_array_identical(s.arrays_rep.a[i], np.array([1, 2, 3], dtype=np.int16))
            assert_array_identical(s.arrays_rep.b[i], np.array([4., 5., 6., 7.], dtype=np.float32))
            assert_array_identical(s.arrays_rep.c[i], np.array([np.complex64(1+2j), np.complex64(7+8j)]))
            assert_array_identical(s.arrays_rep.d[i], np.array(asbytes_nested(["cheese", "bacon", "spam"]), dtype=np.object))

    def test_arrays_replicated_3d(self):

        s = readsav(path.join(DATA_PATH, 'struct_arrays_replicated_3d.sav'), verbose=False)

        # Check column types
        assert_true(s.arrays_rep.a.dtype.type is np.object_)
        assert_true(s.arrays_rep.b.dtype.type is np.object_)
        assert_true(s.arrays_rep.c.dtype.type is np.object_)
        assert_true(s.arrays_rep.d.dtype.type is np.object_)

        # Check column shapes
        assert_equal(s.arrays_rep.a.shape, (4, 3, 2))
        assert_equal(s.arrays_rep.b.shape, (4, 3, 2))
        assert_equal(s.arrays_rep.c.shape, (4, 3, 2))
        assert_equal(s.arrays_rep.d.shape, (4, 3, 2))

        # Check values
        for i in range(4):
            for j in range(3):
                for k in range(2):
                    assert_array_identical(s.arrays_rep.a[i, j, k], np.array([1, 2, 3], dtype=np.int16))
                    assert_array_identical(s.arrays_rep.b[i, j, k], np.array([4., 5., 6., 7.], dtype=np.float32))
                    assert_array_identical(s.arrays_rep.c[i, j, k], np.array([np.complex64(1+2j), np.complex64(7+8j)]))
                    assert_array_identical(s.arrays_rep.d[i, j, k], np.array(asbytes_nested(["cheese", "bacon", "spam"]), dtype=np.object))

    def test_inheritance(self):
        s = readsav(path.join(DATA_PATH, 'struct_inherit.sav'), verbose=False)
        assert_identical(s.fc.x, np.array([0], dtype=np.int16))
        assert_identical(s.fc.y, np.array([0], dtype=np.int16))
        assert_identical(s.fc.r, np.array([0], dtype=np.int16))
        assert_identical(s.fc.c, np.array([4], dtype=np.int16))


class TestPointers:
    '''Check that pointers in .sav files produce references to the same object in Python'''

    def test_pointers(self):
        s = readsav(path.join(DATA_PATH, 'scalar_heap_pointer.sav'), verbose=False)
        assert_identical(s.c64_pointer1, np.complex128(1.1987253647623157e+112-5.1987258887729157e+307j))
        assert_identical(s.c64_pointer2, np.complex128(1.1987253647623157e+112-5.1987258887729157e+307j))
        assert_true(s.c64_pointer1 is s.c64_pointer2)


class TestPointerArray:
    '''Test that pointers in arrays are correctly read in'''

    def test_1d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_pointer_1d.sav'), verbose=False)
        assert_equal(s.array1d.shape, (123, ))
        assert_true(np.all(s.array1d == np.float32(4.)))
        assert_true(np.all(vect_id(s.array1d) == id(s.array1d[0])))

    def test_2d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_pointer_2d.sav'), verbose=False)
        assert_equal(s.array2d.shape, (22, 12))
        assert_true(np.all(s.array2d == np.float32(4.)))
        assert_true(np.all(vect_id(s.array2d) == id(s.array2d[0,0])))

    def test_3d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_pointer_3d.sav'), verbose=False)
        assert_equal(s.array3d.shape, (11, 22, 12))
        assert_true(np.all(s.array3d == np.float32(4.)))
        assert_true(np.all(vect_id(s.array3d) == id(s.array3d[0,0,0])))

    def test_4d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_pointer_4d.sav'), verbose=False)
        assert_equal(s.array4d.shape, (4, 5, 8, 7))
        assert_true(np.all(s.array4d == np.float32(4.)))
        assert_true(np.all(vect_id(s.array4d) == id(s.array4d[0,0,0,0])))

    def test_5d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_pointer_5d.sav'), verbose=False)
        assert_equal(s.array5d.shape, (4, 3, 4, 6, 5))
        assert_true(np.all(s.array5d == np.float32(4.)))
        assert_true(np.all(vect_id(s.array5d) == id(s.array5d[0,0,0,0,0])))

    def test_6d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_pointer_6d.sav'), verbose=False)
        assert_equal(s.array6d.shape, (3, 6, 4, 5, 3, 4))
        assert_true(np.all(s.array6d == np.float32(4.)))
        assert_true(np.all(vect_id(s.array6d) == id(s.array6d[0,0,0,0,0,0])))

    def test_7d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_pointer_7d.sav'), verbose=False)
        assert_equal(s.array7d.shape, (2, 1, 2, 3, 4, 3, 2))
        assert_true(np.all(s.array7d == np.float32(4.)))
        assert_true(np.all(vect_id(s.array7d) == id(s.array7d[0,0,0,0,0,0,0])))

    def test_8d(self):
        s = readsav(path.join(DATA_PATH, 'array_float32_pointer_8d.sav'), verbose=False)
        assert_equal(s.array8d.shape, (4, 3, 2, 1, 2, 3, 5, 4))
        assert_true(np.all(s.array8d == np.float32(4.)))
        assert_true(np.all(vect_id(s.array8d) == id(s.array8d[0,0,0,0,0,0,0,0])))

class TestPointerStructures:
    '''Test that structures are correctly read in'''

    def test_scalars(self):
        s = readsav(path.join(DATA_PATH, 'struct_pointers.sav'), verbose=False)
        assert_identical(s.pointers.g, np.array(np.float32(4.), dtype=np.object_))
        assert_identical(s.pointers.h, np.array(np.float32(4.), dtype=np.object_))
        assert_true(id(s.pointers.g[0]) == id(s.pointers.h[0]))

    def test_pointers_replicated(self):
        s = readsav(path.join(DATA_PATH, 'struct_pointers_replicated.sav'), verbose=False)
        assert_identical(s.pointers_rep.g, np.repeat(np.float32(4.), 5).astype(np.object_))
        assert_identical(s.pointers_rep.h, np.repeat(np.float32(4.), 5).astype(np.object_))
        assert_true(np.all(vect_id(s.pointers_rep.g) == vect_id(s.pointers_rep.h)))

    def test_pointers_replicated_3d(self):
        s = readsav(path.join(DATA_PATH, 'struct_pointers_replicated_3d.sav'), verbose=False)
        assert_identical(s.pointers_rep.g, np.repeat(np.float32(4.), 24).reshape(4, 3, 2).astype(np.object_))
        assert_identical(s.pointers_rep.h, np.repeat(np.float32(4.), 24).reshape(4, 3, 2).astype(np.object_))
        assert_true(np.all(vect_id(s.pointers_rep.g) == vect_id(s.pointers_rep.h)))

    def test_arrays(self):
        s = readsav(path.join(DATA_PATH, 'struct_pointer_arrays.sav'), verbose=False)
        assert_array_identical(s.arrays.g[0], np.repeat(np.float32(4.), 2).astype(np.object_))
        assert_array_identical(s.arrays.h[0], np.repeat(np.float32(4.), 3).astype(np.object_))
        assert_true(np.all(vect_id(s.arrays.g[0]) == id(s.arrays.g[0][0])))
        assert_true(np.all(vect_id(s.arrays.h[0]) == id(s.arrays.h[0][0])))
        assert_true(id(s.arrays.g[0][0]) == id(s.arrays.h[0][0]))

    def test_arrays_replicated(self):

        s = readsav(path.join(DATA_PATH, 'struct_pointer_arrays_replicated.sav'), verbose=False)

        # Check column types
        assert_true(s.arrays_rep.g.dtype.type is np.object_)
        assert_true(s.arrays_rep.h.dtype.type is np.object_)

        # Check column shapes
        assert_equal(s.arrays_rep.g.shape, (5, ))
        assert_equal(s.arrays_rep.h.shape, (5, ))

        # Check values
        for i in range(5):
            assert_array_identical(s.arrays_rep.g[i], np.repeat(np.float32(4.), 2).astype(np.object_))
            assert_array_identical(s.arrays_rep.h[i], np.repeat(np.float32(4.), 3).astype(np.object_))
            assert_true(np.all(vect_id(s.arrays_rep.g[i]) == id(s.arrays_rep.g[0][0])))
            assert_true(np.all(vect_id(s.arrays_rep.h[i]) == id(s.arrays_rep.h[0][0])))

    def test_arrays_replicated_3d(self):
        warn_ctx = WarningManager()
        warn_ctx.__enter__()
        try:
            warnings.filterwarnings('ignore', message="warning: multi-dimensional structures")
            s = readsav(path.join(DATA_PATH, 'struct_pointer_arrays_replicated_3d.sav'), verbose=False)
        finally:
            warn_ctx.__exit__()

        # Check column types
        assert_true(s.arrays_rep.g.dtype.type is np.object_)
        assert_true(s.arrays_rep.h.dtype.type is np.object_)

        # Check column shapes
        assert_equal(s.arrays_rep.g.shape, (4, 3, 2))
        assert_equal(s.arrays_rep.h.shape, (4, 3, 2))

        # Check values
        for i in range(4):
            for j in range(3):
                for k in range(2):
                    assert_array_identical(s.arrays_rep.g[i, j, k], np.repeat(np.float32(4.), 2).astype(np.object_))
                    assert_array_identical(s.arrays_rep.h[i, j, k], np.repeat(np.float32(4.), 3).astype(np.object_))
                    assert_true(np.all(vect_id(s.arrays_rep.g[i, j, k]) == id(s.arrays_rep.g[0, 0, 0][0])))
                    assert_true(np.all(vect_id(s.arrays_rep.h[i, j, k]) == id(s.arrays_rep.h[0, 0, 0][0])))


if __name__ == "__main__":
    run_module_suite()
