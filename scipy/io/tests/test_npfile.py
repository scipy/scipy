import os
from StringIO import StringIO
from tempfile import mkstemp
from numpy.testing import *
import numpy as np

from scipy.io.npfile import npfile, sys_endian_code

class TestNpFile(TestCase):

    def test_init(self):
        fd, fname = mkstemp()
        os.close(fd)
        npf = npfile(fname)
        arr = np.reshape(np.arange(10), (5,2))
        self.assertRaises(IOError, npf.write_array, arr)
        npf.close()
        npf = npfile(fname, 'w')
        npf.write_array(arr)
        npf.rewind()
        self.assertRaises(IOError, npf.read_array,
                          arr.dtype, arr.shape)
        npf.close()
        os.remove(fname)

        npf = npfile(StringIO(), endian='>', order='F')
        assert npf.endian == '>', 'Endian not set correctly'
        assert npf.order == 'F', 'Order not set correctly'
        npf.endian = '<'
        assert npf.endian == '<', 'Endian not set correctly'

    def test_parse_endian(self):
        npf = npfile(StringIO())
        swapped_code = sys_endian_code == '<' and '>' or '<'
        assert npf.parse_endian('native') == sys_endian_code
        assert npf.parse_endian('swapped') == swapped_code
        assert npf.parse_endian('l') == '<'
        assert npf.parse_endian('B') == '>'
        assert npf.parse_endian('dtype') == 'dtype'
        self.assertRaises(ValueError, npf.parse_endian, 'nonsense')

    def test_read_write_raw(self):
        npf = npfile(StringIO())
        str = 'test me with this string'
        npf.write_raw(str)
        npf.rewind()
        assert str == npf.read_raw(len(str))

    def test_remaining_bytes(self):
        npf = npfile(StringIO())
        assert npf.remaining_bytes() == 0
        npf.write_raw('+' * 10)
        assert npf.remaining_bytes() == 0
        npf.rewind()
        assert npf.remaining_bytes() == 10
        npf.seek(5)
        assert npf.remaining_bytes() == 5

    def test_read_write_array(self):
        npf = npfile(StringIO())
        arr = np.reshape(np.arange(10), (5,2))
        # Arr as read in fortran order
        f_arr = arr.reshape((2,5)).T
        # Arr written in fortran order read in C order
        cf_arr = arr.T.reshape((5,2))
        # Byteswapped array
        bo = arr.dtype.byteorder
        swapped_code = sys_endian_code == '<' and '>' or '<'
        if bo in ['=', sys_endian_code]:
            nbo = swapped_code
        else:
            nbo = sys_endian_code
        bs_arr = arr.newbyteorder(nbo)
        adt = arr.dtype
        shp = arr.shape
        npf.write_array(arr)
        npf.rewind()
        assert_array_equal(npf.read_array(adt), arr.flatten())
        npf.rewind()
        assert_array_equal(npf.read_array(adt, shp), arr)
        npf.rewind()
        assert_array_equal(npf.read_array(adt, shp, endian=swapped_code),
                           bs_arr)
        npf.rewind()
        assert_array_equal(npf.read_array(adt, shp, order='F'),
                           f_arr)
        npf.rewind()
        npf.write_array(arr, order='F')
        npf.rewind()
        assert_array_equal(npf.read_array(adt), arr.flatten('F'))
        npf.rewind()
        assert_array_equal(npf.read_array(adt, shp),
                           cf_arr)

        npf = npfile(StringIO(), endian='swapped', order='F')
        npf.write_array(arr)
        npf.rewind()
        assert_array_equal(npf.read_array(adt, shp), arr)
        npf.rewind()
        assert_array_equal(npf.read_array(adt, shp, endian='dtype'), bs_arr)
        npf.rewind()
        assert_array_equal(npf.read_array(adt, shp, order='C'), cf_arr)

if __name__ == "__main__":
    run_module_suite()
