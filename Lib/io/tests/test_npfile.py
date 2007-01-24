from StringIO import StringIO
import os
from tempfile import mkstemp
from numpy.testing import *
import numpy as N

set_package_path()
from io.npfile import npfile, sys_endian_code
restore_path()

class test_npfile(NumpyTestCase):

    def test_init(self):
        fd, fname = mkstemp()
        npf = npfile(fname)
        arr = N.reshape(N.arange(10), (5,2))
        self.assertRaises(IOError, npf.write, arr)
        npf.close()
        npf = npfile(fname, 'w')
        npf.write(arr)
        npf.rewind()
        self.assertRaises(IOError, npf.read,
                          arr.shape,
                          arr.dtype)
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

    def test_raw_read_write(self):
        npf = npfile(StringIO())
        str = 'test me with this string'
        npf.raw_write(str)
        npf.rewind()
        assert str == npf.raw_read(len(str))
        
    def test_read_write(self):
        npf = npfile(StringIO())
        arr = N.reshape(N.arange(10), (5,2))
        f_arr = arr.reshape((2,5)).T
        bo = arr.dtype.byteorder
        swapped_code = sys_endian_code == '<' and '>' or '<'
        if bo in ['=', sys_endian_code]:
            nbo = swapped_code
        else:
            nbo = sys_endian_code
        bs_arr = arr.newbyteorder(nbo)
        adt = arr.dtype
        shp = arr.shape
        npf.write(arr)
        npf.rewind()
        assert_array_equal(npf.read(shp, adt), arr)
        npf.rewind()
        assert_array_equal(npf.read(shp, adt, endian=swapped_code),
                           bs_arr)
        npf.rewind()
        assert_array_equal(npf.read(shp, adt, order='F'),
                           f_arr)
        
        npf = npfile(StringIO(), endian='swapped', order='F')
        npf.write(arr)
        npf.rewind()
        assert_array_equal(npf.read(shp, adt), arr)
        npf.rewind()
        assert_array_equal(npf.read(shp, adt, endian='dtype'), bs_arr)
        npf.rewind()
        # assert_array_equal(npf.read(shp, adt, order='C'), f_arr)
        
