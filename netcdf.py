"""NetCDF file reader.

This is adapted from Roberto De Almeida's Pupynere PUre PYthon NEtcdf REader.

classes changed to underscore_separated instead of CamelCase

TODO:

  Add write capability.
"""

#__author__ = "Roberto De Almeida <rob@pydap.org>"


__all__ = ['netcdf_file', 'netcdf_variable']

import struct
import mmap

from numpy import ndarray, zeros, array


ABSENT       = '\x00' * 8
ZERO         = '\x00' * 4
NC_BYTE      = '\x00\x00\x00\x01'
NC_CHAR      = '\x00\x00\x00\x02'
NC_SHORT     = '\x00\x00\x00\x03'
NC_INT       = '\x00\x00\x00\x04'
NC_FLOAT     = '\x00\x00\x00\x05'
NC_DOUBLE    = '\x00\x00\x00\x06'
NC_DIMENSION = '\x00\x00\x00\n'
NC_VARIABLE  = '\x00\x00\x00\x0b'
NC_ATTRIBUTE = '\x00\x00\x00\x0c'


class netcdf_file(object):
    """A NetCDF file parser."""

    def __init__(self, file, mode):
        mode += 'b'
        self._buffer = open(file, mode)
        if mode in ['rb', 'r+b']:
            self._parse()
        elif mode == 'ab':
            raise NotImplementedError

    def flush(self):
        pass

    def sync(self):
        pass

    def close(self):
        pass

    def create_dimension(self, name, length):
        pass

    def create_variable(self, name, type, dimensions):
        pass

    def read(self, size=-1):
        """Alias for reading the file buffer."""
        return self._buffer.read(size)

    def _parse(self):
        """Initial parsing of the header."""
        # Check magic bytes.
        assert self.read(3) == 'CDF'

        # Read version byte.
        byte = self.read(1)
        self.version_byte = struct.unpack('>b', byte)[0]

        # Read header info.
        self._numrecs()
        self._dim_array()
        self._gatt_array()
        self._var_array()

    def _numrecs(self):
        """Read number of records."""
        self._nrecs = self._unpack_int()

    def _dim_array(self):
        """Read a dict with dimensions names and sizes."""
        assert self.read(4) in [ZERO, NC_DIMENSION]
        count = self._unpack_int()

        self.dimensions = {}
        self._dims = []
        for dim in range(count):
            name = self._read_string()
            length = self._unpack_int()
            if length == 0: length = None # record dimension
            self.dimensions[name] = length
            self._dims.append(name)  # preserve dim order

    def _gatt_array(self):
        """Read global attributes."""
        self.attributes = self._att_array()

        # Update __dict__ for compatibility with S.IO.N
        self.__dict__.update(self.attributes)

    def _att_array(self):
        """Read a dict with attributes."""
        assert self.read(4) in [ZERO, NC_ATTRIBUTE]
        count = self._unpack_int()

        # Read attributes.
        attributes = {}
        for attribute in range(count):
            name = self._read_string()
            nc_type = self._unpack_int()
            n = self._unpack_int()

            # Read value for attributes.
            attributes[name] = self._read_values(n, nc_type)

        return attributes

    def _var_array(self):
        """Read all variables."""
        assert self.read(4) in [ZERO, NC_VARIABLE]

        # Read size of each record, in bytes.
        self._read_recsize()

        # Read variables.
        self.variables = {}
        count = self._unpack_int()
        for variable in range(count):
            name = self._read_string()
            self.variables[name] = self._read_var()

    def _read_recsize(self):
        """Read all variables and compute record bytes."""
        pos = self._buffer.tell()

        recsize = 0
        count = self._unpack_int()
        for variable in range(count):
            name = self._read_string()
            n = self._unpack_int()
            isrec = False
            for i in range(n):
                dimid = self._unpack_int()
                name = self._dims[dimid]
                dim = self.dimensions[name]
                if dim is None and i == 0:
                    isrec = True
            attributes = self._att_array()
            nc_type = self._unpack_int()
            vsize = self._unpack_int()
            begin = [self._unpack_int, self._unpack_int64][self.version_byte-1]()

            if isrec: recsize += vsize

        self._recsize = recsize
        self._buffer.seek(pos)

    def _read_var(self):
        dimensions = []
        shape = []
        n = self._unpack_int()
        isrec = False
        for i in range(n):
            dimid = self._unpack_int()
            name = self._dims[dimid]
            dimensions.append(name)
            dim = self.dimensions[name]
            if dim is None and i == 0:
                dim = self._nrecs
                isrec = True
            shape.append(dim)
        dimensions = tuple(dimensions)
        shape = tuple(shape)

        attributes = self._att_array()
        nc_type = self._unpack_int()
        vsize = self._unpack_int()

        # Read offset.
        begin = [self._unpack_int, self._unpack_int64][self.version_byte-1]()

        return netcdf_variable(self._buffer.fileno(), nc_type, vsize, begin, shape, dimensions, attributes, isrec, self._recsize)

    def _read_values(self, n, nc_type):
        bytes = [1, 1, 2, 4, 4, 8]
        typecodes = ['b', 'c', 'h', 'i', 'f', 'd']

        count = n * bytes[nc_type-1]
        values = self.read(count)
        padding = self.read((4 - (count % 4)) % 4)

        typecode = typecodes[nc_type-1]
        if nc_type != 2:  # not char
            values = struct.unpack('>%s' % (typecode * n), values)
            values = array(values, dtype=typecode)
        else:
            # Remove EOL terminator.
            if values.endswith('\x00'): values = values[:-1]

        return values

    def _unpack_int(self):
        return struct.unpack('>i', self.read(4))[0]
    _unpack_int32 = _unpack_int

    def _unpack_int64(self):
        return struct.unpack('>q', self.read(8))[0]

    def _read_string(self):
        count = struct.unpack('>i', self.read(4))[0]
        s = self.read(count)
        # Remove EOL terminator.
        if s.endswith('\x00'): s = s[:-1]
        padding = self.read((4 - (count % 4)) % 4)
        return s

    def close(self):
        self._buffer.close()


class netcdf_variable(object):
    def __init__(self, fileno, nc_type, vsize, begin, shape, dimensions, attributes, isrec=False, recsize=0):
        self._nc_type = nc_type
        self._vsize = vsize
        self._begin = begin
        self.shape = shape
        self.dimensions = dimensions
        self.attributes = attributes  # for ``dap.plugins.netcdf``
        self.__dict__.update(attributes)
        self._is_record = isrec

        # Number of bytes and type.
        self._bytes = [1, 1, 2, 4, 4, 8][self._nc_type-1]
        type_ = ['i', 'S', 'i', 'i', 'f', 'f'][self._nc_type-1]
        dtype = '>%s%d' % (type_, self._bytes)
        bytes = self._begin + self._vsize

        if isrec:
            # Record variables are not stored contiguosly on disk, so we
            # need to create a separate array for each record.
            #
            # TEO:  This will copy data from the newly-created array
            #  into the __array_data__ region, thus removing any benefit of using
            #  a memory-mapped file.  You might as well just read the data
            #  in directly.
            self.__array_data__ = zeros(shape, dtype)
            bytes += (shape[0] - 1) * recsize
            for n in range(shape[0]):
                offset = self._begin + (n * recsize)
                mm = mmap.mmap(fileno, bytes, access=mmap.ACCESS_READ)
                self.__array_data__[n] = ndarray.__new__(ndarray, shape[1:], dtype=dtype, buffer=mm, offset=offset, order=0)
        else:
            # Create buffer and data.
            mm = mmap.mmap(fileno, bytes, access=mmap.ACCESS_READ)
            self.__array_data__ = ndarray.__new__(ndarray, shape, dtype=dtype, buffer=mm, offset=self._begin, order=0)

        # N-D array interface
        self.__array_interface__ = {'shape'  : shape,
                                    'typestr': dtype,
                                    'data'   : self.__array_data__,
                                    'version': 3,
                                   }

    def __getitem__(self, index):
        return self.__array_data__.__getitem__(index)

    def getValue(self):
        """For scalars."""
        return self.__array_data__.item()

    def assignValue(self, value):
        """For scalars."""
        self.__array_data__.itemset(value)

    def typecode(self):
        return ['b', 'c', 'h', 'i', 'f', 'd'][self._nc_type-1]


def _test():
    import doctest
    doctest.testmod()
