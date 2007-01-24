# Author: Matthew Brett, Travis Oliphant

"""
Class for reading and writing numpy arrays from / to files
"""

import sys

import numpy as N

sys_endian_code = (sys.byteorder == 'little') and '<' or '>'

class npfile(object):
    ''' Class for reading and writing numpy arrays to/from files
    
    Inputs:
      file_name -- The complete path name to the file to open
                   or an open file-like object
      permission -- Open the file with given permissions: ('r', 'H', 'a')
                    for reading, writing, or appending.  This is the same
                    as the mode argument in the builtin open command.
      format -- The byte-ordering of the file:
                (['native', 'n'], ['ieee-le', 'l'], ['ieee-be', 'B']) for
                native, little-endian, or big-endian respectively.

    Attributes
      endian   -- default endian code for reading / writing
      order    -- default order for reading writing ('C' or 'F')
      file -- file object containing read / written data
    '''

    def __init__(self, file_name,
                 permission='rb',
                 endian = 'dtype',
                 order = 'C'):
        if 'b' not in permission: permission += 'b'
        if isinstance(file_name, basestring):
            self.file = file(file_name, permission)
        else:
            try:
                closed = file_name.closed
            except AttributeError:
                raise TypeError, 'Need filename or file object as input'
            if closed:
                raise TypeError, 'File object should be open'
            self.file = file_name
            
        self.endian = endian
        self.order = order

    def get_endian(self):
        return self._endian
    def set_endian(self, endian_code):
        self._endian = self.parse_endian(endian_code)
    endian = property(get_endian, set_endian, None, 'get/set endian code')
                                     
    def parse_endian(self, endian_code):
        ''' Returns valid endian code from wider input options'''
        if endian_code in ['native', 'n', 'N','default', '=']:
            return sys_endian_code
        elif endian_code in ['swapped', 's', 'S']:
            return sys_endian_code == '<' and '>' or '<'
        elif endian_code in ['ieee-le','l','L','little-endian',
                             'little','le','<']:
            return '<'
        elif endian_code in ['ieee-be','B','b','big-endian',
                             'big','be', '>']:
            return '>'
        elif endian_code == 'dtype':
            return 'dtype'
        else:
            raise ValueError, "Unrecognized endian code: " + endian_code
        return

    def __del__(self):
        try:
            self.file.close()
        except:
            pass

    def close(self):
        self.file.close()

    def seek(self, *args):
        self.file.seek(*args)

    def tell(self):
        return self.file.tell()

    def rewind(self,howmany=None):
        """Rewind a file to its beginning or by a specified amount.
        """
        if howmany is None:
            self.seek(0)
        else:
            self.seek(-howmany,1)

    def size(self):
        """Return the size of the file.

        Cached once found
        """
        try:
            sz = self.thesize
        except AttributeError:
            curpos = self.tell()
            self.seek(0,2)
            sz = self.tell()
            self.seek(curpos)
            self.thesize = sz
        return sz

    def raw_read(self, size=-1):
        """Read raw bytes from file as string."""
        return self.file.read(size)

    def raw_write(self, str):
        """Write string to file as raw bytes."""
        return self.file.write(str)

    def _endian_order(self, endian, order):
        ''' Housekeeping function to return endian, order from input args '''
        if endian is None:
            endian = self.endian
        else:
            endian = self.parse_endian(endian)
        if order is None:
            order = self.order
        return endian, order

    def _endian_from_dtype(self, dt):
        dt_endian = dt.byteorder
        if dt_endian == '=':
            dt_endian = sys_endian_code
        return dt_endian
    
    def write(self, data, endian=None, order=None):
        ''' Write to open file object the flattened numpy array data

        Inputs
        data      - numpy array or object convertable to array
        endian    - endianness of written data
                    (can be None, 'dtype', '<', '>')
                    (default from self.endian)
        order     - order of array to write (C, F)
                    (default from self.order)
        '''
        endian, order = self._endian_order(endian, order)
        data = N.asarray(data)
        dt_endian = self._endian_from_dtype(data.dtype)
        if not endian == 'dtype':
            if dt_endian != endian:
                data = data.byteswap()
        self.file.write(data.tostring(order=order))
        
    fwrite = write
    
    def read(self, shape, dt, endian=None, order=None):
        ''' Read data from file and return it in a numpy array
        
        Inputs
        shape     - shape of output array, or number of elements
        dt        - dtype of array to write
        endian    - endianness of written data
                    (can be None, 'dtype', '<', '>')
                    (default from self.endian)
        order     - order of array to write (C, F)
                    (default from self.order)
        '''
        endian, order = self._endian_order(endian, order)
        try:
            shape = tuple(shape)
        except TypeError:
            shape = (shape,)
        dt = N.dtype(dt)
        dt_endian = self._endian_from_dtype(dt)
        if not endian == 'dtype':
            if dt_endian != endian:
                dt = dt.newbyteorder(endian)
        sz = dt.itemsize * N.product(shape)
        buf = self.file.read(sz)
        return N.ndarray(shape=shape,
                         dtype=dt,
                         buffer=buf,
                         order=order).copy()
    

    fread = read
