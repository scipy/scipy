# Author: Matthew Brett

''' File-like interfact for memmapped array '''

from numpy import *

class ByteStream(object):
    ''' Overlays file-like interface on memmapped array

    This may speed up array reading from files
    
    @byte_array         - uint array or string containing bytes
    '''

    def __init__(self, byte_array):
        if isinstance(byte_array, ndarray):
            if not byte_array.dtype == uint8:
                raise ValueError, 'Need uint8 byte array as array input'
            self.bytes = byte_array
        elif isinstance(byte_array, basestring):
            self.bytes = ndarray(
                shape=(len(byte_array)),
                dtype=uint8,
                buffer=byte_array)
        else:
            raise ValueError, "Need string or byte array as input"
        self.array_len = len(byte_array)
        self.seek(0)

    # current file position
    def get_pos(self):
        return self._pos
    def set_pos(self, offset):
        if offset < 0:
            raise IOError, 'Invalid argument'
        self._pos = offset
    pos = property(get_pos,
                   set_pos,
                   None,
                   'get/set current position')
    
    def seek(self, offset, whence=0):
        """ Method emulates seek method of file objects """
        if whence == 0:
            self.pos = offset
        elif whence == 1: # seek relative to the current position
            self.pos += offset
        elif whence == 2: # relative to end
            self.pos = self.array_len + offset
        else:
            raise ValueError, 'Invalid value %d for whence parameter' % whence

    def tell(self):
        return self.pos
    
    def read(self, num_bytes=-1):
        if num_bytes < 0:
            num_bytes = self.array_len
        if self.pos >= self.array_len:
            return array([], dtype=uint8)
        next_pos = min(self.pos + num_bytes, self.array_len)
        res = self.bytes[self.pos:next_pos]
        self.pos = next_pos
        return res

    def write(self, data):
        assert False, 'Not implemented'
        

