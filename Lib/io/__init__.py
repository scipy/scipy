"""Binary File input and output.

  Classes

     fopen -- a class for easily reading and writing binary data.

  Functions

     read_array -- reading ascii streams into Numeric arrays
     write_array -- write an array to an ascii stream
     loadmat -- read a MATLAB (version <= 4) style Mat file
     savemat -- write a MATLAB (version <= 4) style Mat file

     fread -- low-level reading
     fwrite -- low-level writing
     bswap -- in-place byte-swapping
     packbits -- Pack a binary array of 1's and 0's into an array of bytes
     unpackbits -- Unpack an array packed by packbits.
"""
__all__ = []

_moddict = {'numpyio' : ['packbits', 'unpackbits','bswap', 'fread',
                         'fwrite', 'convert_objectarray'],
            'array_import' : ['read_array',
                              'write_array']
            }


import scipy
scipy.names2all(__all__, ['mio'], globals())
scipy.somenames2all(__all__, _moddict, globals())
del scipy
