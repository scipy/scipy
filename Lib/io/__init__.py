"""Binary File input and output.

  Classes

     fopen -- a class for easily reading and writing binary data.

  Functions

     fread -- low-level reading
     fwrite -- low-level writing
     bswap -- in-place byte-swapping
     packbits -- Pack a binary array of 1's and 0's into an array of bytes
     unpackbits -- Unpack an array packed by packbits.

"""
from mio import *
from numpyio import packbits, unpackbits, bswap, fread, fwrite
