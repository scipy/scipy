"""
Data input and output
=====================

  Classes

     fopen -- a class for easily reading and writing binary data.

  Functions

     read_array -- reading ascii streams into Numeric arrays
     write_array -- write an array to an ascii stream
     loadmat -- read a MATLAB (version <= 4) style mat file
     savemat -- write a MATLAB (version <= 4) style mat file

     fread -- low-level reading
     fwrite -- low-level writing
     bswap -- in-place byte-swapping
     packbits -- Pack a binary array of 1's and 0's into an array of bytes
     unpackbits -- Unpack an array packed by packbits.

     save --- simple storing of Python dictionary into module
               that can then be imported and the data accessed as
               attributes of the module.

"""
