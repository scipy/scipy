#
# io - Data input and output
#

from info import __doc__

from numpy import deprecate_with_doc

# These are all deprecated (until the end deprecated tag)
from npfile import npfile
from data_store import save, load, create_module, create_shelf
from array_import import read_array, write_array
from pickler import objload, objsave

from numpyio import packbits, unpackbits, bswap, fread, fwrite, \
     convert_objectarray

fread = deprecate_with_doc("""
scipy.io.fread is can be replaced with raw reading capabilities of NumPy
including fromfile as well as memory-mapping capabilities.
""")(fread)

fwrite = deprecate_with_doc("""
scipy.io.fwrite can be replaced with raw writing capabilities of
NumPy.  Also, remember that files can be directly memory-mapped into NumPy
arrays which is often a better way of reading especially large files.

Look at the tofile methods as well as save and savez for writing arrays into
easily transported files of data.
""")(fwrite)

bswap = deprecate_with_doc("""
scipy.io.bswap is easily replaced with the byteswap method on an array.
out = scipy.io.bswap(arr) --> out = arr.byteswap(True)
""")(bswap)

packbits = deprecate_with_doc("""
The functionality of scipy.io.packbits is now available as numpy.packbits
The calling convention is a bit different as the 2-d case is not specialized.

However, you can simulate scipy.packbits by raveling the last 2 dimensions
of the array and calling numpy.packbits with an axis=-1 keyword:

def scipy_packbits(inp):
    a = np.asarray(inp)
    if a.ndim < 2:
       return np.packbits(a)
    oldshape = a.shape
    newshape = oldshape[:-2] + (oldshape[-2]*oldshape[-1],)
    a = np.reshape(a, newshape)
    return np.packbits(a, axis=-1).ravel()
""")(packbits)

unpackbits = deprecate_with_doc("""
The functionality of scipy.io.unpackbits is now available in numpy.unpackbits
The calling convention is different however as the 2-d case is no longer
specialized.

Thus, the scipy.unpackbits behavior must be simulated using numpy.unpackbits.

def scipy_unpackbits(inp, els_per_slice, out_type=None):
    inp = np.asarray(inp)
    num4els = ((els_per_slice-1) >> 3) + 1
    inp = np.reshape(inp, (-1,num4els))
    res = np.unpackbits(inp, axis=-1)[:,:els_per_slice]
    return res.ravel()
""")(unpackbits)

convert_objectarray = deprecate_with_doc("""
The same functionality can be obtained using NumPy string arrays and the
.astype method (except for the optional missing value feature).
""")(convert_objectarray)

# end deprecated

# matfile read and write
from matlab.mio import loadmat, savemat

# netCDF file support
from netcdf import netcdf_file, netcdf_variable

from recaster import sctype_attributes, Recaster
import matlab.byteordercodes as byteordercodes
from data_store import save_as_module
from mmio import mminfo, mmread, mmwrite

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
