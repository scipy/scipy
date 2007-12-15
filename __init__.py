#
# io - Data input and output
#

from info import __doc__

from numpyio import packbits, unpackbits, bswap, fread, fwrite, \
     convert_objectarray

# matfile read and write
from matlab.mio import loadmat, savemat

# netCDF file support
from netcdf import netcdf_file, netcdf_variable

from npfile import npfile

from recaster import sctype_attributes, Recaster

from array_import import read_array, write_array
from data_store import save as save_as_module
from data_store import load, create_module, create_shelf
from pickler import objload, objsave

from mmio import mminfo, mmread, mmwrite

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import NumpyTest
test = NumpyTest().test
