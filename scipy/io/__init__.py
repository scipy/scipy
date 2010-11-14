#
# io - Data input and output
#

from info import __doc__

from numpy import deprecate


# matfile read and write
from matlab import loadmat, savemat

# netCDF file support
from netcdf import netcdf_file, netcdf_variable

from matlab import byteordercodes
from data_store import save_as_module
from mmio import mminfo, mmread, mmwrite
from idl import readsav

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
