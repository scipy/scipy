"""
Data input and output

SciPy has many modules, classes, and functions available to read data
from and write data to a variety of file formats.

Modules
-------

.. autosummary::
   :toctree: generated/

   arff - Read ARFF files, the standard data format for WEKA
   byteordercodes - System byteorder utilities - NumPy byteorder encoding
   data_store - Load or save values to a file
   dumbdbm_patched - A dumb and slow but simple dbm clone
   matlab - Utilities for dealing with MATLAB(R) files
   mmio - Matrix Market I/O in Python
   netcdf - NetCDF reader/writer module
   wavfile - module to read / write wav files using numpy arrays

Classes
-------

.. autosummary::
   :toctree: generated/

   netcdf_file - A file object for NetCDF data
   netcdf_variable - A data object for the netcdf module

Functions
---------

.. autosummary::
   :toctree: generated/

   loadmat - Read a MATLAB style mat file (version 4 through 7.1)
   savemat - Write a MATLAB style mat file (version 4 through 7.1)
   mminfo - Query matrix info from Matrix Market formatted file
   mmread - Read matrix from Matrix Market formatted file
   mmwrite - Write matrix to Matrix Market formatted file
   save_as_module - Data saved as module, accessed on load as attirbutes

"""
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
