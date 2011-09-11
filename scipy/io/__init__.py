# -*- encoding:utf-8 -*-
"""
==================================
Input and output (:mod:`scipy.io`)
==================================

.. currentmodule:: scipy.io

SciPy has many modules, classes, and functions available to read data
from and write data to a variety of file formats.

.. seealso:: :ref:`numpy-reference.routines.io` (in Numpy)

MATLAB® files
=============

.. autosummary::
   :toctree: generated/

   loadmat - Read a MATLAB style mat file (version 4 through 7.1)
   savemat - Write a MATLAB style mat file (version 4 through 7.1)

IDL® files
==========

.. autosummary::
   :toctree: generated/

   readsav - Read an IDL 'save' file

Matrix Market files
===================

.. autosummary::
   :toctree: generated/

   mminfo - Query matrix info from Matrix Market formatted file
   mmread - Read matrix from Matrix Market formatted file
   mmwrite - Write matrix to Matrix Market formatted file

Other
=====

.. autosummary::
   :toctree: generated/

   save_as_module - Data saved as module, accessed on load as attirbutes

Wav sound files (:mod:`scipy.io.wavfile`)
=========================================

.. module:: scipy.io.wavfile

.. autosummary::
   :toctree: generated/

   read
   write

Arff files (:mod:`scipy.io.arff`)
=================================

.. module:: scipy.io.arff

.. autosummary::
   :toctree: generated/

   loadarff

Netcdf (:mod:`scipy.io.netcdf`)
===============================

.. module:: scipy.io.netcdf

.. autosummary::
   :toctree: generated/

   netcdf_file - A file object for NetCDF data
   netcdf_variable - A data object for the netcdf module

"""
# matfile read and write
from matlab import loadmat, savemat, byteordercodes

# netCDF file support
from netcdf import netcdf_file, netcdf_variable

from data_store import save_as_module
from mmio import mminfo, mmread, mmwrite
from idl import readsav
from harwell_boeing import hb_read, hb_write

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test

