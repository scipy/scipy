"""
Data input and output
=====================

  Functions

     loadmat -- read a MATLAB style mat file (version 4 through 7.1)
     savemat -- write a MATLAB (version through 7.1) style mat file
     netcdf_file -- read NetCDF files (version of ``pupynere`` package)
        This version cannot write NetCDF files.
     save_as_module -- simple storing of Python dictionary into module
               that can then be imported and the data accessed as
               attributes of the module.
     mminfo   -- query matrix info from Matrix Market formatted file
     mmread   -- read matrix from Matrix Market formatted file
     mmwrite  -- write matrix to Matrix Market formatted file

"""
postpone_import = 1
