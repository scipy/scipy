# The netCDF interface is mostly written in C; this module only
# adds some luxury and puts everything in its proper place in the
# package hierarchy.

import sys

if sys.modules.has_key('pythondoc'):

    # Fake code just for the docstrings!

    from Scientific_netcdf import _C_API

    class NetCDFFile:

        """netCDF file

        Constructor: NetCDFFile(|filename|, |mode|='"r"')

        Arguments:

        |filename| -- name of the netCDF file. By convention, netCDF files
                      have the extension ".nc", but this is not enforced.
                      The filename may contain a home directory indication
                      starting with "~".

        |mode| -- access mode. "r" means read-only; no data can be modified.
                  "w" means write; a new file is created, an existing
                  file with the same name is deleted. "a" means append
                  (in analogy with serial files); an existing file is
                  opened for reading and writing, and if the file does
                  not exist it is created. "r+" is similar to "a",
                  but the file must already exist. An "s" can be appended
                  to any of the modes listed above; it indicates that the
                  file will be opened or created in "share" mode, which
                  reduces buffering in order to permit simultaneous read
                  access by other processes to a file that is being written.

        A NetCDFFile object has two standard attributes: 'dimensions' and
        'variables'. The values of both are dictionaries, mapping dimension
        names to their associated lengths and variable names to variables,
        respectively. Application programs should never modify these
        dictionaries.

        All other attributes correspond to global attributes defined in the
        netCDF file. Global file attributes are created by assigning to
        an attribute of the NetCDFFile object. 
        """

        def __init__(self, *args):
            raise ImportError, "this code should never be executed"

        def close(self):
            """Closes the file. Any read or write access to the file
            or one of its variables after closing raises an exception."""
            pass
        
        def createDimension(self, name, length):
            """Creates a new dimension with the given |name| and
            |length|. |length| must be a positive integer or 'None',
            which stands for the unlimited dimension. Note that there can
            be only one unlimited dimension in a file."""
            pass

        def createVariable(self, name, type, dimensions):
            """Creates a new variable with the given |name|, |type|, and
            |dimensions|. The |type| is a one-letter string with the same
            meaning as the typecodes for arrays in module Numeric; in
            practice the predefined type constants from Numeric should
            be used. |dimensions| must be a tuple containing dimension
            names (strings) that have been defined previously.

            The return value is the NetCDFVariable object describing the
            new variable."""
            pass

        def sync(self):
            "Writes all buffered data to the disk file."
        flush = sync

    class NetCDFVariable:

        """Variable in a netCDF file

        NetCDFVariable objects are constructed by calling the method
        'createVariable' on the NetCDFFile object.

        NetCDFVariable objects behave much like array objects defined
        in module Numeric, except that their data resides in a file.
        Data is read by indexing and written by assigning to an
        indexed subset; the entire array can be accessed by the index
        '[:]' or using the methods 'getValue' and
        'assignValue'. NetCDFVariable objects also have attribute
        "shape" with the same meaning as for arrays, but the shape
        cannot be modified. There is another read-only attribute
        "dimensions", whose value is the tuple of dimension names.

        All other attributes correspond to variable attributes defined in the
        netCDF file. Variable attributes are created by assigning to
        an attribute of the NetCDFVariable object. 

        Note:
        If a file open for reading is simultaneously written by another program,
        the size of the unlimited dimension may change. Every time the shape
        of a variable is requested, the current size will be obtained from
        the file. For reading and writing, the size obtained during the last
        shape request is used. This ensures consistency: foo[-1] means the
        same thing no matter how often it is evaluated, as long as the shape
        is not re-evaluated in between.
        """

        def __init__(self, *args):
            raise ImportError, "this code should never be executed"

        def assignValue(self, value):
            """Assigns |value| to the variable. This method allows
            assignment to scalar variables, which cannot be indexed."""
            pass

        def getValue(self):
            """Returns the value of the variable. This method allows
            access to scalar variables, which cannot be indexed."""
            pass

        def typecode(self):
            "Return the variable's type code (a string)."
            pass

else:

    # This is the real code.
    from Scientific_netcdf import *
    from Scientific_netcdf import _C_API
    import os

    _NetCDFFile = NetCDFFile
    def NetCDFFile(filename, mode=None, history=None):
        filename = os.path.expanduser(filename)
        args = (filename,)
        if mode is not None:
            args = args + (mode,)
            if history is not None:
                args = args + (history,)
        return apply(_NetCDFFile, args)

    del sys
