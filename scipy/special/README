# Author:  Travis E. Oliphant
#          oliphant.travis@altavista.net
#
# This can be freely distributed provided this notification remain and it is
# understood that no warranty is expressed or implied.

CEPHES module for Numerical Python


DESCRIPTION: 
=====================
This package contains all of the important functions from the cephes
library with some complex extensions from the amos library 
in a python module.  There are now some hooks for including toms algorithms 
as well (thanks Lorenzo Catucci).

The module is patterned after the umath module and uses the ufunc
object as the interface so that all of the functions are accessible on
multidimensional arrays.  The ufunc methods (reduce, outer, etc.) are
also available to the functions which take two arguments (the need for 
these methods is dubious for these functions, they are just along for the 
ride as part of the ufunc object).  The ufunc
object interface is general enough to allow even functions which take
more than two arguments and return more than one argument to be
implemented. 

In the ./docs directory is the documentation for the all of the cephes
library, and it includes an HTML page showing brief documentation for 
all included functions calls.

INSTALLATION:

Part of SciPy
 

CHANGES:
===========
0.5        Initial Release

0.51       Fixed the name of a subroutine that was conflicting with
	   one from umath so that Numeric and cephes work together.

0.52       Fixed the real reason Numeric and cephes weren't working together:
           I had copied a setup routine from umathmodule.c that wasn't
           necessary and shouldn't have been called.

0.60       Added kolmogorov and kolmogi functions from cephes file.  These 
	   functions were not really functions in the cephes library, but they
	   existed in a standalone program.  They may not be as robust. 
	   Kolmogorov hangs if given 0 as an argument.  Someone could
	   probably fix this.

0.61       Changed so Kolmogorov doesn't hang on 0 argument.  Returns 1 
           Also kolmogi returns 0 on 1.0 as argument.

0.70       Added many more functions from the library.  Wrote wrappers
	   so that functions which took an integer as one of two
	   arguments could be added.  Also added wrappers around
	   functions which returned more than one value:  I made a new
	   function for each returned value.  yn,kn,jn Student T, and 
	   Poisson distribution functions are notable inclusions.

0.71       Added info about mconf.h in the README and created a patchfile
           to patch mconf.h to a generic BIGENDIAN computer from an IBM_PC
           Changed gamma to Gamma in the c-code so that Solaris's gamma
	   which actually computes lgamma is not called.  Removed log2
	   from library as it gave poor results on UNK machines. 

1.0        Added all of the remaining cephes calls to the library and
	   eliminated the wrapper functions after learning how to
	   really use the ufunc object correctly to handle generic
	   functions.

1.1	   Added complex extensions to bessel and airy functions by
	   incorporating the amos libraries into the module with some 
	   wrappers around the fortran code.

1.11       Fixed small bugs.

1.15       Added changes from Lorenzo Catucci: a new wofz function and a
           testing script.
  
1.2        Added a loopfunc method to the cephes module.  This method is a 
	   kind of generalized map that takes a Python function with scalar
           inputs and outputs and a tuple of input arrays or nested sequence
	   objects and returns an array created by applying the Python function
	   in turn to input tuples created from the elements of the input
           arrays.  It uses the same broadcasting rules as are available with
	   ufuncs.  Basically, this turns a regular Python function into
	   the equivalent of a ufunc.   

1.3        Added besselpoly function.  Added empty docstrings to compile with
	   Numeric 15.3 (thanks Warren Focke).  Changed makefiles to compile 
	   shared libraries with libtool.


