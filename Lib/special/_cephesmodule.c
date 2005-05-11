
/* Cephes module version 1.5
 *  This module defines the functions in the cephes and amos libraries as
 *   Numerical python ufunc objects so that they can operate on arbitrary 
 *   NumPy arrays with broadcasting and typecasting rules implemented.
 *  
 *  Copyright 1999  Travis E. Oliphant
 * Revisions 2002 (added functions from cdflib)
 */

#if defined(NUMARRAY)

#include "_na_cephesmodule.c"

#else

#include "_nc_cephesmodule.c"

#endif  
