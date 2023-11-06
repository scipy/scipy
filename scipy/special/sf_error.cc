#include <Python.h>

#include <stdlib.h>
#include <stdarg.h>

#include "error.h"
// #include "sf_error.h"

extern "C" {
#include "sf_error.c"
}


#ifdef SP_SPECFUN_ERROR
void scipy::special::set_error(const char *func_name, sf_error_t code,
			       const char *fmt, ...) {
    /* Definition of error handling for scipy::special C++ library of special
     * functions used in SciPy.
     *
     * The arg "code" must take an int that is a valid value for enum sf_error_t.
     * See sf_error.h for the definition of this enum.
     *
     * Other packages making use of this library can supply their own implementation. 
     */
    va_list ap;
    va_start(ap, fmt);
    // Version of sf_error that takes va_list instead of variadic args.
    sf_error_v(func_name, code, fmt, ap);
    va_end(ap);
}
#endif

