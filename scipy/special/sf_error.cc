#include <Python.h>

#include <stdlib.h>
#include <stdarg.h>

#include "special/error.h"
// #include "sf_error.h"

extern "C" {
#include "sf_error.c"
}


#ifdef SP_SPECFUN_ERROR
void special::set_error(const char *func_name, sf_error_t code,
			       const char *fmt, ...) {
    /* Definition of error handling for special C++ library of special
     * functions used in SciPy.
     *
     * See special/error.h for info on valid codes in enum sf_error_t.
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
