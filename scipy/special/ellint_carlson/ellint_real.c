#ifdef ELLINT_POLY_COMPLEX
#undef ELLINT_POLY_COMPLEX
#endif

#define ELLINT_POLY_REAL
#include "ellint_carlson.h"
#include "ellint_poly.h"
#include "ellint_common.h"
#include "_arithmetic.c"
#include "_rf_impl.c"
#include "_rc_impl.c"
#include "_rd_impl.c"
#include "_rg_impl.c"
#include "_rj_impl.c"
#undef ELLINT_POLY_REAL
