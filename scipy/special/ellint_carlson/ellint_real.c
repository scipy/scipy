#ifdef ELLINT_POLY_COMPLEX
#undef ELLINT_POLY_COMPLEX
#endif

#define ELLINT_POLY_REAL
#include "_rf_impl.c"
#include "_rc_impl.c"
#include "_rd_impl.c"
#include "_rg_impl.c"
#include "_rj_impl.c"
#undef ELLINT_POLY_REAL
