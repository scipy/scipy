from . cimport sf_error

from libc.math cimport NAN, isnan, isinf, isfinite

cdef extern from "cdflib.h" nogil:
    cdef struct TupleDDI:
        double d1
        double d2
        int i1

    cdef struct TupleDID:
        double d1
        int i1
        double d2

    cdef struct TupleDDID:
        double d1
        double d2
        int i1
        double d3

    TupleDID cdfbin_which2(double, double, double, double, double)
    TupleDID cdfbin_which3(double, double, double, double, double)
    TupleDID cdfchi_which3(double, double, double)
    TupleDID cdff_which4(double, double, double, double);
    TupleDID cdffnc_which3(double, double, double, double, double);
    TupleDID cdffnc_which4(double, double, double, double, double);
    TupleDID cdffnc_which5(double, double, double, double, double);
    TupleDID cdfgam_which2(double, double, double, double);
    TupleDID cdfgam_which4(double, double, double, double);
    TupleDID cdfnbn_which2(double, double, double, double, double);
    TupleDID cdfnbn_which3(double, double, double, double, double);
    TupleDID cdfnor_which3(double, double, double, double);
    TupleDID cdfnor_which4(double, double, double, double);
    TupleDID cdfpoi_which2(double, double, double);
    TupleDID cdft_which3(double, double, double);
    TupleDID cdftnc_which3(double, double, double, double);
    TupleDID cdftnc_which4(double, double, double, double);


cdef inline double get_result(
        char *name,
        char **argnames,
        double result,
        int status,
        double bound,
        int return_bound
) noexcept nogil:
    cdef char *arg
    """Get result and perform error handling from cdflib output."""
    if status < 0:
        arg = argnames[-(status + 1)]
        sf_error.error(name, sf_error.ARG,
                       "Input parameter %s is out of range", arg)
        return NAN
    if status == 0:
        return result
    if status == 1:
        sf_error.error(name, sf_error.OTHER,
                       "Answer appears to be lower than lowest search bound (%g)", bound)
        return bound if return_bound else NAN
    if status == 2:
        sf_error.error(name, sf_error.OTHER,
                       "Answer appears to be higher than highest search bound (%g)", bound)
        return bound if return_bound else NAN
    if status == 3 or status == 4:
        sf_error.error(name, sf_error.OTHER,
                       "Two internal parameters that should sum to 1.0 do not.")
        return NAN
    if status == 10:
        sf_error.error(name, sf_error.OTHER, "Computational error")
        return NAN
    sf_error.error(name, sf_error.OTHER, "Unknown error.")
    return NAN


cdef inline double bdtrik(double p, double xn, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status
        char *argnames[5]
        TupleDID ret

    if isnan(p) or not isfinite(xn) or isnan(pr):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "xn"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    ret = cdfbin_which2(p, q, xn, pr, ompr)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("btdtrik", argnames, result, status, bound, 1)


cdef inline double bdtrin(double s, double p, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status = 10
        char *argnames[5]
        TupleDID ret

    if isnan(s) or isnan(p) or isnan(pr):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "s"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    ret = cdfbin_which3(p, q, s, pr, ompr)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("btdtrin", argnames, result, status, bound, 1)


cdef inline double fdtridfd(double dfn, double p, double f) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]
        TupleDID ret

    if isnan(dfn) or isnan(p) or isnan(f):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "f"
    argnames[3] = "dfn"

    ret = cdff_which4(p, q, f, dfn)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("fdtridfd", argnames, result, status, bound, 1)


cdef inline double gdtria(double p, double shp, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]
        TupleDID ret

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "shp"

    if isnan(p) or isnan(shp) or isnan(x):
      return NAN

    ret = cdfgam_which4(p, q, x, shp)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("gdtria", argnames, result, status, bound, 1)


cdef inline double gdtrix(double scl, double shp, double p) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]
        TupleDID ret

    if isnan(scl) or isnan(shp) or isnan(p):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "shp"
    argnames[3] = "scl"

    ret = cdfgam_which2(p, q, shp, scl)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("gdtrix", argnames, result, status, bound, 1)


cdef inline double nbdtrik(double p, double xn, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status = 10
        char *argnames[5]
        TupleDID ret

    if isnan(p) or not isfinite(xn) or isnan(pr):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "xn"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    ret = cdfnbn_which2(p, q, xn, pr, ompr)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("nbdtrik", argnames, result, status, bound, 1)


cdef inline double nbdtrin(double s, double p, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status = 10
        char *argnames[5]
        TupleDID ret

    if isnan(s) or isnan(p) or isnan(pr):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "s"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    ret = cdfnbn_which3(p, q, s, pr, ompr)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("nbdtrin", argnames, result, status, bound, 1)


cdef inline double ncfdtridfd(double dfn, double p, double nc, double f) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[5]
        TupleDID ret

    if isnan(dfn) or isnan(p) or isnan(nc) or isnan(f):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "f"
    argnames[3] = "dfn"
    argnames[4] = "nc"

    ret = cdffnc_which4(p, q, f, dfn, nc)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("ncfdtridfd", argnames, result, status, bound, 1)


cdef inline double ncfdtridfn(double p, double dfd, double nc, double f) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[5]
        TupleDID ret

    if isnan(p) or isnan(dfd) or isnan(nc) or isnan(f):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "f"
    argnames[3] = "dfd"
    argnames[4] = "nc"

    ret = cdffnc_which3(p, q, f, dfd, nc)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("ncfdtridfn", argnames, result, status, bound, 1)


cdef inline double ncfdtrinc(double dfn, double dfd, double p, double f) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[5]
        TupleDID ret

    if isnan(dfn) or isnan(dfd) or isnan(p) or isnan(f):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "f"
    argnames[3] = "dfn"
    argnames[4] = "dfd"

    ret = cdffnc_which5(p, q, f, dfn, dfd)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("ncfdtrinc", argnames, result, status, bound, 1)


cdef inline double nctdtridf(double p, double nc, double t) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]
        TupleDID ret

    if isnan(p) or isnan(nc) or isnan(t):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "t"
    argnames[3] = "nc"

    ret = cdftnc_which3(p, q, t, nc)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("nctdtridf", argnames, result, status, bound, 1)


cdef inline double nctdtrinc(double df, double p, double t) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]
        TupleDID ret

    if isnan(df) or isnan(p) or isnan(t):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "t"
    argnames[3] = "df"

    ret = cdftnc_which4(p, q, t, df)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("nctdtrinc", argnames, result, status, bound, 1)


cdef inline double nrdtrimn(double p, double std, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]
        TupleDID ret

    if isnan(p) or isnan(std) or isnan(x):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "std"

    ret = cdfnor_which3(p, q, x, std)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("nrdtrimn", argnames, result, status, bound, 1)


cdef inline double nrdtrisd(double mn, double p, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]
        TupleDID ret

    if isnan(mn) or isnan(p) or isnan(x):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "mn"

    ret = cdfnor_which4(p, q, x, mn)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("nrdtrisd", argnames, result, status, bound, 1)


cdef inline double pdtrik(double p, double xlam) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[3]
        TupleDID ret

    if isnan(p) or isnan(xlam):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "xlam"

    ret = cdfpoi_which2(p, q, xlam)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("pdtrik", argnames, result, status, bound, 1)


cdef inline double stdtridf(double p, double t) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[3]
        TupleDID ret

    if isnan(p) or isnan(q) or isnan(t):
        return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "t"

    ret = cdft_which3(p, q, t)
    result, status, bound = ret.d1, ret.i1, ret.d2
    return get_result("stdtridf", argnames, result, status, bound, 1)
