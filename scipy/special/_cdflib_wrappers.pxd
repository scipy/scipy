from . cimport sf_error

from libc.math cimport NAN, isnan, isinf, isfinite
cdef extern from "cephes.h" nogil:
    double ndtr(double a)
    double ndtri(double y0)

from ._cdflib cimport (
    cdfbet_which3,
    cdfbet_which4,
    cdfbin_which2,
    cdfbin_which3,
    cdfchi_which3,
    cdfchn_which1,
    cdfchn_which2,
    cdfchn_which3,
    cdfchn_which4,
    cdff_which4,
    cdfgam_which2,
    cdfgam_which3,
    cdfgam_which4,
    cdfnbn_which2,
    cdfnbn_which3,
    cdffnc_which1,
    cdffnc_which2,
    cdffnc_which3,
    cdffnc_which4,
    cdffnc_which5,
    cdftnc_which1,
    cdftnc_which2,
    cdftnc_which3,
    cdftnc_which4,
    cdfnor_which3,
    cdfnor_which4,
    cdfpoi_which2,
    cdft_which1,
    cdft_which2,
    cdft_which3,
)


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


cdef inline double btdtria(double p, double b, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double y = 1.0 - x
        double result, bound
        int status = 10
        char *argnames[5]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "y"
    argnames[4] = "b"
    
    if isnan(p) or isnan(b) or isnan(x):
      return NAN

    result, status, bound = cdfbet_which3(p, q, x, y, b)
    return get_result("btdtria", argnames, result, status, bound, 1)



cdef inline double btdtrib(double a, double p, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double y = 1.0 - x
        double result, bound
        int status = 10
        char *argnames[5]

    if isnan(a) or isnan(p) or isnan(x):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "y"
    argnames[4] = "a"

    result, status, bound = cdfbet_which4(p, q, x, y, a)
    return get_result("btdtrib", argnames, result, status, bound, 1)


cdef inline double bdtrik(double p, double xn, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status
        char *argnames[5]

    if isnan(p) or not isfinite(xn) or isnan(pr):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "xn"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    result, status, bound = cdfbin_which2(p, q, xn, pr, ompr)
    return get_result("btdtrik", argnames, result, status, bound, 1)


cdef inline double bdtrin(double s, double p, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status = 10
        char *argnames[5]

    if isnan(s) or isnan(p) or isnan(pr):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "s"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    result, status, bound = cdfbin_which3(p, q, s, pr, ompr)
    return get_result("btdtrin", argnames, result, status, bound, 1)


cdef inline double chdtriv(double p, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[3]

    if isnan(p) or isnan(x):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"

    result, status, bound = cdfchi_which3(p, q, x)
    return get_result("chdtriv", argnames, result, status, bound, 1)


cdef inline double chndtr(double x, double df, double nc) noexcept nogil:
    cdef:
        double result, _, bound
        int status = 10
        char *argnames[3]

    if isnan(x) or isnan(df) or isnan(nc):
      return NAN

    argnames[0] = "x"
    argnames[1] = "df"
    argnames[2] = "nc"

    result, _, status, bound = cdfchn_which1(x, df, nc)
    return get_result("chndtr", argnames, result, status, bound, 1)


cdef inline double chndtridf(double x, double p, double nc) noexcept nogil:
    cdef:
        double result, bound
        int status = 10
        char *argnames[3]

    if isnan(x) or isnan(p) or isnan(nc):
      return NAN

    argnames[0] = "p"
    argnames[1] = "x"
    argnames[2] = "nc"

    result, status, bound = cdfchn_which3(p, x, nc)
    return get_result("chndtridf", argnames, result, status, bound, 1)


cdef inline double chndtrinc(double x, double df, double p) noexcept nogil:
    cdef:
        double result, bound
        int status = 10
        char *argnames[3]

    if isnan(x) or isnan(df) or isnan(p):
      return NAN

    argnames[0] = "p"
    argnames[1] = "x"
    argnames[2] = "df"

    result, status, bound = cdfchn_which4(p, x, df)
    return get_result("chndtrinc", argnames, result, status, bound, 1)


cdef inline double chndtrix(double p, double df, double nc) noexcept nogil:
    cdef:
        double result, bound
        int status = 10
        char *argnames[3]

    if isnan(p) or isnan(df) or isnan(nc):
      return NAN

    argnames[0] = "p"
    argnames[1] = "df"
    argnames[2] = "nc"

    result, status, bound = cdfchn_which2(p, df, nc)
    return get_result("chndtrix", argnames, result, status, bound, 0)


cdef inline double fdtridfd(double dfn, double p, double f) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    if isnan(dfn) or isnan(p) or isnan(f):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "f"
    argnames[3] = "dfn"

    result, status, bound = cdff_which4(p, q, f, dfn)
    return get_result("fdtridfd", argnames, result, status, bound, 1)


cdef inline double gdtria(double p, double shp, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "shp"

    if isnan(p) or isnan(shp) or isnan(x):
      return NAN

    result, status, bound = cdfgam_which4(p, q, x, shp)
    return get_result("gdtria", argnames, result, status, bound, 1)


cdef inline double gdtrib(double scl, double p, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    if isnan(scl) or isnan(p) or isnan(x):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "scl"

    result, status, bound = cdfgam_which3(p, q, x, scl)
    return get_result("gdtrib", argnames, result, status, bound, 1)


cdef inline double gdtrix(double scl, double shp, double p) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    if isnan(scl) or isnan(shp) or isnan(p):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "shp"
    argnames[3] = "scl"

    result, status, bound = cdfgam_which2(p, q, shp, scl)
    return get_result("gdtrix", argnames, result, status, bound, 1)


cdef inline double nbdtrik(double p, double xn, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status = 10
        char *argnames[5]

    if isnan(p) or not isfinite(xn) or isnan(pr):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "xn"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    result, status, bound = cdfnbn_which2(p, q, xn, pr, ompr)
    return get_result("nbdtrik", argnames, result, status, bound, 1)


cdef inline double nbdtrin(double s, double p, double pr) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double ompr = 1.0 - pr
        double result, bound
        int status = 10
        char *argnames[5]

    if isnan(s) or isnan(p) or isnan(pr):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "s"
    argnames[3] = "pr"
    argnames[4] = "ompr"

    result, status, bound = cdfnbn_which3(p, q, s, pr, ompr)
    return get_result("nbdtrin", argnames, result, status, bound, 1)


cdef inline double ncfdtr(double dfn, double dfd, double nc, double f) noexcept nogil:
    cdef:
        double result, _, bound
        int status = 10
        char *argnames[4]

    if isnan(dfn) or isnan(dfd) or isnan(nc) or isnan(f):
      return NAN

    argnames[0] = "f"
    argnames[1] = "dfn"
    argnames[2] = "dfd"
    argnames[3] = "nc"

    result, _, status, bound = cdffnc_which1(f, dfn, dfd, nc)
    return get_result("ncfdtr", argnames, result, status, bound, 0)


cdef inline double ncfdtri(double dfn, double dfd, double nc, double p) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[5]

    if isnan(dfn) or isnan(dfd) or isnan(nc) or isnan(p):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "dfn"
    argnames[3] = "dfd"
    argnames[4] = "nc"

    result, status, bound = cdffnc_which2(p, q, dfn, dfd, nc)
    return get_result("ncfdtri", argnames, result, status, bound, 1)


cdef inline double ncfdtridfd(double dfn, double p, double nc, double f) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[5]

    if isnan(dfn) or isnan(p) or isnan(nc) or isnan(f):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "f"
    argnames[3] = "dfn"
    argnames[4] = "nc"

    result, status, bound = cdffnc_which4(p, q, f, dfn, nc)
    return get_result("ncfdtridfd", argnames, result, status, bound, 1)


cdef inline double ncfdtridfn(double p, double dfd, double nc, double f) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[5]

    if isnan(p) or isnan(dfd) or isnan(nc) or isnan(f):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "f"
    argnames[3] = "dfd"
    argnames[4] = "nc"

    result, status, bound = cdffnc_which3(p, q, f, dfd, nc)
    return get_result("ncfdtridfn", argnames, result, status, bound, 1)


cdef inline double ncfdtrinc(double dfn, double dfd, double p, double f) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[5]

    if isnan(dfn) or isnan(dfd) or isnan(p) or isnan(f):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "f"
    argnames[3] = "dfn"
    argnames[4] = "dfd"

    result, status, bound = cdffnc_which5(p, q, f, dfn, dfd)
    return get_result("ncfdtrinc", argnames, result, status, bound, 1)


cdef inline double nctdtr(double df, double nc, double t) noexcept nogil:
    cdef:
        double result, _, bound
        int status = 10
        char *argnames[3]

    if isnan(df) or isnan(nc) or isnan(t):
      return NAN

    argnames[0] = "t"
    argnames[1] = "df"
    argnames[2] = "nc"

    result, _, status, bound = cdftnc_which1(t, df, nc)
    return get_result("nctdtr", argnames, result, status, bound, 1)


cdef inline double nctdtridf(double p, double nc, double t) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    if isnan(p) or isnan(nc) or isnan(t):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "t"
    argnames[3] = "nc"

    result, status, bound = cdftnc_which3(p, q, t, nc)
    return get_result("nctdtridf", argnames, result, status, bound, 1)


cdef inline double nctdtrinc(double df, double p, double t) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    if isnan(df) or isnan(p) or isnan(t):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "t"
    argnames[3] = "df"

    result, status, bound = cdftnc_which4(p, q, t, df)
    return get_result("nctdtrinc", argnames, result, status, bound, 1)


cdef inline double nctdtrit(double df, double nc, double p) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    if isnan(df) or isnan(nc) or isnan(p):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "df"
    argnames[3] = "nc"

    result, status, bound = cdftnc_which2(p, q, df, nc)
    return get_result("nctdtrit", argnames, result, status, bound, 1)


cdef inline double nrdtrimn(double p, double std, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    if isnan(p) or isnan(std) or isnan(x):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "std"

    result, status, bound = cdfnor_which3(p, q, x, std)
    return get_result("nrdtrimn", argnames, result, status, bound, 1)


cdef inline double nrdtrisd(double mn, double p, double x) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[4]

    if isnan(mn) or isnan(p) or isnan(x):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "x"
    argnames[3] = "mn"

    result, status, bound = cdfnor_which4(p, q, x, mn)
    return get_result("nrdtrisd", argnames, result, status, bound, 1)


cdef inline double pdtrik(double p, double xlam) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[3]

    if isnan(p) or isnan(xlam):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "xlam"

    result, status, bound = cdfpoi_which2(p, q, xlam)
    return get_result("pdtrik", argnames, result, status, bound, 1)


cdef inline double stdtr(double df, double t) noexcept nogil:
    cdef:
        double result, _, bound
        int status = 10
        char *argnames[2]

    argnames[0] = "t"
    argnames[1] = "df"

    if isinf(df) and df > 0:
        return NAN if isnan(t) else ndtr(t)

    if isnan(df) or isnan(t):
      return NAN

    result, _, status, bound = cdft_which1(t, df)
    return get_result("stdtr", argnames, result, status, bound, 1)


cdef inline double stdtridf(double p, double t) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[3]

    if isnan(p) or isnan(q) or isnan(t):
        return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "t"

    result, status, bound = cdft_which3(p, q, t)
    return get_result("stdtridf", argnames, result, status, bound, 1)


cdef inline double stdtrit(double df, double p) noexcept nogil:
    cdef:
        double q = 1.0 - p
        double result, bound
        int status = 10
        char *argnames[3]

    if isinf(df) and df > 0:
        return NAN if isnan(p) else ndtri(p)

    if isnan(p) or isnan(df):
      return NAN

    argnames[0] = "p"
    argnames[1] = "q"
    argnames[2] = "df"

    result, status, bound = cdft_which2(p, q, df)
    return get_result("stdtrit", argnames, result, status, bound, 1)
