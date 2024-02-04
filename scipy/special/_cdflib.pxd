cdef double algdiv(double a, double b) noexcept nogil
cdef double alngam(double x) noexcept nogil
cdef double alnrel(double a) noexcept nogil
cdef double apser(double a, double b, double x, double eps) noexcept nogil
cdef double basym(double a, double b, double lmbda, double eps) noexcept nogil
cdef double bcorr(double a0, double b0) noexcept nogil
cdef double betaln(double a0, double b0) noexcept nogil
cdef double bfrac(double a, double b, double x, double y, double lmbda, double eps) noexcept nogil
cdef (double, int) bgrat(double a, double b, double x , double y, double w, double eps) noexcept nogil
cdef double bpser(double a, double b, double x, double eps) noexcept nogil
cdef (double, double, int) bratio(double a, double b, double x, double y) noexcept nogil
cdef double brcmp1(int mu, double a, double b, double x, double y) noexcept nogil
cdef double brcomp(double a, double b, double x, double y) noexcept nogil
cdef double bup(double a, double b, double x, double y, int n, double eps) noexcept nogil
cdef (double, double, int, double) cdfbet_which1( double x, double y,double a, double b) noexcept nogil
cdef (double, double, int, double) cdfbet_which2(double p, double q, double a, double b) noexcept nogil
cdef (double, int, double) cdfbet_which3(double p, double q, double x, double y, double b) noexcept nogil
cdef (double, int, double) cdfbet_which4(double p, double q, double x, double y, double a) noexcept nogil
cdef (double, double, int, double) cdfbin_which1(double s, double xn, double pr, double ompr) noexcept nogil
cdef (double, int, double) cdfbin_which2(double p, double q, double xn, double pr, double ompr) noexcept nogil
cdef (double, int, double) cdfbin_which3(double p, double q, double s, double pr, double ompr) noexcept nogil
cdef (double, double, int, double) cdfbin_which4(double p, double q, double s, double xn) noexcept nogil
cdef (double, double, int, double) cdfchi_which1(double x, double df) noexcept nogil
cdef (double, int, double) cdfchi_which2(double p, double q, double df) noexcept nogil
cdef (double, int, double) cdfchi_which3(double p, double q, double x) noexcept nogil
cdef (double, double, int, double) cdfchn_which1(double x, double df, double pnonc) noexcept nogil
cdef (double, int, double) cdfchn_which2(double p, double df, double pnonc) noexcept nogil
cdef (double, int, double) cdfchn_which3(double p, double x, double pnonc) noexcept nogil
cdef (double, int, double) cdfchn_which4(double p, double x, double df) noexcept nogil
cdef (double, double, int, double) cdff_which1(double f, double dfn, double dfd) noexcept nogil
cdef (double, int, double) cdff_which2(double p, double q, double dfn, double dfd) noexcept nogil
cdef (double, int, double) cdff_which3(double p, double q, double f, double dfd) noexcept nogil
cdef (double, int, double) cdff_which4(double p, double q, double f, double dfn) noexcept nogil
cdef (double, double, int, double) cdffnc_which1(double f, double dfn, double dfd, double phonc) noexcept nogil
cdef (double, int, double) cdffnc_which2(double p, double q, double dfn, double dfd, double phonc) noexcept nogil
cdef (double, int, double) cdffnc_which3(double p, double q, double f, double dfd, double phonc) noexcept nogil
cdef (double, int, double) cdffnc_which4(double p, double q, double f, double dfn, double phonc) noexcept nogil
cdef (double, int, double) cdffnc_which5(double p, double q, double f, double dfn, double dfd) noexcept nogil
cdef (double, double, int, double) cdfgam_which1(double x, double shape, double scale) noexcept nogil
cdef (double, int, double) cdfgam_which2(double p, double q, double shape, double scale) noexcept nogil
cdef (double, int, double) cdfgam_which3(double p, double q, double x, double scale) noexcept nogil
cdef (double, int, double) cdfgam_which4(double p, double q, double x, double shape) noexcept nogil
cdef (double, double, int, double) cdfnbn_which1(double s, double xn, double pr, double ompr) noexcept nogil
cdef (double, int, double) cdfnbn_which2(double p, double q, double xn, double pr, double ompr) noexcept nogil
cdef (double, int, double) cdfnbn_which3(double p, double q, double s, double pr, double ompr) noexcept nogil
cdef (double, double, int, double) cdfnbn_which4(double p, double q, double s, double xn) noexcept nogil
cdef (double, double, int, double) cdfnor_which1(double x, double mean, double sd) noexcept nogil
cdef (double, int, double) cdfnor_which2(double p, double q, double mean, double sd) noexcept nogil
cdef (double, int, double) cdfnor_which3(double p, double q, double x, double sd) noexcept nogil
cdef (double, int, double) cdfnor_which4(double p, double q, double x, double mean) noexcept nogil
cdef (double, double, int, double) cdfpoi_which1(double s, double xlam) noexcept nogil
cdef (double, int, double) cdfpoi_which2(double p, double q, double xlam) noexcept nogil
cdef (double, int, double) cdfpoi_which3(double p, double q, double s) noexcept nogil
cdef (double, double, int, double) cdft_which1(double t, double df) noexcept nogil
cdef (double, int, double) cdft_which2(double p, double q, double df) noexcept nogil
cdef (double, int, double) cdft_which3(double p, double q, double t) noexcept nogil
cdef (double, double, int, double) cdftnc_which1(double t, double df, double pnonc) noexcept nogil
cdef (double, int, double) cdftnc_which2(double p, double q, double df, double pnonc) noexcept nogil
cdef (double, int, double) cdftnc_which3(double p, double q, double t, double pnonc) noexcept nogil
cdef (double, int, double) cdftnc_which4(double p, double q, double t, double df) noexcept nogil
cdef (double, double) cumbet(double x, double y, double a, double b) noexcept nogil
cdef (double, double) cumbin(double s, double xn, double pr, double ompr) noexcept nogil
cdef (double, double) cumchi(double x, double df) noexcept nogil
cdef (double, double) cumchn(double x, double df, double pnonc) noexcept nogil
cdef (double, double) cumf(double f, double dfn, double dfd) noexcept nogil
cdef (double, double, int) cumfnc(double f, double dfn, double dfd, double pnonc) noexcept nogil
cdef (double, double) cumgam(double x, double a) noexcept nogil
cdef (double, double) cumnbn(double s, double xn, double pr, double ompr) noexcept nogil
cdef (double, double) cumnor(double x) noexcept nogil
cdef (double, double) cumpoi(double s, double xlam) noexcept nogil
cdef (double, double) cumt(double t, double df) noexcept nogil
cdef (double, double) cumtnc(double t, double df, double pnonc) noexcept nogil
cdef double devlpl(double *a, int n, double x) noexcept nogil
cdef double dinvnr(double p, double q) noexcept nogil
cdef double dt1(double p, double q, double df) noexcept nogil
cdef double erf(double x) noexcept nogil
cdef double erfc1(int ind, double x) noexcept nogil
cdef double esum(int mu, double x) noexcept nogil
cdef double fpser(double a, double b, double x, double eps) noexcept nogil
cdef double gam1(double a) noexcept nogil
cdef (double, int) gaminv(double a, double p, double q, double x0) noexcept nogil
cdef double gaminv_helper_30(double a, double s, double z, double y) noexcept nogil
cdef double gamln(double a) noexcept nogil
cdef double gamln1(double a) noexcept nogil
cdef double gamma(double a) noexcept nogil
cdef (double, double) grat1(double a, double x, double r, double eps) noexcept nogil
cdef (double, double) gratio(double a, double x, int ind) noexcept nogil
cdef double gsumln(double a, double b) noexcept nogil
cdef double psi(double xx) noexcept nogil
cdef double rcomp(double a, double x) noexcept nogil
cdef double rexp(double x) noexcept nogil
cdef double rlog(double x) noexcept nogil
cdef double rlog1(double x) noexcept nogil
cdef double stvaln(double p) noexcept nogil
