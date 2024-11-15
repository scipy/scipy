from libc.stdlib cimport free

cdef extern from "./c_src/levyst.h":
    struct nolan_precanned:
        double (*g)(nolan_precanned *, double)
        double alpha
        double zeta
        double xi
        double zeta_prefactor
        double alpha_exp
        double alpha_xi
        double zeta_offset
        double two_beta_div_pi
        double pi_div_two_beta
        double x0_div_term
        double c1
        double c2
        double c3

    nolan_precanned *nolan_precan(double alpha, double beta, double x0)

cdef class Nolan:
    cdef nolan_precanned * p

    def __init__(self, alpha, beta, x0):
        self.p = nolan_precan(alpha, beta, x0)

    def g(self, theta):
       return self.p.g(self.p, theta)

    @property
    def zeta(self):
        return self.p.zeta

    @property
    def xi(self):
        return self.p.xi

    @property
    def c1(self):
        return self.p.c1

    @property
    def c2(self):
        return self.p.c2

    @property
    def c3(self):
        return self.p.c3

    def __dealloc__(self):
        free(self.p)
