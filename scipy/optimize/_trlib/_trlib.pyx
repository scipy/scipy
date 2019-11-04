from scipy.optimize._trustregion import (_minimize_trust_region, BaseQuadraticSubproblem)
import numpy as np
from . cimport ctrlib
cimport libc.stdio
cimport numpy as np

from scipy._lib.messagestream cimport MessageStream


class TRLIBQuadraticSubproblem(BaseQuadraticSubproblem):

    def __init__(self, x, fun, jac, hess, hessp, tol_rel_i=-2.0, tol_rel_b=-3.0,
                 disp=False):
        super(TRLIBQuadraticSubproblem, self).__init__(x, fun, jac, hess, hessp)
        self.tol_rel_i = tol_rel_i
        self.tol_rel_b = tol_rel_b
        self.disp = disp
        self.itmax = int(min(1e9/self.jac.shape[0], 2*self.jac.shape[0]))
        cdef long itmax, iwork_size, fwork_size, h_pointer
        itmax = self.itmax
        ctrlib.trlib_krylov_memory_size(itmax, &iwork_size, &fwork_size,
                                        &h_pointer)
        self.h_pointer = h_pointer
        self.fwork = np.empty([fwork_size])
        cdef double [:] fwork_view = self.fwork
        cdef double *fwork_ptr = NULL
        if fwork_view.shape[0] > 0:
            fwork_ptr = &fwork_view[0]
        ctrlib.trlib_krylov_prepare_memory(itmax, fwork_ptr)
        self.iwork = np.zeros([iwork_size], dtype=np.int)
        self.s  = np.empty(self.jac.shape)
        self.g  = np.empty(self.jac.shape)
        self.v  = np.empty(self.jac.shape)
        self.gm = np.empty(self.jac.shape)
        self.p  = np.empty(self.jac.shape)
        self.Hp = np.empty(self.jac.shape)
        self.Q  = np.empty([self.itmax+1, self.jac.shape[0]])
        self.timing = np.zeros([ctrlib.trlib_krylov_timing_size()],
                               dtype=np.int)
        self.init = ctrlib._TRLIB_CLS_INIT

    def solve(self, double trust_radius):

        cdef long equality = 0
        cdef long itmax_lanczos = 100
        cdef double tol_r_i = self.tol_rel_i
        cdef double tol_a_i =  0.0 
        cdef double tol_r_b = self.tol_rel_b
        cdef double tol_a_b =  0.0 
        cdef double zero = 2e-16
        cdef double obj_lb = -1e20
        cdef long ctl_invariant = 0
        cdef long convexify = 1
        cdef long earlyterm = 1
        cdef double g_dot_g = 0.0
        cdef double v_dot_g = 0.0
        cdef double p_dot_Hp = 0.0
        cdef long refine = 1
        cdef long verbose = 0
        cdef long unicode = 1
        cdef long ret = 0
        cdef long action = 0
        cdef long it = 0
        cdef long ityp = 0
        cdef long itmax = self.itmax
        cdef long init  = self.init
        cdef double flt1 = 0.0
        cdef double flt2 = 0.0
        cdef double flt3 = 0.0
        prefix = b""
        cdef long   [:] iwork_view  = self.iwork
        cdef double [:] fwork_view  = self.fwork
        cdef long   [:] timing_view = self.timing
        cdef long   *iwork_ptr = NULL
        cdef double *fwork_ptr = NULL
        cdef long   *timing_ptr = NULL

        if self.disp:
            verbose = 2

        if iwork_view.shape[0] > 0:
            iwork_ptr = &iwork_view[0]
        if fwork_view.shape[0] > 0:
            fwork_ptr = &fwork_view[0]
        if timing_view.shape[0] > 0:
            timing_ptr = &timing_view[0]

        cdef MessageStream messages = MessageStream()
        try:
            while True:
                messages.clear()
                ret = ctrlib.trlib_krylov_min(init, trust_radius, equality,
                          itmax, itmax_lanczos, tol_r_i, tol_a_i,
                          tol_r_b, tol_a_b, zero, obj_lb, ctl_invariant,
                          convexify, earlyterm, g_dot_g, v_dot_g, p_dot_Hp,
                          iwork_ptr, fwork_ptr, refine, verbose, unicode,
                          prefix, messages.handle,
                          timing_ptr, &action, &it, &ityp, &flt1, &flt2, &flt3)
                if self.disp:
                    msg = messages.get()
                    if msg:
                        print(msg)

                init = 0
                if action == ctrlib._TRLIB_CLA_INIT:
                    self.s[:]  = 0.0
                    self.gm[:] = 0.0
                    self.g[:]  = self.jac
                    self.v[:]  = self.g
                    g_dot_g = np.dot(self.g, self.g)
                    v_dot_g = np.dot(self.v, self.g)
                    self.p[:]  = - self.v
                    self.Hp[:] = self.hessp(self.p)
                    p_dot_Hp = np.dot(self.p, self.Hp)
                    self.Q[0,:] = self.v/np.sqrt(v_dot_g)
                if action == ctrlib._TRLIB_CLA_RETRANSF:
                    self.s[:] = np.dot(self.fwork[self.h_pointer:self.h_pointer+it+1],
                                       self.Q[:it+1,:])
                if action == ctrlib._TRLIB_CLA_UPDATE_STATIO:
                    if ityp == ctrlib._TRLIB_CLT_CG:
                        self.s += flt1 * self.p
                if action == ctrlib._TRLIB_CLA_UPDATE_GRAD:
                    if ityp == ctrlib._TRLIB_CLT_CG:
                        self.Q[it,:] = flt2*self.v
                        self.gm[:] = self.g
                        self.g    += flt1*self.Hp
                    if ityp == ctrlib._TRLIB_CLT_L:
                        self.s[:]  = self.Hp + flt1*self.g + flt2*self.gm
                        self.gm[:] = flt3*self.g
                        self.g[:]  = self.s
                    self.v[:] = self.g
                    g_dot_g = np.dot(self.g, self.g)
                    v_dot_g = np.dot(self.v, self.g)
                if action == ctrlib._TRLIB_CLA_UPDATE_DIR:
                    self.p[:]  = flt1 * self.v + flt2 * self.p
                    self.Hp[:] = self.hessp(self.p)
                    p_dot_Hp = np.dot(self.p, self.Hp)
                    if ityp == ctrlib._TRLIB_CLT_L:
                        self.Q[it,:] = self.p
                if action == ctrlib._TRLIB_CLA_OBJVAL:
                    g_dot_g = .5*np.dot(self.s, self.hessp(self.s))
                    g_dot_g += np.dot(self.s, self.jac)
                if ret < 10:
                    break
                self.init = ctrlib._TRLIB_CLS_HOTSTART
            self.lam = self.fwork[7]
        finally:
            messages.close()

        return self.s, self.lam > 0.0
