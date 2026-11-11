#ifndef IPX_C_H_
#define IPX_C_H_

#include "ipm/ipx/ipx_config.h"
#include "ipm/ipx/ipx_info.h"
#include "ipm/ipx/ipx_parameters.h"
#include "ipm/ipx/ipx_status.h"

#ifdef __cplusplus
extern "C"{
#endif
    /* These functions call their equivalent method of LpSolver for
       the object pointed to by @self. See highs/lp_solver.h for
       documentation of the methods. */
    ipxint ipx_load_model(void* self, ipxint num_var, const double* obj,
                          const double* lb, const double* ub, ipxint num_constr,
                          const ipxint* Ap, const ipxint* Ai, const double* Ax,
                          const double* rhs, const char* constr_type);
    ipxint ipx_load_ipm_starting_point(void* self, const double* x,
                                       const double* xl, const double* xu,
                                       const double* slack, const double* y,
                                       const double* zl, const double* zu);
    ipxint ipx_solve(void* self);
    struct ipx_info ipx_get_info(void* self);
    ipxint ipx_get_interior_solution(void* self, double* x, double* xl,
                                     double* xu, double* slack, double* y,
                                     double* zl, double* zu);
    ipxint ipx_get_basic_solution(void* self, double* x, double* slack,
                                  double* y, double* z,
                                  ipxint* cbasis, ipxint* vbasis);
    struct ipx_parameters ipx_get_parameters(void* self);
    void ipx_set_parameters(void* self, struct ipx_parameters);
    void ipx_clear_model(void* self);

    /* for debugging */
    ipxint ipx_get_iterate(void* self, double* x, double* y, double* zl,
                           double* zu, double* xl, double* xu);
    ipxint ipx_get_basis(void* self, ipxint* cbasis, ipxint* vbasis);
    ipxint ipx_get_kktmatrix(void* self, ipxint* AIp, ipxint* AIi, double* AIx,
                             double* g);
    ipxint ipx_symbolic_invert(void* self, ipxint* rowcounts,
                               ipxint* colcounts);
#ifdef __cplusplus
}
#endif

#endif  /* IPX_C_H_ */
