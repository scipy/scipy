#ifndef IPX_C_H_
#define IPX_C_H_

#include "ipx_config.h"
#include "ipx_info.h"
#include "ipx_parameters.h"
#include "ipx_status.h"

#ifdef __cplusplus
extern "C"{
#endif
    /* Returns an ipx_parameters struct with default values. */
    struct ipx_parameters ipx_default_parameters();

    /* Allocates a new LpSolver object. On success, *p_self holds a pointer to
       the new object. If the memory allocation fails, *p_self becomes NULL
       (the required memory is tiny). The function does nothing if @p_self is
       NULL. */
    void ipx_new(void** p_self);

    /* Deallocates the LpSolver object pointed to by *p_self and sets *p_self to
       NULL. If either p_self or *p_self is NULL, the function does nothing. */
    void ipx_free(void** p_self);

    /* The remaining functions call their equivalent method of LpSolver for the
       object pointed to by @self. See src/lp_solver.h for documentation of the
       methods. */
    ipxint ipx_solve(void* self, ipxint num_var, const double* obj,
                     const double* lb, const double* ub, ipxint num_constr,
                     const ipxint* Ap, const ipxint* Ai, const double* Ax,
                     const double* rhs, const char* constr_type);
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
