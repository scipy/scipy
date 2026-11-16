#ifndef FACTOR_HIGHS_C_API_H
#define FACTOR_HIGHS_C_API_H

#include "util/HighsType.h"

/*  C API to HiPO linear solver
    It is meant to be used outside of HiGHS as a standalone linear solver.
    Refer to FactorHighs.h for a description of the functions.
    Refer to FactorHighs_c_api_example.c for a small example.
*/

#ifdef __cplusplus
extern "C" {
#endif

// Initialise and create/destroy objects
HighsInt FactorHighs_initialise(HighsInt threads);
void FactorHighs_terminate(void);
void* FactorHighs_symbolic_create(void);
void FactorHighs_symbolic_destroy(void* S);
void* FactorHighs_create(void);
void FactorHighs_destroy(void* FH);

// Main functions
HighsInt FactorHighs_analyse(void* FH, void* S, HighsInt n, HighsInt nz,
                             const HighsInt* rows, const HighsInt* ptr,
                             const HighsInt* signs, const HighsInt* perm);
HighsInt FactorHighs_factorise(void* FH, const void* S, HighsInt n, HighsInt nz,
                               const HighsInt* rows, const HighsInt* ptr,
                               const double* vals);
HighsInt FactorHighs_solve(void* FH, double* x, HighsInt k);

// Partial solves
HighsInt FactorHighs_forwardSolve(void* FH, double* x, HighsInt k);
HighsInt FactorHighs_diagSolve(void* FH, double* x, HighsInt k);
HighsInt FactorHighs_backwardSolve(void* FH, double* x, HighsInt k);

// Reordering
HighsInt FactorHighs_reorderMetis(void* FH, HighsInt n, HighsInt nz,
                                  const HighsInt* rows, const HighsInt* ptr,
                                  HighsInt* perm);
HighsInt FactorHighs_reorderAmd(void* FH, HighsInt n, HighsInt nz,
                                const HighsInt* rows, const HighsInt* ptr,
                                HighsInt* perm);
HighsInt FactorHighs_reorderRcm(void* FH, HighsInt n, HighsInt nz,
                                const HighsInt* rows, const HighsInt* ptr,
                                HighsInt* perm);

// Set options
void FactorHighs_setRegularisation(void* FH, double reg_p, double reg_d);
void FactorHighs_setBlockSize(void* FH, HighsInt nb);
void FactorHighs_setPivoting(void* FH, HighsInt pivoting);
void FactorHighs_setLogging(void* FH, HighsInt display);
void FactorHighs_setOneIndexing(void* FH, HighsInt one_indexing);

// Extract information
void FactorHighs_getRegularisation(void* FH, double* reg);
void FactorHighs_iperm(void* FH, void* S, HighsInt* ip);
void FactorHighs_inertia(void* FH, HighsInt* pos, HighsInt* neg, HighsInt* zero,
                         double tol);

// Other
void FactorHighs_symbolic_print(void* FH, void* S, HighsInt verbose);
void FactorHighs_newIter(void* FH);

#ifdef __cplusplus
}
#endif

#endif