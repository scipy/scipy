#ifndef _BASICLU_H
#define _BASICLU_H

#ifdef __cplusplus
extern "C"{
#endif

/*
 * BASICLU integer type
 *
 * All integers in BASICLU code are of type lu_int, which must be a signed
 * integer type. It is required that all integer values arising in the
 * computation can be stored in double variables and converted back to lu_int
 * without altering their value.
 *
 * LU_INT_MAX must be the maximum value of a variable of type lu_int.
 *
 * The default is 64 bit integers to make the code easily callable from Julia.
 * int64_t is optional in the C99 standard, but available on most systems.
 *
 */
#include <stdint.h>
typedef int64_t lu_int;
#define LU_INT_MAX INT64_MAX
/* #include <limits.h> */
/* typedef long lu_int; */
/* #define LU_INT_MAX LONG_MAX */

/* size of istore */
#define BASICLU_SIZE_ISTORE_1 1024
#define BASICLU_SIZE_ISTORE_M 21

/* size of xstore */
#define BASICLU_SIZE_XSTORE_1 1024
#define BASICLU_SIZE_XSTORE_M 4

/* ------------ */
/* status codes */
/* ------------ */

#define BASICLU_OK 0
#define BASICLU_REALLOCATE 1
#define BASICLU_WARNING_singular_matrix 2
#define BASICLU_ERROR_invalid_store (-1)
#define BASICLU_ERROR_invalid_call (-2)
#define BASICLU_ERROR_argument_missing (-3)
#define BASICLU_ERROR_invalid_argument (-4)
#define BASICLU_ERROR_maximum_updates (-5)
#define BASICLU_ERROR_singular_update (-6)

#define BASICLU_ERROR_invalid_object (-8)
#define BASICLU_ERROR_out_of_memory (-9)

/* ------------------------ */
/* public entries in xstore */
/* ------------------------ */

/* user parameters */
#define BASICLU_MEMORYL 1
#define BASICLU_MEMORYU 2
#define BASICLU_MEMORYW 3
#define BASICLU_DROP_TOLERANCE 4
#define BASICLU_ABS_PIVOT_TOLERANCE 5
#define BASICLU_REL_PIVOT_TOLERANCE 6
#define BASICLU_BIAS_NONZEROS 7
#define BASICLU_MAXN_SEARCH_PIVOT 8
#define BASICLU_PAD 9
#define BASICLU_STRETCH 10
#define BASICLU_COMPRESSION_THRESHOLD 11
#define BASICLU_SPARSE_THRESHOLD 12
#define BASICLU_REMOVE_COLUMNS 13
#define BASICLU_SEARCH_ROWS 14

/* user readable */
#define BASICLU_DIM 64
#define BASICLU_STATUS 65
#define BASICLU_ADD_MEMORYL 66
#define BASICLU_ADD_MEMORYU 67
#define BASICLU_ADD_MEMORYW 68

#define BASICLU_NUPDATE 70
#define BASICLU_NFORREST 71
#define BASICLU_NFACTORIZE 72
#define BASICLU_NUPDATE_TOTAL 73
#define BASICLU_NFORREST_TOTAL 74
#define BASICLU_NSYMPERM_TOTAL 75
#define BASICLU_LNZ 76
#define BASICLU_UNZ 77
#define BASICLU_RNZ 78
#define BASICLU_MIN_PIVOT 79
#define BASICLU_MAX_PIVOT 80
#define BASICLU_UPDATE_COST 81
#define BASICLU_TIME_FACTORIZE 82
#define BASICLU_TIME_SOLVE 83
#define BASICLU_TIME_UPDATE 84
#define BASICLU_TIME_FACTORIZE_TOTAL 85
#define BASICLU_TIME_SOLVE_TOTAL 86
#define BASICLU_TIME_UPDATE_TOTAL 87
#define BASICLU_LFLOPS 88
#define BASICLU_UFLOPS 89
#define BASICLU_RFLOPS 90
#define BASICLU_CONDEST_L 91
#define BASICLU_CONDEST_U 92
#define BASICLU_MAX_ETA 93
#define BASICLU_NORM_L 94
#define BASICLU_NORM_U 95
#define BASICLU_NORMEST_LINV 96
#define BASICLU_NORMEST_UINV 97
#define BASICLU_MATRIX_ONENORM 98
#define BASICLU_MATRIX_INFNORM 99
#define BASICLU_RESIDUAL_TEST 111

#define BASICLU_MATRIX_NZ 100
#define BASICLU_RANK 101
#define BASICLU_BUMP_SIZE 102
#define BASICLU_BUMP_NZ 103
#define BASICLU_NSEARCH_PIVOT 104
#define BASICLU_NEXPAND 105
#define BASICLU_NGARBAGE 106
#define BASICLU_FACTOR_FLOPS 107
#define BASICLU_TIME_SINGLETONS 108
#define BASICLU_TIME_SEARCH_PIVOT 109
#define BASICLU_TIME_ELIM_PIVOT 110

#define BASICLU_PIVOT_ERROR 120

/* ----------------------- */
/* user callable functions */
/* ----------------------- */

#include "basiclu_initialize.h"
#include "basiclu_factorize.h"
#include "basiclu_get_factors.h"
#include "basiclu_solve_dense.h"
#include "basiclu_solve_sparse.h"
#include "basiclu_solve_for_update.h"
#include "basiclu_update.h"

#include "basiclu_object.h"
#include "basiclu_obj_initialize.h"
#include "basiclu_obj_factorize.h"
#include "basiclu_obj_get_factors.h"
#include "basiclu_obj_solve_dense.h"
#include "basiclu_obj_solve_sparse.h"
#include "basiclu_obj_solve_for_update.h"
#include "basiclu_obj_update.h"
#include "basiclu_obj_free.h"

#ifdef __cplusplus
}
#endif

#endif
