#ifndef SCIPY_SLU_CONFIG_H
#define SCIPY_SLU_CONFIG_H

#include <stdlib.h>

/*
 * Support routines: route SuperLU's allocation and abort through Python.
 * These are implemented by the `_superlu` extension module, in
 * `scipy/sparse/linalg/_dsolve/_superlu_utils.c`, and resolved when it links
 * against this library.
 */
void superlu_python_module_abort(char *msg);
void *superlu_python_module_malloc(size_t size);
void superlu_python_module_free(void *ptr);

#define USER_ABORT  superlu_python_module_abort
#define USER_MALLOC superlu_python_module_malloc
#define USER_FREE   superlu_python_module_free

/*
 * BLAS integer type: matches the BLAS library's integer ABI.
 * When HAVE_BLAS_ILP64 is defined (ILP64 build), this is int64_t;
 * otherwise it is int (LP64, the default).
 *
 * This is a macro rather than a typedef: upstream's superlu_config.h also
 * declares `slu_blasint` (as `int`) unless it is already defined as a macro,
 * and it may be included before this header (e.g. via slu_ddefs.h). Defining
 * a macro overrides upstream's typedef for all later uses in either order.
 */
#ifndef slu_blasint
#ifdef HAVE_BLAS_ILP64
#include <stdint.h>
#define slu_blasint int64_t
#else
#define slu_blasint int
#endif
#endif

#endif
