#ifndef SCIPY_SLU_CONFIG_H
#define SCIPY_SLU_CONFIG_H

#include <stdlib.h>

/*
 * Support routines
 */
void superlu_python_module_abort(char *msg);
void *superlu_python_module_malloc(size_t size);
void superlu_python_module_free(void *ptr);

#define USER_ABORT  superlu_python_module_abort
#define USER_MALLOC superlu_python_module_malloc
#define USER_FREE   superlu_python_module_free

#define SCIPY_SPECIFIC_FIX 1

/* 
 * Fortran configuration 
 */
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define UpCase 1
#else
#define NoChange 1
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#error Uppercase and trailing slash in Fortran names not supported
#else
#define Add_ 1
#endif
#endif

#endif
