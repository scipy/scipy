//------------------------------------------------------------------------------
// SuiteSparse_config/SuiteSparse_config.h: common utilites for SuiteSparse
//------------------------------------------------------------------------------

// SuiteSparse_config, Copyright (c) 2012-2023, Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

// Configuration file for SuiteSparse: a Suite of Sparse matrix packages: AMD,
// COLAMD, CCOLAMD, CAMD, CHOLMOD, UMFPACK, CXSparse, SuiteSparseQR, ParU, ...

// The SuiteSparse_config.h file is configured by CMake to be specific to the
// C/C++ compiler and BLAS library being used for SuiteSparse.  The original
// file is SuiteSparse_config/SuiteSparse_config.h.in.  Do not edit the
// SuiteSparse_config.h file directly.

#ifndef SUITESPARSE_CONFIG_H
#define SUITESPARSE_CONFIG_H

//------------------------------------------------------------------------------
// SuiteSparse-wide ANSI C11 #include files
//------------------------------------------------------------------------------

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <stddef.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

void *Highs_SuiteSparse_malloc    // pointer to allocated block of memory
(
    size_t nitems,          // number of items to malloc (>=1 is enforced)
    size_t size_of_item     // sizeof each item
) ;

void *Highs_SuiteSparse_free      // always returns NULL
(
    void *p                 // block to free
) ;

#define SUITESPARSE_DATE "Nov 4, 2025"
#define SUITESPARSE_MAIN_VERSION    7
#define SUITESPARSE_SUB_VERSION     12
#define SUITESPARSE_SUBSUB_VERSION  1

#endif

