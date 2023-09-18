/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2.h                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2022 Case Van Horsen                                   *
 *                                                                         *
 * This file is part of GMPY2.                                             *
 *                                                                         *
 * GMPY2 is free software: you can redistribute it and/or modify it under  *
 * the terms of the GNU Lesser General Public License as published by the  *
 * Free Software Foundation, either version 3 of the License, or (at your  *
 * option) any later version.                                              *
 *                                                                         *
 * GMPY2 is distributed in the hope that it will be useful, but WITHOUT    *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public    *
 * License for more details.                                               *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with GMPY2; if not, see <http://www.gnu.org/licenses/>    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*
  gmpy C API extension header file.
  Part of Python's gmpy module since version 0.4

  Created by Pearu Peterson <pearu@cens.ioc.ee>, November 2000.
  Edited by A. Martelli <aleaxit@yahoo.com>, December 2000.
  Edited by Case Van Horsen <casevh@gmail.com>, 2009, 2010, 2011.

  Version 1.02, February 2007.
  Version 1.03, June 2008
  Version 1.04, June 2008 (no changes)
  Version 1.05, February 2009 (support MPIR)
  Version 1.20, January 2010 (remove obsolete MS hacks) casevh
  Version 2.00, April 2010 (change to gmpy2) casevh
                October 2010 (added Py_hash_t) casevh
                December 2010 (added mpfr, mpc) casevh
                January 2011 (add Pygmpy_context) casevh
                April 2011 (split into multiple files) casevh
  Version 2.10  August 2014 (reflect major rewrite during 2013/2014) casevh
 */

#ifndef Py_GMPYMODULE_H
#define Py_GMPYMODULE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Structure of gmpy2.h
 *
 * Revised 17-APR-2017 casevh
 *
 * 1. Checks for specific Python versions.
 * 2. Include headers for GMP/MPIR, MPFR, and MPC.
 * 3. Define public C-API.
 *    1. Define gmpy2 types.
 *    2. Define the public API.
 *
 */

/* Check for minimum Python version requirements. */

#if PY_VERSION_HEX < 0x02060000
#  error "GMPY2 requires Python 2.6 or later."
#endif

/* Include headers for GMP/MPIR, MPFR, and MPC. */

#ifdef MPIR
#  include <mpir.h>
#else
#  include <gmp.h>
#endif

#include <mpfr.h>
#include <mpc.h>

/* Check MPFR and MPC versions. */

#if (!defined(MPC_VERSION) || (MPC_VERSION < MPC_VERSION_NUM(1,0,3)))
#  error "GMPY2 requires MPC 1.0.3 or later."
#endif

#if (defined(MPC_VERSION) && (MPC_VERSION >= MPC_VERSION_NUM(1,1,0)))
#  define MPC_110
#endif


#if PY_VERSION_HEX < 0x030200A4
typedef long Py_hash_t;
typedef unsigned long Py_uhash_t;
#  define _PyHASH_IMAG 1000003
#endif

/* GMPY2 Public API */

/* Types
 *    MPZ_Object
 *    XMPZ_Object        (mutable version of MPZ_Object)
 *    MPQ_Object
 *    XMPQ_Object        (mutable version of MPQ_Object)
 *    MPFR_Object
 *    XMPFR_Object       (mutable version of MPFR_Object)
 *    MPC_Object
 *    XMPC_Object        (mutable version of MPC_Object)
 *    CTXT_Object
 *    CTXT_Manager_Object
 *    RandomState_Object
 */

typedef struct {
    PyObject_HEAD
    mpz_t z;
    Py_hash_t hash_cache;
} MPZ_Object;

typedef struct {
    PyObject_HEAD
    mpz_t z;
} XMPZ_Object;

typedef struct {
    PyObject_HEAD
    mpq_t q;
    Py_hash_t  hash_cache;
} MPQ_Object;

typedef struct {
    PyObject_HEAD
    mpfr_t f;
    Py_hash_t hash_cache;
    int rc;
} MPFR_Object;

typedef struct {
    PyObject_HEAD
    mpc_t c;
    Py_hash_t hash_cache;
    int rc;
} MPC_Object;

typedef struct {
    PyObject_HEAD
    gmp_randstate_t state;
} RandomState_Object;

typedef struct {
    mpfr_prec_t mpfr_prec;   /* current precision in bits, for MPFR */
    mpfr_rnd_t mpfr_round;   /* current rounding mode for float (MPFR) */
    mpfr_exp_t emax;         /* maximum exponent */
    mpfr_exp_t emin;         /* minimum exponent */
    int subnormalize;        /* if 1, subnormalization is performed */
    int underflow;           /* did an underflow occur? */
    int overflow;            /* did an overflow occur? */
    int inexact;             /* was the result inexact? */
    int invalid;             /* invalid operation (i.e. NaN)? */
    int erange;              /* did a range error occur? */
    int divzero;             /* divided by zero? */
    int traps;               /* if 0, do not trap any exceptions */
                             /* if not 0, then raise traps per bits above  */
    mpfr_prec_t real_prec;   /* current precision in bits, for Re(MPC) */
    mpfr_prec_t imag_prec;   /* current precision in bits, for Im(MPC) */
    mpfr_rnd_t real_round;   /* current rounding mode for Re(MPC) */
    mpfr_rnd_t imag_round;   /* current rounding mode for Im(MPC) */
    int allow_complex;       /* if 1, allow mpfr functions to return an mpc */
    int rational_division;   /* if 1, mpz/mpz returns an mpq result */
    int allow_release_gil;   /* if 1, allow mpz functions to release the GIL */
} gmpy_context;

typedef struct {
    PyObject_HEAD
    gmpy_context ctx;
#ifndef WITHOUT_THREADS
    PyThreadState *tstate;
#endif
} CTXT_Object;

typedef struct {
    PyObject_HEAD
    CTXT_Object *new_context; /* Context that will be returned when
                               * __enter__ is called. */
    CTXT_Object *old_context; /* Context that will restored when
                               * __exit__ is called. */
} CTXT_Manager_Object;

#define MPZ(obj)  (((MPZ_Object*)(obj))->z)
#define MPQ(obj)  (((MPQ_Object*)(obj))->q)
#define MPFR(obj) (((MPFR_Object*)(obj))->f)
#define MPC(obj)  (((MPC_Object*)(obj))->c)

/* Start of the C-API definitions */

#define MPZ_Type_NUM          0
#define XMPZ_Type_NUM         1
#define MPQ_Type_NUM          2
#define XMPQ_Type_NUM         3
#define MPFR_Type_NUM         4
#define XMPFR_Type_NUM        5
#define MPC_Type_NUM          6
#define XMPC_Type_NUM         7
#define CTXT_Type_NUM         8
#define CTXT_Manager_Type_NUM 9
#define RandomState_Type_NUM  10

/* The following functions are found in gmpy2_cache. */

#define GMPy_MPZ_New_NUM            11
#define GMPy_MPZ_New_RETURN         MPZ_Object *
#define GMPy_MPZ_New_PROTO          (CTXT_Object *context)

#define GMPy_MPZ_NewInit_NUM        12
#define GMPy_MPZ_NewInit_RETURN     PyObject *
#define GMPy_MPZ_NewInit_PROTO      (PyTypeObject *type, PyObject *args, PyObject *keywds)

#define GMPy_MPZ_Dealloc_NUM        13
#define GMPy_MPZ_Dealloc_RETURN     void
#define GMPy_MPZ_Dealloc_PROTO      (MPZ_Object *self)

/* The following function is found in gmpy2_convert_gmp. */

#define GMPy_MPZ_ConvertArg_NUM     14
#define GMPy_MPZ_ConvertArg_RETURN  int
#define GMPy_MPZ_ConvertArg_PROTO   (PyObject *arg, PyObject **ptr)

/* The following functions are found in gmpy2_cache. */

#define GMPy_XMPZ_New_NUM           15
#define GMPy_XMPZ_New_RETURN        XMPZ_Object *
#define GMPy_XMPZ_New_PROTO         (CTXT_Object *context)

#define GMPy_XMPZ_NewInit_NUM       16
#define GMPy_XMPZ_NewInit_RETURN    PyObject *
#define GMPy_XMPZ_NewInit_PROTO     (PyTypeObject *type, PyObject *args, PyObject *keywds)

#define GMPy_XMPZ_Dealloc_NUM       17
#define GMPy_XMPZ_Dealloc_RETURN    void
#define GMPy_XMPZ_Dealloc_PROTO     (XMPZ_Object *self)

/* The following functions are found in gmpy2_cache. */

#define GMPy_MPQ_New_NUM            18
#define GMPy_MPQ_New_RETURN         MPQ_Object *
#define GMPy_MPQ_New_PROTO          (CTXT_Object *context)

#define GMPy_MPQ_NewInit_NUM        19
#define GMPy_MPQ_NewInit_RETURN     PyObject *
#define GMPy_MPQ_NewInit_PROTO      (PyTypeObject *type, PyObject *args, PyObject *keywds)

#define GMPy_MPQ_Dealloc_NUM        20
#define GMPy_MPQ_Dealloc_RETURN     void
#define GMPy_MPQ_Dealloc_PROTO      (MPQ_Object *self)

/* The following function is found in gmpy2_convert_gmp. */

#define GMPy_MPQ_ConvertArg_NUM     21
#define GMPy_MPQ_ConvertArg_RETURN  int
#define GMPy_MPQ_ConvertArg_PROTO   (PyObject *arg, PyObject **ptr)

/* The following functions are found in gmpy2_cache. */

#define GMPy_MPFR_New_NUM           22
#define GMPy_MPFR_New_RETURN        MPFR_Object *
#define GMPy_MPFR_New_PROTO         (mpfr_prec_t bits, CTXT_Object *context)

#define GMPy_MPFR_NewInit_NUM       23
#define GMPy_MPFR_NewInit_RETURN    PyObject *
#define GMPy_MPFR_NewInit_PROTO     (PyTypeObject *type, PyObject *args, PyObject *keywds)

#define GMPy_MPFR_Dealloc_NUM       24
#define GMPy_MPFR_Dealloc_RETURN    void
#define GMPy_MPFR_Dealloc_PROTO     (MPFR_Object *self)

/* The following function is found in gmpy2_convert_gmp. */

#define GMPy_MPFR_ConvertArg_NUM    25
#define GMPy_MPFR_ConvertArg_RETURN int
#define GMPy_MPFR_ConvertArg_PROTO  (PyObject *arg, PyObject **ptr)

/* The following functions are found in gmpy2_cache. */

#define GMPy_MPC_New_NUM             26
#define GMPy_MPC_New_RETURN          MPC_Object *
#define GMPy_MPC_New_PROTO           (mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context)

#define GMPy_MPC_NewInit_NUM         27
#define GMPy_MPC_NewInit_RETURN      PyObject *
#define GMPy_MPC_NewInit_PROTO       (PyTypeObject *type, PyObject *args, PyObject *keywds)

#define GMPy_MPC_Dealloc_NUM        28
#define GMPy_MPC_Dealloc_RETURN     void
#define GMPy_MPC_Dealloc_PROTO      (MPC_Object *self)

/* The following function is found in gmpy2_convert_gmp. */

#define GMPy_MPC_ConvertArg_NUM     29
#define GMPy_MPC_ConvertArg_RETURN  int
#define GMPy_MPC_ConvertArg_PROTO   (PyObject *arg, PyObject **ptr)

/* Total number of C-API pointers. */

#define GMPy_API_pointers 30

/* End of C-API definitions. */

#ifdef GMPY2_MODULE

/* Define various macros to deal with differences between Python 2 and 3. */

#if (PY_MAJOR_VERSION == 3)
#define PY3
#define Py2or3String_FromString     PyUnicode_FromString
#define Py2or3String_FromFormat     PyUnicode_FromFormat
#define Py2or3String_Check          PyUnicode_Check
#define Py2or3String_Format         PyUnicode_Format
#define Py2or3String_Type           Py_UCS4
#define Py2or3String_1Char(obj)     (PyUnicode_READY(obj) ? (Py_UCS4)0 : PyUnicode_READ_CHAR(obj, 0))
#define PyStrOrUnicode_Check(op)    (PyBytes_Check(op) || PyUnicode_Check(op))
#define PyIntOrLong_FromLong        PyLong_FromLong
#define PyIntOrLong_Check(op)       (PyLong_Check(op))
#define PyIntOrLong_CheckExact(op)  PyLong_CheckExact(op)
#define PyIntOrLong_FromSize_t      PyLong_FromSize_t
#define PyIntOrLong_FromSsize_t     PyLong_FromSsize_t
#define PyIntOrLong_AsSsize_t       PyLong_AsSsize_t
#define PyIntOrLong_AsLong          PyLong_AsLong
#else
#define PY2
#define Py2or3String_FromString     PyString_FromString
#define Py2or3String_FromFormat     PyString_FromFormat
#define Py2or3String_Check          PyString_Check
#define Py2or3String_Format         PyString_Format
#define Py2or3String_Type           char
#define Py2or3String_1Char(obj)     PyString_AsString(obj)[0]
#define PyStrOrUnicode_Check(op)    (PyString_Check(op) || PyUnicode_Check(op))
#define PyIntOrLong_FromLong        PyInt_FromLong
#define PyIntOrLong_Check(op)       (PyInt_Check(op) || PyLong_Check(op))
#define PyIntOrLong_CheckExact(op)  (PyInt_CheckExact(op) || PyLong_CheckExact(op))
#define PyIntOrLong_FromSize_t      PyInt_FromSize_t
#define PyIntOrLong_FromSsize_t     PyInt_FromSsize_t
#define PyIntOrLong_AsSsize_t       PyInt_AsSsize_t
#define PyIntOrLong_AsLong          PyInt_AsLong
#endif

#ifndef ABS
#  define ABS(a)  (((a) < 0) ? -(a) : (a))
#endif

#if defined(MS_WIN32) && defined(_MSC_VER)
   /* so one won't need to link explicitly to gmp.lib...: */
#  if defined(MPIR)
#    pragma comment(lib,"mpir.lib")
#  else
#    pragma comment(lib,"gmp.lib")
#  endif
#  define USE_ALLOCA 1
#  define inline __inline
#endif

#ifdef __GNUC__
#  define USE_ALLOCA 1
#endif

#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#   ifdef _MSC_VER
#    include <malloc.h>
#    define alloca _alloca
#   else
#    if HAVE_ALLOCA_H
#     include <alloca.h>
#    else
       char *alloca ();
#    endif
#   endif
# endif
#endif

#define ALLOC_THRESHOLD 8192

#define INDEX_ERROR(msg)    PyErr_SetString(PyExc_IndexError, msg)
#define TYPE_ERROR(msg)     PyErr_SetString(PyExc_TypeError, msg)
#define VALUE_ERROR(msg)    PyErr_SetString(PyExc_ValueError, msg)
#define ZERO_ERROR(msg)     PyErr_SetString(PyExc_ZeroDivisionError, msg)
#define SYSTEM_ERROR(msg)   PyErr_SetString(PyExc_SystemError, msg)
#define OVERFLOW_ERROR(msg) PyErr_SetString(PyExc_OverflowError, msg)
#define RUNTIME_ERROR(msg)  PyErr_SetString(PyExc_RuntimeError, msg)

#define GMPY_DEFAULT -1

/* To prevent excessive memory usage, we don't want to save very large
 * numbers in the cache. The default value specified in the options
 * structure is 128 words (512 bytes on 32-bit platforms, 1024 bytes on
 * 64-bit platforms).
 */
#define MAX_CACHE_LIMBS 16384

/* The maximum number of objects that can be saved in a cache is specified
 * here. The default value is 100.*/
#define MAX_CACHE 1000

#ifdef USE_ALLOCA
#  define TEMP_ALLOC(B, S)     \
    if(S < ALLOC_THRESHOLD) {  \
        B = alloca(S);         \
    } else {                   \
        if(!(B = malloc(S))) { \
            PyErr_NoMemory();  \
            return NULL;       \
        }                      \
    }
#  define TEMP_FREE(B, S) if(S >= ALLOC_THRESHOLD) free(B)
#else
#  define TEMP_ALLOC(B, S)     \
    if(!(B = malloc(S)))  {    \
        PyErr_NoMemory();      \
        return NULL;           \
    }
#  define TEMP_FREE(B, S) free(B)
#endif

/* Various defs to mask differences between Python versions. */
#if PY_VERSION_HEX < 0x03050000
#define Py_RETURN_NOTIMPLEMENTED \
    return Py_INCREF(Py_NotImplemented), Py_NotImplemented
#endif

#ifndef Py_SIZE
#  define Py_SIZE(ob)     (((PyVarObject*)(ob))->ob_size)
#endif

#ifndef Py_TYPE
#  define Py_TYPE(ob)     (((PyObject*)(ob))->ob_type)
#endif

/* Import a collection of general purpose macros. */

#include "gmpy2_macros.h"

/* Import the files that complete the definition of the types defined above. */

#include "gmpy2_mpz.h"
#include "gmpy2_xmpz.h"
#include "gmpy2_mpq.h"
#include "gmpy2_mpfr.h"
#include "gmpy2_mpc.h"
#include "gmpy2_context.h"
#include "gmpy2_random.h"

/* Import the header files that provide the various functions. */

/* Support object caching, creation, and deletion. */

#include "gmpy2_cache.h"

/* Suport for miscellaneous functions (ie. version, license, etc.). */

#include "gmpy2_misc.h"

/* Support conversion to/from binary format. */

#include "gmpy2_binary.h"

/* Support for mpz/xmpz specific functions. */

#include "gmpy2_convert.h"
#include "gmpy2_convert_utils.h"
#include "gmpy2_convert_gmp.h"
#include "gmpy2_convert_mpfr.h"
#include "gmpy2_convert_mpc.h"

#include "gmpy2_mpz_divmod.h"
#include "gmpy2_mpz_divmod2exp.h"
#include "gmpy2_mpz_pack.h"
#include "gmpy2_mpz_bitops.h"
#include "gmpy2_mpz_misc.h"

#include "gmpy2_xmpz_inplace.h"
#include "gmpy2_xmpz_misc.h"
#include "gmpy2_xmpz_limbs.h"

/* Support for mpq specific functions. */

#include "gmpy2_mpq_misc.h"

/* Support for mpfr specific functions. */

#include "gmpy2_mpfr_misc.h"

/* Support for mpc specific functions. */

#include "gmpy2_mpc_misc.h"

/* Support Lucas sequences. */

#include "gmpy_mpz_lucas.h"

/* Support probable-prime tests. */

#include "gmpy_mpz_prp.h"

/* Support higher-level Python methods and functions; generally not
 * specific to a single type.
 */

#include "gmpy2_abs.h"
#include "gmpy2_add.h"
#include "gmpy2_divmod.h"
#include "gmpy2_floordiv.h"
#include "gmpy2_minus.h"
#include "gmpy2_mod.h"
#include "gmpy2_mul.h"
#include "gmpy2_plus.h"
#include "gmpy2_pow.h"
#include "gmpy2_sub.h"
#include "gmpy2_truediv.h"
#include "gmpy2_math.h"
#include "gmpy2_const.h"
#include "gmpy2_square.h"
#include "gmpy2_format.h"
#include "gmpy2_hash.h"
#include "gmpy2_fused.h"
#include "gmpy2_muldiv_2exp.h"
#include "gmpy2_predicate.h"
#include "gmpy2_sign.h"
#include "gmpy2_richcompare.h"
#include "gmpy2_cmp.h"

#ifdef VECTOR
#  include "gmpy2_vector.h"
#endif /* defined(VECTOR) */

#else /* defined(GMPY2_MODULE) */

/* This section is used for other C-coded modules that use gmpy2's API. */

static void **GMPy_C_API;

#define MPZ_Check(op)    ((op)->ob_type == (PyTypeObject*)GMPy_C_API[MPZ_Type_NUM])
#define XMPZ_Check(op)   ((op)->ob_type == (PyTypeObject*)GMPy_C_API[XMPZ_Type_NUM])
#define MPQ_Check(op)    ((op)->ob_type == (PyTypeObject*)GMPy_C_API[MPQ_Type_NUM])
#define XMPQ_Check(op)   ((op)->ob_type == (PyTypeObject*)GMPy_C_API[XMPQ_Type_NUM])
#define MPFR_Check(op)   ((op)->ob_type == (PyTypeObject*)GMPy_C_API[MPFR_Type_NUM])
#define XMPFR_Check(op)  ((op)->ob_type == (PyTypeObject*)GMPy_C_API[XMPFR_Type_NUM])
#define MPC_Check(op)    ((op)->ob_type == (PyTypeObject*)GMPy_C_API[MPC_Type_NUM])
#define XMPC_Check(op)   ((op)->ob_type == (PyTypeObject*)GMPy_C_API[XMPC_Type_NUM])

#define GMPy_MPZ_New         (*(GMPy_MPZ_New_RETURN         (*)GMPy_MPZ_New_PROTO)         GMPy_C_API[GMPy_MPZ_New_NUM])
#define GMPy_MPZ_NewInit     (*(GMPy_MPZ_NewInit_RETURN     (*)GMPy_MPZ_NewInit_PROTO)     GMPy_C_API[GMPy_MPZ_NewInit_NUM])
#define GMPy_MPZ_Dealloc     (*(GMPy_MPZ_Dealloc_RETURN     (*)GMPy_MPZ_Dealloc_PROTO)     GMPy_C_API[GMPy_MPZ_Dealloc_NUM])
#define GMPy_MPZ_ConvertArg  (*(GMPy_MPZ_ConvertArg_RETURN  (*)GMPy_MPZ_ConvertArg_PROTO)  GMPy_C_API[GMPy_MPZ_ConvertArg_NUM])

#define GMPy_XMPZ_New        (*(GMPy_XMPZ_New_RETURN        (*)GMPy_XMPZ_New_PROTO)        GMPy_C_API[GMPy_XMPZ_New_NUM])
#define GMPy_XMPZ_NewInit    (*(GMPy_XMPZ_NewInit_RETURN    (*)GMPy_XMPZ_NewInit_PROTO)    GMPy_C_API[GMPy_XMPZ_NewInit_NUM])
#define GMPy_XMPZ_Dealloc    (*(GMPy_XMPZ_Dealloc_RETURN    (*)GMPy_XMPZ_Dealloc_PROTO)    GMPy_C_API[GMPy_XMPZ_Dealloc_NUM])

#define GMPy_MPQ_New         (*(GMPy_MPQ_New_RETURN         (*)GMPy_MPQ_New_PROTO)         GMPy_C_API[GMPy_MPQ_New_NUM])
#define GMPy_MPQ_NewInit     (*(GMPy_MPQ_NewInit_RETURN     (*)GMPy_MPQ_NewInit_PROTO)     GMPy_C_API[GMPy_MPQ_NewInit_NUM])
#define GMPy_MPQ_Dealloc     (*(GMPy_MPQ_Dealloc_RETURN     (*)GMPy_MPQ_Dealloc_PROTO)     GMPy_C_API[GMPy_MPQ_Dealloc_NUM])
#define GMPy_MPQ_ConvertArg  (*(GMPy_MPQ_ConvertArg_RETURN  (*)GMPy_MPQ_ConvertArg_PROTO)  GMPy_C_API[GMPy_MPQ_ConvertArg_NUM])

#define GMPy_MPFR_New        (*(GMPy_MPFR_New_RETURN        (*)GMPy_MPFR_New_PROTO)        GMPy_C_API[GMPy_MPFR_New_NUM])
#define GMPy_MPFR_NewInit    (*(GMPy_MPFR_NewInit_RETURN    (*)GMPy_MPFR_NewInit_PROTO)    GMPy_C_API[GMPy_MPFR_NewInit_NUM])
#define GMPy_MPFR_Dealloc    (*(GMPy_MPFR_Dealloc_RETURN    (*)GMPy_MPFR_Dealloc_PROTO)    GMPy_C_API[GMPy_MPFR_Dealloc_NUM])
#define GMPy_MPFR_ConvertArg (*(GMPy_MPFR_ConvertArg_RETURN (*)GMPy_MPFR_ConvertArg_PROTO) GMPy_C_API[GMPy_MPFR_ConvertArg_NUM])

#define GMPy_MPC_New         (*(GMPy_MPC_New_RETURN         (*)GMPy_MPC_New_PROTO)         GMPy_C_API[GMPy_MPC_New_NUM])
#define GMPy_MPC_NewInit     (*(GMPy_MPC_NewInit_RETURN     (*)GMPy_MPC_NewInit_PROTO)     GMPy_C_API[GMPy_MPC_NewInit_NUM])
#define GMPy_MPC_Dealloc     (*(GMPy_MPC_Dealloc_RETURN     (*)GMPy_MPC_Dealloc_PROTO)     GMPy_C_API[GMPy_MPC_Dealloc_NUM])
#define GMPy_MPC_ConvertArg  (*(GMPy_MPC_ConvertArg_RETURN  (*)GMPy_MPC_ConvertArg_PROTO)  GMPy_C_API[GMPy_MPC_ConvertArg_NUM])

static int
import_gmpy2(void)
{
    GMPy_C_API = (void **)PyCapsule_Import("gmpy2._C_API", 0);
    return (GMPy_C_API != NULL) ? 0 : -1;
}

#endif /* defined(GMPY2_MODULE) */

#ifdef __cplusplus
}
#endif /* defined(__cplusplus */
#endif /* !defined(Py_GMPYMODULE_H */
