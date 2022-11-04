/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! @file slu_Cnames.h
 * \brief Macros defining how C routines will be called
 *
 * <pre>
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 1, 1997
 *
 * These macros define how C routines will be called.  ADD_ assumes that
 * they will be called by fortran, which expects C routines to have an
 * underscore postfixed to the name (Suns, and the Intel expect this).
 * NOCHANGE indicates that fortran will be calling, and that it expects
 * the name called by fortran to be identical to that compiled by the C
 * (RS6K's do this).  UPCASE says it expects C routines called by fortran
 * to be in all upcase (CRAY wants this). 
 * </pre>
 */
#ifndef __SUPERLU_CNAMES /* allow multiple inclusions */
#define __SUPERLU_CNAMES

#include "scipy_slu_config.h"

#define ADD_       0
#define ADD__      1
#define NOCHANGE   2
#define UPCASE     3
#define OLD_CRAY   4
#define C_CALL     5

#ifdef UpCase
#define F77_CALL_C  UPCASE
#endif

#ifdef NoChange
#define F77_CALL_C  NOCHANGE
#endif

#ifdef Add_
#define F77_CALL_C  ADD_
#endif

#ifdef Add__
#define F77_CALL_C  ADD__
#endif

#ifdef _CRAY
#define F77_CALL_C  OLD_CRAY
#endif

/* Default */
#ifndef F77_CALL_C
#define F77_CALL_C  ADD_
#endif


#if (F77_CALL_C == ADD_)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine
 * No redefinition necessary to have following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm_(...)
 *
 * This is the default.
 */

#endif

#if (F77_CALL_C == ADD__)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm__(...)
 */
/* BLAS */
#define sswap_    sswap__
#define saxpy_    saxpy__
#define sasum_    sasum__
#define isamax_   isamax__
#define scopy_    scopy__
#define sscal_    sscal__
#define sger_     sger__
#define snrm2_    snrm2__
#define ssymv_    ssymv__
#define sdot_     sdot__
#define saxpy_    saxpy__
#define ssyr2_    ssyr2__
#define srot_     srot__
#define sgemv_    sgemv__
#define strsv_    strsv__
#define sgemm_    sgemm__
#define strsm_    strsm__

#define dswap_    dswap__
#define daxpy_    daxpy__
#define dasum_    dasum__
#define idamax_   idamax__
#define dcopy_    dcopy__
#define dscal_    dscal__
#define dger_     dger__
#define dnrm2_    dnrm2__
#define dsymv_    dsymv__
#define ddot_     ddot__
#define dsyr2_    dsyr2__
#define drot_     drot__
#define dgemv_    dgemv__
#define dtrsv_    dtrsv__
#define dgemm_    dgemm__
#define dtrsm_    dtrsm__

#define cswap_    cswap__
#define caxpy_    caxpy__
#define scasum_   scasum__
#define icamax_   icamax__
#define ccopy_    ccopy__
#define cscal_    cscal__
#define scnrm2_   scnrm2__
#define caxpy_    caxpy__
#define cgemv_    cgemv__
#define ctrsv_    ctrsv__
#define cgemm_    cgemm__
#define ctrsm_    ctrsm__
#define cgerc_    cgerc__
#define chemv_    chemv__
#define cher2_    cher2__

#define zswap_    zswap__
#define zaxpy_    zaxpy__
#define dzasum_   dzasum__
#define izamax_   izamax__
#define zcopy_    zcopy__
#define zscal_    zscal__
#define dznrm2_   dznrm2__
#define zaxpy_    zaxpy__
#define zgemv_    zgemv__
#define ztrsv_    ztrsv__
#define zgemm_    zgemm__
#define ztrsm_    ztrsm__
#define zgerc_    zgerc__
#define zhemv_    zhemv__
#define zher2_    zher2__

/* LAPACK */
#define dlacon_   dlacon__
#define slacon_   slacon__
#define icmax1_   icmax1__
#define scsum1_   scsum1__
#define clacon_   clacon__
#define dzsum1_   dzsum1__
#define izmax1_   izmax1__
#define zlacon_   zlacon__

/* Fortran interface */
#define c_bridge_dgssv_ c_bridge_dgssv__
#define c_fortran_sgssv_ c_fortran_sgssv__
#define c_fortran_dgssv_ c_fortran_dgssv__
#define c_fortran_cgssv_ c_fortran_cgssv__
#define c_fortran_zgssv_ c_fortran_zgssv__
#endif

#if (F77_CALL_C == UPCASE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void DGEMM(...)
 */
/* BLAS */
#define sswap_    SSWAP
#define saxpy_    SAXPY
#define sasum_    SASUM
#define isamax_   ISAMAX
#define scopy_    SCOPY
#define sscal_    SSCAL
#define sger_     SGER
#define snrm2_    SNRM2
#define ssymv_    SSYMV
#define sdot_     SDOT
#define saxpy_    SAXPY
#define ssyr2_    SSYR2
#define srot_     SROT
#define sgemv_    SGEMV
#define strsv_    STRSV
#define sgemm_    SGEMM
#define strsm_    STRSM

#define dswap_    DSWAP
#define daxpy_    DAXPY
#define dasum_    DASUM
#define idamax_   IDAMAX
#define dcopy_    DCOPY
#define dscal_    DSCAL
#define dger_     DGER
#define dnrm2_    DNRM2
#define dsymv_    DSYMV
#define ddot_     DDOT
#define dsyr2_    DSYR2
#define drot_     DROT
#define dgemv_    DGEMV
#define dtrsv_    DTRSV
#define dgemm_    DGEMM
#define dtrsm_    DTRSM

#define cswap_    CSWAP
#define caxpy_    CAXPY
#define scasum_   SCASUM
#define icamax_   ICAMAX
#define ccopy_    CCOPY
#define cscal_    CSCAL
#define scnrm2_   SCNRM2
#define cgemv_    CGEMV
#define ctrsv_    CTRSV
#define cgemm_    CGEMM
#define ctrsm_    CTRSM
#define cgerc_    CGERC
#define chemv_    CHEMV
#define cher2_    CHER2

#define zswap_    ZSWAP
#define zaxpy_    ZAXPY
#define dzasum_   DZASUM
#define izamax_   IZAMAX
#define zcopy_    ZCOPY
#define zscal_    ZSCAL
#define dznrm2_   DZNRM2
#define zgemv_    ZGEMV
#define ztrsv_    ZTRSV
#define zgemm_    ZGEMM
#define ztrsm_    ZTRSM
#define zgerc_    ZGERC
#define zhemv_    ZHEMV
#define zher2_    ZHER2

/* LAPACK */
#define dlacon_   DLACON
#define slacon_   SLACON
#define icmax1_   ICMAX1
#define scsum1_   SCSUM1
#define clacon_   CLACON
#define dzsum1_   DZSUM1
#define izmax1_   IZMAX1
#define zlacon_   ZLACON

/* Fortran interface */
#define c_bridge_dgssv_ C_BRIDGE_DGSSV
#define c_fortran_sgssv_ C_FORTRAN_SGSSV
#define c_fortran_dgssv_ C_FORTRAN_DGSSV
#define c_fortran_cgssv_ C_FORTRAN_CGSSV
#define c_fortran_zgssv_ C_FORTRAN_ZGSSV
#endif


#if (F77_CALL_C == OLD_CRAY)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void SGEMM(...)
 */
/* BLAS */
#define sswap_    SSWAP
#define saxpy_    SAXPY
#define sasum_    SASUM
#define isamax_   ISAMAX
#define scopy_    SCOPY
#define sscal_    SSCAL
#define sger_     SGER
#define snrm2_    SNRM2
#define ssymv_    SSYMV
#define sdot_     SDOT
#define ssyr2_    SSYR2
#define srot_     SROT
#define sgemv_    SGEMV
#define strsv_    STRSV
#define sgemm_    SGEMM
#define strsm_    STRSM

#define dswap_    SSWAP
#define daxpy_    SAXPY
#define dasum_    SASUM
#define idamax_   ISAMAX
#define dcopy_    SCOPY
#define dscal_    SSCAL
#define dger_     SGER
#define dnrm2_    SNRM2
#define dsymv_    SSYMV
#define ddot_     SDOT
#define dsyr2_    SSYR2
#define drot_     SROT
#define dgemv_    SGEMV
#define dtrsv_    STRSV
#define dgemm_    SGEMM
#define dtrsm_    STRSM

#define cswap_    CSWAP
#define caxpy_    CAXPY
#define scasum_   SCASUM
#define icamax_   ICAMAX
#define ccopy_    CCOPY
#define cscal_    CSCAL
#define scnrm2_   SCNRM2
#define caxpy_    CAXPY
#define cgemv_    CGEMV
#define ctrsv_    CTRSV
#define cgemm_    CGEMM
#define ctrsm_    CTRSM
#define cgerc_    CGERC
#define chemv_    CHEMV
#define cher2_    CHER2

#define zswap_    ZSWAP
#define zaxpy_    ZAXPY
#define dzasum_   DZASUM
#define izamax_   IZAMAX
#define zcopy_    ZCOPY
#define zscal_    ZSCAL
#define dznrm2_   DZNRM2
#define zgemv_    ZGEMV
#define ztrsv_    ZTRSV
#define zgemm_    ZGEMM
#define ztrsm_    ZTRSM
#define zgerc_    ZGERC
#define zhemv_    ZHEMV
#define zher2_    ZHER2

/* LAPACK */
#define dlacon_   DLACON
#define slacon_   SLACON
#define icmax1_   ICMAX1
#define scsum1_   SCSUM1
#define clacon_   CLACON
#define dzsum1_   DZSUM1
#define izmax1_   IZMAX1
#define zlacon_   ZLACON

/* Fortran interface */
#define c_bridge_dgssv_ C_BRIDGE_DGSSV
#define c_fortran_sgssv_ C_FORTRAN_SGSSV
#define c_fortran_dgssv_ C_FORTRAN_DGSSV
#define c_fortran_cgssv_ C_FORTRAN_CGSSV
#define c_fortran_zgssv_ C_FORTRAN_ZGSSV
#endif


#if (F77_CALL_C == NOCHANGE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgemm(...)           void dgemm(...)
 */
/* BLAS */
#define sswap_    sswap
#define saxpy_    saxpy
#define sasum_    sasum
#define isamax_   isamax
#define scopy_    scopy
#define sscal_    sscal
#define sger_     sger
#define snrm2_    snrm2
#define ssymv_    ssymv
#define sdot_     sdot
#define saxpy_    saxpy
#define ssyr2_    ssyr2
#define srot_     srot
#define sgemv_    sgemv
#define strsv_    strsv
#define sgemm_    sgemm
#define strsm_    strsm

#define dswap_    dswap
#define daxpy_    daxpy
#define dasum_    dasum
#define idamax_   idamax
#define dcopy_    dcopy
#define dscal_    dscal
#define dger_     dger
#define dnrm2_    dnrm2
#define dsymv_    dsymv
#define ddot_     ddot
#define dsyr2_    dsyr2
#define drot_     drot
#define dgemv_    dgemv
#define dtrsv_    dtrsv
#define dgemm_    dgemm
#define dtrsm_    dtrsm

#define cswap_    cswap
#define caxpy_    caxpy
#define scasum_   scasum
#define icamax_   icamax
#define ccopy_    ccopy
#define cscal_    cscal
#define scnrm2_   scnrm2
#define cgemv_    cgemv
#define ctrsv_    ctrsv
#define cgemm_    cgemm
#define ctrsm_    ctrsm
#define cgerc_    cgerc
#define chemv_    chemv
#define cher2_    cher2

#define zswap_    zswap
#define zaxpy_    zaxpy
#define dzasum_   dzasum
#define izamax_   izamax
#define zcopy_    zcopy
#define zscal_    zscal
#define dznrm2_   dznrm2
#define zgemv_    zgemv
#define ztrsv_    ztrsv
#define zgemm_    zgemm
#define ztrsm_    ztrsm
#define zgerc_    zgerc
#define zhemv_    zhemv
#define zher2_    zher2

/* LAPACK */
#define dlacon_   dlacon
#define slacon_   slacon
#define icmax1_   icmax1
#define scsum1_   scsum1
#define clacon_   clacon
#define dzsum1_   dzsum1
#define izmax1_   izmax1
#define zlacon_   zlacon

/* Fortran interface */
#define c_bridge_dgssv_ c_bridge_dgssv
#define c_fortran_sgssv_ c_fortran_sgssv
#define c_fortran_dgssv_ c_fortran_dgssv
#define c_fortran_cgssv_ c_fortran_cgssv
#define c_fortran_zgssv_ c_fortran_zgssv
#endif


#endif /* __SUPERLU_CNAMES */
