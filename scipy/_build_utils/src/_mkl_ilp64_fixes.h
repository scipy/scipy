/*
 * Fix missing ILP64 symbols in MKL (symbols without trailing underscore).
 * See https://community.intel.com/t5/Intel-oneAPI-Math-Kernel-Library/Missing-cspr-64-symbol-in-ILP64-MKL/td-p/1703471
 *
 * MKL ILP64 provides these symbols as `foo_64` instead of `foo_64_`
 * (i.e. without the Fortran trailing underscore after the ILP64 suffix).
 *
 * Include this header wherever BLAS/LAPACK ILP64 symbols are referenced.
 */
#ifdef FIX_MKL_2025_ILP64_MISSING_SYMBOL
#define cgetc2_64_ cgetc2_64
#define cgtts2_64_ cgtts2_64
#define clarfx_64_ clarfx_64
#define clartv_64_ clartv_64
#define cptts2_64_ cptts2_64
#define cspr_64_ cspr_64
#define dgetc2_64_ dgetc2_64
#define dgtts2_64_ dgtts2_64
#define dlarfx_64_ dlarfx_64
#define dlartv_64_ dlartv_64
#define dlasd3_64_ dlasd3_64
#define dlasd4_64_ dlasd4_64
#define dptts2_64_ dptts2_64
#define sgetc2_64_ sgetc2_64
#define sgtts2_64_ sgtts2_64
#define slartv_64_ slartv_64
#define slasd3_64_ slasd3_64
#define slasd4_64_ slasd4_64
#define sptts2_64_ sptts2_64
#define zgetc2_64_ zgetc2_64
#define zgtts2_64_ zgtts2_64
#define zlartv_64_ zlartv_64
#define zptts2_64_ zptts2_64
#endif
