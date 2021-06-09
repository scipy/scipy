dnl Process this file with autoconf to produce a configure script.
dnl $Id: configure.ac 5330 2011-04-19 10:50:18Z leydold $

dnl Initialize.
AC_PREREQ(2.63)
AC_INIT([unuran],[1.8.1],[unuran@statmath.wu.ac.at])
AC_CONFIG_AUX_DIR(autoconf)
AM_INIT_AUTOMAKE([unuran],[1.8.1])
AC_CONFIG_SRCDIR([src/unuran_config.h])
AC_CONFIG_HEADER([config.h])

dnl Library versioning (current:revision:age)
dnl See the libtool manual for an explanation of the numbers
dnl 
dnl unuran-0.4.0   0:0:0
dnl unuran-0.4.2   1:1:0
dnl unuran-0.4.4   1:2:0
dnl unuran-0.5.0   2:0:0
dnl unuran-0.6.0   3:0:0
dnl unuran-0.7.0   4:0:0
dnl unuran-0.7.1   4:1:0
dnl unuran-0.7.2   4:2:0
dnl unuran-0.8.1   5:0:0
dnl unuran-0.9.0   6:0:0
dnl unuran-1.0.0   7:0:0
dnl unuran-1.0.1   7:1:0
dnl unuran-1.1.0   8:0:1
dnl unuran-1.2.0   9:0:0
dnl unuran-1.2.1   9:1:0
dnl unuran-1.2.2   9:2:0
dnl unuran-1.2.3   9:3:2
dnl unuran-1.2.4   9:4:2
dnl unuran-1.3.0   10:0:0
dnl unuran-1.3.1   10:0:0
dnl unuran-1.4.0   11:0:1
dnl unuran-1.4.1   11:1:1
dnl unuran-1.5.0   12:0:0
dnl unuran-1.6.0   13:0:0
dnl unuran-1.7.0   14:0:0
dnl unuran-1.8.0   15:0:0
dnl unuran-1.8.1   15:0:0

UNURAN_LT_VERSION="15:0:0"
AC_SUBST([UNURAN_LT_VERSION])

dnl use maintainer mode
AM_MAINTAINER_MODE

dnl Check for which system.
AC_CANONICAL_HOST

dnl Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_CHECK_TOOL(AR, ar, :)
LT_INIT([disable-shared])

dnl Checks for libraries.
dnl libm: the standard math library
AC_CHECK_LIB([m], [sqrt])

dnl libgsl: the GNU Scientific library 
dnl (used for experimental code or 
dnl when the GSL uniform random number generators are used, see below)
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])


dnl Optional support for external uniform random number generators

dnl RngStream: Pierre L'Ecuyer's multiple random stream generator
AX_ADD_URNG_LIB([rngstream],[RNGSTREAM],[Pierre L'Ecuyer's RNGSTREAM library],
                [rngstreams],[RngStream_CreateStream],[RngStream.h])

dnl prng: Otmar Lendl's library for pseudo random number generators
AX_ADD_URNG_LIB([prng],[PRNG],[Otmar Lendl's PRNG library],
                [prng],[prng_new],[prng.h])

dnl GSL: GNU Scientific Library
AX_ADD_URNG_LIB([gsl],[GSL],[random number generators from GNU Scientific Library],
                [gsl],[gsl_rng_alloc],[gsl/gsl_rng.h])

dnl Default URNG
AC_ARG_WITH(urng_default,
	[AS_HELP_STRING([--with-urng-default],
	   [URNG: global default (builtin|rngstream) @<:@default=builtin@:>@])],
	[case "${withval}" in
	   builtin) ;; rngstream) ;; 
	   *) AC_MSG_ERROR(bad value '${withval}' for --with-urng-default) ;;
	esac],
	[with_urng_default=builtin])
AS_IF([test "x$with_urng_default" == xrngstream],
      AS_IF([test "x$with_urng_rngstream" == xno],
	     [with_urng_default=builtin
              AC_MSG_NOTICE([
********************************************************
*
*   --with-urng-default=rngstream
*   but --with-urng-rngstream=no.
*   use --with-urng-default=builtin instead
*
********************************************************
             ]) ], 
             [AC_DEFINE([UNUR_URNG_DEFAULT_RNGSTREAM], [1], 
                [Define to 1 if you use RNGSTREAM as global URNG.])]), 
      [])

dnl Optional support for Rmath standalone library
AX_ADD_RMATH()

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h sys/time.h unistd.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME

dnl Checks for library functions.
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strcasecmp strchr strtol strtoul])

dnl For the following function replacements exist
AH_BOTTOM([/* Substitute UNU.RAN functions for missing system functions */])

AC_CHECK_DECLS(log1p,,,[#include <math.h>])
AH_BOTTOM([#if !HAVE_DECL_LOG1P
#define log1p _unur_log1p
#endif])

dnl Checks for GNU extensions
AC_CHECK_DECLS([getopt, snprintf, vsnprintf])

dnl Check for C99 functions */
AC_CHECK_DECLS(isfinite,,,[#include <math.h>])
AC_CHECK_DECLS(isinf,,,[#include <math.h>])
AC_CHECK_DECLS(isnan,,,[#include <math.h>])
AC_CHECK_DECLS(INFINITY,,,[#include <math.h>])
AC_CHECK_DECLS(HUGE_VAL,,,[#include <math.h>])

dnl Check for alarm and signal handler
AC_CHECK_DECLS(alarm,,,[#include <unistd.h>])
AC_CHECK_DECLS(signal,,,[#include <signal.h>])
AC_CHECK_FUNCS([alarm signal])

dnl Check IEEE comparisons, whether "x != x" is true for NaNs
AX_IEEE_COMP

dnl Some routines are deprecated.
dnl Compiling and linking these into the library can be controlled
dnl by a configure flag. 
dnl Default is disabled.
AC_ARG_ENABLE(deprecated,
	[AS_HELP_STRING([--enable-deprecated],
		[enable support for deprecated UNU.RAN routines @<:@default=no@:>@])],
	[],
	[enable_deprecated=no])
AM_CONDITIONAL(ENABLE_DEPRECATED, test x$enable_deprecated = xyes)
AS_IF([test "x$enable_deprecated" = xyes],
	[AC_DEFINE([USE_DEPRECATED_CODE], [1], 
	           [Define to 1 if you want to use deprecated code.])], 
	[])

dnl UNU.RAN has additional checks against invalid NULL pointers or 
dnl pointers to wrong structures. For the latter cookies are inserted
dnl into structures and checked whenever such structure is accessed.
dnl Both checks can be enabled by the '--enable-check-struct' flag.
dnl Default is disabled.
AC_ARG_ENABLE(check-struct,
	[AS_HELP_STRING([--enable-check-struct],
		[Debug: check validity of pointers to structures @<:@default=no@:>@])],
	[],
	[enable_check_struct=no])
AS_IF([test "x$enable_check_struct" = xyes],
	[AC_DEFINE([UNUR_ENABLE_CHECKNULL], [1], 
	           [Define to 1 if you want to perform additional checks against an invalid NULL pointer])
	 AC_DEFINE([UNUR_COOKIES], [1], 
	           [Define to 1 if you want to use magic cookies to validate type of a pointer])], 
	[])


dnl UNU.RAN allows to write information about distribution and a log file.
dnl The default name of the logfile is set in 'src/unuran_config.h'.
dnl The logfile is enabled by the '--enable-logging' flag.
AC_ARG_ENABLE(logging,
	[AS_HELP_STRING([--enable-logging],
		[Debug: print informations about generator into logfile @<:@default=no@:>@])],
	[],
	[enable_logging=no])
AS_IF([test "x$enable_logging" = xyes],
	[AC_DEFINE([UNUR_ENABLE_LOGGING], [1], 
	           [Define to 1 if you want to use a logfile for logging informations about UNU.RAN objects])], 
	[])

dnl UNU.RAN provides a function that returns a string with some information
dnl about a given generator object.
dnl This function is enabled by the '--enable-info' flag.
AC_ARG_ENABLE(info,
	[AS_HELP_STRING([--enable-info],
		[Info: provide function with information about generator objects @<:@default=yes@:>@])],
	[],
	[enable_info=yes])
AS_IF([test "x$enable_info" = xyes],
	[AC_DEFINE([UNUR_ENABLE_INFO], [1], 
	           [Define to 1 if you want to use the info routine for printing informations about UNU.RAN objects])], 
	[])


dnl UNU.RAN contains some experimental code that can be enabled
dnl by means of the '--enable-experimental' flag.
AC_ARG_ENABLE(experimental,
	[AS_HELP_STRING([--enable-experimental],
		[Experimental code (not documented) @<:@default=no@:>@])],
	[],
	[enable_experimental=no])
AS_IF([test "x$enable_experimental" = xyes],
	[AC_DEFINE([USE_EXPERIMENTAL_CODE], [1], 
	           [Define to 1 if you want to use experimental code])], 
	[])
AM_CONDITIONAL(ENABLE_EXPERIMENTAL, test x$enable_experimental = xyes)


dnl Some tests compare the return values of PDFs, CDFs and of
dnl the function parser with the output of corresponding
dnl Mathematica(R) scripts.
dnl Thus for these tests we need Mathematica(R).
dnl However, simply checking for the existance of programm 'math'
dnl breaks the configure script when 'math' is available but
dnl ask for a password. 
dnl Hence these tests are enabled by means of a configure option.
AC_ARG_WITH(mathematica,
	[AS_HELP_STRING([--with-mathematica],
		[Use Mathematica(R) for some checks @<:@default=no@:>@])],
	[case "${withval}" in
	  yes) mathematica=true ;;
	  no)  mathematica=false ;;
	  *) AC_MSG_ERROR(bad value '${withval}' for --with-mathematica) ;;
	esac],
	[mathematica=false])
AM_CONDITIONAL(HAVE_MATHEMATICA, test x$mathematica = xtrue)

dnl check for Valgrind
dnl (Valgrind is a GPL'd system for debugging and profiling x86-Linux programs)
AC_CHECK_PROG(have_valgrind,valgrind,true,false) 
AM_CONDITIONAL(HAVE_VALGRIND, test x$have_valgrind = xtrue)

dnl Set flags for compiler

dnl GCC
EXTRA_GCC_FLAGS=
EXTRA_GCC_FLAGS="$EXTRA_GCC_FLAGS -Wshadow"
##EXTRA_GCC_FLAGS="$EXTRA_GCC_FLAGS -Wbad-function-cast -Wcast-align -Wmissing-declarations -Wmissing-format-attribute -Wmissing-prototypes -Wnested-externs -Wpointer-arith -Wredundant-decls -Wstrict-prototypes -Wwrite-strings -Wconversion"

## the following flags raise warnings:
##EXTRA_GCC_FLAGS="$EXTRA_GCC_FLAGS -Wcast-qual -Wfloat-equal -Wunreachable-code"
## -Wcast-qual
##     Warn whenever a pointer is cast so as to remove a type qualifier from the target type.  For example, warn if a
##     "const char *" is cast to an ordinary "char *".
## -Wfloat-equal
##     Warn if floating point values are used in equality comparisons.
## -Wunreachable-code
##     Warn if the compiler detects that code will never be executed.

CC_FLAGS=
if test X"$GCC" = Xyes ; then
   CC_FLAGS="$CC_FLAGS -Wall -Wextra $EXTRA_GCC_FLAGS"
fi


case "${host_os}" in
linux*)
	case  "${CC}" in
	## Intel compiler
	*icc)
		## -mp is the preferred form of -mieee-fp
		CC_FLAGS="$CC_FLAGS -mp"
		;;
	esac
	;;
esac

AC_SUBST(AM_CFLAGS,[$CC_FLAGS])


dnl create Makefiles for source
AC_CONFIG_FILES([\
	Makefile \
	src/Makefile \
	src/distr/Makefile \
	src/distributions/Makefile \
	src/methods/Makefile \
	src/parser/Makefile \
	src/specfunct/Makefile \
	src/tests/Makefile \
	src/uniform/Makefile \
	src/urng/Makefile \
	src/utils/Makefile \
	scripts/Makefile \
	doc/Makefile \
	doc/src/Makefile \
	doc/figures/Makefile \
	tests/Makefile \
	tests/testdistributions/Makefile \
	examples/Makefile \
	experiments/Makefile \
])

AC_OUTPUT

dnl -----------------------------------------

AC_MSG_RESULT([
*========================================================
*
* Configuration for ${PACKAGE}-${VERSION}:
*
*  Source directory:           ${srcdir}
*  Installation directory:     ${prefix}
*  C Compiler:                 ${CC}
*
*  Build shared library:       ${enable_shared}
*  Build static library:       ${enable_static}
*
*  Use deprecated code:        ${enable_deprecated}
*  Enable runtime checks:      ${enable_check_struct}
*  Enable logging of data:     ${enable_logging}
*  Enable info routine:        ${enable_info}
*  Experimenal code:           ${enable_experimental}
*  Use Rmath library:          ${with_Rmath}
*
*  Support for external uniform random number generators
*    Rngstream:                ${with_urng_rngstream}
*    PRNG:                     ${with_urng_prng}
*    GSL:                      ${with_urng_gsl}
*
*  Default URNG:               ${with_urng_default}
*
*========================================================
])

dnl -----------------------------------------




