############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT =

#
#  The name of the libraries to be created/linked to
#
#TMGLIB       = tmglib$(PLAT).a
SUPERLULIB   = ../libs/libsuperlu$(PLAT).a
BLASLIB      = ../../libs/libmyblas2$(PLAT).a

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = cc
CFLAGS       = -fpic -O3
LOADER       = cc
LOADOPTS     = 

#
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
CDEFS        = -DAdd_ -DUSER_MALLOC=superlu_python_module_malloc -DUSER_FREE=superlu_python_module_free -DUSER_ABORT=superlu_python_module_abort
#
# The directory in which Matlab is installed
#
MATLAB	     = /usr/sww/matlab
