############################################################################
#
#  Program:         SuperLU
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   October 2, 1995
#
#  Modified:        February 4, 1997  Version 1.0
#
############################################################################

include make.inc

all: install lib testing

lib: superlulib tmglib

clean: cleanlib cleantesting

install:
	( cd INSTALL; $(MAKE) )
#	( cd INSTALL; cp lsame.c ../SRC/; \
#	  cp dlamch.c ../SRC/; cp slamch.c ../SRC/ )

blaslib:
	( cd CBLAS; $(MAKE) )

superlulib:
	( cd SRC; $(MAKE) )

tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

matlabmex:
	( cd MATLAB; $(MAKE) )

testing:
	( cd TESTING ; $(MAKE) )

cleanlib:
	( cd SRC; $(MAKE) clean )
	( cd TESTING/MATGEN; $(MAKE) clean )
	( cd CBLAS; $(MAKE) clean )

cleantesting:
	( cd INSTALL; $(MAKE) clean )
	( cd TESTING; $(MAKE) clean )
	( cd MATLAB;  $(MAKE) clean )
	( cd EXAMPLE; $(MAKE) clean )
