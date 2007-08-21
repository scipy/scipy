############################################################################
#
#  Program:         ARPACK
#
#  Module:          Makefile
#
#  Purpose:         Sources Makefile
#
#  Creation date:   February 22, 1996
#
#  Modified:        September 6, 1996
#
#  Send bug reports, comments or suggestions to arpack.caam.rice.edu
#
############################################################################
#\SCCS Information: @(#)
# FILE: Makefile   SID: 2.1   DATE OF SID: 9/9/96   RELEASE: 2
 
include ../ARmake.inc
 
############################################################################
#  To create or add to the library, enter make followed by one or
#  more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The name of the library is defined by $(ARPACKLIB) in
#  ../ARmake.inc and is created at the next higher directory level.
#

OBJS  = icnteq.o icopy.o iset.o iswap.o ivout.o second.o

SOBJ  = svout.o  smout.o

DOBJ  = dvout.o  dmout.o

COBJ  = cvout.o  cmout.o 

ZOBJ  = zvout.o  zmout.o

.SUFFIXES:      .o .F .f

.f.o:
	$(FC) $(FFLAGS) -c $<
#
#  make the library containing both single and double precision
#
all: single double complex complex16
 
single: $(SOBJ) $(OBJS)
	$(AR) $(ARFLAGS) $(ARPACKLIB) $(SOBJ) $(OBJS)
	$(RANLIB) $(ARPACKLIB)
 
double: $(DOBJ) $(OBJS) $(ZOBJ)
	$(AR) $(ARFLAGS) $(ARPACKLIB) $(DOBJ) $(OBJS)
	$(RANLIB) $(ARPACKLIB)

complex: $(SOBJ) $(OBJS) $(COBJ)
	$(AR) $(ARFLAGS) $(ARPACKLIB) $(SOBJ) $(COBJ) $(OBJS)
	$(RANLIB) $(ARPACKLIB)

complex16: $(DOBJ) $(OBJS) $(ZOBJ)
	$(AR) $(ARFLAGS) $(ARPACKLIB) $(DOBJ) $(ZOBJ) $(OBJS)
	$(RANLIB) $(ARPACKLIB)
#
sdrv:
	
ddrv:
	
cdrv:
	
zdrv:
	
#
#  clean	- remove all object files
#
clean:
	rm -f *.o a.out core

