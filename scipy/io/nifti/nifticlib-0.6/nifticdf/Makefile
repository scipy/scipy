include ../Makefile

PROJNAME = nifticdflib

INCFLAGS = $(NIFTI_INC) $(NIFTICDF_INC)

SRCS=nifticdf.c
OBJS=nifticdf.o

depend:	
	$(RM) -f depend.mk
	$(MAKE) depend.mk

depend.mk:
	$(CC) $(DEPENDFLAGS) $(INCFLAGS) $(SRCS) >> depend.mk

lib: libnifticdf.a

test: $(TESTXFILES)

nifticdf.o: nifticdf.c nifticdf.h
	$(CC) -c $(CFLAGS) $(USEZLIB) $(INCFLAGS) $<

libnifticdf.a: $(OBJS)
	$(AR) -r libnifticdf.a $(OBJS)
	$(RANLIB) $@

include depend.mk
