include ../Makefile


PROJNAME = fslio

INCFLAGS = $(ZLIB_INC) $(ZNZ_INC) $(NIFTI_INC)
LIBS = $(ZLIB_LIBS) $(ZNZ_LIBS) $(NIFTI_LIBS) $(FSL_LIBS)

SRCS=fslio.c
OBJS=fslio.o

SCRIPTS = remove_ext fsloutputtype imtest imglob imcp imln imrm immv



depend:	
	$(RM) -f depend.mk
	$(MAKE) depend.mk

depend.mk:
	$(CC) $(DEPENDFLAGS) $(INCFLAGS) $(SRCS) >> depend.mk

lib: libfslio.a

libfslio.a: ${OBJS}
	${AR} -r libfslio.a ${OBJS}


include depend.mk
