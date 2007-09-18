include ../Makefile

PROJNAME = znzlib

INCFLAGS = $(ZLIB_INC)

SRCS=znzlib.c
OBJS=znzlib.o

TESTXFILES = testprog

depend:	
	$(RM) -f depend.mk
	$(MAKE) depend.mk

depend.mk:
	$(CC) $(DEPENDFLAGS) $(INCFLAGS) $(SRCS) >> depend.mk

lib: libznz.a

test: $(TESTXFILES)

znzlib.o: znzlib.c znzlib.h
	$(CC) -c $(CFLAGS) $(USEZLIB) $(INCFLAGS) $<

libznz.a: $(OBJS)
	$(AR) -r libznz.a $(OBJS)
	$(RANLIB) $@

testprog: libznz.a testprog.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o testprog testprog.c $(ZLIB_LIBS)

include depend.mk
