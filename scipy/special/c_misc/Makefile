LIB=c_misc
ROOTDIR=/usr/local
BUILDIR=$(ROOTDIR)
CFLAGS=-O2
CCC=gcc
CC=$(CCC)
OBJ=besselpoly.o
lib$(LIB).a:    $(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a

install:	lib$(LIB).a
	install -c lib$(LIB).a $(BUILDIR)/lib

clean:
	rm -f *.o *.a *.so

