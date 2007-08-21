LIB=toms
ROOTDIR=/usr/local
BUILDIR=$(ROOTDIR)
FFC=g77
FC=$(FFC) 
FFLAGS=-O
OBJ=wofz.o

lib$(LIB).a:    $(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a 

install:	lib$(LIB).a
	libtool install -c lib$(LIB).a $(BUILDIR)/lib

clean:
	rm -f *.o *.a *.so

