LIB=sparsekit
PREFIX=/usr/local
FFLAGS=-O2 -funroll-all-loops
FC=g77
OBJ = \
	sformats.o \
	sunary.o \
	sblassm.o \
	smatvec.o \
	dformats.o \
	dunary.o \
	dblassm.o \
	dmatvec.o \
	cformats.o \
	cunary.o \
	cblassm.o \
	cmatvec.o \
	zformats.o \
	zunary.o \
	zblassm.o \
	zmatvec.o
lib$(LIB).a:	$(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a

install:	lib$(LIB).a
	cp lib$(LIB).a ../libs/
	rm *.o

clean:
	rm -f *.o *.a *.so
