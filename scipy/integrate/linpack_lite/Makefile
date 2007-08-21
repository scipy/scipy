LIB=linpack_lite
FC=f77
FFLAGS=-O
OBJ = \
        dgtsl.o \
        dgbfa.o \
        dgbsl.o  \
        dgefa.o   \
        dgesl.o  \

lib$(LIB).a:	$(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a

clean:
	rm *.o *.a *.so *.pyc
