LIB=amos
ROOTDIR=/usr/local
BUILDIR=$(ROOTDIR)
FFLAGS=-O2
FFC=g77
FC=$(FFC) 
OBJ = \
	d1mach.o	\
	dgamln.o	\
	fdump.o	\
	i1mach.o	\
	xerror.o	\
	zabs.o	\
	zacai.o	\
	zacon.o	\
	zairy.o	\
	zasyi.o	\
	zbesh.o	\
	zbesi.o	\
	zbesj.o	\
	zbesk.o	\
	zbesy.o	\
	zbinu.o	\
	zbiry.o	\
	zbknu.o	\
	zbuni.o	\
	zbunk.o	\
	zdiv.o	\
	zexp.o	\
	zkscl.o	\
	zlog.o	\
	zmlri.o	\
	zmlt.o	\
	zrati.o	\
	zs1s2.o	\
	zseri.o	\
	zshch.o	\
	zsqrt.o	\
	zuchk.o	\
	zunhj.o	\
	zuni1.o	\
	zuni2.o	\
	zunik.o	\
	zunk1.o	\
	zunk2.o	\
	zuoik.o	\
	zwrsk.o
lib$(LIB).a:	$(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a

install:	lib$(LIB).a
	install -c lib$(LIB).a $(BUILDIR)/lib 

clean:
	rm -fr *.o *.a *.so *.lo *.la
