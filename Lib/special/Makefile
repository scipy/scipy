CCC = gcc
FFC = g77
ROOTDIR=/usr
BUILDIR=$(ROOTDIR)
INCLUDES = -I/usr/local/include/python2.1
OPT=-O2 -march=pentium
CFLAGS = $(OPT) $(INCLUDES) 

FFLAGS=$(OPT)
INSTALLDIR = $(BUILDIR)/lib/python2.1/site-packages

OBJS = ufunc_extras.o amos_wrappers.o toms_wrappers.o cephesmodule.o 
LIBS = amos/libamos.a toms/libtoms.a c_misc/libmiscfuncs.a cephes/libcephes.a


all: 	$(OBJS) $(LIBS)
	$(FFC) -shared $(FFLAGS) -o cephesmodule.so $(OBJS) \
	-L./misc -lmiscfuncs -L./cephes -lcephes \
	-L./amos -lamos -L./toms -ltoms -lm 

cephesmodule.o: cephesmodule.c cephes_doc.h
	$(CCC) -c $(CFLAGS) -fPIC cephesmodule.c

cephes_doc.h: docs/included_functions.html
	python gendoc.py 

amos/libamos.a:
	cd amos && make CCC=$(CCC) CFLAGS="$(OPT)" FFC=$(FFC) \
	FFLAGS="$(FFLAGS)" ROOTDIR="$(ROOTDIR)" BUILDIR="$(BUILDIR)"

cephes/libcephes.a:
	cd cephes && make CCC=$(CCC) CFLAGS="$(OPT)" FFC=$(FFC) \
	FFLAGS="$(FFLAGS)" ROOTDIR="$(ROOTDIR)" BUILDIR="$(BUILDIR)"

toms/libtoms.a:
	cd toms && make CCC=$(CCC) CFLAGS="$(OPT)" FFC=$(FFC) \
	FFLAGS="$(FFLAGS)" ROOTDIR="$(ROOTDIR)" BUILDIR="$(BUILDIR)"

c_misc/libmiscfuncs.a:
	cd misc && make CCC=$(CCC) CFLAGS="$(OPT)" FFC=$(FFC) \
	FFLAGS="$(FFLAGS)" ROOTDIR="$(ROOTDIR)" BUILDIR="$(BUILDIR)"

install:
	install -d $(INSTALLDIR)/SpecialFuncs
	install -c cephesmodule.so $(INSTALLDIR)/SpecialFuncs
	cp SpecialFuncs.py $(INSTALLDIR)/SpecialFuncs
	echo "SpecialFuncs" > $(INSTALLDIR)/SpecialFuncs.pth
	cd misc && make install BUILDIR="$(BUILDIR)"
	cd toms && make install BUILDIR="$(BUILDIR)"
	cd cephes && make install BUILDIR="$(BUILDIR)"
	cd amos && make install BUILDIR="$(BUILDIR)"

testref:
	cd tests && make testref

test:
	cd tests && make test

clean:
	rm -f *.o *.so

distclean:
	rm -f *.o *.so
	rm -fr .libs
	rm -f cephes/*.o cephes/*.a cephes/*.l?
	rm -fr cephes/.libs
	rm -f amos/*.o amos/*.a amos/*.l?
	rm -fr amos/.libs
	rm -f toms/*.o toms/*.a toms/*.l?
	rm -fr toms/.libs
	rm -f misc/*.o misc/*.a misc/*.l? 
	rm -fr misc/.libs
	find . "(" -name '*~' -o -name '*.pyc' -o -name '*.pyo' \
        -o -name 'ref_*.pkl' ")" -exec rm {} \;

