CC = gcc
FC = g77
CFLAGS = -O2
FFLAGS = -O -shared -c
BLASLOC = ../multipack/lib 

all:	libs _sparsekitmodule.so _zsuperlumodule.so _dsuperlumodule.so _ssuperlumodule.so _csuperlumodule.so

.f.o:
	$(FC) -c $<
.c.o:
	$(CC) -c -I/usr/include/python1.5 -I./SuperLU/SRC -fpic $(CFLAGS) $<
_sparsekitmodule.so:	_sparsekitmodule.o
	$(FC) -o  _sparsekitmodule.so -shared _sparsekitmodule.o -L./libs -lsparsekit -L$(BLASLOC) -lblas

_zsuperlumodule.so:	_zsuperlumodule.o _superlu_utils.o
	$(CC) -o _zsuperlumodule.so -shared _zsuperlumodule.o _superlu_utils.o -L./libs -lsuperlu -lmyblas2 -L$(BLASLOC) -lblas

_dsuperlumodule.so:	_dsuperlumodule.o 
	$(CC) -o _dsuperlumodule.so -shared _dsuperlumodule.o  -L./libs -lsuperlu -lmyblas2  -L$(BLASLOC) -lblas

_ssuperlumodule.so:	_ssuperlumodule.o
	$(CC) -o _ssuperlumodule.so -shared _ssuperlumodule.o -L./libs -lsuperlu -lmyblas2 -L$(BLASLOC) -lblas

_csuperlumodule.so:	_csuperlumodule.o
	$(CC) -o _csuperlumodule.so -shared _csuperlumodule.o -L./libs -lsuperlu -lmyblas2 -L$(BLASLOC) -lblas

libs:
	mkdir -p libs
	(cd sparsekit; make; make install)
	(cd SuperLU/SRC; make)
	(cd SuperLU/CBLAS; make)

