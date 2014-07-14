LIB=odepack
FFLAGS=-O
OBJ=\
adjlr.o\
aigbt.o\
ainvg.o\
blkdta000.o\
bnorm.o\
cdrv.o\
cfode.o\
cntnzu.o\
decbt.o\
ewset.o\
fnorm.o\
intdy.o\
iprep.o\
jgroup.o\
lsoda.o\
lsodar.o\
lsode.o\
lsodes.o\
lsodi.o\
lsoibt.o\
md.o\
mdi.o\
mdm.o\
mdp.o\
mdu.o\
nnfc.o\
nnsc.o\
nntc.o\
nroc.o\
nsfc.o\
odrv.o\
pjibt.o\
prep.o\
prepj.o\
prepji.o\
prja.o\
prjs.o\
rchek.o\
roots.o\
slsbt.o\
slss.o\
solbt.o\
solsy.o\
srcar.o\
srcma.o\
srcms.o\
srcom.o\
sro.o\
stoda.o\
stode.o\
stodi.o\
vmnorm.o\
vnorm.o\
xerrwv.o\
xsetf.o\
xsetun.o

lib$(LIB).a:	$(OBJ)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a

install:	lib$(LIB).a
	mv lib$(LIB).a /usr/local/lib
	rm *.o

test: test.o
	f77 test.o -l$(LIB)
	time a.out
