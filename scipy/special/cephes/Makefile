LIB = cephes
ROOTDIR=/usr/local
BUILDIR=$(ROOTDIR)
CCC=gcc
CC = $(CCC) 
OPT=-O2 -march=pentiumpro
CFLAGS=$(OPT)

OBJS = airy.o bdtr.o beta.o \
btdtr.o cbrt.o chbevl.o chdtr.o const.o \
dawsn.o ellie.o ellik.o ellpe.o ellpj.o ellpk.o \
exp10.o exp2.o expn.o fabs.o fdtr.o \
fresnl.o gamma.o gdtr.o hyp2f1.o hyperg.o i0.o i1.o igami.o \
incbet.o incbi.o igam.o iv.o j0.o j1.o jn.o jv.o k0.o k1.o \
kn.o kolmogorov.o nbdtr.o ndtr.o ndtri.o pdtr.o \
polevl.o polmisc.o polyn.o psi.o rgamma.o round.o \
shichi.o sici.o sindg.o spence.o stdtr.o struve.o \
tandg.o unity.o yn.o zeta.o zetac.o \
setprec.o mtherr.o


lib$(LIB).a:	$(OBJS)
	ar ru lib$(LIB).a $?
	ranlib lib$(LIB).a

install:	lib$(LIB).a
	install -c lib$(LIB).a $(BUILDIR)/lib

clean:
	rm -f *.o *.so *.lo

