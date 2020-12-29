#
#     (C) Rasmus Munk Larsen, Stanford University, 2004
#

include ../../make.inc
OBJ = dbdsdc.o  dlamrg.o  dlascl.o  dlasd4.o  dlasda.o  dlasq2.o  dlasr.o  ilaenv.o \
dbdsqr.o  dlanst.o  dlasd0.o  dlasd5.o  dlasdq.o  dlasq3.o  dlasrt.o  lsame.o \
dlacpy.o  dlapy2.o  dlasd1.o  dlasd6.o  dlasdt.o  dlasq4.o  dlassq.o \
dlaed6.o  dlartg.o  dlasd2.o  dlasd7.o  dlaset.o  dlasq5.o  dlasv2.o \
dlas2.o   dlasd3.o  dlasd8.o  dlasq1.o  dlasq6.o  dlamch.o xerbla.o \
ieeeck.o  dlartg.o  dlarnv.o  dlaruv.o  zlarnv.o  zlascl.o

.F.o:
	$(F77) $(FFLAGS) $(IPATH) -c -o $*.o $*.F

.f.o:
	$(F77) $(FFLAGS) $(IPATH) -c -o $*.o $*.f

dlamch.o: dlamch.f
	$(F77)	$(FFLAGS_NOOPT) $(IPATH) -c -o $*.o $*.f

ieeeck.o: ieeeck.f
	$(F77)	$(FFLAGS_NOOPT) $(IPATH) -c -o $*.o $*.f

ilaenv.o: ilaenv.f
	$(F77)	$(FFLAGS_NOOPT) $(IPATH) -c -o $*.o $*.f

dlasd4.o: dlasd4.f
	$(F77)	$(FFLAGS_ACCURATE) $(IPATH) -c -o $*.o $*.f


lib: ../libzlapack_util_$(PLAT).a

../libzlapack_util_$(PLAT).a: $(OBJ)
	$(MAKELIB) libzlapack_util_$(PLAT).a $(OBJ)	
	$(RANLIB)  libzlapack_util_$(PLAT).a
	mv libzlapack_util_$(PLAT).a ..

clean:	
	rm -f  *.o *.il
	rm -rf rii_files	

cleanall:
	rm -f  *.o *.x *.a *~ *.il ifc*
	rm -rf rii_files	

