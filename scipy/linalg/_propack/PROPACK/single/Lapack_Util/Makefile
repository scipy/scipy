#
#     (C) Rasmus Munk Larsen, Stanford University, 2004
#

include ../../make.inc
OBJ = sbdsdc.o  slamrg.o  slascl.o  slasd4.o  slasda.o  slasq2.o  slasr.o  \
ilaenv.o sbdsqr.o  slanst.o  slasd0.o  slasd5.o  slasdq.o  slasq3.o  slasrt.o \
lsame.o  slacpy.o  slapy2.o  slasd1.o  slasd6.o  slasdt.o  slasq4.o  slassq.o \
slaed6.o slartg.o  slasd2.o  slasd7.o  slaset.o  slasq5.o  slasv2.o \
slas2.o  slasd3.o  slasd8.o  slasq1.o  slasq6.o  \
ieeeck.o slartg.o  slarnv.o  slaruv.o  slamch.o xerbla.o 

.F.o:
	$(F77) $(FFLAGS) $(IPATH) -c -o $*.o $*.F

.f.o:
	$(F77) $(FFLAGS) $(IPATH) -c -o $*.o $*.f

slamch.o: slamch.f
	$(F77)	$(FFLAGS_NOOPT) $(IPATH) -c -o $*.o $*.f

ieeeck.o: ieeeck.f
	$(F77)	$(FFLAGS_NOOPT) $(IPATH) -c -o $*.o $*.f

ilaenv.o: ilaenv.f
	$(F77)	$(FFLAGS_NOOPT) $(IPATH) -c -o $*.o $*.f

slasd4.o: slasd4.f
	$(F77)	$(FFLAGS_ACCURATE) $(IPATH) -c -o $*.o $*.f


lib: ../libslapack_util_$(PLAT).a

../libslapack_util_$(PLAT).a: $(OBJ)
	$(MAKELIB) libslapack_util_$(PLAT).a $(OBJ)	
	$(RANLIB)  libslapack_util_$(PLAT).a
	mv libslapack_util_$(PLAT).a ..

clean:
	rm -f  *.o *.il 
	rm -rf rii_files	

cleanall:
	rm -f  *.o *.x *.a *~ *.il ifc*
	rm -rf rii_files	

