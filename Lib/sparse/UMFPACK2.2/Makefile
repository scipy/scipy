
UMFD	= umd2fb.o umd2co.o umd2fa.o umd2f0.o umd2f1.o umd2f2.o \
	  umd2fg.o umd2in.o umd2of.o umd2s2.o umd2sl.o umd2so.o umd2su.o \
	  umd2er.o umd2p1.o umd2p2.o umd2lt.o umd2ut.o \
	  umd2rf.o umd2ra.o umd2r0.o umd2r2.o umd2rg.o umd21i.o

UMFS	= ums2fb.o ums2co.o ums2fa.o ums2f0.o ums2f1.o ums2f2.o \
	  ums2fg.o ums2in.o ums2of.o ums2s2.o ums2sl.o ums2so.o ums2su.o \
	  ums2er.o ums2p1.o ums2p2.o ums2lt.o ums2ut.o \
	  ums2rf.o ums2ra.o ums2r0.o ums2r2.o ums2rg.o ums21i.o

UMFC	= umc2fb.o umc2co.o umc2fa.o umc2f0.o umc2f1.o umc2f2.o \
	  umc2fg.o umc2in.o umc2of.o umc2s2.o umc2sl.o umc2so.o umc2su.o \
	  umc2er.o umc2p1.o umc2p2.o \
	  umc2rf.o umc2ra.o umc2r0.o umc2r2.o umc2rg.o umc21i.o

UMFZ	= umz2fb.o umz2co.o umz2fa.o umz2f0.o umz2f1.o umz2f2.o \
	  umz2fg.o umz2in.o umz2of.o umz2s2.o umz2sl.o umz2so.o umz2su.o \
	  umz2er.o umz2p1.o umz2p2.o \
	  umz2rf.o umz2ra.o umz2r0.o umz2r2.o umz2rg.o umz21i.o

HARWELL	= mc21b.o mc13e.o

BLASLIB	= libblas.a
UMFPACKLIB = libumfpack.a

FFLAGS	= -O
F77	= f77

ALL = $(UMFPACKLIB) dmain smain cmain zmain dmain.out smain.out cmain.out zmain.out

all:  $(ALL)

libumfpack.a: $(UMFD) $(UMFS) $(UMFC) $(UMFZ) $(HARWELL)
	ar -rv libumfpack.a $(UMFD) $(UMFS) $(UMFC) $(UMFZ) $(HARWELL)

dmain.out: dmain in
	dmain < in > dmain.out

smain.out: smain in
	smain < in > smain.out

cmain.out: cmain cin
	cmain < cin > cmain.out

zmain.out: zmain cin
	zmain < cin > zmain.out

clean:
	- rm *.o

purge: clean
	- rm $(ALL)

dmain: dmain.o $(UMFPACKLIB) $(BLASLIB)
	$(F77) $(FFLAGS) -o dmain dmain.o $(UMFPACKLIB) $(BLASLIB)

smain: smain.o $(UMFPACKLIB) $(BLASLIB)
	$(F77) $(FFLAGS) -o smain smain.o $(UMFPACKLIB) $(BLASLIB)

cmain: cmain.o $(UMFPACKLIB) $(BLASLIB)
	$(F77) $(FFLAGS) -o cmain cmain.o $(UMFPACKLIB) $(BLASLIB) 

zmain: zmain.o $(UMFPACKLIB) $(BLASLIB)
	$(F77) $(FFLAGS) -o zmain zmain.o $(UMFPACKLIB) $(BLASLIB) 

.SUFFIXES:
.SUFFIXES: .o .f

.f.o:
	$(F77) $(FFLAGS) -c $<
