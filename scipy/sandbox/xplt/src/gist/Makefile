SHELL=/bin/sh
MAKE=make
include ../Make.cfg

CFLAGS=$(COPTIONS) -I. -I../play
X11LIB=$(XLIB) -lX11
CCLOAD=$(CC) $(LDOPTIONS) $(LDFLAGS)
SYS_LIBS=$(X11LIB) $(MATHLIB) $(FPELIB)

D_GISTPATH='-DGISTPATH="~/gist:~/Gist:$(Y_SITE)/g"'

OBJS=gist.o tick.o tick60.o engine.o gtext.o draw.o draw0.o clip.o \
  gread.o gcntr.o hlevel.o ps.o cgm.o xfancy.o xbasic.o
BOBJS=browser.o cgmin.o eps.o

all: gist

# always check to see whether ../play/libplay.a needs to be rebuilt
libgist:
	@cd ../play; $(MAKE)
	@$(MAKE) libgist.a

# libyor puts the object files in ../yorick/libyor.a instead of libgist.a
libyor:
	@cd ../play; $(MAKE)
	@$(MAKE) make.done

# libyor puts the object files in ../libpyg.a instead of libgist.a
libpyg:
	@cd ../play; $(MAKE)
	@$(MAKE) LIBYOR_A=../libpyg.a make.done

# make.done clobbers libyor.a if its rule actually runs
LIBYOR_A=../yorick/libyor.a
make.done: $(OBJS) ../play/libplay.a
	cp ../play/libplay.a $(LIBYOR_A)
	$(AR) r $(LIBYOR_A) $(OBJS)
	$(RANLIB) $(LIBYOR_A)
	@touch $@

libgist.a: $(OBJS) ../play/libplay.a
	cp ../play/libplay.a $@
	$(AR) r $@ $(OBJS)
	$(RANLIB) $@

gist: libgist $(BOBJS) main.o
	$(CCLOAD) -o $@ main.o $(BOBJS) -L. -lgist $(SYS_LIBS)

bench: libgist bench.o main.o
	$(CCLOAD) -o $@ main.o bench.o -L. -lgist $(SYS_LIBS)

main.c: ../play/unix/pmain.c
	cp ../play/unix/pmain.c main.c

clean::
	rm -f *~ '#'* *.o *.a *.done core a.out main.c gist bench

distclean:: clean
	rm -f .gdbinit

config:
	@:

# ------------- dependencies -------------

# cgm.h: gist.h engine.h
# draw.h: gist.h
# engine.h: gist.h
# gtext.h: gist.h
# hlevel.h: gist.h
# ps.h: gist.h engine.h
# xbasic.h: gist.h engine.h play.h
# xfancy.h: xbasic.h

PLAY_H=../play/play.h
PSTDLIB_H=../play/pstdlib.h
PSTDIO_H=../play/pstdio.h
PLAY_ALL=$(PLAY_H) $(PSTDLIB_H) $(PSTDIO_H)

cgm.o: cgm.h gtext.h $(PLAY_ALL)  gist.h engine.h
clip.o: clip.h $(PSTDLIB_H)
draw.o: draw.c draw.h gtext.h $(PSTDLIB_H)  gist.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(NO_EXP10) -c draw.c
draw0.o: draw.h gtext.h $(PLAY_H) $(PSTDLIB_H)  gist.h
engine.o: engine.h gist.h draw.h $(PSTDLIB_H)
gcntr.o: gist.h
gist.o: gist.h engine.h clip.h $(PLAY_H) $(PSTDLIB_H)
gread.o: gread.c gist.h $(PLAY_ALL)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(D_GISTPATH) -c gread.c
gtext.o: gtext.h  gist.h
hlevel.o: hlevel.h xbasic.h engine.h $(PSTDIO_H)  gist.h $(PLAY_H)
ps.o: ps.h gtext.h $(PLAY_H) $(PSTDIO_H)  gist.h engine.h
tick.o: tick.c gist.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(NO_EXP10) -c tick.c
tick60.o: gist.h
xbasic.o: xbasic.h gtext.h $(PSTDIO_H)  gist.h engine.h $(PLAY_H)
xfancy.o: xfancy.c xfancy.h draw.h  gist.h xbasic.h engine.h $(PLAY_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(NO_EXP10) -c xfancy.c

browser.o: xbasic.h cgm.h cgmin.h eps.h   gist.h engine.h $(PLAY_ALL)
cgmin.o: cgmin.h engine.h $(PLAY_ALL)   gist.h
eps.o: eps.h ps.h gtext.h $(PSTDIO_H) $(PSTDLIB_H)   gist.h engine.h

bench.o: hlevel.h xfancy.h $(PLAY_ALL)  gist.h xbasic.h engine.h

# Always compile main.c with -g so new packages can be debugged.
main.o: main.c
	$(CC) $(CPPFLAGS) -I../play -g -c main.c
