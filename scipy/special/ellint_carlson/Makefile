CC = gcc-9
CFLAGS ?= -std=c11 -fstack-protector -D_FORTIFY_SOURCE=1 -DNDEBUG
#CFLAGS = -std=c11 -UNDEBUG -g -ggdb -glldb
WFLAGS = -W -Wall -Wextra -pedantic -Wshadow -Wduplicated-branches -Wduplicated-cond -Winit-self -Wunused -fstrict-aliasing -Wstrict-aliasing -Wtrampolines -Wunsafe-loop-optimizations -Wexpansion-to-defined -Wbad-function-cast -Wcast-qual -Wcast-align -Wconversion -Wmissing-prototypes -Wpadded -Wredundant-decls -Wnested-externs -Winline -Wdisabled-optimization 
OPTFLAGS ?= -march=native -mtune=native -O3 -fno-unsafe-math-optimizations -frounding-math -flto -fgcse-sm -fgcse-las -ftree-vectorize
#OPTFLAGS = -O0 -fno-unsafe-math-optimizations -fno-tree-vectorize
export CC CFLAGS WFLAGS OPTFLAGS

.PHONY: clean test

CTEMPLATES = _rf_impl.c _rc_impl.c _rd_impl.c _rj_impl.c
HEADERS = ellint_carlson.h ellint_common.h ellint_poly.h ellint_argcheck.h

libellint_carlson.so: ellint_real.o ellint_complex.o
	$(CC) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) -fPIC -shared $+ -o $@ -lm

ellint_real.o: ellint_real.c $(CTEMPLATES) $(HEADERS)
	$(CC) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) -fPIC -c $< -I.

ellint_complex.o: ellint_complex.c $(CTEMPLATES) $(HEADERS)
	$(CC) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) -fPIC -c $< -I.

test: libellint_carlson.so
	$(MAKE) -C tests test

clean:
	$(MAKE) -C tests clean
	rm -rf *.o *.so *.dylib *.dSYM
