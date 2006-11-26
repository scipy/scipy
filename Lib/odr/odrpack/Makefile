# Make the ODRPACK library

FC = g77
FOPT = -O
LIB = libodrpack.a
OBJS = d_odr.o d_mprec.o d_lpkbls.o dlunoc.o

.f.o:
	$(FC) $(FOPT) -c $<

$(LIB): $(OBJS)
	ar cru $(LIB) $(OBJS)

all: $(LIB)

clean:
	rm -f $(LIB) $(OBJS)

