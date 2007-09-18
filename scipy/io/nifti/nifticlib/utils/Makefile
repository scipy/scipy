CC	=	gcc
CFLAGS	=	-ansi -pedantic

FSLIO_INCS	=	-I../include
NIFTI_INCS	=	-I../include
ZNZ_INCS	=	-I/usr/include

FSLIO_LIBS	=	-L../lib -lfslio
NIFTI_LIBS	=	-L../lib -lniftiio
NIFTICDF_LIBS	=	-L../lib -lnifticdf
ZNZ_LIBS	=	-L/usr/lib -L../lib -lznz -lm -lz



##############################################################
# platform specific redefines (to use, set ARCH appropriately)

## ARCH = X86_64

ifeq ($(ARCH),SGI) ## SGI 32bit
ZNZ_INCS	=	-I/usr/freeware/include
ZNZ_LIBS	=	-L/usr/freeware/lib32 -L../lib -lznz -lm -lz
else
ifeq ($(ARCH),I386) ## 32-bit Linux
ZNZ_INCS	=	-I/usr/include
ZNZ_LIBS	=	-L/usr/lib -L../lib -lznz -lm -lz
else
ifeq ($(ARCH),X86_64) ## 64-bit Linux
ZNZ_INCS	=
ZNZ_LIBS	=	-L../lib -lznz -lm -lz
endif
endif
endif

####################################################


UTILS_PROGS     =       nifti_stats nifti_tool nifti1_test



all:	$(UTILS_PROGS)

clean:
	rm -f $(UTILS_PROGS)


nifti_stats: nifti_stats.c
	$(CC) $(CFLAGS) -o nifti_stats nifti_stats.c $(NIFTI_INCS) $(ZNZ_INCS) $(NIFTICDF_LIBS) $(NIFTI_LIBS) $(ZNZ_LIBS)

nifti_tool: nifti_tool.c nifti_tool.h
	$(CC) $(CFLAGS) -Wall -o nifti_tool nifti_tool.c $(NIFTI_INCS) $(ZNZ_INCS) $(NIFTI_LIBS) $(ZNZ_LIBS)

nifti1_test: nifti1_test.c
	$(CC) $(CFLAGS) -o nifti1_test nifti1_test.c $(NIFTI_INCS) $(ZNZ_INCS) $(NIFTI_LIBS) $(ZNZ_LIBS)
