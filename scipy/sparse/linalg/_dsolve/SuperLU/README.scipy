This directory contains a bundled version of SuperLU from
https://github.com/xiaoyeli/superlu at
commit a8d2543611ba38b129b619751ca2a7284e624a07.
We use all .c and .h files except `mc64ad.c` which is license
incompatible from the SRC folder.

SciPy-specific configuration (allocators, BLAS integer type and BLAS symbol
renaming) is injected via `../scipy_superlu_config.h`, which meson copies to
`superlu_config.h` in the build directory in place of the header SuperLU would
normally generate from `superlu_config.h.in`.
