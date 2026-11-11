The `superlu` directory is populated by `tools/vendoring/vendor_superlu.sh`.

It contains the sources of [SuperLU](https://github.com/xiaoyeli/superlu),
unmodified: all `.c` and `.h` files from upstream's `SRC` folder, except
`mc64ad.c` which is license-incompatible, plus `superlu_config.h.in`.

SciPy-specific configuration (allocators, BLAS integer type and BLAS symbol
renaming) lives in `scipy_superlu_config.h`, which `meson.build` copies to
`superlu_config.h` in the build directory, in place of the header SuperLU would
normally generate from `superlu_config.h.in`.

The allocator and abort hooks declared in `scipy_superlu_config.h` are
implemented by the `_superlu` extension module in
`scipy/sparse/linalg/_dsolve`.
