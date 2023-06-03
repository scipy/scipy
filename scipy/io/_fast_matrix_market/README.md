# fast_matrix_market

Fast and full-featured Matrix Market I/O library for C++ and Python.

Originally written by Adam Lugowski and published at https://github.com/alugowski/fast_matrix_market

The Python bindings are interchangeable with scipy.io._mmio methods to read and write Matrix Market files.

The core implementation is written in C++17 and exposed to Python with pybind11.
The C++ code directly reads and writes numpy arrays for a zero-copy interface.

Parallelism is via a thread pool built on C++ threads.

# Dependencies

The core speedup over a pure Python implementation comes from two sources: fast sequential parsing/formatting and parallelism.

C++17 offers `<charconv>` with very fast sequential methods. Integer overloads work very well, but the floating-point overloads are only available in very recent versions of only some compilers.

Note that the standard library methods like `strtod` internally lock on the system locale.
This locking severely restricts usable parallelism with these methods.

The dependencies below fill the gap until floating-point `<charconv>` support becomes dependable.

* [fast_float](fast_matrix_market/dependencies/fast_float) is a fast implementation of floating-point `std::from_chars`. This algorithm has been incorporated into GCC 12, LLVM as of 2021, and with ongoing work to update MSVC's existing floating-point `std::from_chars` to fast_float.
* [Ryu](fast_matrix_market/dependencies/ryu) is used for floating-point formatting. Basis for many of the new `to_chars` implementations.

For future maintainers:
Once floating-point `<charconv>` is available then add the following definitions to enable it:
* `'-DFMM_FROM_CHARS_DOUBLE_SUPPORTED'`
* `'-DFMM_FROM_CHARS_LONG_DOUBLE_SUPPORTED'`
* `'-DFMM_TO_CHARS_DOUBLE_SUPPORTED'`
* `'-DFMM_TO_CHARS_LONG_DOUBLE_SUPPORTED'`

and the dependencies may be dropped.