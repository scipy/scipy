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

# Code Outline

The C++ library is in [fast_matrix_market/include/fast_matrix_market/](fast_matrix_market/include/fast_matrix_market),
and the pybind11-based Python bindings are in [src/](src).

## C++ Library

Matrix Market files are read and written in two phases: header and body.

The header routines are in `header.hpp` and `types.hpp`.

The entrypoints for the body methods are in `read_body.hpp` and `write_body.hpp`, with `*_threads.hpp` versions that
parallelize the core method.

`parse_handlers.hpp` and `formatters.hpp` contain methods to define I/O from file to datastructure and datastructure to file, respectively.
These include the logic on how to read or populate dense arrays, coordinate matrices, compressed CSR/CSR, etc.
Both parse handlers and formatters allow breaking down the datastructure into independent chunks for processing in parallel.

`thirdparty/task_thread_pool` is a simple thread pool class similar to Python's `ThreadPoolExecutor`.
This works well and avoids a dependency on a heavy package like TBB.

`field_conv.hpp` contains abstractions to read/write integer, floating point, and complex numbers. This exists primarily
to use the best available method from:

* fast_float (reads) or ryu (writes)
* `<charconv>` methods `std::from_chars` and `std::to_chars`
* stdlib methods like `strtod()`

`chunking.hpp` enables reading from a stream in large chunks with newline boundaries.

## Python bindings

The Python bindings build an extension called `_core` from `src/_core*.cpp`.
[src/pystreambuf.h](src/pystreambuf.h) adapts Python stream objects to C++ streams.
This enables using `BytesIO`, `GZipFile`, etc, from C++.
