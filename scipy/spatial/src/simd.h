/*
 * SIMD definitions. This file defines types for SIMD (SSE, Neon, AVX)
 * operations, using the GCC vector extensions. These are supported by GCC
 * and Clang, are architecture-independent, and work even when no SIMD
 * instructions are available (e.g., ia32 with default GCC options).
 *
 * This header defines the macro SCIPY_SIMD if SIMD types are enabled.
 *
 * TODO move to scipy/_lib and reuse in other modules!
 */

#ifndef SCIPY_SIMD_H
#define SCIPY_SIMD_H

/* XXX earlier versions of GCC and Clang might work */
#if defined(__GNUC__) && __GNUC__ > 4 \
    || defined(__clang__) && __clang_major__ >= 3 && __clang_minor >= 4
#define SCIPY_SIMD

typedef double Double2 __attribute__((vector_size(2 * sizeof(double))));
typedef double Double4 __attribute__((vector_size(4 * sizeof(double))));

#endif

#endif
