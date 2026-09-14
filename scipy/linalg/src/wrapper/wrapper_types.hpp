/**
 * @file
 * @brief The type vocabulary shared by the BLAS/LAPACK wrappers.
 *
 * This header defines the four types for the BLAS/LAPACK wrappers and some helper
 * template specializations to provide convenience.
 *
 */
#pragma once

#include <complex>
#include <type_traits>   /* std::false_type / std::true_type, for is_complex */

namespace wrapper {

    /* numpy-style width aliases (float32, ..., complex128); they keep the s/d/c/z flavor
     * columns of the declaration tables short and aligned. */
    using f32  = float;
    using f64  = double;
    using c64  = std::complex<float>;
    using c128 = std::complex<double>;

    /**
     * @brief Real counterpart of a flavor: f32 -> f32, c64 -> f32, c128 -> f64.
     *
     * The scalar type of the real arguments in the mixed real/complex families (the pyf
     * `<ftypereal>`): `c`/`s` in `csrot`/`zdrot`, the scale factor in `csscal`/`zdscal`,
     * LAPACK's real `rwork` buffers alongside a complex flavor.
     */
    template <class T> struct real_of                  { using type = T; };
    template <class T> struct real_of<std::complex<T>> { using type = T; };
    template <class T> using real_of_t = typename real_of<T>::type;

    /**
     * @brief Whether a flavor is complex, for the `if constexpr` branches where the real and
     *        complex routines take genuinely different argument lists (`gees`'s split
     *        `wr`/`wi` versus a single complex `w`).
     */
    template <class T> struct is_complex                  : std::false_type {};
    template <class T> struct is_complex<std::complex<T>> : std::true_type  {};
    template <class T> inline constexpr bool is_complex_v = is_complex<T>::value;

    /**
     * @brief Scalar type -> routine-name prefix, regular naming (s/d/c/z): `saxpy`, `zgemv`, ...
     *
     * The irregular prefixes are library-specific and live with their library: BLAS's
     * two-letter Level-1 families in `blas_helpers.hpp`, LAPACK's in `lapack_helpers.hpp`.
     */
    template <class T> constexpr const char *flavor();
    template <> inline constexpr const char *flavor<f32>()  { return "s"; }
    template <> inline constexpr const char *flavor<f64>()  { return "d"; }
    template <> inline constexpr const char *flavor<c64>()  { return "c"; }
    template <> inline constexpr const char *flavor<c128>() { return "z"; }

}  // namespace wrapper
