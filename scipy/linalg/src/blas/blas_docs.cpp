/**
 * @file
 * @brief Docstrings for the `_fblas` / `_fblas_64` wrappers, built on demand.
 *
 * Each wrapper's docstring is assembled at runtime from a template that is parameterised by the
 * element type.  The single entry point `blas::capi::build_doc(name)` is called lazily by the
 * `BlasFunc.__doc__` getter (see `blas_module.cpp`) the first time a routine's docstring is
 * requested, and the result is cached on the callable.  Nothing is built at import time, and
 * only the shared templates -- not ~130 fully expanded strings -- live in the compiled module.
 *
 * Adding a routine: give its family a `doc_<family>(name, Dtype)` builder and list its typed
 * names in `doc_table` -- `DOC_FAMILY(fam)` for the four regular `s`/`d`/`c`/`z` prefixes,
 * `DOC_FAMILY_R(fam)` / `DOC_FAMILY_C(fam)` for a real-only (`s`/`d`) or complex-only (`c`/`z`)
 * pair, or `DOC_FAMILY2(...)` for families whose names are irregular (`scnrm2`, `isamax`, ...).
 * Routines that genuinely need their own text (mixed-type, tuple-returning) get an explicit row.
 */
#include <Python.h>
#include <cstring>
#include <string>

namespace blas {
    namespace capi {

        /**
         * @brief The element-type-dependent facts a docstring template needs.
         *
         * `scalar` is the numpydoc type a real-or-complex scalar argument takes (`float` or
         * `complex`); `is_complex` selects the wording where the real and complex routines
         * genuinely differ (e.g. what `asum` sums).  Single- vs double-precision does not appear
         * in the prose -- it is carried by the routine name -- so `S`/`D` and `C`/`Z` coincide.
         */
        struct Dtype {
            const char *scalar;
            bool is_complex;
        };

        static constexpr Dtype S{"float", false};
        static constexpr Dtype D{"float", false};
        static constexpr Dtype C{"complex", true};
        static constexpr Dtype Z{"complex", true};

        /* Complete numpydoc `Parameters` entries for the offset/stride/count arguments common to
         * the Level-1 routines, written once so every routine documents them identically.  Drop
         * one into a docstring with `s += P_INCX;`. */
        static constexpr const char *P_N           = "n : int, optional\n    Number of elements to process. Default is the whole vector.\n";
        static constexpr const char *P_OFFX        = "offx : int, optional\n    Index of the first element of `x` to use. Default is 0.\n";
        static constexpr const char *P_INCX        = "incx : int, optional\n    Stride between successive elements of `x`. Default is 1.\n";
        static constexpr const char *P_OFFY        = "offy : int, optional\n    Index of the first element of `y` to use. Default is 0.\n";
        static constexpr const char *P_INCY        = "incy : int, optional\n    Stride between successive elements of `y`. Default is 1.\n";
        static constexpr const char *P_OVERWRITE_X = "overwrite_x : int, optional\n    If nonzero, `x` may be overwritten in place. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_Y = "overwrite_y : int, optional\n    If nonzero, `y` may be overwritten in place. Default is 0.\n";

        /* Shared Level-2 fragments.  The `trans`/`lower`/`diag` flags are integers in every L2
         * wrapper (never characters), so their encodings are documented once here. */
        static constexpr const char *P_TRANS        = "trans : int, optional\n    Operation on the matrix: 0 none, 1 transpose, 2 conjugate transpose. Default is 0.\n";
        static constexpr const char *P_LOWER        = "lower : int, optional\n    If 0, the upper triangle is used; if 1, the lower. Default is 0.\n";
        static constexpr const char *P_DIAG         = "diag : int, optional\n    If 0, the matrix has a non-unit diagonal; if 1, it is unit-triangular. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_A  = "overwrite_a : int, optional\n    If nonzero, `a` may be overwritten in place. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_AP = "overwrite_ap : int, optional\n    If nonzero, `ap` may be overwritten in place. Default is 0.\n";

        /* Shared Level-3 fragments (all flags are integers here too). */
        static constexpr const char *P_TRANS_A      = "trans_a : int, optional\n    Operation on `a`: 0 none, 1 transpose, 2 conjugate transpose. Default is 0.\n";
        static constexpr const char *P_TRANS_B      = "trans_b : int, optional\n    Operation on `b`: 0 none, 1 transpose, 2 conjugate transpose. Default is 0.\n";
        static constexpr const char *P_SIDE         = "side : int, optional\n    If 0, `a` multiplies from the left; if 1, from the right. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_C  = "overwrite_c : int, optional\n    If nonzero, `c` may be overwritten in place. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_B  = "overwrite_b : int, optional\n    If nonzero, `b` may be overwritten in place. Default is 0.\n";

        static std::string
        doc_axpy(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, y, n=None, a=1.0, offx=0, incx=1, offy=0, incy=1)\n\n";
            s += "Compute ``y = a*x + y`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input vector.\n";
            s += "y : ndarray\n    Input/output vector; on exit it holds ``a*x + y``.\n";
            s += P_N;
            s += "a : " + std::string(t.scalar) + ", optional\n    Scalar multiplier. Default is 1.\n";
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;

            s += "\nReturns\n-------\n";
            s += "y : ndarray\n    The updated vector ``a*x + y``.\n";
            return s;
        }

        static std::string
        doc_nrm2(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, n=None, offx=0, incx=1)\n\n";
            s += "Compute the Euclidean norm ``sqrt(sum(abs(x)**2))`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input vector.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;

            s += "\nReturns\n-------\n";
            s += "nrm : float\n    The Euclidean norm of `x`.\n";
            return s;
        }

        static std::string
        doc_asum(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, n=None, offx=0, incx=1)\n\n";
            if (t.is_complex) {
                s += "Compute ``sum(abs(x.real) + abs(x.imag))`` (BLAS ``" + std::string(name) + "``).\n\n";
            } else {
                s += "Compute the sum of absolute values ``sum(abs(x))`` (BLAS ``" + std::string(name) + "``).\n\n";
            }

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input vector.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;

            s += "\nReturns\n-------\n";
            s += "s : float\n    The sum of absolute values.\n";
            return s;
        }

        static std::string
        doc_iamax(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, n=None, offx=0, incx=1)\n\n";
            s += "Return the index of the element of `x` with the largest absolute value "
                 "(BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input vector.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;

            s += "\nReturns\n-------\n";
            s += "k : int\n";
            if (t.is_complex) {
                s += "    Index of the element maximizing ``abs(x.real) + abs(x.imag)``. This is a\n";
            } else {
                s += "    Index of the element with the largest absolute value. This is a\n";
            }
            s += "    0-based index into `x`; the underlying BLAS routine returns a 1-based index,\n";
            s += "    which SciPy shifts by one.\n";
            return s;
        }

        static std::string
        doc_copy(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, y, n=None, offx=0, incx=1, offy=0, incy=1)\n\n";
            s += "Copy `x` into `y` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input vector.\n";
            s += "y : ndarray\n    Output vector; on exit it holds a copy of `x`.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;

            s += "\nReturns\n-------\n";
            s += "y : ndarray\n    The vector `y` holding the copied elements of `x`.\n";
            return s;
        }

        static std::string
        doc_swap(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, y, n=None, offx=0, incx=1, offy=0, incy=1)\n\n";
            s += "Swap the contents of `x` and `y` in place (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input/output vector; on exit it holds the original `y`.\n";
            s += "y : ndarray\n    Input/output vector; on exit it holds the original `x`.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    The vector `x` after the swap.\n";
            s += "y : ndarray\n    The vector `y` after the swap.\n";
            return s;
        }

        static std::string
        doc_scal(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, x, n=None, offx=0, incx=1)\n\n";
            s += "Scale `x` in place by `a`: ``x = a*x`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : " + std::string(t.scalar) + "\n    Scalar multiplier.\n";
            s += "x : ndarray\n    Input/output vector; overwritten with ``a*x``.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    The scaled vector ``a*x``.\n";
            return s;
        }

        /* csscal/zdscal: the scale factor is real even though `x` is complex, and the input is
         * copied unless overwrite_x is set (unlike the regular scal, which always works in place). */
        static std::string
        doc_scal_real(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, x, n=None, offx=0, incx=1, overwrite_x=0)\n\n";
            s += "Scale complex `x` by the real scalar `a`: ``x = a*x`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : float\n    Real scalar multiplier.\n";
            s += "x : ndarray\n    Complex input vector.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;
            s += P_OVERWRITE_X;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    The scaled vector ``a*x``.\n";
            return s;
        }

        static std::string
        doc_dot(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, y, n=None, offx=0, incx=1, offy=0, incy=1)\n\n";
            s += "Compute the dot product ``sum(x*y)`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input vector.\n";
            s += "y : ndarray\n    Input vector.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;

            s += "\nReturns\n-------\n";
            s += "xy : float\n    The dot product of `x` and `y`.\n";
            return s;
        }

        static std::string
        doc_dotu(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, y, n=None, offx=0, incx=1, offy=0, incy=1)\n\n";
            s += "Compute the unconjugated dot product ``sum(x*y)`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input vector.\n";
            s += "y : ndarray\n    Input vector.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;

            s += "\nReturns\n-------\n";
            s += "xy : complex\n    The dot product ``x . y``.\n";
            return s;
        }

        static std::string
        doc_dotc(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, y, n=None, offx=0, incx=1, offy=0, incy=1)\n\n";
            s += "Compute the dot product with `x` conjugated: ``sum(conj(x)*y)`` "
                 "(BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input vector (conjugated).\n";
            s += "y : ndarray\n    Input vector.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;

            s += "\nReturns\n-------\n";
            s += "xy : complex\n    The dot product ``conj(x) . y``.\n";
            return s;
        }

        static std::string
        doc_rotg(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b)\n\n";
            s += "Construct a Givens plane rotation ``(c, s)`` that eliminates `b` "
                 "(BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : " + std::string(t.scalar) + "\n    First component of the vector to rotate.\n";
            s += "b : " + std::string(t.scalar) + "\n    Second component of the vector to rotate.\n";

            s += "\nReturns\n-------\n";
            s += "c : " + std::string(t.scalar) + "\n    Cosine of the rotation.";
            if (t.is_complex) { s += "  Its imaginary part is zero."; }
            s += "\n";
            s += "s : " + std::string(t.scalar) + "\n    Sine of the rotation.\n";
            return s;
        }

        static std::string
        doc_rotmg(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d1, d2, x1, y1)\n\n";
            s += "Construct a modified Givens rotation (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "d1 : float\n    Scaling factor for the first coordinate.\n";
            s += "d2 : float\n    Scaling factor for the second coordinate.\n";
            s += "x1 : float\n    First coordinate of the vector to rotate.\n";
            s += "y1 : float\n    Second coordinate of the vector to rotate.\n";

            s += "\nReturns\n-------\n";
            s += "param : ndarray\n    Length-5 array describing the rotation: ``param[0]`` is a flag\n";
            s += "    selecting the matrix form and ``param[1:5]`` are its entries.\n";
            return s;
        }

        static std::string
        doc_rot(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, y, c, s, n=None, offx=0, incx=1, offy=0, incy=1, "
                 "overwrite_x=0, overwrite_y=0)\n\n";
            s += "Apply a Givens plane rotation to `x` and `y`: "
                 "``x, y = c*x + s*y, c*y - s*x`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input/output vector; on exit it holds ``c*x + s*y``.\n";
            s += "y : ndarray\n    Input/output vector; on exit it holds ``c*y - s*x``.\n";
            if (t.is_complex) {
                s += "c : float\n    Cosine of the rotation (real, even though `x` and `y` are complex).\n";
                s += "s : float\n    Sine of the rotation (real, even though `x` and `y` are complex).\n";
            } else {
                s += "c : float\n    Cosine of the rotation.\n";
                s += "s : float\n    Sine of the rotation.\n";
            }
            s += P_N;
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;
            s += P_OVERWRITE_X;
            s += P_OVERWRITE_Y;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    The rotated vector `x`.\n";
            s += "y : ndarray\n    The rotated vector `y`.\n";
            return s;
        }

        static std::string
        doc_rotm(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, y, param, n=None, offx=0, incx=1, offy=0, incy=1, "
                 "overwrite_x=0, overwrite_y=0)\n\n";
            s += "Apply a modified Givens rotation to `x` and `y` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    Input/output vector.\n";
            s += "y : ndarray\n    Input/output vector.\n";
            s += "param : ndarray\n    Length-5 array describing the rotation, as returned by ``rotmg``.\n";
            s += P_N;
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;
            s += P_OVERWRITE_X;
            s += P_OVERWRITE_Y;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    The transformed vector `x`.\n";
            s += "y : ndarray\n    The transformed vector `y`.\n";
            return s;
        }

        /* ---- Level 2 ---------------------------------------------------------------------- */

        static std::string
        doc_gemv(const char *name, const Dtype &t) noexcept
        {
            const std::string sc = t.scalar;
            std::string s;
            s += std::string(name) + "(alpha, a, x, beta=0.0, y=None, offx=0, incx=1, offy=0, incy=1, trans=0, overwrite_y=0)\n\n";
            s += "Compute the matrix-vector product ``y = alpha*op(a)@x + beta*y`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + sc + "\n    Scalar multiplier for the product.\n";
            s += "a : ndarray\n    Input matrix.\n";
            s += "x : ndarray\n    Input vector.\n";
            s += "beta : " + sc + ", optional\n    Scalar multiplier for `y`. Default is 0.\n";
            s += "y : ndarray, optional\n    Input/output vector; if omitted a new zero vector is allocated. Default is None.\n";
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;
            s += P_TRANS;
            s += P_OVERWRITE_Y;

            s += "\nReturns\n-------\n";
            s += "y : ndarray\n    The vector ``alpha*op(a)@x + beta*y``, where ``op(a)`` is `a`, ``a.T``,\n";
            s += "    or ``a.conj().T`` according to `trans`.\n";
            return s;
        }

        static std::string
        doc_gbmv(const char *name, const Dtype &t) noexcept
        {
            const std::string sc = t.scalar;
            std::string s;
            s += std::string(name) + "(m, n, kl, ku, alpha, a, x, incx=1, offx=0, beta=0.0, y=None, incy=1, offy=0, trans=0, overwrite_y=0)\n\n";
            s += "Compute ``y = alpha*op(a)@x + beta*y`` for a banded matrix `a` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "m : int\n    Number of rows of the full matrix.\n";
            s += "n : int\n    Number of columns of the full matrix.\n";
            s += "kl : int\n    Number of sub-diagonals.\n";
            s += "ku : int\n    Number of super-diagonals.\n";
            s += "alpha : " + sc + "\n    Scalar multiplier for the product.\n";
            s += "a : ndarray\n    Banded storage of the matrix, shape ``(kl + ku + 1, n)``.\n";
            s += "x : ndarray\n    Input vector.\n";
            s += P_INCX;
            s += P_OFFX;
            s += "beta : " + sc + ", optional\n    Scalar multiplier for `y`. Default is 0.\n";
            s += "y : ndarray, optional\n    Input/output vector; if omitted a new zero vector is allocated. Default is None.\n";
            s += P_INCY;
            s += P_OFFY;
            s += P_TRANS;
            s += P_OVERWRITE_Y;

            s += "\nReturns\n-------\n";
            s += "y : ndarray\n    The vector ``alpha*op(a)@x + beta*y``.\n";
            return s;
        }

        /* sbmv (symmetric) / hbmv (Hermitian) banded matrix-vector product. */
        static std::string
        doc_sbmv_hbmv(const char *name, const Dtype &t, const char *kind) noexcept
        {
            const std::string sc = t.scalar;
            std::string s;
            s += std::string(name) + "(k, alpha, a, x, incx=1, offx=0, beta=0.0, y=None, incy=1, offy=0, lower=0, overwrite_y=0)\n\n";
            s += "Compute ``y = alpha*a@x + beta*y`` for a banded " + std::string(kind) + " matrix `a` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "k : int\n    Number of super-diagonals (or sub-diagonals) of the band.\n";
            s += "alpha : " + sc + "\n    Scalar multiplier for the product.\n";
            s += "a : ndarray\n    Banded storage of the " + std::string(kind) + " matrix, shape ``(k + 1, n)``.\n";
            s += "x : ndarray\n    Input vector.\n";
            s += P_INCX;
            s += P_OFFX;
            s += "beta : " + sc + ", optional\n    Scalar multiplier for `y`. Default is 0.\n";
            s += "y : ndarray, optional\n    Input/output vector; if omitted a new zero vector is allocated. Default is None.\n";
            s += P_INCY;
            s += P_OFFY;
            s += P_LOWER;
            s += P_OVERWRITE_Y;

            s += "\nReturns\n-------\n";
            s += "y : ndarray\n    The vector ``alpha*a@x + beta*y``.\n";
            return s;
        }
        static std::string doc_sbmv(const char *name, const Dtype &t) { return doc_sbmv_hbmv(name, t, "symmetric"); }
        static std::string doc_hbmv(const char *name, const Dtype &t) { return doc_sbmv_hbmv(name, t, "Hermitian"); }

        /* symv (symmetric) / hemv (Hermitian) matrix-vector product. */
        static std::string
        doc_symv_hemv(const char *name, const Dtype &t, const char *kind) noexcept
        {
            const std::string sc = t.scalar;
            std::string s;
            s += std::string(name) + "(alpha, a, x, beta=0.0, y=None, offx=0, incx=1, offy=0, incy=1, lower=0, overwrite_y=0)\n\n";
            s += "Compute ``y = alpha*a@x + beta*y`` for a " + std::string(kind) + " matrix `a` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + sc + "\n    Scalar multiplier for the product.\n";
            s += "a : ndarray\n    " + std::string(kind) + " matrix; only the triangle selected by `lower` is referenced.\n";
            s += "x : ndarray\n    Input vector.\n";
            s += "beta : " + sc + ", optional\n    Scalar multiplier for `y`. Default is 0.\n";
            s += "y : ndarray, optional\n    Input/output vector; if omitted a new zero vector is allocated. Default is None.\n";
            s += P_OFFX;
            s += P_INCX;
            s += P_OFFY;
            s += P_INCY;
            s += P_LOWER;
            s += P_OVERWRITE_Y;

            s += "\nReturns\n-------\n";
            s += "y : ndarray\n    The vector ``alpha*a@x + beta*y``.\n";
            return s;
        }
        static std::string doc_symv(const char *name, const Dtype &t) { return doc_symv_hemv(name, t, "symmetric"); }
        static std::string doc_hemv(const char *name, const Dtype &t) { return doc_symv_hemv(name, t, "Hermitian"); }

        /* spmv (symmetric) / hpmv (Hermitian) packed matrix-vector product. */
        static std::string
        doc_spmv_hpmv(const char *name, const Dtype &t, const char *kind) noexcept
        {
            const std::string sc = t.scalar;
            std::string s;
            s += std::string(name) + "(n, alpha, ap, x, incx=1, offx=0, beta=0.0, y=None, incy=1, offy=0, lower=0, overwrite_y=0)\n\n";
            s += "Compute ``y = alpha*a@x + beta*y`` for a packed " + std::string(kind) + " matrix (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "n : int\n    Order of the matrix.\n";
            s += "alpha : " + sc + "\n    Scalar multiplier for the product.\n";
            s += "ap : ndarray\n    The " + std::string(kind) + " matrix in packed storage.\n";
            s += "x : ndarray\n    Input vector.\n";
            s += P_INCX;
            s += P_OFFX;
            s += "beta : " + sc + ", optional\n    Scalar multiplier for `y`. Default is 0.\n";
            s += "y : ndarray, optional\n    Input/output vector; if omitted a new zero vector is allocated. Default is None.\n";
            s += P_INCY;
            s += P_OFFY;
            s += P_LOWER;
            s += P_OVERWRITE_Y;

            s += "\nReturns\n-------\n";
            s += "y : ndarray\n    The vector ``alpha*a@x + beta*y``.\n";
            return s;
        }
        static std::string doc_spmv(const char *name, const Dtype &t) { return doc_spmv_hpmv(name, t, "symmetric"); }
        static std::string doc_hpmv(const char *name, const Dtype &t) { return doc_spmv_hpmv(name, t, "Hermitian"); }

        /* spr / hpr: packed rank-1 update.  spr's alpha has the matrix type; hpr's alpha is real. */
        static std::string
        doc_spr(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, alpha, x, ap, incx=1, offx=0, lower=0, overwrite_ap=0)\n\n";
            s += "Packed symmetric rank-1 update ``ap = alpha*outer(x, x) + ap`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "n : int\n    Order of the matrix.\n";
            s += "alpha : " + std::string(t.scalar) + "\n    Scalar multiplier.\n";
            s += "x : ndarray\n    Input vector.\n";
            s += "ap : ndarray\n    Symmetric matrix in packed storage, updated in place.\n";
            s += P_INCX;
            s += P_OFFX;
            s += P_LOWER;
            s += P_OVERWRITE_AP;

            s += "\nReturns\n-------\n";
            s += "ap : ndarray\n    The updated packed matrix.\n";
            return s;
        }

        static std::string
        doc_hpr(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, alpha, x, ap, incx=1, offx=0, lower=0, overwrite_ap=0)\n\n";
            s += "Packed Hermitian rank-1 update ``ap = alpha*outer(x, conj(x)) + ap`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "n : int\n    Order of the matrix.\n";
            s += "alpha : float\n    Real scalar multiplier.\n";
            s += "x : ndarray\n    Input vector.\n";
            s += "ap : ndarray\n    Hermitian matrix in packed storage, updated in place.\n";
            s += P_INCX;
            s += P_OFFX;
            s += P_LOWER;
            s += P_OVERWRITE_AP;

            s += "\nReturns\n-------\n";
            s += "ap : ndarray\n    The updated packed matrix.\n";
            return s;
        }

        /* spr2 / hpr2: packed rank-2 update. */
        static std::string
        doc_spr2_hpr2(const char *name, const Dtype &t, const char *kind, const char *formula) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, alpha, x, y, ap, incx=1, offx=0, incy=1, offy=0, lower=0, overwrite_ap=0)\n\n";
            s += "Packed " + std::string(kind) + " rank-2 update ``ap = " + std::string(formula) + " + ap`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "n : int\n    Order of the matrix.\n";
            s += "alpha : " + std::string(t.scalar) + "\n    Scalar multiplier.\n";
            s += "x : ndarray\n    First input vector.\n";
            s += "y : ndarray\n    Second input vector.\n";
            s += "ap : ndarray\n    " + std::string(kind) + " matrix in packed storage, updated in place.\n";
            s += P_INCX;
            s += P_OFFX;
            s += P_INCY;
            s += P_OFFY;
            s += P_LOWER;
            s += P_OVERWRITE_AP;

            s += "\nReturns\n-------\n";
            s += "ap : ndarray\n    The updated packed matrix.\n";
            return s;
        }
        static std::string doc_spr2(const char *name, const Dtype &t) { return doc_spr2_hpr2(name, t, "symmetric", "alpha*(outer(x, y) + outer(y, x))"); }
        static std::string doc_hpr2(const char *name, const Dtype &t) { return doc_spr2_hpr2(name, t, "Hermitian", "alpha*outer(x, conj(y)) + conj(alpha)*outer(y, conj(x))"); }

        /* syr: symmetric rank-1 update, returning the full matrix `a`. */
        static std::string
        doc_syr(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(alpha, x, lower=0, incx=1, offx=0, n=None, a=None, overwrite_a=0)\n\n";
            s += "Symmetric rank-1 update ``a = alpha*outer(x, x) + a`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + std::string(t.scalar) + "\n    Scalar multiplier.\n";
            s += "x : ndarray\n    Input vector.\n";
            s += P_LOWER;
            s += P_INCX;
            s += P_OFFX;
            s += "n : int, optional\n    Order of the matrix. Default is inferred from `x`.\n";
            s += "a : ndarray, optional\n    Symmetric matrix to update; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The updated matrix.\n";
            return s;
        }

        /* her: Hermitian rank-1 update.  alpha is accepted as complex but only its real part is used. */
        static std::string
        doc_her(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(alpha, x, lower=0, incx=1, offx=0, n=None, a=None, overwrite_a=0)\n\n";
            s += "Hermitian rank-1 update ``a = alpha*outer(x, conj(x)) + a`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : complex\n    Scalar multiplier. Only its real part is used.\n";
            s += "x : ndarray\n    Input vector.\n";
            s += P_LOWER;
            s += P_INCX;
            s += P_OFFX;
            s += "n : int, optional\n    Order of the matrix. Default is inferred from `x`.\n";
            s += "a : ndarray, optional\n    Hermitian matrix to update; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The updated matrix.\n";
            return s;
        }

        /* syr2 / her2: rank-2 update returning the full matrix `a`. */
        static std::string
        doc_syr2_her2(const char *name, const Dtype &t, const char *kind, const char *formula) noexcept
        {
            std::string s;
            s += std::string(name) + "(alpha, x, y, lower=0, incx=1, offx=0, incy=1, offy=0, n=None, a=None, overwrite_a=0)\n\n";
            s += std::string(kind) + " rank-2 update ``a = " + std::string(formula) + " + a`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + std::string(t.scalar) + "\n    Scalar multiplier.\n";
            s += "x : ndarray\n    First input vector.\n";
            s += "y : ndarray\n    Second input vector.\n";
            s += P_LOWER;
            s += P_INCX;
            s += P_OFFX;
            s += P_INCY;
            s += P_OFFY;
            s += "n : int, optional\n    Order of the matrix. Default is inferred from `x`.\n";
            s += "a : ndarray, optional\n    Matrix to update; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The updated matrix.\n";
            return s;
        }
        static std::string doc_syr2(const char *name, const Dtype &t) { return doc_syr2_her2(name, t, "Symmetric", "alpha*(outer(x, y) + outer(y, x))"); }
        static std::string doc_her2(const char *name, const Dtype &t) { return doc_syr2_her2(name, t, "Hermitian", "alpha*outer(x, conj(y)) + conj(alpha)*outer(y, conj(x))"); }

        /* ger / geru / gerc: general rank-1 update ``a = alpha*outer(x, op(y)) + a``. */
        static std::string
        doc_ger_family(const char *name, const Dtype &t, const char *outer) noexcept
        {
            std::string s;
            s += std::string(name) + "(alpha, x, y, incx=1, incy=1, a=None, overwrite_x=0, overwrite_y=0, overwrite_a=0)\n\n";
            s += "General rank-1 update ``a = alpha*" + std::string(outer) + " + a`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + std::string(t.scalar) + "\n    Scalar multiplier.\n";
            s += "x : ndarray\n    Input vector of length `m`.\n";
            s += "y : ndarray\n    Input vector of length `n`.\n";
            s += P_INCX;
            s += P_INCY;
            s += "a : ndarray, optional\n    ``(m, n)`` matrix to update; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_OVERWRITE_X;
            s += P_OVERWRITE_Y;
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The updated ``(m, n)`` matrix.\n";
            return s;
        }
        static std::string doc_ger(const char *name, const Dtype &t)  { return doc_ger_family(name, t, "outer(x, y)"); }
        static std::string doc_geru(const char *name, const Dtype &t) { return doc_ger_family(name, t, "outer(x, y)"); }
        static std::string doc_gerc(const char *name, const Dtype &t) { return doc_ger_family(name, t, "outer(x, conj(y))"); }

        /* Triangular matrix-vector product (trmv) or solve (trsv). */
        static std::string
        doc_trmv_trsv(const char *name, const char *action, const char *result) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, x, offx=0, incx=1, lower=0, trans=0, diag=0, overwrite_x=0)\n\n";
            s += std::string(action) + " (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Triangular matrix; only the triangle selected by `lower` is referenced.\n";
            s += "x : ndarray\n    Input/output vector.\n";
            s += P_OFFX;
            s += P_INCX;
            s += P_LOWER;
            s += P_TRANS;
            s += P_DIAG;
            s += P_OVERWRITE_X;

            s += "\nReturns\n-------\n";
            s += std::string("x : ndarray\n    ") + result + "\n";
            return s;
        }
        static std::string doc_trmv(const char *name, const Dtype &) { return doc_trmv_trsv(name, "Compute the triangular matrix-vector product ``x = op(a)@x``", "The product ``op(a)@x``."); }
        static std::string doc_trsv(const char *name, const Dtype &) { return doc_trmv_trsv(name, "Solve the triangular system ``op(a)@x = b`` in place", "The solution ``x``."); }

        /* Banded triangular matrix-vector product (tbmv) or solve (tbsv). */
        static std::string
        doc_tbmv_tbsv(const char *name, const char *action, const char *result) noexcept
        {
            std::string s;
            s += std::string(name) + "(k, a, x, incx=1, offx=0, lower=0, trans=0, diag=0, overwrite_x=0)\n\n";
            s += std::string(action) + " for a banded triangular matrix (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "k : int\n    Number of super-diagonals (or sub-diagonals) of the band.\n";
            s += "a : ndarray\n    Banded storage of the triangular matrix, shape ``(k + 1, n)``.\n";
            s += "x : ndarray\n    Input/output vector.\n";
            s += P_INCX;
            s += P_OFFX;
            s += P_LOWER;
            s += P_TRANS;
            s += P_DIAG;
            s += P_OVERWRITE_X;

            s += "\nReturns\n-------\n";
            s += std::string("x : ndarray\n    ") + result + "\n";
            return s;
        }
        static std::string doc_tbmv(const char *name, const Dtype &) { return doc_tbmv_tbsv(name, "Compute the triangular matrix-vector product ``x = op(a)@x``", "The product ``op(a)@x``."); }
        static std::string doc_tbsv(const char *name, const Dtype &) { return doc_tbmv_tbsv(name, "Solve the triangular system ``op(a)@x = b`` in place", "The solution ``x``."); }

        /* Packed triangular matrix-vector product (tpmv) or solve (tpsv). */
        static std::string
        doc_tpmv_tpsv(const char *name, const char *action, const char *result) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, ap, x, incx=1, offx=0, lower=0, trans=0, diag=0, overwrite_x=0)\n\n";
            s += std::string(action) + " for a packed triangular matrix (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "n : int\n    Order of the matrix.\n";
            s += "ap : ndarray\n    Triangular matrix in packed storage.\n";
            s += "x : ndarray\n    Input/output vector.\n";
            s += P_INCX;
            s += P_OFFX;
            s += P_LOWER;
            s += P_TRANS;
            s += P_DIAG;
            s += P_OVERWRITE_X;

            s += "\nReturns\n-------\n";
            s += std::string("x : ndarray\n    ") + result + "\n";
            return s;
        }
        static std::string doc_tpmv(const char *name, const Dtype &) { return doc_tpmv_tpsv(name, "Compute the triangular matrix-vector product ``x = op(a)@x``", "The product ``op(a)@x``."); }
        static std::string doc_tpsv(const char *name, const Dtype &) { return doc_tpmv_tpsv(name, "Solve the triangular system ``op(a)@x = b`` in place", "The solution ``x``."); }

        /* ---- Level 3 ---------------------------------------------------------------------- */

        static std::string
        doc_gemm(const char *name, const Dtype &t) noexcept
        {
            const std::string sc = t.scalar;
            std::string s;
            s += std::string(name) + "(alpha, a, b, beta=0.0, c=None, trans_a=0, trans_b=0, overwrite_c=0)\n\n";
            s += "Compute the matrix-matrix product ``c = alpha*op(a)@op(b) + beta*c`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + sc + "\n    Scalar multiplier for the product.\n";
            s += "a : ndarray\n    First input matrix.\n";
            s += "b : ndarray\n    Second input matrix.\n";
            s += "beta : " + sc + ", optional\n    Scalar multiplier for `c`. Default is 0.\n";
            s += "c : ndarray, optional\n    Input/output matrix; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_TRANS_A;
            s += P_TRANS_B;
            s += P_OVERWRITE_C;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The matrix ``alpha*op(a)@op(b) + beta*c``, where ``op(x)`` is `x`, ``x.T``,\n";
            s += "    or ``x.conj().T`` according to `trans_a` / `trans_b`.\n";
            return s;
        }

        /* symm (symmetric) / hemm (Hermitian) matrix-matrix product. */
        static std::string
        doc_symm_hemm(const char *name, const Dtype &t, const char *kind) noexcept
        {
            const std::string sc = t.scalar;
            std::string s;
            s += std::string(name) + "(alpha, a, b, beta=0.0, c=None, side=0, lower=0, overwrite_c=0)\n\n";
            s += "Compute ``c = alpha*a@b + beta*c`` for a " + std::string(kind) + " matrix `a` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + sc + "\n    Scalar multiplier for the product.\n";
            s += "a : ndarray\n    " + std::string(kind) + " matrix; only the triangle selected by `lower` is referenced.\n";
            s += "b : ndarray\n    General matrix.\n";
            s += "beta : " + sc + ", optional\n    Scalar multiplier for `c`. Default is 0.\n";
            s += "c : ndarray, optional\n    Input/output matrix; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_SIDE;
            s += P_LOWER;
            s += P_OVERWRITE_C;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The matrix ``alpha*a@b + beta*c`` (or ``alpha*b@a + beta*c`` if ``side=1``).\n";
            return s;
        }
        static std::string doc_symm(const char *name, const Dtype &t) { return doc_symm_hemm(name, t, "symmetric"); }
        static std::string doc_hemm(const char *name, const Dtype &t) { return doc_symm_hemm(name, t, "Hermitian"); }

        static std::string
        doc_syrk(const char *name, const Dtype &t) noexcept
        {
            const std::string sc = t.scalar;
            std::string s;
            s += std::string(name) + "(alpha, a, beta=0.0, c=None, trans=0, lower=0, overwrite_c=0)\n\n";
            s += "Symmetric rank-k update ``c = alpha*op(a)@op(a).T + beta*c`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + sc + "\n    Scalar multiplier for the product.\n";
            s += "a : ndarray\n    Input matrix.\n";
            s += "beta : " + sc + ", optional\n    Scalar multiplier for `c`. Default is 0.\n";
            s += "c : ndarray, optional\n    Symmetric matrix to update; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_TRANS;
            s += P_LOWER;
            s += P_OVERWRITE_C;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The updated matrix, where ``op(a)`` is `a` (``trans=0``) or ``a.T`` (``trans=1``).\n";
            return s;
        }

        /* herk: Hermitian rank-k.  alpha and beta are accepted as complex but only their real
         * parts are used (the result is Hermitian). */
        static std::string
        doc_herk(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(alpha, a, beta=0.0, c=None, trans=0, lower=0, overwrite_c=0)\n\n";
            s += "Hermitian rank-k update ``c = alpha*op(a)@op(a).conj().T + beta*c`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : complex\n    Scalar multiplier. Only its real part is used.\n";
            s += "a : ndarray\n    Input matrix.\n";
            s += "beta : complex, optional\n    Scalar multiplier for `c`. Only its real part is used. Default is 0.\n";
            s += "c : ndarray, optional\n    Hermitian matrix to update; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_TRANS;
            s += P_LOWER;
            s += P_OVERWRITE_C;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The updated matrix, where ``op(a)`` is `a` (``trans=0``) or ``a.conj().T`` (``trans=2``).\n";
            return s;
        }

        static std::string
        doc_syr2k(const char *name, const Dtype &t) noexcept
        {
            const std::string sc = t.scalar;
            std::string s;
            s += std::string(name) + "(alpha, a, b, beta=0.0, c=None, trans=0, lower=0, overwrite_c=0)\n\n";
            s += "Symmetric rank-2k update ``c = alpha*(op(a)@op(b).T + op(b)@op(a).T) + beta*c`` (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + sc + "\n    Scalar multiplier for the product.\n";
            s += "a : ndarray\n    First input matrix.\n";
            s += "b : ndarray\n    Second input matrix.\n";
            s += "beta : " + sc + ", optional\n    Scalar multiplier for `c`. Default is 0.\n";
            s += "c : ndarray, optional\n    Symmetric matrix to update; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_TRANS;
            s += P_LOWER;
            s += P_OVERWRITE_C;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The updated matrix.\n";
            return s;
        }

        /* her2k: Hermitian rank-2k.  alpha is used in full; beta is accepted as complex but only
         * its real part is used. */
        static std::string
        doc_her2k(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(alpha, a, b, beta=0.0, c=None, trans=0, lower=0, overwrite_c=0)\n\n";
            s += "Hermitian rank-2k update "
                 "``c = alpha*op(a)@op(b).conj().T + conj(alpha)*op(b)@op(a).conj().T + beta*c`` "
                 "(BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : complex\n    Scalar multiplier for the product.\n";
            s += "a : ndarray\n    First input matrix.\n";
            s += "b : ndarray\n    Second input matrix.\n";
            s += "beta : complex, optional\n    Scalar multiplier for `c`. Only its real part is used. Default is 0.\n";
            s += "c : ndarray, optional\n    Hermitian matrix to update; if omitted a new zero matrix is allocated. Default is None.\n";
            s += P_TRANS;
            s += P_LOWER;
            s += P_OVERWRITE_C;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The updated matrix.\n";
            return s;
        }

        /* trmm (product) / trsm (solve): triangular matrix-matrix operations. */
        static std::string
        doc_trmm_trsm(const char *name, const Dtype &t, const char *action, const char *result) noexcept
        {
            std::string s;
            s += std::string(name) + "(alpha, a, b, side=0, lower=0, trans_a=0, diag=0, overwrite_b=0)\n\n";
            s += std::string(action) + " (BLAS ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "alpha : " + std::string(t.scalar) + "\n    Scalar multiplier.\n";
            s += "a : ndarray\n    Triangular matrix; only the triangle selected by `lower` is referenced.\n";
            s += "b : ndarray\n    Input/output matrix.\n";
            s += P_SIDE;
            s += P_LOWER;
            s += P_TRANS_A;
            s += P_DIAG;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += std::string("b : ndarray\n    ") + result + "\n";
            return s;
        }
        static std::string doc_trmm(const char *name, const Dtype &t) { return doc_trmm_trsm(name, t, "Compute the triangular matrix-matrix product ``b = alpha*op(a)@b``", "The matrix ``alpha*op(a)@b`` (or ``alpha*b@op(a)`` if ``side=1``)."); }
        static std::string doc_trsm(const char *name, const Dtype &t) { return doc_trmm_trsm(name, t, "Solve the triangular system ``op(a)@x = alpha*b`` for `x`", "The solution `x` (or the solution of ``x@op(a) = alpha*b`` if ``side=1``)."); }

        typedef std::string (*DocFn)(const char *name, const Dtype &t);

        struct DocEntry {
            const char *name;
            DocFn fn;
            Dtype t;
        };

        /** @brief Expand a regular family (s/d/c/z prefixes) to its four docstring-table rows. */
        #define DOC_FAMILY(fam)     \
            {"s" #fam, doc_##fam, S}, \
            {"d" #fam, doc_##fam, D}, \
            {"c" #fam, doc_##fam, C}, \
            {"z" #fam, doc_##fam, Z}

        /** @brief Same, for families whose four names are irregular (`scnrm2`, `isamax`, ...). */
        #define DOC_FAMILY2(sname, dname, cname, zname, builder) \
            {#sname, builder, S}, \
            {#dname, builder, D}, \
            {#cname, builder, C}, \
            {#zname, builder, Z}

        /** @brief A real-only family (`s`/`d` prefixes). */
        #define DOC_FAMILY_R(fam) {"s" #fam, doc_##fam, S}, {"d" #fam, doc_##fam, D}

        /** @brief A complex-only family (`c`/`z` prefixes). */
        #define DOC_FAMILY_C(fam) {"c" #fam, doc_##fam, C}, {"z" #fam, doc_##fam, Z}

        static const DocEntry doc_table[] = {
            DOC_FAMILY(axpy),
            DOC_FAMILY(copy),
            DOC_FAMILY(swap),
            DOC_FAMILY(scal),
            DOC_FAMILY(rotg),
            {"csscal", doc_scal_real, C},
            {"zdscal", doc_scal_real, Z},
            DOC_FAMILY_R(dot),
            DOC_FAMILY_C(dotu),
            DOC_FAMILY_C(dotc),
            DOC_FAMILY2(srot, drot, csrot, zdrot, doc_rot),
            DOC_FAMILY_R(rotm),
            DOC_FAMILY_R(rotmg),
            DOC_FAMILY2(snrm2,  dnrm2,  scnrm2, dznrm2, doc_nrm2),
            DOC_FAMILY2(sasum,  dasum,  scasum, dzasum, doc_asum),
            DOC_FAMILY2(isamax, idamax, icamax, izamax, doc_iamax),

            /* Level 2 */
            DOC_FAMILY(gemv),
            DOC_FAMILY(gbmv),
            DOC_FAMILY_R(sbmv),
            DOC_FAMILY_C(hbmv),
            DOC_FAMILY(symv),
            DOC_FAMILY_C(hemv),
            DOC_FAMILY(spmv),
            DOC_FAMILY_C(hpmv),
            DOC_FAMILY(spr),
            DOC_FAMILY_C(hpr),
            DOC_FAMILY_R(spr2),
            DOC_FAMILY_C(hpr2),
            DOC_FAMILY(syr),
            DOC_FAMILY_C(her),
            DOC_FAMILY_R(syr2),
            DOC_FAMILY_C(her2),
            DOC_FAMILY_R(ger),
            DOC_FAMILY_C(geru),
            DOC_FAMILY_C(gerc),
            DOC_FAMILY(trmv),
            DOC_FAMILY(trsv),
            DOC_FAMILY(tbmv),
            DOC_FAMILY(tbsv),
            DOC_FAMILY(tpmv),
            DOC_FAMILY(tpsv),

            /* Level 3 */
            DOC_FAMILY(gemm),
            DOC_FAMILY(symm),
            DOC_FAMILY_C(hemm),
            DOC_FAMILY(syrk),
            DOC_FAMILY_C(herk),
            DOC_FAMILY(syr2k),
            DOC_FAMILY_C(her2k),
            DOC_FAMILY(trmm),
            DOC_FAMILY(trsm),
        };

        /**
         * @brief Build the docstring for `name`, or return nullptr if none is registered.
         *
         * Returns a new reference to a `str`.  A nullptr return with no exception set means the
         * routine has no docstring yet (the caller surfaces `None`); nullptr with an exception set
         * is a real failure to propagate.
         */
        PyObject *
        build_doc(const char *name) noexcept
        {
            for (const DocEntry &e : doc_table) {
                if (std::strcmp(e.name, name) == 0) {
                    const std::string d = e.fn(e.name, e.t);
                    return PyUnicode_FromStringAndSize(d.data(), static_cast<Py_ssize_t>(d.size()));
                }
            }
            return nullptr;
        }

    } // namespace blas::capi
}  // namespace blas
