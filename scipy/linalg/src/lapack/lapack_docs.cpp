/**
 * @file
 * @brief Docstrings for the `_flapack_cpp` wrappers, built on demand.
 *
 * Each wrapper's docstring is assembled at runtime from a template parameterised by the element
 * type.  `lapack::capi::build_doc(name)` is called lazily by the `LapackFunc.__doc__` getter (see
 * `lapack_module.cpp`) the first time a routine's docstring is requested, and the result is
 * cached on the callable.  Nothing is built at import time, and only the shared templates -- not
 * a fully expanded string per routine -- live in the compiled module.
 *
 * This mirrors `blas_docs.cpp` deliberately, down to the `P_*` fragment constants and the
 * `doc_<family>` / `doc_table` / `DOC_FAMILY` shape, so the two files are one idiom.
 *
 * Adding a routine: give its family a `doc_<family>(name, Dtype)` builder and list its typed
 * names in `doc_table` with `DOC_FAMILY(fam)`.  Reuse the `P_*` / `R_*` constants wherever the
 * prose is genuinely the same; write it out where it is not.
 */
#include <Python.h>
#include <cstring>
#include <string>

namespace lapack {
    namespace capi {

        /**
         * @brief The element-type-dependent facts a docstring template needs.
         *
         * `scalar` is the numpydoc type of a value that follows the routine's flavor;
         * `is_complex` selects the wording where the real and complex routines genuinely differ
         * (`gels`'s `trans` letter).  Single- vs double-precision never appears in the prose --
         * the routine name carries it -- so `S`/`D` and `C`/`Z` coincide.
         */
        struct Dtype {
            const char *scalar;
            bool is_complex;
        };

        static constexpr Dtype S{"float", false};
        static constexpr Dtype D{"float", false};
        static constexpr Dtype C{"complex", true};
        static constexpr Dtype Z{"complex", true};

        /* Complete numpydoc `Parameters` entries for the arguments that recur across families,
         * written once so every routine documents them identically.  Drop one into a docstring
         * with `s += P_N;`. */
        static constexpr const char *P_M        = "m : int\n    Number of rows of the matrix.\n";
        static constexpr const char *P_N        = "n : int\n    Number of columns of the matrix.\n";
        static constexpr const char *P_N_ORDER  = "n : int\n    Order of the square matrix.\n";
        static constexpr const char *P_NRHS     = "nrhs : int\n    Number of right-hand sides.\n";
        static constexpr const char *P_COND     =
            "cond : float, optional\n"
            "    Singular values below ``cond`` times the largest are treated as zero.\n"
            "    Default is -1.0, which uses machine precision.\n";
        static constexpr const char *P_LWORK_Q  = "lwork : int, optional\n    Passed through to LAPACK; -1 requests the query. Default is -1.\n";
        static constexpr const char *P_COMPUTE_UV    = "compute_uv : int, optional\n    If nonzero, singular vectors are computed. Default is 1.\n";
        static constexpr const char *P_FULL_MATRICES =
            "full_matrices : int, optional\n"
            "    If nonzero, full-sized `u` and `vt` are computed. Default is 1.\n";
        static constexpr const char *P_COMPUTE_VL    = "compute_vl : int, optional\n    If nonzero, left eigenvectors are computed. Default is 1.\n";
        static constexpr const char *P_COMPUTE_VR    = "compute_vr : int, optional\n    If nonzero, right eigenvectors are computed. Default is 1.\n";
        static constexpr const char *P_LO       = "lo : int, optional\n    Index of the first row/column to reduce, 0-based. Default is 0.\n";
        static constexpr const char *P_HI       = "hi : int, optional\n    Index of the last row/column to reduce, 0-based. Default is ``n - 1``.\n";

        static constexpr const char *P_A_SQUARE  = "a : ndarray\n    Square matrix of shape ``(n, n)``.\n";
        static constexpr const char *P_A_GENERAL = "a : ndarray\n    Matrix of shape ``(m, n)``.\n";
        static constexpr const char *P_B_RHS     = "b : ndarray\n    Right-hand side(s) of shape ``(n, nrhs)``.\n";
        static constexpr const char *P_LU_IN     = "lu : ndarray\n    Combined ``L`` and ``U`` factors, as returned by ``getrf``.\n";
        static constexpr const char *P_PIV_IN    = "piv : ndarray\n    Pivot indices from ``getrf``, 0-based.\n";
        static constexpr const char *P_TRANS_INT =
            "trans : int, optional\n"
            "    System to solve: 0 for ``a @ x = b``, 1 for ``a.T @ x = b``,\n"
            "    2 for ``a.conj().T @ x = b``. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_A   =
            "overwrite_a : int, optional\n"
            "    If nonzero, `a` may be overwritten in place. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_B   =
            "overwrite_b : int, optional\n"
            "    If nonzero, `b` may be overwritten in place. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_LU  =
            "overwrite_lu : int, optional\n"
            "    If nonzero, `lu` may be overwritten in place. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_RHS =
            "overwrite_rhs : int, optional\n"
            "    If nonzero, `rhs` may be overwritten in place. Default is 0.\n";

        /* The tridiagonal families.  A tridiagonal matrix arrives as its three diagonals, and
         * `gttrf` hands its factorization on to `gttrs`, `gtcon` and `gtsvx` as five arrays. */
        static constexpr const char *P_DL     = "dl : ndarray\n    Subdiagonal elements, length ``n - 1``.\n";
        static constexpr const char *P_D_TRI  = "d : ndarray\n    Diagonal elements, length ``n``.\n";
        static constexpr const char *P_DU     = "du : ndarray\n    Superdiagonal elements, length ``n - 1``.\n";
        static constexpr const char *P_GT_FACTORS =
            "dl : ndarray\n    Multipliers defining ``L``, as returned by ``gttrf``.\n"
            "d : ndarray\n    Diagonal of ``U``, as returned by ``gttrf``.\n"
            "du : ndarray\n    First superdiagonal of ``U``, as returned by ``gttrf``.\n"
            "du2 : ndarray\n    Second superdiagonal of ``U``, as returned by ``gttrf``.\n"
            "ipiv : ndarray\n    Pivot indices from ``gttrf``, 1-based.\n";
        static constexpr const char *P_TRANS_STR =
            "trans : str, optional\n"
            "    ``'N'``, ``'T'`` or ``'C'`` for the system, its transpose, or its conjugate\n"
            "    transpose. Default is ``'N'``.\n";

        /* The symmetric tridiagonal eigensolvers. */
        static constexpr const char *P_D_SYM =
            "d : ndarray\n    Diagonal elements of the symmetric tridiagonal matrix, length ``n``.\n";
        static constexpr const char *P_E_SYM = "e : ndarray\n    Off-diagonal elements, length ``n - 1``.\n";
        static constexpr const char *P_COMPUTE_V =
            "compute_v : int, optional\n    If nonzero, eigenvectors are computed. Default is 1.\n";
        static constexpr const char *P_RANGE_SELECT =
            "range : int\n"
            "    Which eigenvalues to compute: 0 for all, 1 for those in ``(vl, vu]``, 2 for\n"
            "    those with indices `il` through `iu`.\n"
            "vl : float\n    Lower bound of the half-open interval; used only when `range` is 1.\n"
            "vu : float\n    Upper bound of the half-open interval; used only when `range` is 1.\n"
            "il : int\n"
            "    Index of the smallest eigenvalue to return, 1-based; used only when\n"
            "    `range` is 2.\n"
            "iu : int\n"
            "    Index of the largest eigenvalue to return, 1-based; used only when\n"
            "    `range` is 2.\n";

        static constexpr const char *P_OVERWRITE_DL =
            "overwrite_dl : int, optional\n"
            "    If nonzero, `dl` may be overwritten in place. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_D  =
            "overwrite_d : int, optional\n"
            "    If nonzero, `d` may be overwritten in place. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_DU =
            "overwrite_du : int, optional\n"
            "    If nonzero, `du` may be overwritten in place. Default is 0.\n";
        static constexpr const char *P_OVERWRITE_E  =
            "overwrite_e : int, optional\n"
            "    If nonzero, `e` may be overwritten in place. Default is 0.\n";

        static constexpr const char *P_JPTV    =
            "jptv : ndarray\n"
            "    Column pivots: a nonzero entry pins that column of `a` to the front.\n"
            "    Modified in place and returned. Use zeros to let the routine choose freely.\n";

        /* Shared `Returns` entries. */
        static constexpr const char *R_INFO    = "info : int\n    0 on success; if negative, the ``-info``-th argument had an illegal value.\n";
        static constexpr const char *R_INFO_LU =
            "info : int\n"
            "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
            "    if positive, ``u[info-1, info-1]`` is exactly zero and the factor is singular.\n";
        static constexpr const char *R_INFO_CV =
            "info : int\n"
            "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
            "    if positive, the algorithm failed to converge.\n";
        static constexpr const char *R_LU_OUT  = "lu : ndarray\n    Combined ``L`` and ``U`` factors; the unit diagonal of ``L`` is not stored.\n";
        static constexpr const char *R_PIV_OUT = "piv : ndarray\n    Pivot indices, 0-based: row ``i`` was interchanged with row ``piv[i]``.\n";
        static constexpr const char *R_X_OUT   = "x : ndarray\n    Solution of the system.\n";
        static constexpr const char *R_WORK_OUT = "work : ndarray\n    Workspace actually used; ``work[0]`` is the optimal `lwork`.\n";
        static constexpr const char *R_TAU     = "tau : ndarray\n    Scalar factors of the elementary reflectors that define ``Q``.\n";
        static constexpr const char *R_RANK    = "rank : int\n    Effective rank of `a`, determined using `cond`.\n";
        static constexpr const char *R_S_SING  = "s : ndarray\n    Singular values of `a`, in descending order.\n";
        static constexpr const char *R_VL      = "vl : ndarray\n    Left eigenvectors, or a ``(1, n)`` placeholder when `compute_vl` is 0.\n";
        static constexpr const char *R_VR      = "vr : ndarray\n    Right eigenvectors, or a ``(1, n)`` placeholder when `compute_vr` is 0.\n";

        /** @brief The `lwork` entry, naming the companion that computes the optimal value. */
        static std::string
        p_lwork(const char *deflt, const char *query) noexcept
        {
            return "lwork : int, optional\n    Size of the workspace. Default is " + std::string(deflt)
                 + ".\n    Use ``" + std::string(query) + "`` for the optimal value.\n";
        }

        /** @brief `gees` and `gges` take their selector under a flavor-prefixed keyword --
         *         `sselect`, `dselect`, `cselect`, `zselect` -- so the signature is built from
         *         the wrapper's own first letter rather than a fixed name. */
        static std::string
        doc_gees(const char *name, const Dtype &t) noexcept
        {
            const std::string sel = std::string(1, name[0]) + "select";
            std::string s;
            s += std::string(name) + "(" + sel + ", a, compute_v=1, sort_t=0, lwork=3*n,\n";
            s += std::string(std::strlen(name) + 1, ' ') + sel + "_extra_args=(), overwrite_a=0)\n\n";
            s += "Compute the Schur factorization of a general matrix (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            if (t.is_complex) {
                s += sel + " : callable\n"
                     "    Eigenvalue selector, called as ``" + sel + "(w)``. Used only when `sort_t` is\n"
                     "    nonzero; return true to move that eigenvalue to the leading block.\n";
            }
            else {
                s += sel + " : callable\n"
                     "    Eigenvalue selector, called as ``" + sel + "(wr, wi)`` with the real and\n"
                     "    imaginary parts. Used only when `sort_t` is nonzero; return true to move that\n"
                     "    eigenvalue to the leading block.\n";
            }
            s += P_A_SQUARE;
            s += "compute_v : int, optional\n    If nonzero, the Schur vectors are computed. Default is 1.\n";
            s += "sort_t : int, optional\n    If nonzero, eigenvalues are sorted using the selector. Default is 0.\n";
            s += p_lwork("``3 * n``", "the value reported in work[0]");
            s += sel + "_extra_args : tuple, optional\n    Extra positional arguments appended to every selector call. Default is ``()``.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "t : ndarray\n    Schur form of `a`.\n";
            s += "sdim : int\n    Number of eigenvalues the selector accepted; 0 when `sort_t` is 0.\n";
            if (t.is_complex) {
                s += "w : ndarray\n    Eigenvalues, in the order they appear on the diagonal of `t`.\n";
            }
            else {
                s += "wr : ndarray\n    Real parts of the eigenvalues.\n";
                s += "wi : ndarray\n    Imaginary parts of the eigenvalues; conjugate pairs appear adjacently.\n";
            }
            s += "vs : ndarray\n    Schur vectors, or a ``(1, n)`` placeholder when `compute_v` is 0.\n";
            s += R_WORK_OUT;
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, the QR algorithm failed or the selected block could not be\n"
                 "    reordered.\n";
            return s;
        }

        static std::string
        doc_gges(const char *name, const Dtype &t) noexcept
        {
            const std::string sel = std::string(1, name[0]) + "select";
            std::string s;
            s += std::string(name) + "(" + sel + ", a, b, jobvsl=1, jobvsr=1, sort_t=0, ldvsl=n, ldvsr=n,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "lwork=8*n+16, " + sel + "_extra_args=(), overwrite_a=0, overwrite_b=0)\n\n";
            s += "Compute the generalized Schur factorization of a matrix pair (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            if (t.is_complex) {
                s += sel + " : callable\n"
                     "    Generalized-eigenvalue selector, called as ``" + sel + "(alpha, beta)``. Used\n"
                     "    only when `sort_t` is nonzero.\n";
            }
            else {
                s += sel + " : callable\n"
                     "    Generalized-eigenvalue selector, called as ``" + sel + "(alphar, alphai, beta)``.\n"
                     "    Used only when `sort_t` is nonzero.\n";
            }
            s += P_A_SQUARE;
            s += "b : ndarray\n    Square matrix of shape ``(n, n)``, the second member of the pair.\n";
            s += "jobvsl : int, optional\n    If nonzero, the left Schur vectors are computed. Default is 1.\n";
            s += "jobvsr : int, optional\n    If nonzero, the right Schur vectors are computed. Default is 1.\n";
            s += "sort_t : int, optional\n    If nonzero, eigenvalues are sorted using the selector. Default is 0.\n";
            s += "ldvsl : int, optional\n    Leading dimension of `vsl`. Default is ``n`` when `jobvsl` is nonzero, else 1.\n";
            s += "ldvsr : int, optional\n    Leading dimension of `vsr`. Default is ``n`` when `jobvsr` is nonzero, else 1.\n";
            s += p_lwork(t.is_complex ? "``2 * n``" : "``8 * n + 16``", "the value reported in work[0]");
            s += sel + "_extra_args : tuple, optional\n    Extra positional arguments appended to every selector call. Default is ``()``.\n";
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    Generalized Schur form of `a`.\n";
            s += "b : ndarray\n    Generalized Schur form of `b`.\n";
            s += "sdim : int\n    Number of eigenvalues the selector accepted; 0 when `sort_t` is 0.\n";
            if (t.is_complex) {
                s += "alpha : ndarray\n    Numerators of the generalized eigenvalues ``alpha / beta``.\n";
            }
            else {
                s += "alphar : ndarray\n    Real parts of the numerators of ``alpha / beta``.\n";
                s += "alphai : ndarray\n    Imaginary parts of the numerators of ``alpha / beta``.\n";
            }
            s += "beta : ndarray\n    Denominators of the generalized eigenvalues.\n";
            s += "vsl : ndarray\n    Left Schur vectors, or a placeholder when `jobvsl` is 0.\n";
            s += "vsr : ndarray\n    Right Schur vectors, or a placeholder when `jobvsr` is 0.\n";
            s += R_WORK_OUT;
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, the QZ algorithm failed or the selected block could not be\n"
                 "    reordered.\n";
            return s;
        }

        static std::string
        doc_gebal(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, scale=0, permute=0, overwrite_a=0)\n\n";
            s += "Balance a general matrix to improve eigenvalue accuracy (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_SQUARE;
            s += "scale : int, optional\n    If nonzero, rows and columns are scaled to make their norms comparable.\n    Default is 0.\n";
            s += "permute : int, optional\n"
                 "    If nonzero, the matrix is permuted to isolate eigenvalues already in the\n"
                 "    corners. Default is 0.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "ba : ndarray\n    The balanced matrix.\n";
            s += "lo : int\n    First index of the balanced middle block, 0-based.\n";
            s += "hi : int\n    Last index of the balanced middle block, 0-based.\n";
            s += "pivscale : ndarray\n    Permutation and scaling applied, packed as LAPACK returns them.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gehrd(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, lo=0, hi=n-1, lwork=max(n,1), overwrite_a=0)\n\n";
            s += "Reduce a matrix to upper Hessenberg form (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_SQUARE;
            s += P_LO;
            s += P_HI;
            s += p_lwork("``max(n, 1)``", "gehrd_lwork");
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "ht : ndarray\n"
                 "    Upper Hessenberg form in the upper triangle and first subdiagonal; below\n"
                 "    that, the elementary reflectors that define ``Q``.\n";
            s += R_TAU;
            s += R_INFO;
            return s;
        }


        /* Builders follow `gen_methods` order (see `lapack_gen.cpp`), which keeps every
         * `_lwork` query immediately after the routine it queries. */

        static std::string
        doc_gehrd_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, lo=0, hi=n-1)\n\n";
            s += "Query the optimal `lwork` for ``" + std::string(1, name[0]) + "gehrd``.\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_LO;
            s += P_HI;

            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gesv(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve ``a @ x = b`` by LU factorization with partial pivoting (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_SQUARE;
            s += P_B_RHS;
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_LU_OUT;
            s += R_PIV_OUT;
            s += R_X_OUT;
            s += R_INFO_LU;
            return s;
        }

        static std::string
        doc_gecon(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, anorm, norm='1')\n\n";
            s += "Estimate the reciprocal condition number of a factorized matrix (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    LU factorization of a square matrix, as returned by ``getrf``.\n";
            s += "anorm : float\n    Norm of the original matrix, in the norm selected by `norm`.\n";
            s += "norm : str, optional\n    ``'1'`` or ``'O'`` for the 1-norm, ``'I'`` for the infinity norm.\n    Default is ``'1'``.\n";

            s += "\nReturns\n-------\n";
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_getrf(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, overwrite_a=0)\n\n";
            s += "Compute the LU factorization ``a = p @ l @ u`` with partial pivoting (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_GENERAL;
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += R_LU_OUT;
            s += R_PIV_OUT;
            s += R_INFO_LU;
            return s;
        }

        static std::string
        doc_getrs(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(lu, piv, b, trans=0, overwrite_b=0)\n\n";
            s += "Solve a system already factorized by ``getrf`` (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_LU_IN;
            s += P_PIV_IN;
            s += P_B_RHS;
            s += P_TRANS_INT;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_X_OUT;
            s += R_INFO;
            return s;
        }

        static std::string
        doc_getc2(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, overwrite_a=0)\n\n";
            s += "Compute the LU factorization with complete pivoting (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_SQUARE;
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += R_LU_OUT;
            s += "ipiv : ndarray\n    Row pivot indices, 0-based.\n";
            s += "jpiv : ndarray\n    Column pivot indices, 0-based.\n";
            s += "info : int\n"
                 "    0 on success; if positive, ``u[info-1, info-1]`` is exactly zero. The factor is\n"
                 "    perturbed rather than left singular, so a following ``gesc2`` still solves.\n";
            return s;
        }

        static std::string
        doc_gesc2(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(lu, rhs, ipiv, jpiv, overwrite_rhs=0)\n\n";
            s += "Solve a system factorized by ``getc2``, scaling to avoid overflow (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_LU_IN;
            s += "rhs : ndarray\n    Right-hand side vector of length ``n``.\n";
            s += "ipiv : ndarray\n    Row pivot indices from ``getc2``, 0-based.\n";
            s += "jpiv : ndarray\n    Column pivot indices from ``getc2``, 0-based.\n";
            s += P_OVERWRITE_RHS;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    Solution, scaled by `scale` to avoid overflow.\n";
            s += "scale : float\n    Scale factor applied to the right-hand side, ``0 <= scale <= 1``.\n";
            return s;
        }

        static std::string
        doc_getri(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(lu, piv, lwork=3*n, overwrite_lu=0)\n\n";
            s += "Invert a matrix from its ``getrf`` factorization (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_LU_IN;
            s += P_PIV_IN;
            s += "lwork : int, optional\n"
                 "    Size of the workspace; must be at least ``n``. Default is ``3 * n``.\n"
                 "    Use ``getri_lwork`` for the optimal value.\n";
            s += P_OVERWRITE_LU;

            s += "\nReturns\n-------\n";
            s += "inv_a : ndarray\n    Inverse of the matrix `lu` was factorized from.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, ``u[info-1, info-1]`` is exactly zero and the matrix is singular.\n";
            return s;
        }

        static std::string
        doc_getri_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(n)\n\n";
            s += "Query the optimal `lwork` for ``" + std::string(1, name[0]) + "getri``.\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;

            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        /** @brief `gesdd` and `gesvd` present the same interface; they differ in algorithm --
         *         divide-and-conquer versus the QR iteration. */
        static std::string
        doc_gesxd_family(const char *name, const char *how, const char *query) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, compute_uv=1, full_matrices=1, lwork=..., overwrite_a=0)\n\n";
            s += "Compute the singular value decomposition ``a = u @ diag(s) @ vt`` " + std::string(how) + "\n"
                 "(LAPACK ``" + std::string(name) + "``).\n"
                 "\n";

            s += "Parameters\n----------\n";
            s += P_A_GENERAL;
            s += P_COMPUTE_UV;
            s += P_FULL_MATRICES;
            s += p_lwork("routine-specific", query);
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "u : ndarray\n    Left singular vectors; a ``(1, 1)`` placeholder when `compute_uv` is 0.\n";
            s += R_S_SING;
            s += "vt : ndarray\n    Right singular vectors, transposed; a placeholder when `compute_uv` is 0.\n";
            s += R_INFO_CV;
            return s;
        }

        static std::string doc_gesdd(const char *name, const Dtype &) noexcept
            { return doc_gesxd_family(name, "by divide-and-conquer", "gesdd_lwork"); }
        static std::string doc_gesvd(const char *name, const Dtype &) noexcept { return doc_gesxd_family(name, "by QR iteration", "gesvd_lwork"); }

        /** @brief `gesdd_lwork` and `gesvd_lwork` take the same arguments. */
        static std::string
        doc_gesdd_family_lwork(const char *name, const Dtype &t, const char *base) noexcept
        {
            std::string s;
            s += std::string(name) + "(m, n, compute_uv=1, full_matrices=1)\n\n";
            s += "Query the optimal `lwork` for ``" + std::string(1, name[0]) + base + "``.\n\n";

            s += "Parameters\n----------\n";
            s += P_M;
            s += P_N;
            s += P_COMPUTE_UV;
            s += P_FULL_MATRICES;

            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        static std::string doc_gesdd_lwork(const char *name, const Dtype &t) noexcept { return doc_gesdd_family_lwork(name, t, "gesdd"); }
        static std::string doc_gesvd_lwork(const char *name, const Dtype &t) noexcept { return doc_gesdd_family_lwork(name, t, "gesvd"); }

        static std::string
        doc_gels(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, trans='N', lwork=..., overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve an overdetermined or underdetermined system by QR or LQ factorization\n(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_GENERAL;
            s += "b : ndarray\n    Right-hand side(s) with ``max(m, n)`` rows.\n";
            if (t.is_complex) {
                s += "trans : str, optional\n    ``'N'`` for the system as given, ``'C'`` for its conjugate transpose.\n    Default is ``'N'``.\n";
            }
            else {
                s += "trans : str, optional\n    ``'N'`` for the system as given, ``'T'`` for its transpose.\n    Default is ``'N'``.\n";
            }
            s += p_lwork("``min(m, n) + max(min(m, n), nrhs)``", "gels_lwork");
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "lqr : ndarray\n    The QR or LQ factorization of `a`, as LAPACK leaves it.\n";
            s += "x : ndarray\n"
                 "    Solution in the first ``n`` rows; the remaining rows hold the residual\n"
                 "    sum of squares for each column when ``m > n``.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, `a` is rank-deficient and no least-squares solution was computed.\n";
            return s;
        }

        static std::string
        doc_gels_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(m, n, nrhs, trans='N')\n\n";
            s += "Query the optimal `lwork` for ``" + std::string(1, name[0]) + "gels``.\n\n";

            s += "Parameters\n----------\n";
            s += P_M;
            s += P_N;
            s += P_NRHS;
            /* The transpose letter is a real difference: the complex routines take 'C'. */
            if (t.is_complex) {
                s += "trans : str, optional\n    ``'N'`` for the system as given, ``'C'`` for its conjugate transpose.\n    Default is ``'N'``.\n";
            }
            else {
                s += "trans : str, optional\n    ``'N'`` for the system as given, ``'T'`` for its transpose.\n    Default is ``'N'``.\n";
            }

            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gelss(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, cond=-1.0, lwork=..., overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve a least-squares problem by singular value decomposition (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_GENERAL;
            s += "b : ndarray\n    Right-hand side(s) with ``max(m, n)`` rows.\n";
            s += P_COND;
            s += p_lwork("routine-specific", "gelss_lwork");
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "v : ndarray\n    The first ``min(m, n)`` right singular vectors of `a`, as LAPACK leaves them.\n";
            s += "x : ndarray\n    Minimum-norm solution in the first ``n`` rows.\n";
            s += R_S_SING;
            s += R_RANK;
            s += R_WORK_OUT;
            s += R_INFO_CV;
            return s;
        }

        static std::string
        doc_gelsy(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, jptv, cond, lwork, overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve a least-squares problem by complete orthogonal factorization\n(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_GENERAL;
            s += "b : ndarray\n    Right-hand side(s) with ``max(m, n)`` rows.\n";
            s += P_JPTV;
            s += "cond : float\n"
                 "    Columns are treated as dependent when the estimated condition number of the\n"
                 "    leading block would exceed ``1 / cond``.\n";
            s += "lwork : int\n    Size of the workspace. Use ``gelsy_lwork`` for the optimal value.\n";
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "v : ndarray\n    The complete orthogonal factorization of `a`, as LAPACK leaves it.\n";
            s += "x : ndarray\n    Minimum-norm solution in the first ``n`` rows.\n";
            s += "jptv : ndarray\n    Column pivots actually used, 1-based as LAPACK numbers them.\n";
            s += R_RANK;
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gelsd(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            if (t.is_complex) {
                s += std::string(name) + "(a, b, lwork, size_rwork, size_iwork, cond=-1.0,\n";
                s += std::string(std::strlen(name) + 1, ' ') + "overwrite_a=0, overwrite_b=0)\n\n";
            }
            else {
                s += std::string(name) + "(a, b, lwork, size_iwork, cond=-1.0, overwrite_a=0, overwrite_b=0)\n\n";
            }
            s += "Solve a least-squares problem by divide-and-conquer SVD (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_GENERAL;
            s += "b : ndarray\n    Right-hand side(s) with ``max(m, n)`` rows.\n";
            s += "lwork : int\n    Size of the workspace. Use ``gelsd_lwork`` for the required value.\n";
            if (t.is_complex) {
                s += "size_rwork : int\n    Size of the real workspace, as reported by ``gelsd_lwork``.\n";
            }
            s += "size_iwork : int\n    Size of the integer workspace, as reported by ``gelsd_lwork``.\n";
            s += P_COND;
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    Minimum-norm solution in the first ``n`` rows.\n";
            s += R_S_SING;
            s += R_RANK;
            s += R_INFO_CV;
            return s;
        }

        /** @brief `gelss_lwork`, `gelsy_lwork` and `gelsd_lwork` share their argument list;
         *         only `gelsy` makes `cond` required. */
        static std::string
        doc_gelsx_lwork(const char *name, const Dtype &t, const char *base, const char *args) noexcept
        {
            std::string s;
            s += std::string(name) + args + "\n\n";
            s += "Query the optimal `lwork` for ``" + std::string(1, name[0]) + base + "``.\n\n";

            s += "Parameters\n----------\n";
            s += P_M;
            s += P_N;
            s += P_NRHS;
            s += P_COND;
            s += P_LWORK_Q;

            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        static std::string doc_gelss_lwork(const char *name, const Dtype &t) noexcept
            { return doc_gelsx_lwork(name, t, "gelss", "(m, n, nrhs, cond=-1.0, lwork=-1)"); }
        static std::string doc_gelsy_lwork(const char *name, const Dtype &t) noexcept
            { return doc_gelsx_lwork(name, t, "gelsy", "(m, n, nrhs, cond, lwork=-1)"); }
        static std::string doc_gelsd_lwork(const char *name, const Dtype &t) noexcept
            { return doc_gelsx_lwork(name, t, "gelsd", "(m, n, nrhs, cond=-1.0, lwork=-1)"); }

        static std::string
        doc_geqp3(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, lwork=3*(n+1), overwrite_a=0)\n\n";
            s += "Compute a QR factorization with column pivoting (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_GENERAL;
            s += p_lwork("``3 * (n + 1)``", "the value reported in work[0]");
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "qr : ndarray\n    ``R`` in the upper triangle; below it, the reflectors that define ``Q``.\n";
            s += "jpvt : ndarray\n"
                 "    Column pivots, 1-based as LAPACK numbers them: column ``jpvt[i] - 1`` of `a`\n"
                 "    became column ``i`` of ``a @ p``.\n";
            s += R_TAU;
            s += R_WORK_OUT;
            s += R_INFO;
            return s;
        }

        /** @brief `geqrf`, `geqrfp` and `gerqf` are the unpivoted factorizations; they differ in
         *         the sign convention on ``R`` and in which side ``Q`` acts from. */
        static std::string
        doc_geqrf_family(const char *name, const char *summary, const char *out, const char *deflt,
                         const char *query, bool returns_work) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, lwork=" + std::string(deflt) + ", overwrite_a=0)\n\n";
            s += std::string(summary) + " (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_GENERAL;
            s += p_lwork(deflt, query);
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += out;
            s += R_TAU;
            if (returns_work) { s += R_WORK_OUT; }
            s += R_INFO;
            return s;
        }

        static std::string doc_geqrf(const char *name, const Dtype &) noexcept
        {
            return doc_geqrf_family(name, "Compute a QR factorization",
                "qr : ndarray\n    ``R`` in the upper triangle; below it, the reflectors that define ``Q``.\n",
                "``3 * n``", "geqrf_lwork", true);
        }

        static std::string doc_geqrfp(const char *name, const Dtype &) noexcept
        {
            return doc_geqrf_family(name, "Compute a QR factorization whose ``R`` has a nonnegative diagonal",
                "qr : ndarray\n    ``R`` in the upper triangle, with ``diag(R) >= 0``; below it, the reflectors\n    that define ``Q``.\n",
                "``max(1, n)``", "geqrfp_lwork", false);
        }

        static std::string doc_gerqf(const char *name, const Dtype &) noexcept
        {
            return doc_geqrf_family(name, "Compute an RQ factorization ``a = r @ q``",
                "qr : ndarray\n    ``R`` on and above the ``(m - n)``-th subdiagonal; elsewhere, the reflectors\n    that define ``Q``.\n",
                "``3 * m``", "the value reported in work[0]", true);
        }

        /** @brief `geqrf_lwork` and `geqrfp_lwork` differ only in which routine they query. */
        static std::string
        doc_geqrf_family_lwork(const char *name, const Dtype &t, const char *base) noexcept
        {
            std::string s;
            s += std::string(name) + "(m, n)\n\n";
            s += "Query the optimal `lwork` for ``" + std::string(1, name[0]) + base + "``.\n\n";

            s += "Parameters\n----------\n";
            s += P_M;
            s += P_N;

            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        static std::string doc_geqrf_lwork(const char *name, const Dtype &t) noexcept  { return doc_geqrf_family_lwork(name, t, "geqrf"); }
        static std::string doc_geqrfp_lwork(const char *name, const Dtype &t) noexcept { return doc_geqrf_family_lwork(name, t, "geqrfp"); }

        /** @brief `geequ` and `geequb` differ only in that `geequb` rounds its factors. */
        static std::string
        doc_geequ_family(const char *name, const Dtype &, const char *what) noexcept
        {
            std::string s;
            s += std::string(name) + "(a)\n\n";
            s += std::string(what) + " (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_GENERAL;

            s += "\nReturns\n-------\n";
            s += "r : ndarray\n    Row scale factors, length ``m``.\n";
            s += "c : ndarray\n    Column scale factors, length ``n``.\n";
            s += "rowcnd : float\n    Ratio of the smallest to the largest row scale factor.\n";
            s += "colcnd : float\n    Ratio of the smallest to the largest column scale factor.\n";
            s += "amax : float\n    Largest element of `a` in absolute value.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if ``0 < info <= m``, row ``info-1`` is exactly zero; if ``info > m``, column\n"
                 "    ``info-m-1`` is exactly zero.\n";
            return s;
        }

        static std::string doc_geequ(const char *name, const Dtype &t) noexcept
            { return doc_geequ_family(name, t, "Compute row and column scale factors that equilibrate a matrix"); }
        static std::string doc_geequb(const char *name, const Dtype &t) noexcept
            { return doc_geequ_family(name, t, "Compute equilibration factors restricted to powers of the radix"); }

        static std::string
        doc_geev(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, compute_vl=1, compute_vr=1, lwork=..., overwrite_a=0)\n\n";
            s += "Compute the eigenvalues and eigenvectors of a general matrix (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_SQUARE;
            s += P_COMPUTE_VL;
            s += P_COMPUTE_VR;
            s += p_lwork(t.is_complex ? "``2 * n``" : "``4 * n``", "geev_lwork");
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            if (t.is_complex) {
                s += "w : ndarray\n    Eigenvalues, in no particular order.\n";
            }
            else {
                s += "wr : ndarray\n    Real parts of the eigenvalues.\n";
                s += "wi : ndarray\n"
                     "    Imaginary parts; a conjugate pair occupies two adjacent entries, and its\n"
                     "    eigenvectors occupy the matching pair of columns as real and imaginary parts.\n";
            }
            s += R_VL;
            s += R_VR;
            s += R_INFO_CV;
            return s;
        }

        static std::string
        doc_geev_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, compute_vl=1, compute_vr=1)\n\n";
            s += "Query the optimal `lwork` for ``" + std::string(1, name[0]) + "geev``.\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_COMPUTE_VL;
            s += P_COMPUTE_VR;

            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_ggev(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, compute_vl=1, compute_vr=1, lwork=..., overwrite_a=0,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_b=0)\n\n";
            s += "Solve the generalized eigenproblem ``a @ v = (alpha / beta) * b @ v``\n(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_SQUARE;
            s += "b : ndarray\n    Square matrix of shape ``(n, n)``, the second member of the pair.\n";
            s += P_COMPUTE_VL;
            s += P_COMPUTE_VR;
            s += p_lwork(t.is_complex ? "``2 * n``" : "``8 * n``", "the value reported in work[0]");
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            if (t.is_complex) {
                s += "alpha : ndarray\n    Numerators of the generalized eigenvalues ``alpha / beta``.\n";
            }
            else {
                s += "alphar : ndarray\n    Real parts of the numerators of ``alpha / beta``.\n";
                s += "alphai : ndarray\n    Imaginary parts of the numerators of ``alpha / beta``.\n";
            }
            s += "beta : ndarray\n    Denominators; an eigenvalue is infinite where ``beta`` is zero.\n";
            s += R_VL;
            s += R_VR;
            s += R_WORK_OUT;
            s += R_INFO_CV;
            return s;
        }

        static std::string
        doc_gesvx(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, fact='E', trans='N', af=None, ipiv=None, equed='B', r=None, c=None,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve ``a @ x = b`` with equilibration, condition estimation and error bounds\n(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_SQUARE;
            s += P_B_RHS;
            s += "fact : str, optional\n"
                 "    ``'E'`` to equilibrate then factorize, ``'N'`` to factorize as given, ``'F'``\n"
                 "    to reuse the `af`, `ipiv`, `equed`, `r` and `c` supplied. Default is ``'E'``.\n";
            s += "trans : str, optional\n"
                 "    ``'N'``, ``'T'`` or ``'C'`` for the system, its transpose, or its conjugate\n"
                 "    transpose. Default is ``'N'``.\n";
            s += "af : ndarray, optional\n    Factorization to reuse when ``fact='F'``; otherwise it is computed.\n";
            s += "ipiv : ndarray, optional\n    Pivot indices to reuse when ``fact='F'``, 0-based; otherwise computed.\n";
            s += "equed : str, optional\n"
                 "    Equilibration already applied when ``fact='F'``: ``'N'``, ``'R'``, ``'C'`` or\n"
                 "    ``'B'``. Otherwise it is an output. Default is ``'B'``.\n";
            s += "r : ndarray, optional\n    Row scale factors, used when ``fact='F'`` and `equed` is ``'R'`` or ``'B'``.\n";
            s += "c : ndarray, optional\n    Column scale factors, used when ``fact='F'`` and `equed` is ``'C'`` or ``'B'``.\n";
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "as_ : ndarray\n    `a`, equilibrated if `equed` is not ``b'N'``.\n";
            s += "lu : ndarray\n    LU factorization of the equilibrated matrix.\n";
            s += R_PIV_OUT;
            s += "equed : bytes\n    Equilibration actually applied: ``b'N'``, ``b'R'``, ``b'C'`` or ``b'B'``.\n";
            s += "rs : ndarray\n    Row scale factors.\n";
            s += "cs : ndarray\n    Column scale factors.\n";
            s += "bs : ndarray\n    `b`, scaled to match the equilibrated system.\n";
            s += R_X_OUT;
            s += "rcond : float\n    Estimate of the reciprocal condition number of the equilibrated matrix.\n";
            s += "ferr : ndarray\n    Estimated forward error bound for each solution vector.\n";
            s += "berr : ndarray\n    Componentwise relative backward error of each solution vector.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if ``0 < info <= n``, ``u[info-1, info-1]`` is exactly zero; if ``info = n+1``,\n"
                 "    the matrix is singular to working precision and `x` may be inaccurate.\n";
            return s;
        }

        /* ============================ flapack_gen_tri.pyf.src ====================== */

        static std::string
        doc_gtsv(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(dl, d, du, b, overwrite_dl=0, overwrite_d=0, overwrite_du=0,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_b=0)\n\n";
            s += "Solve ``a @ x = b`` for a tridiagonal ``a`` by Gaussian elimination with partial\n"
                 "pivoting (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_DL;
            s += P_D_TRI;
            s += P_DU;
            s += P_B_RHS;
            s += P_OVERWRITE_DL;
            s += P_OVERWRITE_D;
            s += P_OVERWRITE_DU;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "du2 : ndarray\n"
                 "    Second superdiagonal of ``U`` in the leading ``n - 2`` entries, written over\n"
                 "    `dl`; the final entry is not part of the factorization.\n";
            s += "d : ndarray\n    Diagonal of ``U``.\n";
            s += "du : ndarray\n    First superdiagonal of ``U``.\n";
            s += R_X_OUT;
            s += R_INFO_LU;
            return s;
        }

        static std::string
        doc_gttrf(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(dl, d, du, overwrite_dl=0, overwrite_d=0, overwrite_du=0)\n\n";
            s += "Compute the LU factorization of a tridiagonal matrix with partial pivoting\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_DL;
            s += P_D_TRI;
            s += P_DU;
            s += P_OVERWRITE_DL;
            s += P_OVERWRITE_D;
            s += P_OVERWRITE_DU;

            s += "\nReturns\n-------\n";
            s += "dl : ndarray\n    Multipliers defining ``L``.\n";
            s += "d : ndarray\n    Diagonal of ``U``.\n";
            s += "du : ndarray\n    First superdiagonal of ``U``.\n";
            s += "du2 : ndarray\n    Second superdiagonal of ``U``, length ``n - 2``.\n";
            s += "ipiv : ndarray\n"
                 "    Pivot indices, 1-based -- unlike ``getrf``, which returns them 0-based.\n"
                 "    ``gttrs``, ``gtcon`` and ``gtsvx`` take them exactly as returned.\n";
            s += R_INFO_LU;
            return s;
        }

        static std::string
        doc_gttrs(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(dl, d, du, du2, ipiv, b, trans='N', overwrite_b=0)\n\n";
            s += "Solve a tridiagonal system using the factorization from ``gttrf``\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_GT_FACTORS;
            s += P_B_RHS;
            s += P_TRANS_STR;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_X_OUT;
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gtcon(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(dl, d, du, du2, ipiv, anorm, norm='1')\n\n";
            s += "Estimate the reciprocal condition number of a tridiagonal matrix from its\n"
                 "``gttrf`` factorization (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_GT_FACTORS;
            s += "anorm : float\n    Norm of the original matrix, in the norm selected by `norm`.\n";
            s += "norm : str, optional\n    ``'1'`` or ``'O'`` for the 1-norm, ``'I'`` for the infinity norm.\n    Default is ``'1'``.\n";

            s += "\nReturns\n-------\n";
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gtsvx(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(dl, d, du, b, fact='N', trans='N', dlf=None, df=None, duf=None,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "du2=None, ipiv=None)\n\n";
            s += "Solve a tridiagonal system with condition estimation and error bounds\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_DL;
            s += P_D_TRI;
            s += P_DU;
            s += P_B_RHS;
            s += "fact : str, optional\n"
                 "    ``'N'`` to factorize the matrix, ``'F'`` to reuse the `dlf`, `df`, `duf`,\n"
                 "    `du2` and `ipiv` supplied. Default is ``'N'``.\n";
            s += P_TRANS_STR;
            s += "dlf : ndarray, optional\n    Multipliers of ``L`` to reuse when ``fact='F'``; otherwise computed.\n";
            s += "df : ndarray, optional\n    Diagonal of ``U`` to reuse when ``fact='F'``; otherwise computed.\n";
            s += "duf : ndarray, optional\n    First superdiagonal of ``U`` to reuse when ``fact='F'``; otherwise computed.\n";
            s += "du2 : ndarray, optional\n    Second superdiagonal of ``U`` to reuse when ``fact='F'``; otherwise computed.\n";
            s += "ipiv : ndarray, optional\n"
                 "    Pivot indices to reuse when ``fact='F'``, 1-based as ``gttrf`` returns them;\n"
                 "    otherwise computed.\n";

            s += "\nReturns\n-------\n";
            s += "dlf : ndarray\n    Multipliers defining ``L``.\n";
            s += "df : ndarray\n    Diagonal of ``U``.\n";
            s += "duf : ndarray\n    First superdiagonal of ``U``.\n";
            s += "du2 : ndarray\n    Second superdiagonal of ``U``.\n";
            s += "ipiv : ndarray\n    Pivot indices, 1-based.\n";
            s += R_X_OUT;
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += "ferr : ndarray\n    Estimated forward error bound for each solution vector.\n";
            s += "berr : ndarray\n    Componentwise relative backward error of each solution vector.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if ``0 < info <= n``, ``u[info-1, info-1]`` is exactly zero; if ``info = n+1``,\n"
                 "    the matrix is singular to working precision and `x` may be inaccurate.\n";
            return s;
        }

        static std::string
        doc_stev(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, compute_v=1, overwrite_d=0, overwrite_e=0)\n\n";
            s += "Compute the eigenvalues and optionally the eigenvectors of a symmetric\n"
                 "tridiagonal matrix (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_SYM;
            s += "e : ndarray\n    Off-diagonal elements, length ``max(n - 1, 1)``.\n";
            s += P_COMPUTE_V;
            s += P_OVERWRITE_D;
            s += P_OVERWRITE_E;

            s += "\nReturns\n-------\n";
            s += "vals : ndarray\n    Eigenvalues in ascending order, written over `d`.\n";
            s += "z : ndarray\n"
                 "    Orthonormal eigenvectors as columns, shape ``(n, n)``, or a ``(1, 1)``\n"
                 "    placeholder when `compute_v` is 0.\n";
            s += R_INFO_CV;
            return s;
        }

        static std::string
        doc_stevd(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, compute_v=1, lwork=1+4*n+n*n, liwork=3+5*n, overwrite_d=0,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_e=0)\n\n";
            s += "Compute the eigenvalues and optionally the eigenvectors of a symmetric\n"
                 "tridiagonal matrix by divide-and-conquer (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_SYM;
            s += "e : ndarray\n    Off-diagonal elements, length ``max(n - 1, 1)``.\n";
            s += P_COMPUTE_V;
            s += "lwork : int, optional\n"
                 "    Size of the workspace. Default is ``1 + 4 * n + n ** 2`` when `compute_v` is\n"
                 "    nonzero and ``1`` otherwise, which are the exact requirements.\n";
            s += "liwork : int, optional\n"
                 "    Size of the integer workspace. Default is ``3 + 5 * n`` when `compute_v` is\n"
                 "    nonzero and ``1`` otherwise.\n";
            s += P_OVERWRITE_D;
            s += P_OVERWRITE_E;

            s += "\nReturns\n-------\n";
            s += "vals : ndarray\n    Eigenvalues in ascending order, written over `d`.\n";
            s += "z : ndarray\n"
                 "    Orthonormal eigenvectors as columns, shape ``(n, n)``, or a ``(1, 1)``\n"
                 "    placeholder when `compute_v` is 0.\n";
            s += R_INFO_CV;
            return s;
        }

        static std::string
        doc_sterf(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, overwrite_d=0, overwrite_e=0)\n\n";
            s += "Compute all eigenvalues of a symmetric tridiagonal matrix by the Pal-Walker-\n"
                 "Kahan QL/QR variant, without eigenvectors (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_SYM;
            s += P_E_SYM;
            s += P_OVERWRITE_D;
            s += P_OVERWRITE_E;

            s += "\nReturns\n-------\n";
            s += "vals : ndarray\n    Eigenvalues in ascending order, written over `d`.\n";
            s += R_INFO_CV;
            return s;
        }

        static std::string
        doc_stebz(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, range, vl, vu, il, iu, tol, order)\n\n";
            s += "Compute selected eigenvalues of a symmetric tridiagonal matrix by bisection\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_SYM;
            s += P_E_SYM;
            s += P_RANGE_SELECT;
            s += "tol : float\n"
                 "    Absolute tolerance for the eigenvalues. If not positive, machine precision\n"
                 "    times the matrix norm is used.\n";
            s += "order : str\n"
                 "    ``'B'`` to order eigenvalues within each block, ``'E'`` to order the whole\n"
                 "    set in ascending order.\n";

            s += "\nReturns\n-------\n";
            s += "m : int\n    Number of eigenvalues found; only ``w[:m]`` is meaningful.\n";
            s += "w : ndarray\n    The computed eigenvalues, in the order selected by `order`.\n";
            s += "iblock : ndarray\n"
                 "    Block index of each eigenvalue: ``w[i]`` belongs to submatrix ``iblock[i]``,\n"
                 "    1-based. A nonpositive entry marks an eigenvalue that failed to converge.\n";
            s += "isplit : ndarray\n"
                 "    Splitting points: submatrix ``i`` spans rows ``isplit[i-1]`` through\n"
                 "    ``isplit[i] - 1``, 1-based.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, some eigenvalues failed to converge or were not counted\n"
                 "    reliably.\n";
            return s;
        }

        static std::string
        doc_stein(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, w, iblock, isplit)\n\n";
            s += "Compute eigenvectors of a symmetric tridiagonal matrix by inverse iteration,\n"
                 "for eigenvalues supplied by ``stebz`` (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_SYM;
            s += P_E_SYM;
            s += "w : ndarray\n    The ``m`` eigenvalues to compute vectors for, as returned by ``stebz``.\n";
            s += "iblock : ndarray\n    Block indices from ``stebz``.\n";
            s += "isplit : ndarray\n    Splitting points from ``stebz``.\n";

            s += "\nReturns\n-------\n";
            s += "z : ndarray\n    Orthonormal eigenvectors as columns, shape ``(n, m)``.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, that many eigenvectors failed to converge.\n";
            return s;
        }

        static std::string
        doc_stemr(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, range, vl, vu, il, iu, compute_v=1, lwork=18*n,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "liwork=10*n, overwrite_d=0, overwrite_e=0)\n\n";
            s += "Compute selected eigenvalues and optionally eigenvectors of a symmetric\n"
                 "tridiagonal matrix by the MRRR algorithm (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_SYM;
            s += "e : ndarray\n"
                 "    Off-diagonal elements in the leading ``n - 1`` entries, length ``n``. The\n"
                 "    final entry need not be set; LAPACK uses it as workspace.\n";
            s += P_RANGE_SELECT;
            s += P_COMPUTE_V;
            /* Written out rather than via p_lwork: the default depends on `compute_v`, so it
             * wraps, and it reads better beside the `liwork` entry that mirrors it. */
            s += "lwork : int, optional\n"
                 "    Size of the workspace. Default is ``18 * n`` when `compute_v` is nonzero\n"
                 "    and ``12 * n`` otherwise. Use ``stemr_lwork`` for the optimal value.\n";
            s += "liwork : int, optional\n"
                 "    Size of the integer workspace. Default is ``10 * n`` when `compute_v` is\n"
                 "    nonzero and ``8 * n`` otherwise. Use ``stemr_lwork`` for the optimal value.\n";
            s += P_OVERWRITE_D;
            s += P_OVERWRITE_E;

            s += "\nReturns\n-------\n";
            s += "m : int\n    Number of eigenvalues found; only ``w[:m]`` and ``z[:, :m]`` are meaningful.\n";
            s += "w : ndarray\n    The computed eigenvalues in ascending order, length ``n``.\n";
            s += "z : ndarray\n"
                 "    Orthonormal eigenvectors as columns, shape ``(n, n)``. Not referenced when\n"
                 "    `compute_v` is 0.\n";
            s += R_INFO_CV;
            return s;
        }

        static std::string
        doc_stemr_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, range, vl, vu, il, iu, compute_v=1, overwrite_d=0,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_e=0)\n\n";
            s += "Query the optimal `lwork` and `liwork` for ``" + std::string(1, name[0]) + "stemr``.\n\n";

            s += "Parameters\n----------\n";
            s += P_D_SYM;
            s += "e : ndarray\n    Off-diagonal elements in the leading ``n - 1`` entries, length ``n``.\n";
            s += P_RANGE_SELECT;
            s += P_COMPUTE_V;
            s += P_OVERWRITE_D;
            s += P_OVERWRITE_E;

            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += "iwork : int\n    Optimal size of the integer workspace, to pass as `liwork`.\n";
            s += R_INFO;
            return s;
        }

        /* ================================ registration ============================= */

        typedef std::string (*DocFn)(const char *, const Dtype &);

        struct DocEntry {
            const char *name;
            DocFn fn;
            Dtype t;
        };

        /** @brief Expand a regular family (s/d/c/z prefixes) to its four docstring-table rows. */
        #define DOC_FAMILY(fam)       \
            {"s" #fam, doc_##fam, S}, \
            {"d" #fam, doc_##fam, D}, \
            {"c" #fam, doc_##fam, C}, \
            {"z" #fam, doc_##fam, Z}

        /* Same order as `gen_methods` in `lapack_gen.cpp`: each `_lwork` query sits
         * directly after the routine it queries. */
        static const DocEntry doc_table[] = {
            DOC_FAMILY(gees),
            DOC_FAMILY(gges),
            DOC_FAMILY(gebal),
            DOC_FAMILY(gehrd),
            DOC_FAMILY(gehrd_lwork),
            DOC_FAMILY(gesv),
            DOC_FAMILY(gecon),
            DOC_FAMILY(getrf),
            DOC_FAMILY(getrs),
            DOC_FAMILY(getc2),
            DOC_FAMILY(gesc2),
            DOC_FAMILY(getri),
            DOC_FAMILY(getri_lwork),
            DOC_FAMILY(gesdd),
            DOC_FAMILY(gesdd_lwork),
            DOC_FAMILY(gesvd),
            DOC_FAMILY(gesvd_lwork),
            DOC_FAMILY(gels),
            DOC_FAMILY(gels_lwork),
            DOC_FAMILY(gelss),
            DOC_FAMILY(gelss_lwork),
            DOC_FAMILY(gelsy),
            DOC_FAMILY(gelsy_lwork),
            DOC_FAMILY(gelsd),
            DOC_FAMILY(gelsd_lwork),
            DOC_FAMILY(geqp3),
            DOC_FAMILY(geqrf),
            DOC_FAMILY(geqrf_lwork),
            DOC_FAMILY(geqrfp),
            DOC_FAMILY(geqrfp_lwork),
            DOC_FAMILY(gerqf),
            DOC_FAMILY(geequ),
            DOC_FAMILY(geequb),
            DOC_FAMILY(geev),
            DOC_FAMILY(geev_lwork),
            DOC_FAMILY(ggev),
            DOC_FAMILY(gesvx),

            /* Same order as `gen_tri_methods` in `lapack_gen_tri.cpp`.  The `ste*` eigensolvers
             * are real-only, so they are spelled out as `s`/`d` pairs the way that table spells
             * them out with ROW. */
            DOC_FAMILY(gtsv),
            DOC_FAMILY(gttrf),
            DOC_FAMILY(gttrs),
            DOC_FAMILY(gtcon),
            DOC_FAMILY(gtsvx),
            {"sstev", doc_stev, S},
            {"dstev", doc_stev, D},
            {"sstebz", doc_stebz, S},
            {"dstebz", doc_stebz, D},
            {"ssterf", doc_sterf, S},
            {"dsterf", doc_sterf, D},
            {"sstein", doc_stein, S},
            {"dstein", doc_stein, D},
            {"sstemr", doc_stemr, S},
            {"dstemr", doc_stemr, D},
            {"sstemr_lwork", doc_stemr_lwork, S},
            {"dstemr_lwork", doc_stemr_lwork, D},
            {"sstevd", doc_stevd, S},
            {"dstevd", doc_stevd, D},
        };

        /** @brief The docstring for @p name, or nullptr when none is registered. */
        PyObject *build_doc(const char *name) noexcept
        {
            for (const DocEntry &e : doc_table) {
                if (std::strcmp(e.name, name) == 0) {
                    const std::string d = e.fn(e.name, e.t);
                    return PyUnicode_FromStringAndSize(d.data(), static_cast<Py_ssize_t>(d.size()));
                }
            }
            return nullptr;
        }

    }  // namespace capi
}  // namespace lapack
