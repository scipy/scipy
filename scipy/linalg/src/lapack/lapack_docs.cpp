/**
 * @file
 * @brief Docstrings for the `_flapack_cpp` wrappers, built on demand.
 *
 * Each wrapper's docstring is assembled at runtime from a template parameterised by the element
 * type.  `lapack::capi::build_doc(name)` is called lazily by the `LapackFunc.__doc__` getter (see
 * `lapack_module.cpp`) the first time a routine's docstring is requested, and the result is
 * cached on the callable.
 *
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

        /* The Cholesky family selects its triangle with an integer flag, not a letter --
         * `pocon` alone takes `uplo`, and says so in its own text. */
        static constexpr const char *P_LOWER =
            "lower : int, optional\n"
            "    If nonzero, the lower triangle of `a` is referenced and the factor is lower\n"
            "    triangular; otherwise the upper. Default is 0.\n";

        /* The positive definite tridiagonal group: `d` is real for every flavor, being the
         * diagonal of D in the L*D*L**H factorization of a Hermitian matrix. */
        static constexpr const char *P_D_PT =
            "d : ndarray\n    Diagonal of the matrix, length ``n``. Real for every flavor.\n";
        static constexpr const char *P_E_PT =
            "e : ndarray\n    Off-diagonal, length ``max(0, n - 1)``.\n";

        /* The symmetric/Hermitian indefinite group.  `lower` selects the stored triangle
         * everywhere in it, and the pivots stay 1-based from end to end. */
        static constexpr const char *P_LOWER_SH =
            "lower : int, optional\n"
            "    If nonzero, the lower triangle of `a` is referenced; otherwise the upper.\n"
            "    Default is 0.\n";
        static constexpr const char *P_IPIV_SH =
            "ipiv : ndarray\n"
            "    Pivot indices from the Bunch-Kaufman factorization, 1-based exactly as\n"
            "    ``sytrf``/``hetrf`` returned them -- this group shifts at neither end.\n";
        static constexpr const char *R_IPIV_SH =
            "ipiv : ndarray\n"
            "    Pivot indices describing the interchanges and the 1x1 / 2x2 block structure,\n"
            "    1-based. Pass them on to ``sytrs``/``sytri``/``sycon`` unchanged.\n";

        /* The symmetric/Hermitian eigensolvers.  Two flag conventions run through them and
         * neither can be normalised away: `syev`/`syevd`/`syevr`/`syevx` take integer
         * `compute_v`/`lower`, while the `sygv*` family takes `jobz`/`uplo` letters. */
        static constexpr const char *P_COMPUTE_V_SH =
            "compute_v : int, optional\n"
            "    If nonzero, eigenvectors are computed as well as eigenvalues. Default is 1.\n";
        static constexpr const char *P_JOBZ =
            "jobz : str, optional\n"
            "    ``'V'`` to compute eigenvectors, ``'N'`` for eigenvalues only. Default is\n"
            "    ``'V'``. This family spells it as a letter where ``syev`` uses `compute_v`.\n";
        static constexpr const char *P_UPLO_SH =
            "uplo : str, optional\n"
            "    ``'U'`` or ``'L'`` for the triangle of `a` to reference. Default is ``'L'``\n"
            "    -- note that is the *opposite* of the `lower=0` default elsewhere in this\n"
            "    group, so the same matrix can give different answers between families.\n";
        static constexpr const char *P_RANGE_SH =
            "range : str, optional\n"
            "    ``'A'`` for all eigenvalues, ``'V'`` for those in ``(vl, vu]``, ``'I'`` for\n"
            "    those indexed `il` through `iu`. Default is ``'A'``.\n";
        static constexpr const char *P_SELECT_BOUNDS =
            "vl : float, optional\n    Lower bound of the interval; used only when `range` is ``'V'``. Default is 0.0.\n"
            "vu : float, optional\n    Upper bound of the interval; used only when `range` is ``'V'``. Default is 1.0.\n"
            "il : int, optional\n    Index of the smallest eigenvalue to return, 1-based; used only when `range`\n    is ``'I'``. Default is 1.\n"
            "iu : int, optional\n    Index of the largest eigenvalue to return, 1-based; used only when `range`\n    is ``'I'``. Default is ``n``.\n"
            "abstol : float, optional\n    Absolute error tolerance for the eigenvalues. Default is 0.0.\n";
        static constexpr const char *P_ITYPE =
            "itype : int, optional\n"
            "    Which generalized problem to solve: 1 for ``a @ x = w * b @ x``, 2 for\n"
            "    ``a @ b @ x = w * x``, 3 for ``b @ a @ x = w * x``. Default is 1.\n";
        static constexpr const char *R_W_EIG =
            "w : ndarray\n    Eigenvalues in ascending order, length ``n``. Real for every flavor.\n";

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
            s += std::string(std::strlen(name) + 1, ' ')
               + (t.is_complex ? "lwork=2*n, " : "lwork=8*n+16, ")
               + sel + "_extra_args=(), overwrite_a=0, overwrite_b=0)\n\n";
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
         *         only `gelsy` makes `cond` required.  `gelsd` alone sizes the integer workspace
         *         as well -- and the real one in the complex flavors -- so only it reports those
         *         alongside `work`. */
        static std::string
        doc_gelsx_lwork(const char *name, const Dtype &t, const char *base, const char *args,
                        bool sizes_workspaces = false) noexcept
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
            if (sizes_workspaces) {
                if (t.is_complex) {
                    s += "rwork : float\n"
                         "    Optimal size of the real workspace, to pass as `size_rwork`. The real\n"
                         "    flavors have no such output, so they return one value fewer.\n";
                }
                s += "iwork : int\n"
                     "    Optimal size of the integer workspace, to pass as `size_iwork`.\n";
            }
            s += R_INFO;
            return s;
        }

        static std::string doc_gelss_lwork(const char *name, const Dtype &t) noexcept
            { return doc_gelsx_lwork(name, t, "gelss", "(m, n, nrhs, cond=-1.0, lwork=-1)"); }
        static std::string doc_gelsy_lwork(const char *name, const Dtype &t) noexcept
            { return doc_gelsx_lwork(name, t, "gelsy", "(m, n, nrhs, cond, lwork=-1)"); }
        static std::string doc_gelsd_lwork(const char *name, const Dtype &t) noexcept
            { return doc_gelsx_lwork(name, t, "gelsd", "(m, n, nrhs, cond=-1.0, lwork=-1)", true); }

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

        /* =========================== flapack_gen_banded.pyf.src ==================== */

        static constexpr const char *P_KL = "kl : int\n    Number of subdiagonals.\n";
        static constexpr const char *P_KU = "ku : int\n    Number of superdiagonals.\n";

        static std::string
        p_ab(const char *rows) noexcept
        {
            return "ab : ndarray\n    Band storage of shape ``(" + std::string(rows) + ", n)``:"
                   " ``a[i, j]`` is held at\n    ``ab[ku + i - j, j]``.\n";
        }

        static std::string
        doc_gbsv(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(kl, ku, ab, b, overwrite_ab=0, overwrite_b=0)\n\n";
            s += "Solve ``a @ x = b`` for a banded ``a`` by LU factorization with partial\n"
                 "pivoting (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_KL;
            s += P_KU;
            s += p_ab("2 * kl + ku + 1");
            s += P_B_RHS;
            s += "overwrite_ab : int, optional\n    If nonzero, `ab` may be overwritten in place. Default is 0.\n";
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "lub : ndarray\n    Combined ``L`` and ``U`` factors in band storage.\n";
            s += R_PIV_OUT;
            s += R_X_OUT;
            s += R_INFO_LU;
            return s;
        }

        static std::string
        doc_gbtrf(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(ab, kl, ku, m=shape(ab, 1), n=shape(ab, 1),\n";
            s += std::string(std::strlen(name) + 1, ' ') + "ldab=shape(ab, 0), overwrite_ab=0)\n\n";
            s += "Compute the LU factorization of a banded matrix with partial pivoting\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += p_ab("2 * kl + ku + 1");
            s += P_KL;
            s += P_KU;
            s += "m : int, optional\n"
                 "    Number of rows to factor. Default is ``shape(ab, 1)``, a square matrix;\n"
                 "    a smaller value factors the leading `m` rows.\n";
            s += "n : int, optional\n    Number of columns; must equal ``shape(ab, 1)``.\n";
            s += "ldab : int, optional\n    Leading dimension; must equal ``shape(ab, 0)``.\n";
            s += "overwrite_ab : int, optional\n    If nonzero, `ab` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "lu : ndarray\n    Combined ``L`` and ``U`` factors in band storage.\n";
            s += "ipiv : ndarray\n    Pivot indices, 0-based, length ``min(m, n)``.\n";
            s += R_INFO_LU;
            return s;
        }

        static std::string
        doc_gbtrs(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(ab, kl, ku, b, ipiv, trans=0, n=shape(ab, 1),\n";
            s += std::string(std::strlen(name) + 1, ' ')
               + "ldab=shape(ab, 0), ldb=shape(b, 0), overwrite_b=0)\n\n";
            s += "Solve a banded system using the factorization from ``gbtrf``\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += p_ab("2 * kl + ku + 1");
            s += P_KL;
            s += P_KU;
            s += P_B_RHS;
            s += "ipiv : ndarray\n    Pivot indices from ``gbtrf``, 0-based.\n";
            s += P_TRANS_INT;
            s += "n : int, optional\n    Order of the system; must equal ``shape(ab, 1)``.\n";
            s += "ldab : int, optional\n    Leading dimension of `ab`; must equal ``shape(ab, 0)``.\n";
            s += "ldb : int, optional\n    Leading dimension of `b`; must equal ``shape(b, 0)``.\n";
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_X_OUT;
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gbcon(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(kl, ku, ab, ipiv, anorm, norm='1', ldab=2*kl+ku+1)\n\n";
            s += "Estimate the reciprocal condition number of a banded matrix from its ``gbtrf``\n"
                 "factorization (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_KL;
            s += P_KU;
            s += p_ab("2 * kl + ku + 1");
            s += "ipiv : ndarray\n    Pivot indices from ``gbtrf``, 0-based.\n";
            s += "anorm : float\n    Norm of the original matrix, in the norm selected by `norm`.\n";
            s += "norm : str, optional\n    ``'1'`` or ``'O'`` for the 1-norm, ``'I'`` for the infinity norm.\n    Default is ``'1'``.\n";
            s += "ldab : int, optional\n"
                 "    Leading dimension of `ab`, at least ``2 * kl + ku + 1``, which is also the\n"
                 "    default.\n";

            s += "\nReturns\n-------\n";
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_langb(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(norm, kl, ku, ab, ldab=kl+ku+1)\n\n";
            s += "Compute a norm of a banded matrix (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "norm : str\n"
                 "    ``'M'`` for the largest absolute value, ``'1'`` or ``'O'`` for the 1-norm,\n"
                 "    ``'I'`` for the infinity norm, ``'F'`` or ``'E'`` for the Frobenius norm.\n"
                 "    Lower case is accepted. Note that ``'M'`` is not a consistent matrix norm.\n";
            s += P_KL;
            s += P_KU;
            /* Read-only, so no fill-in rows: this is the one `ab` in the group without them. */
            s += p_ab("kl + ku + 1");
            s += "ldab : int, optional\n"
                 "    Leading dimension of `ab`, at least ``kl + ku + 1``, which is also the\n"
                 "    default.\n";

            s += "\nReturns\n-------\n";
            s += "n2 : " + std::string(t.is_complex ? "float" : t.scalar)
               + "\n    The requested norm, always real.\n";
            return s;
        }

        /* ============================= flapack_pos_def.pyf.src ===================== */

        /** @brief `pstrf` and `pstf2` are the blocked and unblocked pivoted Cholesky; one
         *         interface, one docstring, differing only in how they are described. */
        static std::string
        doc_pstrf_family(const char *name, const char *how) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, tol=-1.0, lower=0, overwrite_a=0)\n\n";
            s += "Compute the pivoted Cholesky factorization of a positive semidefinite matrix,\n"
                 + std::string(how) + " (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Symmetric or Hermitian positive semidefinite matrix of shape ``(n, n)``.\n";
            s += "tol : float, optional\n"
                 "    Pivots at or below this value stop the factorization. Default is -1.0,\n"
                 "    which uses ``n * eps * max(diag(a))``.\n";
            s += P_LOWER;
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n"
                 "    The factor, in the selected triangle, of ``p.T @ a @ p``. The opposite\n"
                 "    triangle is left as it was in `a`.\n";
            s += "piv : ndarray\n"
                 "    Pivot indices, **1-based**: column ``i`` of ``p`` is column ``piv[i] - 1``\n"
                 "    of the identity. Unlike ``getrf``, this family does not shift them down.\n";
            s += "rank_c : int\n    Computed rank of `a`, as determined by `tol`.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if 1, `a` is not positive definite but the rank-deficient factorization\n"
                 "    completed.\n";
            return s;
        }

        static std::string doc_pstrf(const char *name, const Dtype &) noexcept
            { return doc_pstrf_family(name, "blocked"); }
        static std::string doc_pstf2(const char *name, const Dtype &) noexcept
            { return doc_pstrf_family(name, "unblocked"); }

        static std::string
        doc_posv(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, lower=0, overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve ``a @ x = b`` for a positive definite ``a`` by Cholesky factorization\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Symmetric or Hermitian positive definite matrix of shape ``(n, n)``.\n";
            s += P_B_RHS;
            s += P_LOWER;
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    Cholesky factor of `a`, in the selected triangle.\n";
            s += R_X_OUT;
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, the leading minor of order `info` is not positive definite.\n";
            return s;
        }

        static std::string
        doc_posvx(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, fact='E', af=None, equed='Y', s=None, lower=0,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve ``a @ x = b`` for a positive definite ``a`` with equilibration, condition\n"
                 "estimation and error bounds (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Symmetric or Hermitian positive definite matrix of shape ``(n, n)``.\n";
            s += P_B_RHS;
            s += "fact : str, optional\n"
                 "    ``'E'`` to equilibrate then factorize, ``'N'`` to factorize as given,\n"
                 "    ``'F'`` to reuse the `af`, `equed` and `s` supplied. Default is ``'E'``.\n";
            s += "af : ndarray, optional\n    Cholesky factor to reuse when ``fact='F'``; otherwise it is computed.\n";
            s += "equed : str, optional\n"
                 "    Equilibration already applied when ``fact='F'``: ``'N'`` or ``'Y'``.\n"
                 "    Otherwise it is an output. Default is ``'Y'``.\n";
            s += "s : ndarray, optional\n    Scale factors, used when ``fact='F'`` and `equed` is ``'Y'``.\n";
            s += P_LOWER;
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "a_s : ndarray\n    `a`, equilibrated if `equed` is ``b'Y'``.\n";
            s += "lu : ndarray\n    Cholesky factor of the equilibrated matrix.\n";
            s += "equed : bytes\n    Equilibration actually applied: ``b'N'`` or ``b'Y'``.\n";
            s += "s : ndarray\n    Scale factors.\n";
            s += "b_s : ndarray\n    `b`, scaled to match the equilibrated system.\n";
            s += R_X_OUT;
            s += "rcond : float\n    Estimate of the reciprocal condition number of the equilibrated matrix.\n";
            s += "ferr : ndarray\n    Estimated forward error bound for each solution vector.\n";
            s += "berr : ndarray\n    Componentwise relative backward error of each solution vector.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if ``0 < info <= n``, the leading minor of order `info` is not positive\n"
                 "    definite; if ``info = n+1``, the matrix is singular to working precision\n"
                 "    and `x` may be inaccurate.\n";
            return s;
        }

        static std::string
        doc_pocon(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, anorm, uplo='U')\n\n";
            s += "Estimate the reciprocal condition number of a positive definite matrix from its\n"
                 "Cholesky factorization (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Cholesky factorization of a positive definite matrix, as returned by ``potrf``.\n";
            s += "anorm : float\n    1-norm of the original matrix, which equals its infinity norm here.\n";
            /* The one routine in this group naming its triangle with a letter, not `lower`. */
            s += "uplo : str, optional\n"
                 "    ``'U'`` if `a` holds the upper triangular factor, ``'L'`` if lower.\n"
                 "    Default is ``'U'``. The rest of this group uses `lower` instead.\n";

            s += "\nReturns\n-------\n";
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_potrf(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, lower=0, clean=1, overwrite_a=0)\n\n";
            s += "Compute the Cholesky factorization of a positive definite matrix\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Symmetric or Hermitian positive definite matrix of shape ``(n, n)``.\n";
            s += P_LOWER;
            s += "clean : int, optional\n"
                 "    If nonzero, the triangle LAPACK does not write is zeroed, so `c` is\n"
                 "    genuinely triangular. If 0, that triangle still holds whatever `a` had\n"
                 "    there. Default is 1.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n"
                 "    Cholesky factor: ``c.conj().T @ c == a`` for the upper triangle,\n"
                 "    ``c @ c.conj().T == a`` for the lower.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, the leading minor of order `info` is not positive definite.\n";
            return s;
        }

        static std::string
        doc_potrs(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(c, b, lower=0, overwrite_b=0)\n\n";
            s += "Solve a positive definite system using the factorization from ``potrf``\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "c : ndarray\n    Cholesky factor, as returned by ``potrf``.\n";
            s += P_B_RHS;
            s += "lower : int, optional\n"
                 "    If nonzero, `c` is the lower triangular factor; otherwise the upper. Must\n"
                 "    match what was passed to ``potrf``. Default is 0.\n";
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_X_OUT;
            s += R_INFO;
            return s;
        }

        static std::string
        doc_potri(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(c, lower=0, overwrite_c=0)\n\n";
            s += "Invert a positive definite matrix from the factorization computed by ``potrf``\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "c : ndarray\n    Cholesky factor, as returned by ``potrf``.\n";
            s += "lower : int, optional\n"
                 "    If nonzero, `c` is the lower triangular factor; otherwise the upper. Must\n"
                 "    match what was passed to ``potrf``. Default is 0.\n";
            s += "overwrite_c : int, optional\n    If nonzero, `c` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "inv_a : ndarray\n"
                 "    The inverse, stored in the selected triangle only; the opposite triangle\n"
                 "    is left as it was in `c`.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, ``c[info-1, info-1]`` is zero and the inverse does not exist.\n";
            return s;
        }

        /* =========================== flapack_pos_def_tri.pyf.src =================== */

        static std::string
        doc_ptsv(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, b, overwrite_d=0, overwrite_e=0, overwrite_b=0)\n\n";
            s += "Solve ``a @ x = b`` for a positive definite tridiagonal ``a`` by its\n"
                 "``L @ D @ L.conj().T`` factorization (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_PT;
            s += P_E_PT;
            s += P_B_RHS;
            s += "overwrite_d : int, optional\n    If nonzero, `d` may be overwritten in place. Default is 0.\n";
            s += "overwrite_e : int, optional\n    If nonzero, `e` may be overwritten in place. Default is 0.\n";
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "d : ndarray\n    Diagonal of ``D`` from the factorization.\n";
            s += "du : ndarray\n    Off-diagonal of ``L``, written over `e`.\n";
            s += R_X_OUT;
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, the leading minor of order `info` is not positive definite.\n";
            return s;
        }

        static std::string
        doc_pttrf(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, overwrite_d=0, overwrite_e=0)\n\n";
            s += "Compute the ``L @ D @ L.conj().T`` factorization of a positive definite\n"
                 "tridiagonal matrix (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_PT;
            s += P_E_PT;
            s += "overwrite_d : int, optional\n    If nonzero, `d` may be overwritten in place. Default is 0.\n";
            s += "overwrite_e : int, optional\n    If nonzero, `e` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "d : ndarray\n    Diagonal of ``D``.\n";
            s += "e : ndarray\n    Off-diagonal of the unit bidiagonal ``L``.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, the leading minor of order `info` is not positive definite.\n";
            return s;
        }

        /** @brief The complex flavors take a `lower` the real ones do not have, so the
         *         signature and one parameter entry differ between them. */
        static std::string
        doc_pttrs(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + (t.is_complex ? "(d, e, b, lower=0, overwrite_b=0)\n\n"
                                                   : "(d, e, b, overwrite_b=0)\n\n");
            s += "Solve a positive definite tridiagonal system using the factorization from\n"
                 "``pttrf`` (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "d : ndarray\n    Diagonal of ``D`` from ``pttrf``. Real for every flavor.\n";
            s += "e : ndarray\n    Off-diagonal of ``L`` from ``pttrf``.\n";
            s += P_B_RHS;
            if (t.is_complex) {
                s += "lower : int, optional\n"
                     "    If nonzero, `e` is read as the subdiagonal of ``L``; otherwise as the\n"
                     "    superdiagonal of ``U``, which is its conjugate. Default is 0. Note\n"
                     "    ``pttrf`` writes the ``L`` form, so reusing its output here needs\n"
                     "    ``lower=1``; the default silently solves a different system. The real\n"
                     "    flavors have no such argument, drawing no distinction between the two.\n";
            }
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_X_OUT;
            s += R_INFO;
            return s;
        }

        static std::string
        doc_pteqr(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, z, compute_z=0, overwrite_d=0, overwrite_e=0,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_z=0)\n\n";
            s += "Compute the eigenvalues and optionally the eigenvectors of a positive definite\n"
                 "tridiagonal matrix, to high relative accuracy (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_PT;
            s += "e : ndarray\n    Off-diagonal, length ``max(0, n - 1)``. Real for every flavor.\n";
            s += std::string("z : ndarray\n") +
                 "    Required even when unused. With ``compute_z=1`` it is the matrix whose\n"
                 "    eigenvectors are wanted, of shape ``(max(1, n), n)``" +
                 (t.is_complex ? " and complex;\n" : ";\n") +
                 "    with ``compute_z=2`` its contents are ignored and the eigenvectors of the\n"
                 "    tridiagonal matrix itself are returned; with ``compute_z=0`` it is not\n"
                 "    referenced at all and any shape is accepted.\n";
            s += "compute_z : int, optional\n"
                 "    0 for eigenvalues only, 1 to compute eigenvectors of the original matrix\n"
                 "    from the `z` supplied, 2 for eigenvectors of the tridiagonal matrix.\n"
                 "    Default is 0.\n";
            s += "overwrite_d : int, optional\n    If nonzero, `d` may be overwritten in place. Default is 0.\n";
            s += "overwrite_e : int, optional\n    If nonzero, `e` may be overwritten in place. Default is 0.\n";
            s += "overwrite_z : int, optional\n    If nonzero, `z` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "d : ndarray\n    Eigenvalues in descending order.\n";
            s += "e : ndarray\n    Destroyed by the computation; the contents are not meaningful.\n";
            s += "z : ndarray\n    Orthonormal eigenvectors as columns, or `z` unchanged when `compute_z` is 0.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if ``0 < info <= n``, the leading minor of order `info` is not positive\n"
                 "    definite; if ``info > n``, that many off-diagonal elements failed to\n"
                 "    converge to zero.\n";
            return s;
        }

        static std::string
        doc_ptsvx(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(d, e, b, fact='N', df=None, ef=None)\n\n";
            s += "Solve a positive definite tridiagonal system with condition estimation and\n"
                 "error bounds (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_D_PT;
            s += P_E_PT;
            s += P_B_RHS;
            s += "fact : str, optional\n"
                 "    ``'N'`` to factorize the matrix, ``'F'`` to reuse the `df` and `ef`\n"
                 "    supplied. Default is ``'N'``.\n";
            s += "df : ndarray, optional\n"
                 "    Diagonal of ``D`` to reuse when ``fact='F'``; otherwise computed. Real for\n"
                 "    every flavor.\n";
            s += "ef : ndarray, optional\n    Off-diagonal of ``L`` to reuse when ``fact='F'``; otherwise computed.\n";

            s += "\nReturns\n-------\n";
            s += "df : ndarray\n    Diagonal of ``D``.\n";
            s += "ef : ndarray\n    Off-diagonal of ``L``.\n";
            s += R_X_OUT;
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += "ferr : ndarray\n    Estimated forward error bound for each solution vector.\n";
            s += "berr : ndarray\n    Componentwise relative backward error of each solution vector.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if ``0 < info <= n``, the leading minor of order `info` is not positive\n"
                 "    definite; if ``info = n+1``, the matrix is singular to working precision\n"
                 "    and `x` may be inaccurate.\n";
            return s;
        }

        /* ============================ flapack_sym_herm.pyf.src ===================== */

        /** @brief `sytrf`/`hetrf` and `sytf2` -- blocked and unblocked Bunch-Kaufman. @p what
         *         names the symmetry, @p lwork_line is empty for the unblocked `sytf2`. */
        static std::string
        doc_trf_family(const char *name, const char *what, bool blocked) noexcept
        {
            std::string s;
            s += std::string(name) + (blocked ? "(a, lower=0, lwork=max(n, 1), overwrite_a=0)\n\n"
                                              : "(a, lower=0, overwrite_a=0)\n\n");
            s += "Compute the " + std::string(blocked ? "blocked" : "unblocked")
               + " Bunch-Kaufman factorization of a " + std::string(what) + "\nmatrix (LAPACK ``"
               + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(what) + " matrix of shape ``(n, n)``.\n";
            s += P_LOWER_SH;
            if (blocked) {
                s += p_lwork("``max(n, 1)``", std::string(std::string(1, name[0])
                             + (std::strncmp(name + 1, "he", 2) == 0 ? "hetrf_lwork" : "sytrf_lwork")).c_str());
            }
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "ldu : ndarray\n"
                 "    The block diagonal ``D`` and the multipliers of ``U`` or ``L``, packed\n"
                 "    into the selected triangle.\n";
            s += R_IPIV_SH;
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, ``d[info-1, info-1]`` is exactly zero and the factor is\n"
                 "    singular.\n";
            return s;
        }

        static std::string doc_sytrf(const char *name, const Dtype &) noexcept
            { return doc_trf_family(name, "symmetric", true); }
        static std::string doc_hetrf(const char *name, const Dtype &) noexcept
            { return doc_trf_family(name, "Hermitian", true); }
        static std::string doc_sytf2(const char *name, const Dtype &) noexcept
            { return doc_trf_family(name, "symmetric", false); }

        /** @brief The `(n, lower)` workspace queries; @p parent is the routine queried. */
        static std::string
        doc_sh_lwork(const char *name, const Dtype &t, const char *parent) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, lower=0)\n\n";
            s += "Query the optimal `lwork` for ``" + std::string(1, name[0]) + std::string(parent)
               + "``.\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_LOWER_SH;

            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        static std::string doc_sytrf_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_lwork(name, t, "sytrf"); }
        static std::string doc_hetrf_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_lwork(name, t, "hetrf"); }
        static std::string doc_sysv_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_lwork(name, t, "sysv"); }
        static std::string doc_hesv_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_lwork(name, t, "hesv"); }
        static std::string doc_sysvx_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_lwork(name, t, "sysvx"); }
        static std::string doc_hesvx_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_lwork(name, t, "hesvx"); }

        /** @brief `sytrs`/`hetrs` -- solve using the factorization. */
        static std::string
        doc_trs_family(const char *name, const char *what) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, ipiv, b, lower=0, overwrite_b=0)\n\n";
            s += "Solve a " + std::string(what) + " indefinite system using the factorization\n"
                 "from ``" + std::string(std::strncmp(name + 1, "he", 2) == 0 ? "hetrf" : "sytrf")
               + "`` (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n"
                 "    The factorization, of shape ``(lda, n)``. A taller-than-wide buffer is\n"
                 "    accepted as long as ``lda >= n``.\n";
            s += P_IPIV_SH;
            s += P_B_RHS;
            s += P_LOWER_SH;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_X_OUT;
            s += R_INFO;
            return s;
        }

        static std::string doc_sytrs(const char *name, const Dtype &) noexcept
            { return doc_trs_family(name, "symmetric"); }
        static std::string doc_hetrs(const char *name, const Dtype &) noexcept
            { return doc_trs_family(name, "Hermitian"); }

        /** @brief `sytri`/`hetri` -- invert from the factorization. */
        static std::string
        doc_tri_family(const char *name, const char *what) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, ipiv, lower=0, overwrite_a=0)\n\n";
            s += "Invert a " + std::string(what) + " indefinite matrix from the factorization\n"
                 "computed by ``" + std::string(std::strncmp(name + 1, "he", 2) == 0 ? "hetrf" : "sytrf")
               + "`` (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    The factorization, as returned by the corresponding ``*trf``.\n";
            s += P_IPIV_SH;
            s += P_LOWER_SH;
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "inv_a : ndarray\n"
                 "    The inverse, stored in the selected triangle only; the opposite triangle\n"
                 "    is left as it was.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, ``d[info-1, info-1]`` is zero and the inverse does not exist.\n";
            return s;
        }

        static std::string doc_sytri(const char *name, const Dtype &) noexcept
            { return doc_tri_family(name, "symmetric"); }
        static std::string doc_hetri(const char *name, const Dtype &) noexcept
            { return doc_tri_family(name, "Hermitian"); }

        static std::string
        doc_syconv(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, ipiv, lower=0, way=0, overwrite_a=0)\n\n";
            s += "Convert the Bunch-Kaufman factorization between its packed and its expanded\n"
                 "storage (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    The factorization, as returned by ``sytrf``.\n";
            s += P_IPIV_SH;
            s += P_LOWER_SH;
            s += "way : int, optional\n"
                 "    0 to convert the packed form into `a` plus `e`, 1 to revert. An integer\n"
                 "    here rather than the ``'C'``/``'R'`` letter LAPACK takes. Default is 0.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The converted factorization.\n";
            s += "e : ndarray\n    Off-diagonal entries of the 2x2 blocks of ``D``, length ``n``.\n";
            s += R_INFO;
            return s;
        }

        /** @brief `syequb`/`heequb` -- scale factors that equilibrate the matrix. */
        static std::string
        doc_equb_family(const char *name, const char *what) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, lower=0)\n\n";
            s += "Compute scale factors that equilibrate a " + std::string(what) + " matrix\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(what) + " matrix of shape ``(n, n)``.\n";
            s += P_LOWER_SH;

            s += "\nReturns\n-------\n";
            s += "s : ndarray\n    Scale factors, length ``n``. Real for every flavor.\n";
            s += "scond : float\n    Ratio of the smallest to the largest scale factor.\n";
            s += "amax : float\n    Largest element of `a` in absolute value.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, the `info`-th diagonal entry is nonpositive.\n";
            return s;
        }

        static std::string doc_syequb(const char *name, const Dtype &) noexcept
            { return doc_equb_family(name, "symmetric"); }
        static std::string doc_heequb(const char *name, const Dtype &) noexcept
            { return doc_equb_family(name, "Hermitian"); }

        /** @brief `sycon`/`hecon` -- condition estimate from the factorization. */
        static std::string
        doc_con_family(const char *name, const char *what) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, ipiv, anorm, lower=0)\n\n";
            s += "Estimate the reciprocal condition number of a " + std::string(what)
               + " indefinite\nmatrix from its Bunch-Kaufman factorization (LAPACK ``"
               + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    The factorization, as returned by the corresponding ``*trf``.\n";
            s += P_IPIV_SH;
            s += "anorm : float\n    1-norm of the original matrix.\n";
            s += P_LOWER_SH;

            s += "\nReturns\n-------\n";
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += R_INFO;
            return s;
        }

        static std::string doc_sycon(const char *name, const Dtype &) noexcept
            { return doc_con_family(name, "symmetric"); }
        static std::string doc_hecon(const char *name, const Dtype &) noexcept
            { return doc_con_family(name, "Hermitian"); }

        /** @brief `sysv`/`hesv` -- the simple indefinite driver. */
        static std::string
        doc_sv_family(const char *name, const char *what, const char *out) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, lwork=max(n, 1), lower=0, overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve ``a @ x = b`` for a " + std::string(what)
               + " indefinite ``a`` by Bunch-Kaufman\nfactorization (LAPACK ``"
               + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(what) + " matrix of shape ``(n, n)``.\n";
            s += P_B_RHS;
            s += p_lwork("``max(n, 1)``", std::string(std::string(1, name[0])
                         + (std::strncmp(name + 1, "he", 2) == 0 ? "hesv_lwork" : "sysv_lwork")).c_str());
            s += P_LOWER_SH;
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += std::string(out) + " : ndarray\n    The block diagonal ``D`` and the multipliers, packed into the triangle.\n";
            s += R_IPIV_SH;
            s += R_X_OUT;
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if positive, ``d[info-1, info-1]`` is exactly zero and the solve failed.\n";
            return s;
        }

        static std::string doc_sysv(const char *name, const Dtype &) noexcept
            { return doc_sv_family(name, "symmetric", "udut"); }
        static std::string doc_hesv(const char *name, const Dtype &) noexcept
            { return doc_sv_family(name, "Hermitian", "uduh"); }

        static std::string
        doc_sysvx(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, af=None, ipiv=None, lwork=max(3*n, 1), factored=0,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "lower=0, overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve a symmetric indefinite system with condition estimation and error\n"
                 "bounds (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Symmetric matrix of shape ``(n, n)``.\n";
            s += P_B_RHS;
            s += "af : ndarray, optional\n    Factorization to reuse when ``factored=1``; otherwise it is computed.\n";
            s += "ipiv : ndarray, optional\n    Pivot indices to reuse when ``factored=1``, 1-based; otherwise computed.\n";
            s += "lwork : int, optional\n"
                 "    Size of the workspace. Default is ``max(3 * n, 1)``; the minimum is\n"
                 "    ``3 * n`` for the real flavors and ``2 * n`` for the complex ones. Use\n"
                 "    ``sysvx_lwork`` for the optimal value.\n";
            s += "factored : int, optional\n"
                 "    If nonzero, `af` and `ipiv` already hold the factorization and are reused.\n"
                 "    An integer here rather than the ``'F'``/``'N'`` letter ``gesvx`` takes.\n"
                 "    Default is 0.\n";
            s += P_LOWER_SH;
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "a_s : ndarray\n    `a`, unchanged; returned so the caller can chain calls.\n";
            s += "udut : ndarray\n    The block diagonal ``D`` and the multipliers.\n";
            s += R_IPIV_SH;
            s += "b_s : ndarray\n    `b`, unchanged.\n";
            s += R_X_OUT;
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += "ferr : ndarray\n    Estimated forward error bound for each solution vector.\n";
            s += "berr : ndarray\n    Componentwise relative backward error of each solution vector.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if ``0 < info <= n``, ``d[info-1, info-1]`` is exactly zero; if\n"
                 "    ``info = n+1``, the matrix is singular to working precision.\n";
            return s;
        }

        static std::string
        doc_hesvx(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, af=None, ipiv=None, lwork=max(2*n, 1), factored=0,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "lower=0, overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve a Hermitian indefinite system with condition estimation and error\n"
                 "bounds (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Hermitian matrix of shape ``(n, n)``.\n";
            s += P_B_RHS;
            s += "af : ndarray, optional\n    Factorization to reuse when ``factored=1``; otherwise it is computed.\n";
            s += "ipiv : ndarray, optional\n    Pivot indices to reuse when ``factored=1``, 1-based; otherwise computed.\n";
            s += p_lwork("``max(2 * n, 1)``", std::string(std::string(1, name[0]) + "hesvx_lwork").c_str());
            s += "factored : int, optional\n"
                 "    If nonzero, `af` and `ipiv` already hold the factorization and are reused.\n"
                 "    Default is 0.\n";
            s += P_LOWER_SH;
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            /* Seven outputs, not the nine `sysvx` returns: the Hermitian driver keeps `a` and
             * `b` to itself. */
            s += "uduh : ndarray\n    The block diagonal ``D`` and the multipliers.\n";
            s += R_IPIV_SH;
            s += R_X_OUT;
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += "ferr : ndarray\n    Estimated forward error bound for each solution vector.\n";
            s += "berr : ndarray\n    Componentwise relative backward error of each solution vector.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    if ``0 < info <= n``, ``d[info-1, info-1]`` is exactly zero; if\n"
                 "    ``info = n+1``, the matrix is singular to working precision.\n";
            return s;
        }

        /* ======================= flapack_sym_herm.pyf.src, eigen ==================== */

        /** @brief The `sy`/`he` name of a merged pair, from the wrapper's own first letter. */
        static std::string
        sh_name(const char *name, const char *stem) noexcept
        {
            return std::string(1, name[0]) + (name[1] == 'h' ? "he" : "sy") + std::string(stem);
        }

        /** @brief `syev`/`heev` -- the simple QR-iteration driver. */
        static std::string
        doc_ev(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, compute_v=1, lower=0, lwork=..., overwrite_a=0)\n\n";
            s += "Compute the eigenvalues and optionally the eigenvectors of a "
               + std::string(t.is_complex ? "Hermitian" : "real symmetric") + "\nmatrix (LAPACK ``"
               + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(t.is_complex ? "Hermitian" : "Symmetric")
               + " matrix of shape ``(n, n)``.\n";
            s += P_COMPUTE_V_SH;
            s += P_LOWER_SH;
            s += p_lwork(t.is_complex ? "``max(2 * n - 1, 1)``" : "``max(3 * n - 1, 1)``",
                         sh_name(name, "ev_lwork").c_str());
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += R_W_EIG;
            s += "v : ndarray\n"
                 "    Orthonormal eigenvectors as columns when `compute_v` is nonzero;\n"
                 "    otherwise the contents are not meaningful.\n";
            s += R_INFO_CV;
            return s;
        }

        /** @brief The `(n, lower)` eigen workspace queries. */
        static std::string
        doc_sh_ev_lwork(const char *name, const Dtype &t, const char *stem) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, lower=0)\n\n";
            s += "Query the optimal `lwork` for ``" + sh_name(name, stem) + "``.\n\n";
            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_LOWER_SH;
            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        static std::string doc_ev_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_ev_lwork(name, t, "ev"); }
        static std::string doc_evx_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_ev_lwork(name, t, "evx"); }

        /** @brief `syevd`/`heevd` -- divide-and-conquer.  Split because the complex half takes
         *         an extra `lrwork`; @p extra adds its entry. */
        static std::string
        doc_evd_family(const char *name, const char *what, bool complex_side) noexcept
        {
            std::string s;
            s += std::string(name) + (complex_side
                     ? "(a, compute_v=1, lower=0, lwork=..., liwork=..., lrwork=...,\n"
                     : "(a, compute_v=1, lower=0, lwork=..., liwork=...,\n");
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_a=0)\n\n";
            s += "Compute the eigenvalues and optionally the eigenvectors of a "
               + std::string(what) + "\nmatrix by divide-and-conquer (LAPACK ``"
               + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(what) + " matrix of shape ``(n, n)``.\n";
            s += P_COMPUTE_V_SH;
            s += P_LOWER_SH;
            if (complex_side) {
                s += "lwork : int, optional\n    Default is ``2 * n + n ** 2`` with eigenvectors, ``n + 1`` without.\n";
                s += "liwork : int, optional\n    Default is ``3 + 5 * n`` with eigenvectors, ``1`` without.\n";
                s += "lrwork : int, optional\n"
                     "    Size of the real workspace, which the real flavors do not have. Default\n"
                     "    is ``1 + 5 * n + 2 * n ** 2`` with eigenvectors, ``n`` without.\n";
            }
            else {
                s += "lwork : int, optional\n"
                     "    Default is ``1 + 6 * n + 2 * n ** 2`` with eigenvectors, ``2 * n + 1``\n"
                     "    without.\n";
                s += "liwork : int, optional\n    Default is ``3 + 5 * n`` with eigenvectors, ``1`` without.\n";
            }
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += R_W_EIG;
            s += "v : ndarray\n    Orthonormal eigenvectors as columns when `compute_v` is nonzero.\n";
            s += R_INFO_CV;
            return s;
        }

        static std::string doc_evd(const char *name, const Dtype &t) noexcept
            { return doc_evd_family(name, t.is_complex ? "Hermitian" : "real symmetric", t.is_complex); }

        static std::string
        doc_evd_lwork_family(const char *name, const Dtype &t, bool complex_side) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, compute_v=1, lower=0)\n\n";
            s += "Query the optimal workspace sizes for ``" + sh_name(name, "evd") + "``.\n\n";
            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_COMPUTE_V_SH;
            s += P_LOWER_SH;
            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += "iwork : int\n    Optimal size of the integer workspace, to pass as `liwork`.\n";
            if (complex_side) {
                s += "rwork : float\n"
                     "    Optimal size of the real workspace, to pass as `lrwork`. The real\n"
                     "    flavors have no such output, so they return one value fewer.\n";
            }
            s += R_INFO;
            return s;
        }

        static std::string doc_evd_lwork(const char *name, const Dtype &t) noexcept
            { return doc_evd_lwork_family(name, t, t.is_complex); }

        /** @brief `syevr`/`heevr` -- the MRRR driver. */
        static std::string
        doc_evr_family(const char *name, const char *what, bool complex_side) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, compute_v=1, range='A', lower=0, vl=0.0, vu=1.0, il=1,\n";
            s += std::string(std::strlen(name) + 1, ' ')
               + (complex_side ? "iu=n, abstol=0.0, lwork=..., lrwork=..., liwork=...,\n"
                               : "iu=n, abstol=0.0, lwork=..., liwork=...,\n");
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_a=0)\n\n";
            s += "Compute selected eigenvalues and optionally eigenvectors of a "
               + std::string(what) + "\nmatrix by the MRRR algorithm (LAPACK ``"
               + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(what) + " matrix of shape ``(n, n)``.\n";
            s += P_COMPUTE_V_SH;
            s += P_RANGE_SH;
            s += P_LOWER_SH;
            s += P_SELECT_BOUNDS;
            s += "lwork : int, optional\n"
                 "    Default is ``max(26 * n, 1)`` for the real flavors, ``max(2 * n, 1)``\n"
                 "    for the complex.\n";
            if (complex_side) {
                s += "lrwork : int, optional\n    Size of the real workspace. Default is ``max(24 * n, 1)``.\n";
            }
            s += "liwork : int, optional\n    Size of the integer workspace. Default is ``max(10 * n, 1)``.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += R_W_EIG;
            s += "z : ndarray\n    The `m` computed eigenvectors as columns, or empty when `compute_v` is 0.\n";
            s += "m : int\n    Number of eigenvalues found.\n";
            s += "isuppz : ndarray\n"
                 "    Support of each eigenvector: rows ``isuppz[2*i]`` through ``isuppz[2*i+1]``\n"
                 "    are nonzero, 1-based.\n";
            s += "    Length ``2 * n`` when `range` is ``'A'``, when `range` is ``'I'`` spanning\n"
                 "    the whole index range, or when ``n == 1``; empty otherwise, those being the\n"
                 "    only cases LAPACK writes it.\n";
            s += R_INFO_CV;
            return s;
        }

        static std::string doc_evr(const char *name, const Dtype &t) noexcept
            { return doc_evr_family(name, t.is_complex ? "Hermitian" : "real symmetric", t.is_complex); }

        static std::string
        doc_evr_lwork_family(const char *name, const Dtype &t, bool complex_side) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, lower=0)\n\n";
            s += "Query the optimal workspace sizes for ``" + sh_name(name, "evr") + "``.\n\n";
            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_LOWER_SH;
            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            if (complex_side) {
                /* `rwork` sits before `iwork` here, the other way round from `heevd_lwork`. */
                s += "rwork : float\n    Optimal size of the real workspace, to pass as `lrwork`.\n";
            }
            s += "iwork : int\n    Optimal size of the integer workspace, to pass as `liwork`.\n";
            s += R_INFO;
            return s;
        }

        static std::string doc_evr_lwork(const char *name, const Dtype &t) noexcept
            { return doc_evr_lwork_family(name, t, t.is_complex); }

        /** @brief `syevx`/`heevx` -- bisection plus inverse iteration. */
        static std::string
        doc_evx(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, compute_v=1, range='A', lower=0, vl=0.0, vu=1.0, il=1,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "iu=n, abstol=0.0, lwork=..., overwrite_a=0)\n\n";
            s += "Compute selected eigenvalues and optionally eigenvectors of a "
               + std::string(t.is_complex ? "Hermitian" : "real symmetric")
               + "\nmatrix (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(t.is_complex ? "Hermitian" : "Symmetric")
               + " matrix of shape ``(n, n)``.\n";
            s += P_COMPUTE_V_SH;
            s += P_RANGE_SH;
            s += P_LOWER_SH;
            s += P_SELECT_BOUNDS;
            s += p_lwork(t.is_complex ? "``max(2 * n, 1)``" : "``max(8 * n, 1)``",
                         sh_name(name, "evx_lwork").c_str());
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += R_W_EIG;
            s += "z : ndarray\n    The `m` computed eigenvectors as columns, or empty when `compute_v` is 0.\n";
            s += "m : int\n    Number of eigenvalues found.\n";
            s += "ifail : ndarray\n    Indices of eigenvectors that failed to converge, 1-based.\n";
            s += R_INFO_CV;
            return s;
        }

        /** @brief `sygv`/`hegv` -- the generalized definite problem. */
        static std::string
        doc_gv(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, itype=1, jobz='V', uplo='L', lwork=...,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve the generalized "
               + std::string(t.is_complex ? "Hermitian" : "symmetric")
               + "-definite eigenproblem (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(t.is_complex ? "Hermitian" : "Symmetric")
               + " matrix of shape ``(n, n)``.\n";
            s += "b : ndarray\n    Positive definite matrix of shape ``(n, n)``.\n";
            s += P_ITYPE;
            s += P_JOBZ;
            s += P_UPLO_SH;
            s += p_lwork(t.is_complex ? "``max(2 * n - 1, 1)``" : "``max(3 * n - 1, 1)``",
                         sh_name(name, "gv_lwork").c_str());
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_W_EIG;
            s += "v : ndarray\n    Eigenvectors as columns when `jobz` is ``'V'``, normalized per `itype`.\n";
            s += R_INFO_CV;
            return s;
        }

        /** @brief The `(n, uplo)` generalized workspace queries -- `uplo`, not `lower`. */
        static std::string
        doc_sh_gv_lwork(const char *name, const Dtype &t, const char *stem) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, uplo='L')\n\n";
            s += "Query the optimal `lwork` for ``" + sh_name(name, stem) + "``.\n\n";
            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_UPLO_SH;
            s += "\nReturns\n-------\n";
            s += "work : " + std::string(t.scalar) + "\n    Optimal size of the `work` array, as a scalar of the routine's dtype.\n";
            s += R_INFO;
            return s;
        }

        static std::string doc_gv_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_gv_lwork(name, t, "gv"); }
        static std::string doc_gvx_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_gv_lwork(name, t, "gvx"); }

        /** @brief `sygvd`/`hegvd` -- generalized divide-and-conquer, split for `lrwork`. */
        static std::string
        doc_gvd_family(const char *name, const char *what, bool complex_side) noexcept
        {
            std::string s;
            s += std::string(name) + (complex_side
                     ? "(a, b, itype=1, jobz='V', uplo='L', lwork=..., lrwork=...,\n"
                     : "(a, b, itype=1, jobz='V', uplo='L', lwork=...,\n");
            s += std::string(std::strlen(name) + 1, ' ') + "liwork=..., overwrite_a=0, overwrite_b=0)\n\n";
            s += "Solve the generalized " + std::string(what)
               + "-definite eigenproblem by\ndivide-and-conquer (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(complex_side ? "Hermitian" : "Symmetric")
               + " matrix of shape ``(n, n)``.\n";
            s += "b : ndarray\n    Positive definite matrix of shape ``(n, n)``.\n";
            s += P_ITYPE;
            s += P_JOBZ;
            s += P_UPLO_SH;
            if (complex_side) {
                s += "lwork : int, optional\n    Default is ``n * (n + 2)`` when `jobz` is ``'V'``, ``n + 1`` otherwise.\n";
                s += "lrwork : int, optional\n"
                     "    Size of the real workspace, which the real flavors do not have. Default\n"
                     "    is ``2 * n ** 2 + 5 * n + 1`` when `jobz` is ``'V'``, ``n`` otherwise.\n";
            }
            else {
                s += "lwork : int, optional\n"
                     "    Default is ``1 + 6 * n + 2 * n ** 2`` when `jobz` is ``'V'``,\n"
                     "    ``2 * n + 1`` otherwise.\n";
            }
            s += "liwork : int, optional\n    Default is ``5 * n + 3`` when `jobz` is ``'V'``, ``1`` otherwise.\n";
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_W_EIG;
            s += "v : ndarray\n    Eigenvectors as columns when `jobz` is ``'V'``.\n";
            s += R_INFO_CV;
            return s;
        }

        static std::string doc_gvd(const char *name, const Dtype &t) noexcept
            { return doc_gvd_family(name, t.is_complex ? "Hermitian" : "symmetric", t.is_complex); }

        /** @brief `sygvx`/`hegvx` -- selected generalized eigenvalues. */
        static std::string
        doc_gvx(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, itype=1, jobz='V', range='A', uplo='L', vl=0.0,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "vu=1.0, il=1, iu=n, abstol=0.0, lwork=...,\n";
            s += std::string(std::strlen(name) + 1, ' ') + "overwrite_a=0, overwrite_b=0)\n\n";
            s += "Compute selected eigenvalues of the generalized "
               + std::string(t.is_complex ? "Hermitian" : "symmetric")
               + "-definite\nproblem (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(t.is_complex ? "Hermitian" : "Symmetric")
               + " matrix of shape ``(n, n)``.\n";
            s += "b : ndarray\n    Positive definite matrix of shape ``(n, n)``.\n";
            s += P_ITYPE;
            s += P_JOBZ;
            s += P_RANGE_SH;
            s += P_UPLO_SH;
            s += P_SELECT_BOUNDS;
            s += p_lwork(t.is_complex ? "``max(2 * n, 1)``" : "``max(8 * n, 1)``",
                         sh_name(name, "gvx_lwork").c_str());
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += R_W_EIG;
            s += "z : ndarray\n    The `m` computed eigenvectors as columns, or empty when `jobz` is ``'N'``.\n";
            s += "m : int\n    Number of eigenvalues found.\n";
            s += "ifail : ndarray\n    Indices of eigenvectors that failed to converge, 1-based.\n";
            s += R_INFO_CV;
            return s;
        }

        /** @brief `sytrd`/`hetrd` -- reduction to tridiagonal form. */
        static std::string
        doc_trd(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, lower=0, lwork=max(n, 1), overwrite_a=0)\n\n";
            s += "Reduce a " + std::string(t.is_complex ? "Hermitian" : "real symmetric")
               + " matrix to real symmetric tridiagonal form\nby "
               + (t.is_complex ? "a unitary" : "an orthogonal")
               + " similarity transformation (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(t.is_complex ? "Hermitian" : "Symmetric")
               + " matrix of shape ``(n, n)``.\n";
            s += P_LOWER_SH;
            s += p_lwork("``max(n, 1)``", sh_name(name, "trd_lwork").c_str());
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n"
                 "    The tridiagonal form together with the elementary reflectors that produced\n"
                 "    it, packed into `a`.\n";
            s += "d : ndarray\n    Diagonal of the tridiagonal form, length ``n``. Real for every flavor.\n";
            s += "e : ndarray\n    Off-diagonal, length ``n - 1``. Real for every flavor.\n";
            s += R_TAU;
            s += R_INFO;
            return s;
        }

        static std::string doc_trd_lwork(const char *name, const Dtype &t) noexcept
            { return doc_sh_ev_lwork(name, t, "trd"); }

        /** @brief `sygst`/`hegst` -- reduce the generalized problem to standard form. */
        static std::string
        doc_gst(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, itype=1, lower=0, overwrite_a=0)\n\n";
            s += "Reduce a generalized " + std::string(t.is_complex ? "Hermitian" : "symmetric")
               + "-definite problem to standard\nform (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    " + std::string(t.is_complex ? "Hermitian" : "Symmetric")
               + " matrix of shape ``(n, n)``.\n";
            s += "b : ndarray\n    The Cholesky factor of the positive definite `b`, as returned by ``potrf``.\n";
            s += P_ITYPE;
            s += P_LOWER_SH;
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The reduced matrix, overwriting `a`.\n";
            s += R_INFO;
            return s;
        }

        /* ------------------------- triangular storage conversions ------------------------- */

        /* The three layouts a triangle can live in, written once because all six conversions
         * document the same operands from opposite ends. */
        static constexpr const char *P_AP_IN  =
            "ap : ndarray\n"
            "    The triangle in packed storage (TP), length ``n * (n + 1) / 2``.\n";
        static constexpr const char *P_ARF_IN =
            "arf : ndarray\n"
            "    The triangle in rectangular full packed storage (RFP), length\n"
            "    ``n * (n + 1) / 2``.\n";
        static constexpr const char *P_A_TRI_IN =
            "a : ndarray\n"
            "    Square matrix of shape ``(n, n)``; only the triangle named by `uplo` is read.\n";
        static constexpr const char *R_AP_OUT  = "ap : ndarray\n    The triangle in packed storage (TP).\n";
        static constexpr const char *R_ARF_OUT = "arf : ndarray\n    The triangle in rectangular full packed storage (RFP).\n";
        static constexpr const char *R_A_TRI_OUT =
            "a : ndarray\n"
            "    The triangle in full storage, shape ``(n, n)``; the other triangle is zero.\n";
        static constexpr const char *P_UPLO_TRI =
            "uplo : str, optional\n"
            "    ``'U'`` if the matrix is upper triangular, ``'L'`` if lower. Default is ``'U'``.\n";

        /** @brief `transr` selects the RFP layout, and the transposed one is spelled with the
         *         flavor's own transpose letter -- ``'T'`` for real, ``'C'`` for complex. */
        static std::string p_transr(const Dtype &t) noexcept
        {
            return std::string("transr : str, optional\n"
                   "    ``'N'`` for the normal RFP layout, ``'") + (t.is_complex ? "C" : "T") +
                   "'`` for the transposed one.\n    Default is ``'N'``.\n";
        }

        /** @brief The four conversions that take `n` and one packed operand and hand back the
         *         same triangle in another layout. */
        static std::string
        doc_convert_from_packed(const char *name, const Dtype &t, const char *summary,
                                const char *operand, const char *in, const char *out,
                                bool has_transr) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, " + std::string(operand) + ", "
                 + (has_transr ? "transr='N', " : "") + "uplo='U')\n\n";
            /* The layout names make these summaries long, so the LAPACK reference goes on
             * its own line rather than running past the docstring's width. */
            s += std::string(summary) + "\n(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += in;
            if (has_transr) { s += p_transr(t); }
            s += P_UPLO_TRI;

            s += "\nReturns\n-------\n";
            s += out;
            s += R_INFO;
            return s;
        }

        static std::string doc_tpttf(const char *name, const Dtype &t) noexcept
        {
            return doc_convert_from_packed(name, t,
                "Convert a triangle from packed (TP) to rectangular full packed (RFP) storage",
                "ap", P_AP_IN, R_ARF_OUT, true);
        }

        static std::string doc_tpttr(const char *name, const Dtype &t) noexcept
        {
            return doc_convert_from_packed(name, t,
                "Convert a triangle from packed (TP) to full (TR) storage",
                "ap", P_AP_IN, R_A_TRI_OUT, false);
        }

        static std::string doc_tfttp(const char *name, const Dtype &t) noexcept
        {
            return doc_convert_from_packed(name, t,
                "Convert a triangle from rectangular full packed (RFP) to packed (TP) storage",
                "arf", P_ARF_IN, R_AP_OUT, true);
        }

        static std::string doc_tfttr(const char *name, const Dtype &t) noexcept
        {
            return doc_convert_from_packed(name, t,
                "Convert a triangle from rectangular full packed (RFP) to full (TR) storage",
                "arf", P_ARF_IN, R_A_TRI_OUT, true);
        }

        /** @brief The two conversions out of full storage; they take `a` and derive `n` from it,
         *         so they share no signature line with the packed-input four. */
        static std::string
        doc_convert_from_full(const char *name, const Dtype &t, const char *summary,
                              const char *out, bool has_transr) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, " + (has_transr ? "transr='N', " : "") + "uplo='U')\n\n";
            s += std::string(summary) + "\n(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_A_TRI_IN;
            if (has_transr) { s += p_transr(t); }
            s += P_UPLO_TRI;

            s += "\nReturns\n-------\n";
            s += out;
            s += R_INFO;
            return s;
        }

        static std::string doc_trttf(const char *name, const Dtype &t) noexcept
        {
            return doc_convert_from_full(name, t,
                "Convert a triangle from full (TR) to rectangular full packed (RFP) storage",
                R_ARF_OUT, true);
        }

        static std::string doc_trttp(const char *name, const Dtype &t) noexcept
        {
            return doc_convert_from_full(name, t,
                "Convert a triangle from full (TR) to packed (TP) storage",
                R_AP_OUT, false);
        }

        static std::string doc_tfsm(const char *name, const Dtype &t) noexcept
        {
            const char *letter = t.is_complex ? "C" : "T";
            std::string s;
            s += std::string(name) +
                 "(alpha, a, b, transr='N', side='L', uplo='U', trans='N', diag='N',\n"
                 "      overwrite_b=0)\n\n";
            s += "Solve a triangular system whose matrix is in rectangular full packed (RFP)\n"
                 "storage (LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Solves ``op(A) @ x = alpha * b`` when `side` is ``'L'``, or\n"
                 "``x @ op(A) = alpha * b`` when it is ``'R'``.\n\n";

            s += "Parameters\n----------\n";
            s += std::string("alpha : ") + t.scalar + "\n    Scalar the right-hand side is multiplied by.\n";
            s += "a : ndarray\n"
                 "    The triangle of ``A`` in RFP storage. Its order is the number of rows of\n"
                 "    `b` when `side` is ``'L'``, its number of columns when ``'R'``, and its\n"
                 "    length is ``order * (order + 1) / 2``.\n";
            s += "b : ndarray\n    Right-hand side of shape ``(m, n)``.\n";
            s += p_transr(t);
            s += "side : str, optional\n"
                 "    ``'L'`` for ``op(A) @ x``, ``'R'`` for ``x @ op(A)``. Default is ``'L'``.\n";
            s += P_UPLO_TRI;
            s += std::string("trans : str, optional\n    ``'N'`` for ``op(A) = A``, ``'") + letter +
                 "'`` for its " + (t.is_complex ? "conjugate transpose" : "transpose") + ".\n    Default is ``'N'``.\n";
            s += "diag : str, optional\n"
                 "    ``'U'`` if ``A`` has a unit diagonal, ``'N'`` otherwise. Default is ``'N'``.\n";
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    The solution, overwriting `b`.\n";
            return s;
        }

        /* ---------------------- packed positive definite, RFP rank-k ---------------------- */

        /* The packed group takes `n` as its own argument instead of deriving it, so every
         * routine documents the same pair of leading parameters. */
        static constexpr const char *P_AP_PACKED =
            "ap : ndarray\n"
            "    The triangle named by `lower`, in packed storage: the ``n * (n + 1) / 2``\n"
            "    elements of one triangle, column by column. A longer array is accepted and\n"
            "    its tail ignored.\n";
        static constexpr const char *P_AP_FACTOR =
            "ap : ndarray\n"
            "    Cholesky factorization in packed storage, as returned by ``pptrf``.\n";

        static std::string
        doc_ppcon(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, ap, anorm, lower=0)\n\n";
            s += "Estimate the reciprocal condition number of a positive definite matrix in packed\n"
                 "storage, from its Cholesky factorization (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_AP_FACTOR;
            s += "anorm : float\n    1-norm of the original matrix, which equals its infinity norm here.\n";
            s += P_LOWER;

            s += "\nReturns\n-------\n";
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_ppsv(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, ap, b, lower=0, overwrite_b=0)\n\n";
            s += "Solve a positive definite system whose matrix is in packed storage, by Cholesky\n"
                 "factorization (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_AP_PACKED;
            s += "    Note that LAPACK factorizes in place, so `ap` is overwritten with the\n"
                 "    Cholesky factor even though it is declared an input.\n";
            s += P_B_RHS;
            s += P_LOWER;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    Solution of the system, overwriting `b`.\n";
            s += R_INFO;
            return s;
        }

        /** @brief `pptrf` and `pptri` share a signature; only the summary and the result differ. */
        static std::string
        doc_pp_in_place(const char *name, const char *summary, const char *in, const char *out) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, ap, lower=0, overwrite_ap=0)\n\n";
            s += std::string(summary) + " (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += in;
            s += P_LOWER;
            s += "overwrite_ap : int, optional\n"
                 "    If nonzero, `ap` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += out;
            s += R_INFO;
            return s;
        }

        static std::string doc_pptrf(const char *name, const Dtype &) noexcept
        {
            return doc_pp_in_place(name,
                "Compute the Cholesky factorization of a positive definite matrix in packed\nstorage",
                P_AP_PACKED,
                "ul : ndarray\n    The Cholesky factor in packed storage, overwriting `ap`.\n");
        }

        static std::string doc_pptri(const char *name, const Dtype &) noexcept
        {
            return doc_pp_in_place(name,
                "Invert a positive definite matrix in packed storage from its Cholesky\nfactorization",
                P_AP_FACTOR,
                "uli : ndarray\n    The inverse in packed storage, overwriting `ap`.\n");
        }

        static std::string
        doc_pptrs(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, ap, b, lower=0, overwrite_b=0)\n\n";
            s += "Solve a positive definite system in packed storage using a Cholesky factorization\n"
                 "already computed by ``pptrf`` (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_AP_FACTOR;
            s += P_B_RHS;
            s += P_LOWER;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    Solution of the system, overwriting `b`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_sfrk(const char *name, const Dtype &t) noexcept
        {
            const char *letter = t.is_complex ? "C" : "T";
            const char *adj = t.is_complex ? "conjugate transpose" : "transpose";
            std::string s;
            s += std::string(name) + "(n, k, alpha, a, beta, c, transr='N', uplo='U', trans='N',\n"
                 "      overwrite_c=0)\n\n";
            s += std::string(t.is_complex ? "Hermitian" : "Symmetric") +
                 " rank-k update of a matrix held in rectangular full packed (RFP)\n"
                 "storage (LAPACK ``" + std::string(name) + "``).\n\n";
            s += std::string("Computes ``c := alpha * a @ a.") + (t.is_complex ? "conj().T" : "T") +
                 " + beta * c`` when `trans` is ``'N'``,\nor ``c := alpha * a." +
                 (t.is_complex ? "conj().T" : "T") + " @ a + beta * c`` when it is ``'" + letter + "'``.\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += "k : int\n    Inner dimension of the update.\n";
            s += "alpha : float\n"
                 "    Scalar the rank-k term is multiplied by. Real for every flavor.\n";
            s += "a : ndarray\n"
                 "    Shape ``(max(n, 1), k)`` when `trans` is ``'N'``, ``(max(k, 1), n)``\n"
                 "    otherwise. The leading dimension follows `n` and `k` rather than `a`, so\n"
                 "    the row count must match exactly.\n";
            s += "beta : float\n"
                 "    Scalar `c` is scaled by before the update. Real for every flavor.\n";
            s += "c : ndarray\n"
                 "    The matrix to update, in RFP storage, length ``n * (n + 1) / 2``.\n";
            s += p_transr(t);
            s += "uplo : str, optional\n"
                 "    ``'U'`` if `c` holds the upper triangle, ``'L'`` if lower. Default is ``'U'``.\n";
            s += std::string("trans : str, optional\n    ``'N'`` for ``a @ a.") +
                 (t.is_complex ? "conj().T" : "T") + "``, ``'" + letter + "'`` for the " + adj +
                 " product.\n    Default is ``'N'``.\n";
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "cout : ndarray\n    The updated matrix in RFP storage, overwriting `c`.\n";
            return s;
        }

        /* ------------------- rectangular full packed and banded Cholesky ------------------- */

        static constexpr const char *P_AB_BAND =
            "ab : ndarray\n"
            "    The band of the matrix, shape ``(kd + 1, n)``: the order `n` and the bandwidth\n"
            "    `kd` are both read off `ab` itself.\n";
        static constexpr const char *P_AB_FACTOR =
            "ab : ndarray\n"
            "    Banded Cholesky factor, as returned by ``pbtrf``.\n";
        static constexpr const char *P_LDAB =
            "ldab : int, optional\n"
            "    Leading dimension of `ab`. It must equal ``ab.shape[0]``; the argument exists\n"
            "    only because the ``.pyf`` exposed it. Default is ``ab.shape[0]``.\n";

        /** @brief `pftrf` and `pftri` share a signature; only the summary and the result differ. */
        static std::string
        doc_pf_in_place(const char *name, const Dtype &t, const char *summary,
                        const char *in, const char *out) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, a, transr='N', uplo='U', overwrite_a=0)\n\n";
            s += std::string(summary) + "\n(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += in;
            s += p_transr(t);
            s += "uplo : str, optional\n"
                 "    ``'U'`` if `a` holds the upper triangle, ``'L'`` if lower. Default is ``'U'``.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += out;
            s += R_INFO;
            return s;
        }

        static std::string doc_pftrf(const char *name, const Dtype &t) noexcept
        {
            return doc_pf_in_place(name, t,
                "Compute the Cholesky factorization of a positive definite matrix in\n"
                "rectangular full packed (RFP) storage",
                "a : ndarray\n"
                "    The triangle named by `uplo`, in RFP storage, length ``n * (n + 1) / 2``.\n",
                "achol : ndarray\n    The Cholesky factor in RFP storage, overwriting `a`.\n");
        }

        static std::string doc_pftri(const char *name, const Dtype &t) noexcept
        {
            return doc_pf_in_place(name, t,
                "Invert a positive definite matrix in rectangular full packed (RFP) storage\n"
                "from its Cholesky factorization",
                "a : ndarray\n    Cholesky factorization in RFP storage, as returned by ``pftrf``.\n",
                "ainv : ndarray\n    The inverse in RFP storage, overwriting `a`.\n");
        }

        static std::string
        doc_pftrs(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, a, b, transr='N', uplo='U', overwrite_b=0)\n\n";
            s += "Solve a positive definite system in rectangular full packed (RFP) storage using\n"
                 "a Cholesky factorization already computed by ``pftrf``\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += "a : ndarray\n    Cholesky factorization in RFP storage, as returned by ``pftrf``.\n";
            s += "b : ndarray\n"
                 "    Right-hand side(s). It may have more than `n` rows; only the leading `n`\n"
                 "    are read.\n";
            s += p_transr(t);
            s += "uplo : str, optional\n"
                 "    ``'U'`` if `a` holds the upper triangle, ``'L'`` if lower. Default is ``'U'``.\n";
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    Solution of the system, overwriting `b`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_pbtrf(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(ab, lower=0, ldab=ab.shape[0], overwrite_ab=0)\n\n";
            s += "Compute the Cholesky factorization of a positive definite band matrix\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_AB_BAND;
            s += P_LOWER;
            s += P_LDAB;
            s += "overwrite_ab : int, optional\n"
                 "    If nonzero, `ab` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The banded Cholesky factor, overwriting `ab`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_pbtrs(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(ab, b, lower=0, ldab=ab.shape[0], overwrite_b=0)\n\n";
            s += "Solve a positive definite band system using a Cholesky factorization already\n"
                 "computed by ``pbtrf`` (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_AB_FACTOR;
            s += P_B_RHS;
            s += P_LOWER;
            s += P_LDAB;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    Solution of the system, overwriting `b`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_pbsv(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(ab, b, lower=0, ldab=ab.shape[0], overwrite_ab=0,\n"
                 "      overwrite_b=0)\n\n";
            s += "Solve a positive definite band system by Cholesky factorization\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_AB_BAND;
            s += P_B_RHS;
            s += P_LOWER;
            s += P_LDAB;
            s += "overwrite_ab : int, optional\n"
                 "    If nonzero, `ab` may be overwritten in place. Default is 0.\n";
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The banded Cholesky factor, overwriting `ab`.\n";
            s += "x : ndarray\n    Solution of the system, overwriting `b`.\n";
            s += R_INFO;
            return s;
        }

        /* --------------- triangular solves, inverse, condition, LU auxiliaries -------------- */

        static constexpr const char *P_UNITDIAG =
            "unitdiag : int, optional\n"
            "    If nonzero, the matrix is treated as having a unit diagonal and the stored\n"
            "    diagonal is not referenced. Default is 0.\n";
        static constexpr const char *P_DIAG_LETTER =
            "diag : str, optional\n"
            "    ``'U'`` if the matrix has a unit diagonal, ``'N'`` otherwise. Default is\n"
            "    ``'N'``.\n";

        static std::string
        doc_trtrs(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, lower=0, trans=0, unitdiag=0, lda=a.shape[0],\n"
                 "      overwrite_b=0)\n\n";
            s += "Solve a triangular system (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Triangular matrix; only the triangle named by `lower` is read.\n";
            s += P_B_RHS;
            s += "lower : int, optional\n"
                 "    If nonzero, `a` is lower triangular; otherwise upper. Default is 0.\n";
            s += P_TRANS_INT;
            s += P_UNITDIAG;
            s += "lda : int, optional\n"
                 "    Leading dimension of `a`. It must equal ``a.shape[0]``; the argument exists\n"
                 "    only because the ``.pyf`` exposed it. Default is ``a.shape[0]``.\n";
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    Solution of the system, overwriting `b`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_trcon(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, norm='1', uplo='U', diag='N')\n\n";
            s += "Estimate the reciprocal condition number of a triangular matrix\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n"
                 "    Square triangular matrix; only the triangle named by `uplo` is read.\n";
            s += "norm : str, optional\n"
                 "    ``'1'`` or ``'O'`` for the 1-norm, ``'I'`` for the infinity norm.\n"
                 "    Default is ``'1'``.\n";
            s += "uplo : str, optional\n"
                 "    ``'U'`` if `a` is upper triangular, ``'L'`` if lower. Default is ``'U'``.\n";
            s += P_DIAG_LETTER;

            s += "\nReturns\n-------\n";
            s += "rcond : float\n    Estimate of the reciprocal condition number.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_tbtrs(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(ab, b, uplo='U', trans='N', diag='N', overwrite_b=0)\n\n";
            s += "Solve a triangular band system (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "ab : ndarray\n"
                 "    The band of the triangular matrix, shape ``(kd + 1, n)``: the order `n` and\n"
                 "    the bandwidth `kd` are both read off `ab` itself.\n";
            s += P_B_RHS;
            s += "uplo : str, optional\n"
                 "    ``'U'`` if the matrix is upper triangular, ``'L'`` if lower. Default is\n"
                 "    ``'U'``.\n";
            s += P_TRANS_STR;
            s += P_DIAG_LETTER;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    Solution of the system, overwriting `b`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_trtri(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(c, lower=0, unitdiag=0, overwrite_c=0)\n\n";
            s += "Invert a triangular matrix (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "c : ndarray\n"
                 "    Square triangular matrix; only the triangle named by `lower` is read.\n";
            s += "lower : int, optional\n"
                 "    If nonzero, `c` is lower triangular; otherwise upper. Default is 0.\n";
            s += P_UNITDIAG;
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "inv_c : ndarray\n    The inverse, overwriting `c`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_lauum(const char *name, const Dtype &t) noexcept
        {
            const char *adjoint = t.is_complex ? "conj().T" : "T";
            std::string s;
            s += std::string(name) + "(c, lower=0, overwrite_c=0)\n\n";
            s += "Multiply a triangular matrix by its own adjoint\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            /* The `.pyf`'s own comment has these the wrong way round; LAPACK's is authoritative. */
            s += std::string("Computes ``u @ u.") + adjoint + "`` for an upper triangular `c`, or "
                 "``l." + adjoint + " @ l``\nfor a lower triangular one, writing the result into "
                 "the same triangle.\n\n";

            s += "Parameters\n----------\n";
            s += "c : ndarray\n"
                 "    Square triangular matrix; only the triangle named by `lower` is read.\n";
            if (t.is_complex) {
                /* The complex routines take `REAL(A(I,I))`: this exists to square a Cholesky
                 * factor, whose diagonal is real, and a complex diagonal is silently ignored. */
                s += "    The routine assumes a Cholesky factor and reads only the real part of\n"
                     "    the diagonal; any imaginary part is ignored rather than rejected.\n";
            }
            s += "lower : int, optional\n"
                 "    If nonzero, `c` is lower triangular; otherwise upper. Default is 0.\n";
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The product, in the same triangle of `c`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_laswp(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, piv, k1=0, k2=len(piv)-1, off=0, inc=1,\n"
                 "      overwrite_a=0)\n\n";
            s += "Apply a sequence of row interchanges to a matrix\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Row ``i`` is swapped with row ``piv[i]`` for each ``i`` from `k1` to `k2`.\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Matrix whose rows are interchanged.\n";
            s += "piv : ndarray\n"
                 "    Pivot indices, 0-based, at most ``a.shape[0]`` of them. LAPACK numbers\n"
                 "    these from 1; the shift happens inside the wrapper and the array passed in\n"
                 "    is not modified.\n";
            s += "k1 : int, optional\n    First index of `piv` to apply, 0-based. Default is 0.\n";
            s += "k2 : int, optional\n"
                 "    Last index of `piv` to apply, 0-based. Default is ``len(piv) - 1``.\n";
            s += "off : int, optional\n"
                 "    Number of leading entries of `piv` to skip. Default is 0.\n";
            s += "inc : int, optional\n"
                 "    Stride through `piv`; negative applies the interchanges in reverse order.\n"
                 "    ``piv`` must hold ``k1 + 1 + (k2 - k1) * abs(inc) + off`` entries, which\n"
                 "    for ``abs(inc) > 1`` is more than `k2` alone requires. Default is 1.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The matrix with its rows interchanged.\n";
            return s;
        }

        /* --------------------- Schur exchange and the Sylvester equations ------------------- */

        static constexpr const char *P_IFST_ILST =
            "ifst : int\n"
            "    Current position of the diagonal block to move, 1-based.\n"
            "ilst : int\n"
            "    Position to move it to, 1-based. LAPACK adjusts both when a block turns out to\n"
            "    be 2x2; the adjusted values are not returned.\n";

        static std::string
        doc_trexc(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, q, ifst, ilst, wantq=1, overwrite_a=0, overwrite_q=0)\n\n";
            s += "Reorder a Schur decomposition by moving one diagonal block\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Upper (quasi-)triangular Schur form of shape ``(n, n)``.\n";
            s += "q : ndarray\n    Matrix of Schur vectors; updated when `wantq` is nonzero.\n";
            s += P_IFST_ILST;
            s += "wantq : int, optional\n"
                 "    If nonzero, the transformation is accumulated into `q`. Default is 1.\n";
            s += P_OVERWRITE_A;
            s += "overwrite_q : int, optional\n"
                 "    If nonzero, `q` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The reordered Schur form.\n";
            s += "q : ndarray\n    The updated Schur vectors.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_tgexc(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, q, z, ifst, ilst, wantq=1, wantz=1, ";
            /* Only the real flavors carry a workspace, so only they take `lwork`. */
            s += t.is_complex ? "overwrite_a=0,\n      "
                              : "lwork=4*n+16,\n      overwrite_a=0, ";
            s += "overwrite_b=0, overwrite_q=0, overwrite_z=0)\n\n";
            s += "Reorder a generalized Schur decomposition by moving one diagonal block\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Upper (quasi-)triangular factor of the pencil, shape ``(n, n)``.\n";
            s += "b : ndarray\n    Upper triangular factor of the pencil, shape ``(n, n)``.\n";
            s += "q : ndarray\n    Left Schur vectors; updated when `wantq` is nonzero.\n";
            s += "z : ndarray\n    Right Schur vectors; updated when `wantz` is nonzero.\n";
            s += P_IFST_ILST;
            s += "wantq : int, optional\n"
                 "    If nonzero, the left transformation is accumulated into `q`. Default is 1.\n";
            s += "wantz : int, optional\n"
                 "    If nonzero, the right transformation is accumulated into `z`. Default is 1.\n";
            if (!t.is_complex) {
                s += "lwork : int, optional\n"
                     "    Workspace size, at least ``4 * n + 16``; -1 requests the query.\n"
                     "    Default is ``4 * n + 16``. The complex flavors need no workspace and\n"
                     "    take no `lwork`.\n";
            }
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;
            s += "overwrite_q : int, optional\n"
                 "    If nonzero, `q` may be overwritten in place. Default is 0.\n";
            s += "overwrite_z : int, optional\n"
                 "    If nonzero, `z` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The reordered first factor.\n";
            s += "b : ndarray\n    The reordered second factor.\n";
            s += "q : ndarray\n    The updated left Schur vectors.\n";
            s += "z : ndarray\n    The updated right Schur vectors.\n";
            if (!t.is_complex) {
                s += "work : ndarray\n"
                     "    Workspace; ``work[0]`` holds the optimal `lwork`. The complex flavors\n"
                     "    do not return it.\n";
            }
            s += R_INFO;
            return s;
        }

        static std::string
        doc_trsyl(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, c, trana='N', tranb='N', isgn=1, overwrite_c=0)\n\n";
            s += "Solve the Sylvester equation ``op(a) @ x + isgn * x @ op(b) = scale * c``\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n"
                 "    Upper (quasi-)triangular matrix of shape ``(m, m)``, in Schur canonical\n"
                 "    form.\n";
            s += "b : ndarray\n"
                 "    Upper (quasi-)triangular matrix of shape ``(n, n)``, in Schur canonical\n"
                 "    form.\n";
            s += "c : ndarray\n    Right-hand side of shape ``(m, n)``.\n";
            s += "trana : str, optional\n"
                 "    ``'N'``, ``'T'`` or ``'C'`` for ``op(a)``. Default is ``'N'``.\n";
            s += "tranb : str, optional\n"
                 "    ``'N'``, ``'T'`` or ``'C'`` for ``op(b)``. Default is ``'N'``.\n";
            s += "isgn : int, optional\n"
                 "    Sign of the ``x @ op(b)`` term: 1 or -1. Default is 1.\n";
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    The solution, overwriting `c`.\n";
            s += "scale : float\n"
                 "    Scale factor, at most 1, chosen to avoid overflow in the solution.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    1 if `a` and `b` have common or very close eigenvalues, in which case\n"
                 "    perturbed values were used.\n";
            return s;
        }

        /* ---------------------------- Schur cluster reordering ----------------------------- */

        static constexpr const char *P_SELECT_ARRAY =
            "select : ndarray\n"
            "    Boolean mask of length ``n``: the eigenvalues it marks are gathered into the\n"
            "    leading block. A bool array works; LAPACK reads it as a Fortran LOGICAL.\n";
        static constexpr const char *P_JOB_SEN =
            "job : str, optional\n"
            "    What to compute besides the reordering: ``'N'`` nothing, ``'E'`` the cluster\n"
            "    condition number `s`, ``'V'`` the separation `sep`, ``'B'`` both.\n"
            "    Default is ``'B'``.\n";
        static constexpr const char *P_IJOB_SEN =
            "ijob : int, optional\n"
            "    What to compute besides the reordering, 0 through 5; 0 reorders only.\n"
            "    Default is 4.\n";
        static constexpr const char *R_M_CLUSTER =
            "m : int\n    Size of the selected cluster.\n";
        /* The whole point of the `_lwork` siblings: the workspace these routines need depends on
         * `m`, which is not known until `select` has been scanned. */
        static constexpr const char *P_LWORK_SEN_NOTE =
            "    The default is a bare minimum for the cheapest `job` only -- what is actually\n"
            "    needed depends on the cluster size, so query the matching ``_lwork`` routine.\n";

        static std::string
        doc_trsen(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(select, t, q, job='B', wantq=1, lwork=max(1,n), ";
            s += t.is_complex ? "\n      " : "liwork=1,\n      ";
            s += "overwrite_t=0, overwrite_q=0)\n\n";
            s += "Reorder a Schur decomposition to gather selected eigenvalues into a leading\n"
                 "cluster (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_SELECT_ARRAY;
            s += "t : ndarray\n    Upper (quasi-)triangular Schur form of shape ``(n, n)``.\n";
            s += "q : ndarray\n    Matrix of Schur vectors; updated when `wantq` is nonzero.\n";
            s += P_JOB_SEN;
            s += "wantq : int, optional\n"
                 "    If nonzero, the transformation is accumulated into `q`. Default is 1.\n";
            s += "lwork : int, optional\n"
                 "    Workspace size; -1 requests the query. Default is ``max(1, n)``.\n";
            s += P_LWORK_SEN_NOTE;
            if (!t.is_complex) {
                s += "liwork : int, optional\n"
                     "    Integer workspace size; -1 requests the query. Default is 1. The\n"
                     "    complex flavors have no integer workspace and take no `liwork`.\n";
            }
            s += "overwrite_t : int, optional\n"
                 "    If nonzero, `t` may be overwritten in place. Default is 0.\n";
            s += "overwrite_q : int, optional\n"
                 "    If nonzero, `q` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "ts : ndarray\n    The reordered Schur form.\n";
            s += "qs : ndarray\n    The updated Schur vectors.\n";
            if (t.is_complex) {
                s += "w : ndarray\n    Reordered eigenvalues, length ``n``.\n";
            }
            else {
                s += "wr : ndarray\n    Real parts of the reordered eigenvalues, length ``n``.\n";
                s += "wi : ndarray\n    Imaginary parts of the reordered eigenvalues, length ``n``.\n";
            }
            s += R_M_CLUSTER;
            s += "s : float\n    Condition number of the selected cluster; 0 unless `job` asked for it.\n";
            s += "sep : float\n    Separation of the cluster; 0 unless `job` asked for it.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_trsen_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(select, t, job='B')\n\n";
            s += "Query the workspace ``trsen`` needs for this `select` and `job`\n"
                 "(LAPACK ``" + std::string(name).substr(0, std::strlen(name) - 6) + "`` with "
                 "``lwork = -1``).\n\n";
            s += "`select` and `t` are read for real: LAPACK has to scan them to work out the\n"
                 "cluster size before it can size the workspace.\n\n";

            s += "Parameters\n----------\n";
            s += P_SELECT_ARRAY;
            s += "t : ndarray\n    Upper (quasi-)triangular Schur form of shape ``(n, n)``.\n";
            s += P_JOB_SEN;

            s += "\nReturns\n-------\n";
            s += std::string("work : ") + t.scalar + "\n    Optimal `lwork` for the matching ``trsen`` call.\n";
            if (!t.is_complex) {
                s += "iwork : int\n"
                     "    Optimal `liwork`. The complex flavors have no integer workspace and do\n"
                     "    not return it.\n";
            }
            s += R_INFO;
            return s;
        }

        static std::string
        doc_tgsen(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(select, a, b, q, z, ijob=4, wantq=1, wantz=1,\n"
                 "      lwork=..., liwork=..., overwrite_a=0, overwrite_b=0, overwrite_q=0,\n"
                 "      overwrite_z=0)\n\n";
            s += "Reorder a generalized Schur decomposition to gather selected eigenvalues into a\n"
                 "leading cluster (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_SELECT_ARRAY;
            s += "a : ndarray\n    Upper (quasi-)triangular factor of the pencil, shape ``(n, n)``.\n";
            s += "b : ndarray\n    Upper triangular factor of the pencil, shape ``(n, n)``.\n";
            s += "q : ndarray\n    Left Schur vectors; updated when `wantq` is nonzero.\n";
            s += "z : ndarray\n    Right Schur vectors; updated when `wantz` is nonzero.\n";
            s += P_IJOB_SEN;
            s += "wantq : int, optional\n"
                 "    If nonzero, the left transformation is accumulated into `q`. Default is 1.\n";
            s += "wantz : int, optional\n"
                 "    If nonzero, the right transformation is accumulated into `z`. Default is 1.\n";
            s += "lwork : int, optional\n    Workspace size; -1 requests the query.\n";
            s += t.is_complex ? "    Default is ``1`` when `ijob` is 0 and ``n + 2`` otherwise.\n"
                              : "    Default is ``4 * n + 16``.\n";
            s += P_LWORK_SEN_NOTE;
            s += "liwork : int, optional\n    Integer workspace size; -1 requests the query.\n";
            s += t.is_complex ? "    Default is ``1`` when `ijob` is 0 and ``n + 2`` otherwise.\n"
                              : "    Default is ``n + 6``.\n";
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;
            s += "overwrite_q : int, optional\n"
                 "    If nonzero, `q` may be overwritten in place. Default is 0.\n";
            s += "overwrite_z : int, optional\n"
                 "    If nonzero, `z` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "as : ndarray\n    The reordered first factor.\n";
            s += "bs : ndarray\n    The reordered second factor.\n";
            if (t.is_complex) {
                s += "alpha : ndarray\n    Numerators of the reordered eigenvalues, length ``n``.\n";
            }
            else {
                s += "alphar : ndarray\n    Real parts of the eigenvalue numerators, length ``n``.\n";
                s += "alphai : ndarray\n    Imaginary parts of the eigenvalue numerators, length ``n``.\n";
            }
            s += "beta : ndarray\n    Denominators of the reordered eigenvalues, length ``n``.\n";
            s += "qs : ndarray\n    The updated left Schur vectors.\n";
            s += "zs : ndarray\n    The updated right Schur vectors.\n";
            s += R_M_CLUSTER;
            s += "pl : float\n    Lower bound on the reciprocal condition number of the left cluster.\n";
            s += "pr : float\n    Lower bound on the reciprocal condition number of the right cluster.\n";
            s += "dif : ndarray\n    Two estimates of the Difu and Difl separations.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_tgsen_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + (t.is_complex ? "(select, a, b, ijob=4)\n\n"
                                                   : "(select, a, ijob=4)\n\n");
            s += "Query the workspace ``tgsen`` needs for this `select` and `ijob`\n"
                 "(LAPACK ``" + std::string(name).substr(0, std::strlen(name) - 6) + "`` with "
                 "``lwork = -1``).\n\n";

            s += "Parameters\n----------\n";
            s += P_SELECT_ARRAY;
            s += "a : ndarray\n    Upper (quasi-)triangular factor of the pencil, shape ``(n, n)``.\n";
            if (t.is_complex) {
                s += "b : ndarray\n"
                     "    Upper triangular factor of the pencil. The complex query reads both\n"
                     "    factors where the real one needs only `a`.\n";
            }
            s += P_IJOB_SEN;

            s += "\nReturns\n-------\n";
            s += std::string("work : ") + t.scalar + "\n    Optimal `lwork` for the matching ``tgsen`` call.\n";
            s += "iwork : int\n    Optimal `liwork`.\n";
            s += R_INFO;
            return s;
        }

        /* --------------- generating and applying the orthogonal/unitary factor -------------- */

        static constexpr const char *P_TAU_IN =
            "tau : ndarray\n    Reflector scalars from the matching factorization.\n";

        /** @brief The transpose letter these routines take: real flavors ``'T'``, complex
         *         ``'C'``, and no flavor takes the other's. */
        static std::string p_trans_qr(const Dtype &t, bool optional) noexcept
        {
            std::string s = "trans : str";
            s += optional ? ", optional\n" : "\n";
            s += std::string("    ``'N'`` to apply ``Q``, ``'") + (t.is_complex ? "C" : "T") +
                 "'`` to apply its " + (t.is_complex ? "conjugate transpose" : "transpose") + ".\n";
            if (optional) { s += "    Default is ``'N'``.\n"; }
            return s;
        }

        static std::string
        doc_orghr(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, tau, lo=0, hi=n-1, lwork=max(hi-lo,1), overwrite_a=0)\n\n";
            s += "Generate the orthogonal/unitary matrix of a Hessenberg reduction\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    The reduced matrix ``gehrd`` returned, shape ``(n, n)``.\n";
            s += "tau : ndarray\n"
                 "    Reflector scalars from ``gehrd``, length ``n - 1`` -- empty when ``n`` is 1,\n"
                 "    which is what ``gehrd`` returns there.\n";
            s += P_LO;
            s += P_HI;
            s += "lwork : int, optional\n"
                 "    Workspace size, at least ``hi - lo``. Default is ``max(hi - lo, 1)``.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "ht : ndarray\n    The orthogonal/unitary factor ``Q``.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_orghr_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, lo=0, hi=n-1)\n\n";
            s += "Query the optimal workspace for ``" + std::string(name).substr(0, std::strlen(name) - 6) + "``\n"
                 "(LAPACK ``" + std::string(name).substr(0, std::strlen(name) - 6) + "`` with ``lwork = -1``).\n\n";

            s += "Parameters\n----------\n";
            s += P_N_ORDER;
            s += P_LO;
            s += P_HI;

            s += "\nReturns\n-------\n";
            s += std::string("work : ") + t.scalar + "\n    Optimal `lwork`.\n";
            s += R_INFO;
            return s;
        }

        /** @brief `orgqr` and `orgrq` differ only in which side the reflectors act from, which
         *         shows up as the dimension that floors `lwork`. */
        static std::string
        doc_org_factor(const char *name, const char *which, const char *producer,
                       const char *bound) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, tau, lwork=3*" + std::string(bound) + ", overwrite_a=0)\n\n";
            s += "Generate the orthogonal/unitary factor ``Q`` of a " + std::string(which) +
                 " factorization\n(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    The factorization the matching ``" + std::string(producer) +
                 "`` returned.\n";
            s += P_TAU_IN;
            s += "lwork : int, optional\n    Workspace size, at least ``" + std::string(bound) +
                 "``; -1 requests the query.\n    Default is ``3 * " + std::string(bound) + "``.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "q : ndarray\n    The factor ``Q``, overwriting `a`.\n";
            s += R_WORK_OUT;
            s += R_INFO;
            return s;
        }

        static std::string doc_orgqr(const char *name, const Dtype &) noexcept
        { return doc_org_factor(name, "QR", "geqrf", "n"); }

        static std::string doc_orgrq(const char *name, const Dtype &) noexcept
        { return doc_org_factor(name, "RQ", "gerqf", "m"); }

        static std::string
        doc_ormqr(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(side, trans, a, tau, c, lwork, overwrite_c=0)\n\n";
            s += "Multiply a matrix by the ``Q`` of a QR factorization\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Computes ``Q @ c`` or ``c @ Q`` according to `side`, with `trans` selecting\n"
                 "``Q`` or its adjoint.\n\n";

            s += "Parameters\n----------\n";
            s += "side : str\n    ``'L'`` for ``Q @ c``, ``'R'`` for ``c @ Q``.\n";
            s += p_trans_qr(t, false);
            s += "a : ndarray\n"
                 "    The factorization ``geqrf`` returned. Its column count is taken as the\n"
                 "    reflector count, so `a` must be tall and `tau` must match that count.\n";
            s += P_TAU_IN;
            s += "c : ndarray\n    Matrix to multiply.\n";
            s += "lwork : int\n    Workspace size; -1 requests the query. Required, with no default.\n";
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "cq : ndarray\n    The product, overwriting `c`.\n";
            s += R_WORK_OUT;
            s += R_INFO;
            return s;
        }

        static std::string
        doc_ormrz(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, tau, c, side='L', trans='N', lwork=..., overwrite_c=0)\n\n";
            s += "Multiply a matrix by the ``Q`` of an RZ factorization\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "``Q`` is the factor ``tzrzf`` produced when reducing a trapezoidal matrix to\n"
                 "triangular form.\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n"
                 "    The factorization ``tzrzf`` returned, shape ``(k, nt)`` with ``nt >= k``.\n"
                 "    Its trailing ``nt - k`` columns carry the ``Z`` part.\n";
            s += "tau : ndarray\n    Reflector scalars from ``tzrzf``, length ``k``.\n";
            s += "c : ndarray\n    Matrix to multiply, shape ``(m, n)``.\n";
            s += "side : str, optional\n"
                 "    ``'L'`` for ``Q @ c``, ``'R'`` for ``c @ Q``. Default is ``'L'``.\n";
            s += p_trans_qr(t, true);
            s += "lwork : int, optional\n"
                 "    Workspace size, at least `n` when `side` is ``'L'`` and `m` when ``'R'``;\n"
                 "    -1 requests the query. Default is that minimum.\n";
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "cq : ndarray\n    The product, overwriting `c`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_ormrz_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(m, n, side='L', trans='N')\n\n";
            s += "Query the optimal workspace for ``" + std::string(name).substr(0, std::strlen(name) - 6) + "``\n"
                 "(LAPACK ``" + std::string(name).substr(0, std::strlen(name) - 6) + "`` with ``lwork = -1``).\n\n";

            s += "Parameters\n----------\n";
            s += P_M;
            s += P_N;
            s += "side : str, optional\n"
                 "    ``'L'`` for ``Q @ c``, ``'R'`` for ``c @ Q``. Default is ``'L'``.\n";
            s += p_trans_qr(t, true);

            s += "\nReturns\n-------\n";
            s += std::string("work : ") + t.scalar + "\n    Optimal `lwork`.\n";
            s += R_INFO;
            return s;
        }

        /* --------------------------- compact-WY and blocked QR ----------------------------- */

        static constexpr const char *P_NB_BLOCK =
            "nb : int\n"
            "    Block size, at least 1 and no larger than the factorization allows.\n";
        static constexpr const char *P_T_REFLECTORS =
            "t : ndarray\n"
            "    Upper triangular block reflectors, shape ``(nb, k)``. Its row count is taken as\n"
            "    the block size and its column count as the reflector count.\n";
        static constexpr const char *R_T_REFLECTORS =
            "t : ndarray\n    The upper triangular block reflectors.\n";

        static std::string
        doc_geqrt(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(nb, a, overwrite_a=0)\n\n";
            s += "Blocked QR factorization in the compact WY representation\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "``Q`` is described by the reflectors left in `a` together with the block\n"
                 "reflectors `t`, and applied with ``gemqrt``.\n\n";

            s += "Parameters\n----------\n";
            s += P_NB_BLOCK;
            s += P_A_GENERAL;
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    ``R`` in the upper triangle; below it, the reflectors.\n";
            s += R_T_REFLECTORS;
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gemqrt(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(v, t, c, side='L', trans='N', overwrite_c=0)\n\n";
            s += "Multiply a matrix by the ``Q`` of a compact-WY QR factorization\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "v : ndarray\n    The reflectors ``geqrt`` left in its `a`.\n";
            s += P_T_REFLECTORS;
            s += "c : ndarray\n    Matrix to multiply, shape ``(m, n)``.\n";
            s += "side : str, optional\n"
                 "    ``'L'`` for ``Q @ c``, ``'R'`` for ``c @ Q``. Default is ``'L'``.\n";
            s += p_trans_qr(t, true);
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The product, overwriting `c`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_tpqrt(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(l, nb, a, b, overwrite_a=0, overwrite_b=0)\n\n";
            s += "Blocked QR factorization of a triangular-pentagonal matrix\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "The operand is the square upper triangular `a` stacked on the pentagonal `b`.\n\n";

            s += "Parameters\n----------\n";
            s += "l : int\n"
                 "    Order of the trapezoidal part of `b`; 0 makes `b` rectangular.\n";
            s += P_NB_BLOCK;
            s += "a : ndarray\n    Square upper triangular block, shape ``(n, n)``.\n";
            s += "b : ndarray\n    Pentagonal block, shape ``(m, n)``.\n";
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The triangular block of the factorization.\n";
            s += "b : ndarray\n    The reflectors, overwriting `b`.\n";
            s += R_T_REFLECTORS;
            s += R_INFO;
            return s;
        }

        static std::string
        doc_tpmqrt(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(l, v, t, a, b, side='L', trans='N', overwrite_a=0,\n"
                 "      overwrite_b=0)\n\n";
            s += "Multiply a matrix by the ``Q`` of a triangular-pentagonal QR factorization\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "The matrix multiplied is the pair ``(a, b)``, the same split ``tpqrt`` uses.\n\n";

            s += "Parameters\n----------\n";
            s += "l : int\n    Order of the trapezoidal part of `v`; 0 makes it rectangular.\n";
            s += "v : ndarray\n    The reflectors ``tpqrt`` left in its `b`.\n";
            s += P_T_REFLECTORS;
            s += "a : ndarray\n"
                 "    Leading block; shape ``(k, n)`` when `side` is ``'L'`` and ``(m, k)``\n"
                 "    when ``'R'``.\n";
            s += "b : ndarray\n    Trailing block, shape ``(m, n)``.\n";
            s += "side : str, optional\n"
                 "    ``'L'`` for ``Q @ c``, ``'R'`` for ``c @ Q``. Default is ``'L'``.\n";
            s += p_trans_qr(t, true);
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;

            s += "\nReturns\n-------\n";
            s += "a : ndarray\n    The leading block of the product.\n";
            s += "b : ndarray\n    The trailing block of the product.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_tzrzf(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, lwork=max(m,1), overwrite_a=0)\n\n";
            s += "Reduce an upper trapezoidal matrix to upper triangular form\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Factors ``a`` as ``(R 0) @ Z`` with ``Z`` orthogonal/unitary; ``ormrz`` applies\n"
                 "the ``Z`` this produces.\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n"
                 "    Upper trapezoidal matrix of shape ``(m, n)`` with ``n >= m``.\n";
            s += "lwork : int, optional\n"
                 "    Workspace size, at least `m`. Default is ``max(m, 1)``.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "rz : ndarray\n    ``R`` in the leading triangle; the rest describes ``Z``.\n";
            s += "tau : ndarray\n    Reflector scalars, length ``m``.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_tzrzf_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(m, n)\n\n";
            s += "Query the optimal workspace for ``" + std::string(name).substr(0, std::strlen(name) - 6) + "``\n"
                 "(LAPACK ``" + std::string(name).substr(0, std::strlen(name) - 6) + "`` with ``lwork = -1``).\n\n";

            s += "Parameters\n----------\n";
            s += P_M;
            s += P_N;

            s += "\nReturns\n-------\n";
            s += std::string("work : ") + t.scalar + "\n    Optimal `lwork`.\n";
            s += R_INFO;
            return s;
        }

        /* -------------------- CS decomposition and the one-off solvers --------------------- */

        static std::string
        doc_orcsd(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(x11, x12, x21, x22, compute_u1=1, compute_u2=1,\n"
                 "      compute_v1t=1, compute_v2t=1, trans=0, signs=0, lwork=...,";
            s += t.is_complex ? " lrwork=...,\n      " : "\n      ";
            s += "overwrite_x11=0, overwrite_x12=0, overwrite_x21=0, overwrite_x22=0)\n\n";
            s += std::string("CS decomposition of a partitioned ") +
                 (t.is_complex ? "unitary" : "orthogonal") + " matrix\n(LAPACK ``" +
                 std::string(name) + "``).\n\n";
            s += "The matrix arrives as its four blocks and comes back transformed, alongside\n"
                 "the angles `theta` and the four factors.\n\n";

            s += "Parameters\n----------\n";
            s += "x11 : ndarray\n    Leading block, shape ``(p, q)``; it sets `p` and `q`.\n";
            s += "x12 : ndarray\n    Upper right block, shape ``(p, m - q)``.\n";
            s += "x21 : ndarray\n    Lower left block, shape ``(m - p, q)``.\n";
            s += "x22 : ndarray\n    Trailing block, shape ``(m - p, m - q)``; with `x11` it sets `m`.\n";
            s += "compute_u1 : int, optional\n    If nonzero, `u1` is computed. Default is 1.\n";
            s += "compute_u2 : int, optional\n    If nonzero, `u2` is computed. Default is 1.\n";
            s += "compute_v1t : int, optional\n    If nonzero, `v1t` is computed. Default is 1.\n";
            s += "compute_v2t : int, optional\n    If nonzero, `v2t` is computed. Default is 1.\n";
            s += "trans : int, optional\n"
                 "    If nonzero, the blocks are read row-major rather than column-major.\n"
                 "    Default is 0.\n";
            s += "signs : int, optional\n"
                 "    If nonzero, the alternate sign convention is used. Default is 0.\n";
            s += "lwork : int, optional\n"
                 "    Workspace size; -1 requests the query. Query it with ``" +
                 std::string(name) + "_lwork``.\n";
            if (t.is_complex) {
                s += "lrwork : int, optional\n"
                     "    Real workspace size; -1 requests the query. The real flavors have no\n"
                     "    second workspace and take no `lrwork`.\n";
            }
            s += "overwrite_x11, overwrite_x12, overwrite_x21, overwrite_x22 : int, optional\n"
                 "    If nonzero, the matching block may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "cs11, cs12, cs21, cs22 : ndarray\n    The transformed blocks.\n";
            s += "theta : ndarray\n"
                 "    The CS angles, length ``min(p, m - p, q, m - q)``. Real for every flavor.\n";
            s += "u1, u2, v1t, v2t : ndarray\n"
                 "    The factors, each empty unless its ``compute_*`` flag was set.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_orcsd_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(m, p, q)\n\n";
            s += "Query the workspace ``" + std::string(name).substr(0, std::strlen(name) - 6) + "`` needs\n"
                 "(LAPACK ``" + std::string(name).substr(0, std::strlen(name) - 6) + "`` with ``lwork = -1``).\n\n";
            s += "The query runs with every factor requested, no transpose and the default sign\n"
                 "convention, which is what the ``.pyf`` hard-codes.\n\n";

            s += "Parameters\n----------\n";
            s += "m : int\n    Order of the partitioned matrix.\n";
            s += "p : int\n    Row count of the leading block.\n";
            s += "q : int\n    Column count of the leading block.\n";

            s += "\nReturns\n-------\n";
            s += std::string("work : ") + t.scalar + "\n    Optimal `lwork`.\n";
            if (t.is_complex) {
                s += "rwork : float\n"
                     "    Optimal `lrwork`. The real flavors have no second workspace and do not\n"
                     "    return it.\n";
            }
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gejsv(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, joba=4, jobu=0, jobv=0, jobr=1, jobt=0, jobp=1,\n"
                 "      lwork=..., overwrite_a=0)\n\n";
            s += "Singular value decomposition by the Jacobi method, for high relative accuracy\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "LAPACK ships no complex counterpart, so only the real flavors exist.\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Matrix of shape ``(m, n)`` with ``m >= n``.\n";
            s += "joba : int, optional\n"
                 "    Accuracy level, an index into ``'CEFGAR'``: 0 through 5. Default is 4.\n";
            s += "jobu : int, optional\n"
                 "    Left singular vectors, an index into ``'UFWN'``: 0 through 3. Default is 0.\n";
            s += "jobv : int, optional\n"
                 "    Right singular vectors, an index into ``'VJWN'``: 0 through 3. `jobv` of 1\n"
                 "    requires `jobu` below 2. Default is 0.\n";
            s += "jobr : int, optional\n"
                 "    If nonzero, extreme columns are pruned. Default is 1.\n";
            s += "jobt : int, optional\n"
                 "    If nonzero, the transpose may be factored instead. Default is 0.\n";
            s += "jobp : int, optional\n"
                 "    If nonzero, denormal values are flushed. Default is 1.\n";
            s += "lwork : int, optional\n"
                 "    Workspace size, at least 7. Default is the largest of the several bounds\n"
                 "    the routine documents.\n";
            s += P_OVERWRITE_A;

            s += "\nReturns\n-------\n";
            s += "sva : ndarray\n    Singular values, length ``n``.\n";
            s += "u : ndarray\n    Left singular vectors; empty when `jobu` and `jobt` ask for none.\n";
            s += "v : ndarray\n    Right singular vectors; empty when `jobv` and `jobt` ask for none.\n";
            s += "workout : ndarray\n"
                 "    The first seven workspace entries, which carry the scaling diagnostics.\n";
            s += "iworkout : ndarray\n"
                 "    The first three integer workspace entries: numerical rank, nullity and the\n"
                 "    count of accurate singular values.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gglse(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, c, d, lwork=m+n+p, overwrite_a=0, overwrite_b=0,\n"
                 "      overwrite_c=0, overwrite_d=0)\n\n";
            s += "Solve an equality-constrained least squares problem\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Minimizes ``norm(c - a @ x)`` subject to ``b @ x == d``.\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Least squares matrix, shape ``(m, n)``.\n";
            s += "b : ndarray\n    Constraint matrix, shape ``(p, n)`` with ``n - m <= p <= n``.\n";
            s += "c : ndarray\n    Least squares right-hand side, length ``m``.\n";
            s += "d : ndarray\n"
                 "    Constraint right-hand side, length ``p``; empty when there are no\n"
                 "    constraints.\n";
            s += "lwork : int, optional\n"
                 "    Workspace size; -1 requests the query. Default is ``m + n + p``.\n";
            s += P_OVERWRITE_A;
            s += P_OVERWRITE_B;
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";
            s += "overwrite_d : int, optional\n"
                 "    If nonzero, `d` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "t : ndarray\n    The generalized RQ factorization left in `a`.\n";
            s += "r : ndarray\n    The triangular factor left in `b`.\n";
            s += "res : ndarray\n    The residual information left in `c`.\n";
            s += "x : ndarray\n    The solution, length ``n``.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_gglse_lwork(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(m, n, p)\n\n";
            s += "Query the optimal workspace for ``" + std::string(name).substr(0, std::strlen(name) - 6) + "``\n"
                 "(LAPACK ``" + std::string(name).substr(0, std::strlen(name) - 6) + "`` with ``lwork = -1``).\n\n";

            s += "Parameters\n----------\n";
            s += P_M;
            s += P_N;
            s += "p : int\n    Row count of the constraint matrix.\n";

            s += "\nReturns\n-------\n";
            s += std::string("work : ") + t.scalar + "\n    Optimal `lwork`.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_lasd4(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(i, d, z, rho=1.0)\n\n";
            s += "Compute one root of the secular equation of a rank-one modified diagonal\n"
                 "matrix (LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Used by the divide-and-conquer SVD update. LAPACK ships no complex\n"
                 "counterpart, so only the real flavors exist.\n\n";

            s += "Parameters\n----------\n";
            s += "i : int\n    Which root to compute, 0-based.\n";
            s += "d : ndarray\n    The diagonal, in increasing order.\n";
            s += "z : ndarray\n    The rank-one update vector, the same length as `d`.\n";
            s += std::string("rho : ") + t.scalar + ", optional\n    Scalar the update is scaled by. Default is 1.0.\n";

            s += "\nReturns\n-------\n";
            s += "delta : ndarray\n    The differences ``d[j] - sigma``.\n";
            s += std::string("sigma : ") + t.scalar + "\n    The computed root.\n";
            s += "work : ndarray\n    The sums ``d[j] + sigma``.\n";
            s += "info : int\n    0 on success; 1 if the iteration failed to converge.\n";
            return s;
        }

        /* ------------------ symmetric and Hermitian banded eigensolvers -------------------- */

        static constexpr const char *P_AB_HERM_BAND =
            "ab : ndarray\n"
            "    The band of the symmetric/Hermitian matrix, shape ``(kd + 1, n)``: the order\n"
            "    `n` and the bandwidth `kd` are both read off `ab` itself. It is consumed rather\n"
            "    than returned.\n";
        static constexpr const char *P_COMPUTE_V_BAND =
            "compute_v : int, optional\n"
            "    If nonzero, eigenvectors are computed and `z` is sized for them; otherwise `z`\n"
            "    comes back as a 1x1 placeholder. Default is 1.\n";
        static constexpr const char *P_LOWER_BAND =
            "lower : int, optional\n"
            "    If nonzero, `ab` holds the lower band; otherwise the upper. Default is 0.\n";
        static constexpr const char *R_W_BAND =
            "w : ndarray\n    Eigenvalues in ascending order, length ``n``. Real for every flavor.\n";

        static std::string
        doc_sbev(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(ab, compute_v=1, lower=0, ldab=ab.shape[0],\n"
                 "      overwrite_ab=0)\n\n";
            s += "Eigenvalues and eigenvectors of a symmetric band matrix\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "``sbevd`` is the divide-and-conquer version and is generally preferred.\n"
                 "LAPACK ships no Hermitian counterpart under this name, so only the real\n"
                 "flavors exist.\n\n";

            s += "Parameters\n----------\n";
            s += P_AB_HERM_BAND;
            s += P_COMPUTE_V_BAND;
            s += P_LOWER_BAND;
            s += P_LDAB;
            /* These three alone default the flag to 1: the `.pyf` declares `ab`
             * `intent(in, overwrite)`, so reuse is the norm and copying is opt-in. */
            s += "overwrite_ab : int, optional\n"
                 "    If nonzero, `ab` may be overwritten in place -- and it is consumed either\n"
                 "    way, never returned. Unlike the rest of this module, the default is 1;\n"
                 "    pass 0 to keep your array intact.\n";

            s += "\nReturns\n-------\n";
            s += R_W_BAND;
            s += "z : ndarray\n    Eigenvectors as columns, when `compute_v` asked for them.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_sbevd(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(ab, compute_v=1, lower=0, ldab=ab.shape[0], ";
            s += t.is_complex ? "lrwork=...,\n      liwork=..., " : "liwork=...,\n      ";
            s += "overwrite_ab=0)\n\n";
            s += std::string("Eigenvalues and eigenvectors of a ") +
                 (t.is_complex ? "Hermitian" : "symmetric") + " band matrix, by divide and\n"
                 "conquer (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_AB_HERM_BAND;
            if (t.is_complex) {
                s += "    The order must be at least 1: the workspace formulas this routine uses\n"
                     "    do not cover an empty matrix, and the ``.pyf`` rejects it outright.\n";
            }
            s += P_COMPUTE_V_BAND;
            s += P_LOWER_BAND;
            s += P_LDAB;
            if (t.is_complex) {
                s += "lrwork : int, optional\n"
                     "    Real workspace size, at least ``1 + 5n + 2n**2`` with eigenvectors and\n"
                     "    ``n`` without. The real flavors have no second workspace and take no\n"
                     "    `lrwork`.\n";
            }
            s += "liwork : int, optional\n"
                 "    Integer workspace size, at least ``3 + 5n`` with eigenvectors and 1\n"
                 "    without.\n";
            /* These three alone default the flag to 1: the `.pyf` declares `ab`
             * `intent(in, overwrite)`, so reuse is the norm and copying is opt-in. */
            s += "overwrite_ab : int, optional\n"
                 "    If nonzero, `ab` may be overwritten in place -- and it is consumed either\n"
                 "    way, never returned. Unlike the rest of this module, the default is 1;\n"
                 "    pass 0 to keep your array intact.\n";

            s += "\nReturns\n-------\n";
            s += R_W_BAND;
            s += "z : ndarray\n    Eigenvectors as columns, when `compute_v` asked for them.\n";
            s += R_INFO;
            return s;
        }

        static std::string
        doc_sbevx(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(ab, vl, vu, il, iu, ldab=ab.shape[0], compute_v=1,\n"
                 "      range=0, lower=0, abstol=0.0, mmax=..., overwrite_ab=0)\n\n";
            s += std::string("Selected eigenvalues and eigenvectors of a ") +
                 (t.is_complex ? "Hermitian" : "symmetric") + " band matrix\n(LAPACK ``" +
                 std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_AB_HERM_BAND;
            s += "vl : float\n    Lower bound of the half-open interval; used only when `range` is 1.\n";
            s += "vu : float\n    Upper bound of the half-open interval; used only when `range` is 1.\n";
            s += "il : int\n"
                 "    Index of the smallest eigenvalue to return, 1-based; used only when\n"
                 "    `range` is 2. This group keeps Fortran indexing, deliberately.\n";
            s += "iu : int\n"
                 "    Index of the largest eigenvalue to return, 1-based; used only when `range`\n"
                 "    is 2.\n";
            s += P_LDAB;
            s += P_COMPUTE_V_BAND;
            s += "range : int, optional\n"
                 "    Which eigenvalues to compute: 0 for all, 1 for those in ``(vl, vu]``, 2 for\n"
                 "    those with indices `il` through `iu`. Default is 0.\n";
            s += P_LOWER_BAND;
            s += "abstol : float, optional\n"
                 "    Absolute error tolerance. Accuracy is best at twice the underflow\n"
                 "    threshold, ``2 * lamch('S')``, rather than at 0. Default is 0.0.\n";
            s += "mmax : int, optional\n"
                 "    How many eigenvectors `z` has room for. With ``range == 1`` the count is\n"
                 "    not known until the call returns, so the default is the whole spectrum;\n"
                 "    pass a tighter bound to avoid allocating for eigenvectors you will not get.\n";
            /* These three alone default the flag to 1: the `.pyf` declares `ab`
             * `intent(in, overwrite)`, so reuse is the norm and copying is opt-in. */
            s += "overwrite_ab : int, optional\n"
                 "    If nonzero, `ab` may be overwritten in place -- and it is consumed either\n"
                 "    way, never returned. Unlike the rest of this module, the default is 1;\n"
                 "    pass 0 to keep your array intact.\n";

            s += "\nReturns\n-------\n";
            s += R_W_BAND;
            s += "z : ndarray\n    Eigenvectors as columns, ``(n, mmax)`` when `compute_v` is set.\n";
            s += "m : int\n    Number of eigenvalues actually found.\n";
            s += "ifail : ndarray\n"
                 "    Indices of the eigenvectors that failed to converge; meaningful only when\n"
                 "    `info` is positive.\n";
            s += R_INFO;
            return s;
        }

        /* ------------------- auxiliaries, norms and machine parameters --------------------- */

        static constexpr const char *P_NORM_LETTER =
            "norm : str\n"
            "    Which norm: ``'M'`` for the largest absolute entry, ``'1'`` or ``'O'`` for the\n"
            "    1-norm, ``'I'`` for the infinity norm, ``'F'`` or ``'E'`` for the Frobenius\n"
            "    norm. Lowercase is accepted.\n";

        static std::string
        doc_lamch(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(cmach)\n\n";
            s += "Machine parameters for this precision (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += "cmach : str\n"
                 "    Which parameter: ``'E'`` relative machine epsilon, ``'S'`` safe minimum,\n"
                 "    ``'B'`` base, ``'P'`` ``eps * base``, ``'N'`` mantissa digits, ``'R'``\n"
                 "    rounding, ``'M'`` minimum exponent, ``'U'`` underflow threshold, ``'L'``\n"
                 "    maximum exponent, ``'O'`` overflow threshold. Lowercase is accepted, and\n"
                 "    anything else returns 0.\n";

            s += "\nReturns\n-------\n";
            s += std::string("x : ") + t.scalar + "\n    The requested parameter.\n";
            return s;
        }

        static std::string
        doc_lange(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(norm, a)\n\n";
            s += "Norm of a general matrix (LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_NORM_LETTER;
            s += P_A_GENERAL;

            s += "\nReturns\n-------\n";
            s += "n2 : float\n    The norm. Real for every flavor.\n";
            return s;
        }

        static std::string
        doc_lantr(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(norm, a, uplo='U', diag='N')\n\n";
            s += "Norm of a triangular or trapezoidal matrix\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";

            s += "Parameters\n----------\n";
            s += P_NORM_LETTER;
            s += "a : ndarray\n"
                 "    Matrix of shape ``(m, n)``; only the triangle named by `uplo` is read.\n";
            s += "uplo : str, optional\n"
                 "    ``'U'`` if `a` is upper triangular, ``'L'`` if lower. Default is ``'U'``.\n";
            s += P_DIAG_LETTER;

            s += "\nReturns\n-------\n";
            s += "n2 : float\n    The norm. Real for every flavor.\n";
            return s;
        }

        static std::string
        doc_larfg(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(n, alpha, x, incx=1, overwrite_x=0)\n\n";
            s += "Generate an elementary reflector (LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Produces the ``H`` with ``H @ [alpha, x] == [beta, 0]``, described by the\n"
                 "returned `alpha` (which becomes ``beta``), `x` and `tau`.\n\n";

            s += "Parameters\n----------\n";
            s += "n : int\n    Order of the reflector, at least 1.\n";
            s += std::string("alpha : ") + t.scalar + "\n    The leading entry of the vector to reflect.\n";
            s += "x : ndarray\n"
                 "    The remaining entries, ``1 + (n - 2) * abs(incx)`` of them.\n";
            s += "incx : int, optional\n    Stride through `x`. Default is 1.\n";
            s += "overwrite_x : int, optional\n"
                 "    If nonzero, `x` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += std::string("alpha : ") + t.scalar + "\n    The resulting ``beta``.\n";
            s += "x : ndarray\n    The reflector vector.\n";
            s += std::string("tau : ") + t.scalar + "\n    The reflector scalar.\n";
            return s;
        }

        static std::string
        doc_larf(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(v, tau, c, work, side='L', incv=1, overwrite_c=0)\n\n";
            s += "Apply an elementary reflector to a matrix\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Computes ``H @ c`` or ``c @ H`` for the reflector ``H`` described by `v` and\n"
                 "`tau`.\n\n";

            s += "Parameters\n----------\n";
            s += "v : ndarray\n"
                 "    The reflector vector, ``1 + (order - 1) * abs(incv)`` entries, where the\n"
                 "    order is the row count of `c` when `side` is ``'L'`` and its column count\n"
                 "    when ``'R'``.\n";
            s += std::string("tau : ") + t.scalar + "\n    The reflector scalar.\n";
            s += "c : ndarray\n    Matrix to apply the reflector to.\n";
            /* The `.pyf` calls this out as a mistake it kept; the wrapper keeps it too. */
            s += "work : ndarray\n"
                 "    Workspace, length `n` when `side` is ``'L'`` and `m` when ``'R'``. This is\n"
                 "    a required argument only for backwards compatibility -- the ``.pyf`` says\n"
                 "    so itself -- and would otherwise be allocated internally.\n";
            s += "side : str, optional\n"
                 "    ``'L'`` for ``H @ c``, ``'R'`` for ``c @ H``. Default is ``'L'``.\n";
            s += "incv : int, optional\n    Stride through `v`. Default is 1.\n";
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "c : ndarray\n    The result, overwriting `c`.\n";
            return s;
        }

        static std::string
        doc_lartg(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(f, g)\n\n";
            s += "Generate a Givens rotation (LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Produces ``cs`` and ``sn`` with ``[[cs, sn], [-conj(sn), cs]] @ [f, g] ==\n"
                 "[r, 0]``.\n\n";

            s += "Parameters\n----------\n";
            s += std::string("f : ") + t.scalar + "\n    First component.\n";
            s += std::string("g : ") + t.scalar + "\n    Second component, the one rotated to zero.\n";

            s += "\nReturns\n-------\n";
            s += "cs : float\n    The cosine. Real for every flavor.\n";
            s += std::string("sn : ") + t.scalar + "\n    The sine.\n";
            s += std::string("r : ") + t.scalar + "\n    The resulting first component.\n";
            return s;
        }

        static std::string
        doc_rot(const char *name, const Dtype &t) noexcept
        {
            std::string s;
            s += std::string(name) + "(x, y, c, s, n=..., offx=0, incx=1, offy=0, incy=1,\n"
                 "      overwrite_x=0, overwrite_y=0)\n\n";
            s += "Apply a plane rotation to a pair of complex vectors\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Only the complex flavors exist here; for real vectors this is BLAS ``rot``.\n\n";

            s += "Parameters\n----------\n";
            s += "x : ndarray\n    First vector.\n";
            s += "y : ndarray\n    Second vector.\n";
            s += "c : float\n    The cosine. Real even though the vectors are complex.\n";
            s += std::string("s : ") + t.scalar + "\n    The sine.\n";
            s += "n : int, optional\n"
                 "    How many elements to rotate. Default is as many as `x` allows from `offx`\n"
                 "    at stride `incx`.\n";
            s += "offx : int, optional\n    Index in `x` to start from. Default is 0.\n";
            s += "incx : int, optional\n    Stride through `x`. Default is 1.\n";
            s += "offy : int, optional\n    Index in `y` to start from. Default is 0.\n";
            s += "incy : int, optional\n    Stride through `y`. Default is 1.\n";
            s += "overwrite_x : int, optional\n"
                 "    If nonzero, `x` may be overwritten in place. Default is 0.\n";
            s += "overwrite_y : int, optional\n"
                 "    If nonzero, `y` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "x : ndarray\n    The rotated first vector.\n";
            s += "y : ndarray\n    The rotated second vector.\n";
            return s;
        }

        static std::string
        doc_ilaver(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "()\n\n";
            s += "Version of the LAPACK library this module is linked against\n"
                 "(LAPACK ``ilaver``).\n\n";
            s += "The one routine here with no flavor at all, so it takes no arguments and has\n"
                 "no ``s``/``d``/``c``/``z`` prefix.\n\n";

            s += "Parameters\n----------\n";
            s += "None\n\n";

            s += "Returns\n-------\n";
            s += "major : int\n    Major version number.\n";
            s += "minor : int\n    Minor version number.\n";
            s += "patch : int\n    Patch level.\n";
            return s;
        }

        static std::string
        doc_tgsyl(const char *name, const Dtype &) noexcept
        {
            std::string s;
            s += std::string(name) + "(a, b, c, d, e, f, trans='N', ijob=0, lwork=2*m*n,\n"
                 "      overwrite_c=0, overwrite_f=0)\n\n";
            s += "Solve the generalized Sylvester equations\n"
                 "(LAPACK ``" + std::string(name) + "``).\n\n";
            s += "Solves ``a @ r - l @ b == scale * c`` together with\n"
                 "``d @ r - l @ e == scale * f`` for the pair ``(r, l)``, which overwrites\n"
                 "``(c, f)``. The ``.pyf`` declares only the real flavors.\n\n";

            s += "Parameters\n----------\n";
            s += "a : ndarray\n    Upper quasi-triangular matrix of shape ``(m, m)``.\n";
            s += "b : ndarray\n    Upper quasi-triangular matrix of shape ``(n, n)``.\n";
            s += "c : ndarray\n    Right-hand side of the first equation, shape ``(m, n)``.\n";
            s += "d : ndarray\n    Upper triangular matrix of shape ``(m, m)``.\n";
            s += "e : ndarray\n    Upper triangular matrix of shape ``(n, n)``.\n";
            s += "f : ndarray\n    Right-hand side of the second equation, shape ``(m, n)``.\n";
            s += "trans : str, optional\n"
                 "    ``'N'`` for the equations above, ``'T'`` for the transposed pair.\n"
                 "    Default is ``'N'``.\n";
            s += "ijob : int, optional\n"
                 "    0 solves only; 1 through 4 also estimate `dif`, by different methods.\n"
                 "    Default is 0.\n";
            s += "lwork : int, optional\n"
                 "    Workspace size. Only ``ijob`` of 1 or 2 with ``trans == 'N'`` needs the\n"
                 "    full ``2 * m * n``; the rest get by with `n`. Default is ``2 * m * n``.\n";
            s += "overwrite_c : int, optional\n"
                 "    If nonzero, `c` may be overwritten in place. Default is 0.\n";
            s += "overwrite_f : int, optional\n"
                 "    If nonzero, `f` may be overwritten in place. Default is 0.\n";

            s += "\nReturns\n-------\n";
            s += "r : ndarray\n    First component of the solution, overwriting `c`.\n";
            s += "l : ndarray\n    Second component of the solution, overwriting `f`.\n";
            s += "scale : float\n"
                 "    Scale factor, at most 1, chosen to avoid overflow in the solution.\n";
            s += "dif : float\n"
                 "    Estimate of the separation of the two pencils; 0 unless `ijob` asked\n"
                 "    for it.\n";
            s += "info : int\n"
                 "    0 on success; if negative, the ``-info``-th argument had an illegal value;\n"
                 "    1 if the pencils have common or very close eigenvalues.\n";
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

        /** @brief A symmetric/Hermitian family: `s/dsy<fam>`, `c/zhe<fam>` */
        #define DOC_FAMILY_SYHE(fam)     \
            {"ssy" #fam, doc_##fam, S},  \
            {"dsy" #fam, doc_##fam, D},  \
            {"che" #fam, doc_##fam, C},  \
            {"zhe" #fam, doc_##fam, Z}

        /** @brief An orthogonal/unitary family: `s/dor<fam>`, `c/zun<fam>`.
         *
         * The builder keeps the real spelling (`doc_orghr`, not `doc_ghr`) because that is a
         * routine name someone can grep for, where the bare stem is not.  That is the one way
         * this differs from DOC_FAMILY_SYHE above. */
        #define DOC_FAMILY_ORUN(fam)        \
            {"sor" #fam, doc_or##fam, S},   \
            {"dor" #fam, doc_or##fam, D},   \
            {"cun" #fam, doc_or##fam, C},   \
            {"zun" #fam, doc_or##fam, Z}

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
            DOC_FAMILY(gtsv),
            DOC_FAMILY(gttrf),
            DOC_FAMILY(gttrs),
            DOC_FAMILY(gtcon),
            DOC_FAMILY(gtsvx),
            DOC_FAMILY(gbsv),
            DOC_FAMILY(gbtrf),
            DOC_FAMILY(gbtrs),
            DOC_FAMILY(gbcon),
            DOC_FAMILY(langb),
            DOC_FAMILY(pstrf),
            DOC_FAMILY(pstf2),
            DOC_FAMILY(posv),
            DOC_FAMILY(posvx),
            DOC_FAMILY(pocon),
            DOC_FAMILY(potrf),
            DOC_FAMILY(potrs),
            DOC_FAMILY(potri),
            DOC_FAMILY(ptsv),
            DOC_FAMILY(pttrf),
            DOC_FAMILY(pttrs),
            DOC_FAMILY(pteqr),
            DOC_FAMILY(ptsvx),
            DOC_FAMILY(sytrf),
            DOC_FAMILY(sytrf_lwork),
            DOC_FAMILY(sytf2),
            DOC_FAMILY(sytrs),
            DOC_FAMILY(sytri),
            DOC_FAMILY(syconv),
            DOC_FAMILY(syequb),
            DOC_FAMILY(sycon),
            DOC_FAMILY(sysv),
            DOC_FAMILY(sysv_lwork),
            DOC_FAMILY(sysvx),
            DOC_FAMILY(sysvx_lwork),
            DOC_FAMILY_SYHE(ev),
            DOC_FAMILY_SYHE(evx),
            DOC_FAMILY_SYHE(evd),
            DOC_FAMILY_SYHE(evr),
            DOC_FAMILY_SYHE(gvd),
            DOC_FAMILY_SYHE(gv),
            DOC_FAMILY_SYHE(gvx),
            DOC_FAMILY_SYHE(trd),
            DOC_FAMILY_SYHE(gst),
            DOC_FAMILY_SYHE(ev_lwork),
            DOC_FAMILY_SYHE(evx_lwork),
            DOC_FAMILY_SYHE(gv_lwork),
            DOC_FAMILY_SYHE(gvx_lwork),
            DOC_FAMILY_SYHE(trd_lwork),
            DOC_FAMILY_SYHE(evd_lwork),
            DOC_FAMILY_SYHE(evr_lwork),

            // Irregular names that don't fit the DOC_FAMILY pattern.
            {"chetrf", doc_hetrf, C},
            {"chetrf_lwork", doc_hetrf_lwork, C},
            {"zhetrf", doc_hetrf, Z},
            {"zhetrf_lwork", doc_hetrf_lwork, Z},
            {"chetrs", doc_hetrs, C},
            {"zhetrs", doc_hetrs, Z},
            {"chetri", doc_hetri, C},
            {"zhetri", doc_hetri, Z},
            {"cheequb", doc_heequb, C},
            {"zheequb", doc_heequb, Z},
            {"checon", doc_hecon, C},
            {"zhecon", doc_hecon, Z},
            {"chesv", doc_hesv, C},
            {"zhesv", doc_hesv, Z},
            {"chesvx", doc_hesvx, C},
            {"zhesvx", doc_hesvx, Z},
            {"sstev", doc_stev, S},
            {"dstev", doc_stev, D},
            {"sstebz", doc_stebz, S},
            {"dstebz", doc_stebz, D},
            {"ssterf", doc_sterf, S},
            {"dsterf", doc_sterf, D},
            {"sstein", doc_stein, S},
            {"dstein", doc_stein, D},
            {"sstemr", doc_stemr, S},
            {"sstemr_lwork", doc_stemr_lwork, S},
            {"dstemr", doc_stemr, D},
            {"dstemr_lwork", doc_stemr_lwork, D},
            {"sstevd", doc_stevd, S},
            {"dstevd", doc_stevd, D},
            {"chesv_lwork", doc_hesv_lwork, C},
            {"zhesv_lwork", doc_hesv_lwork, Z},
            {"chesvx_lwork", doc_hesvx_lwork, C},
            {"zhesvx_lwork", doc_hesvx_lwork, Z},

            /* flapack_other: triangular storage conversions */
            DOC_FAMILY(tpttf),
            DOC_FAMILY(tpttr),
            DOC_FAMILY(tfttp),
            DOC_FAMILY(tfttr),
            DOC_FAMILY(trttf),
            DOC_FAMILY(trttp),
            DOC_FAMILY(tfsm),

            /* flapack_other: packed positive definite and the RFP rank-k update */
            DOC_FAMILY(ppcon),
            DOC_FAMILY(ppsv),
            DOC_FAMILY(pptrf),
            DOC_FAMILY(pptri),
            DOC_FAMILY(pptrs),
            {"ssfrk", doc_sfrk, S},
            {"dsfrk", doc_sfrk, D},
            {"chfrk", doc_sfrk, C},
            {"zhfrk", doc_sfrk, Z},

            /* flapack_other: RFP and banded Cholesky */
            DOC_FAMILY(pftrf),
            DOC_FAMILY(pftri),
            DOC_FAMILY(pftrs),
            DOC_FAMILY(pbtrf),
            DOC_FAMILY(pbtrs),
            DOC_FAMILY(pbsv),

            /* flapack_other: triangular solves and the LU auxiliaries */
            DOC_FAMILY(trtrs),
            DOC_FAMILY(trcon),
            DOC_FAMILY(tbtrs),
            DOC_FAMILY(trtri),
            DOC_FAMILY(lauum),
            DOC_FAMILY(laswp),

            /* flapack_other: Schur exchange and Sylvester */
            DOC_FAMILY(trexc),
            DOC_FAMILY(tgexc),
            DOC_FAMILY(trsyl),
            {"stgsyl", doc_tgsyl, S},  {"dtgsyl", doc_tgsyl, D},

            /* flapack_other: Schur cluster reordering */
            DOC_FAMILY(trsen),
            DOC_FAMILY(trsen_lwork),
            DOC_FAMILY(tgsen),
            DOC_FAMILY(tgsen_lwork),

            /* flapack_other: the `or`/`un` spelling pairs */
            DOC_FAMILY_ORUN(ghr),
            DOC_FAMILY_ORUN(ghr_lwork),
            DOC_FAMILY_ORUN(gqr),
            DOC_FAMILY_ORUN(grq),
            DOC_FAMILY_ORUN(mqr),
            DOC_FAMILY_ORUN(mrz),
            DOC_FAMILY_ORUN(mrz_lwork),

            /* flapack_other: compact-WY and blocked QR */
            DOC_FAMILY(geqrt),
            DOC_FAMILY(gemqrt),
            DOC_FAMILY(tpqrt),
            DOC_FAMILY(tpmqrt),
            DOC_FAMILY(tzrzf),
            DOC_FAMILY(tzrzf_lwork),

            /* flapack_other: CS decomposition and the one-off solvers */
            DOC_FAMILY_ORUN(csd),
            DOC_FAMILY_ORUN(csd_lwork),
            {"sgejsv", doc_gejsv, S},   {"dgejsv", doc_gejsv, D},
            {"slasd4", doc_lasd4, S},   {"dlasd4", doc_lasd4, D},
            DOC_FAMILY(gglse),
            DOC_FAMILY(gglse_lwork),

            /* flapack_other: symmetric and Hermitian banded eigensolvers */
            {"ssbev", doc_sbev, S},    {"dsbev", doc_sbev, D},
            {"ssbevd", doc_sbevd, S},  {"dsbevd", doc_sbevd, D},
            {"chbevd", doc_sbevd, C},  {"zhbevd", doc_sbevd, Z},
            {"ssbevx", doc_sbevx, S},  {"dsbevx", doc_sbevx, D},
            {"chbevx", doc_sbevx, C},  {"zhbevx", doc_sbevx, Z},

            /* flapack_other: auxiliaries, norms and machine parameters */
            {"slamch", doc_lamch, S},  {"dlamch", doc_lamch, D},
            DOC_FAMILY(lange),
            DOC_FAMILY(lantr),
            DOC_FAMILY(larfg),
            DOC_FAMILY(larf),
            DOC_FAMILY(lartg),
            {"crot", doc_rot, C},      {"zrot", doc_rot, Z},
            {"ilaver", doc_ilaver, D},
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
