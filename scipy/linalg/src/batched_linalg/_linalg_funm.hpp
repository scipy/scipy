/*
 * Templated loops for `linalg.funm`: the conversion of a real Schur form to the complex one,
 * and the Parlett recurrence for a function of an upper triangular matrix.
 *
 * Both work in place on n-by-n F-contiguous complex arrays and leave the vector operations to
 * BLAS and LAPACK.
 */
#pragma once

namespace sp_linalg {


/**
 * @brief Convert a real Schur form to the complex Schur form, in place.
 *
 * On entry T and Z hold the real Schur factors of A = Z T Z^T, cast to complex. T is quasi upper
 * triangular: a nonzero T[k+1, k] marks a 2x2 diagonal block [[a, b], [c, d]], as ?gees returns
 * them. On exit T is upper triangular and Z unitary, with A = Z T Z^H.
 *
 * Each block is triangularized by the rotation G of ?lartg that takes the eigenvector
 * (lambda - d, c) of the block to (r, 0), where lambda is the eigenvalue from ?lanv2 with the
 * positive imaginary part; ?rot applies T <- G T G^H and Z <- Z G^H. This is the algorithm of
 * `scipy.linalg.rsf2csf` and of the real path of `sqrtm`, except that the cosine of the ?lartg
 * rotation is real: the result differs from theirs by a diagonal unitary similarity, which leaves
 * the diagonal of T alone. The blocks do not overlap, so they can be taken in any order.
 *
 * @param T  n-by-n, F-contiguous.
 * @param Z  n-by-n, F-contiguous.
 * @param n  Order of the matrices.
 */
template<typename R>
void rsf2csf_inplace(std::complex<R> *T, std::complex<R> *Z, const CBLAS_INT n)
{
    using C = std::complex<R>;
    const npy_intp ld = n;
    for (CBLAS_INT k = 0; k + 1 < n; k++) {
        C *Tkk = T + k + k*ld;          // T[k, k]
        const R c = Tkk[1].real();      // T[k+1, k]
        if (c == R(0)) { continue; }

        // ?lanv2 overwrites its arguments with the standardized block; only the eigenvalue is used.
        const R d = Tkk[1 + ld].real();
        R aa = Tkk[0].real(), bb = Tkk[ld].real(), cc = c, dd = d;
        R rt1r, rt1i, rt2r, rt2i, cs_std, sn_std;
        lapack::lanv2(&aa, &bb, &cc, &dd, &rt1r, &rt1i, &rt2r, &rt2i, &cs_std, &sn_std);

        // G = [[cs, sn], [-conj(sn), cs]] with G (lambda - d, c) = (r, 0)
        R cs;
        C sn, r;
        lapack::lartg(C(rt1r - d, rt1i), C(c), &cs, &sn, &r);

        // G T on rows k, k+1 from column k on; (.) G^H on columns k, k+1 of T, down to row k+1,
        // and of Z. A row vector times G^H is ?rot with the sine conjugated.
        lapack::rot(n - k, Tkk, n, Tkk + 1, n, cs, sn);
        lapack::rot(k + 2, T + k*ld, 1, T + (k + 1)*ld, 1, cs, std::conj(sn));
        lapack::rot(n, Z + k*ld, 1, Z + (k + 1)*ld, 1, cs, std::conj(sn));
        Tkk[1] = C(0);
        k++;  // the second column of the block is done
    }
}


/**
 * @brief The Parlett recurrence for F = f(T), T upper triangular; Algorithm 11.1.1 in [1].
 *
 * On entry the diagonal of F holds f(T[i, i]); on exit its strict upper triangle holds that of
 * f(T). Where two diagonal entries of T coincide, the denominator T[j, j] - T[i, i] is zero and
 * the entry is left undivided, as in the original algorithm; the caller turns the smallest
 * denominator into an error estimate.
 *
 * The recurrence runs column by column, and upwards within a column, which needs the same
 * entries as the classical order by superdiagonals. That way both of its sums are BLAS calls:
 *   v[i] = sum_{i<k<j} F[i, k] T[k, j] is the product of the upper triangle of F[0:j-1, 1:j]
 *          with T[1:j, j], one ?trmv per column;
 *   u[i] = sum_{i<k<j} T[i, k] F[k, j] receives T[0:k, k] F[k, j], one ?axpy, as soon as F[k, j]
 *          is known.
 * A ?trsv cannot replace the u updates, since it would divide by the zero denominators.
 *
 * [1] G. H. Golub, C. F. Van Loan, "Matrix Computations", 4th ed.
 *
 * @param F       n-by-n, F-contiguous.
 * @param T       n-by-n, F-contiguous; only its upper triangle is read, and it is not modified.
 * @param n       Order of the matrices.
 * @param minden  Starting value of the minimum of the |T[j, j] - T[i, i]|.
 * @param work    Workspace of 2*n elements.
 * @return min(minden, |T[j, j] - T[i, i]|) over all i < j.
 */
template<typename R>
R funm_parlett(std::complex<R> *F, std::complex<R> *T, const CBLAS_INT n, R minden,
               std::complex<R> *work)
{
    using C = std::complex<R>;
    const npy_intp ld = n;
    C *u = work, *v = work + n;
    for (CBLAS_INT j = 1; j < n; j++) {
        C *Tj = T + j*ld, *Fj = F + j*ld;

        // v[0:j-1] = triu(F[0:j-1, 1:j]) @ T[1:j, j], and v[j-1] = 0 is the empty sum.
        std::copy(Tj + 1, Tj + j, v);
        if (j > 1) { blas::trmv('U', 'N', 'N', j - 1, F + ld, n, v, 1); }
        v[j - 1] = C(0);
        std::fill(u, u + j, C(0));

        const C fjj = Fj[j], tjj = Tj[j];
        for (CBLAS_INT i = j - 1; i >= 0; i--) {
            C s = Tj[i] * (fjj - F[i + i*ld]);
            s += u[i] - v[i];
            const C den = tjj - T[i + i*ld];
            if (den != C(0)) { s /= den; }
            Fj[i] = s;
            minden = std::min(minden, std::abs(den));
            // u[0:i] += T[0:i, i] F[i, j]
            if (i > 0) { blas::axpy(i, s, T + i*ld, 1, u, 1); }
        }
    }
    return minden;
}


}  // namespace sp_linalg
