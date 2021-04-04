/**
 * Author: Damian Eads
 * Date:   September 22, 2007 (moved to new file on June 8, 2008)
 *
 * Copyright (c) 2007, 2008, Damian Eads. All rights reserved.
 * Adapted for incorporation into Scipy, April 9, 2008.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   - Redistributions of source code must retain the above
 *     copyright notice, this list of conditions and the
 *     following disclaimer.
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer
 *     in the documentation and/or other materials provided with the
 *     distribution.
 *   - Neither the name of the author nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "_c99compat.h"


static NPY_INLINE void
_row_norms(const double *X, npy_intp num_rows, const npy_intp num_cols, double *norms_buff){
    /* Compute the row norms. */
    npy_intp i, j;
    for (i = 0; i < num_rows; ++i) {
        for (j = 0; j < num_cols; ++j, ++X) {
            const double curr_val = *X;
            norms_buff[i] += curr_val * curr_val;
        }
        norms_buff[i] = sqrt(norms_buff[i]);
    }
}

static NPY_INLINE double
sqeuclidean_distance_double(const double *u, const double *v, const npy_intp n)
{
    double s = 0.0;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        const double d = u[i] - v[i];
        s += d * d;
    }
    return s;
}

static NPY_INLINE double
euclidean_distance_double(const double *u, const double *v, const npy_intp n)
{
    return sqrt(sqeuclidean_distance_double(u, v, n));
}

#if 0   /* XXX unused */
static NPY_INLINE double
ess_distance(const double *u, const double *v, const npy_intp n)
{
    double s = 0.0, d;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        d = fabs(u[i] - v[i]);
        s += d * d;
    }
    return s;
}
#endif

static NPY_INLINE double
chebyshev_distance_double(const double *u, const double *v, const npy_intp n)
{
    double maxv = 0.0;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        const double d = fabs(u[i] - v[i]);
        if (d > maxv) {
            maxv = d;
        }
    }
    return maxv;
}

static NPY_INLINE double
weighted_chebyshev_distance_double(const double *u, const double *v,
                                   const npy_intp n, const double *w)
{
    npy_intp i;
    double maxv = 0.0;
    for (i = 0; i < n; ++i) {
        if (w[i] == 0.0) continue;
        const double d = fabs(u[i] - v[i]);
        if (d > maxv) {
            maxv = d;
        }
    }
    return maxv;
}

static NPY_INLINE double
canberra_distance_double(const double *u, const double *v, const npy_intp n)
{
    double tot = 0.;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        const double x = u[i], y = v[i];
        const double snum = fabs(x - y);
        const double sdenom = fabs(x) + fabs(y);
        if (sdenom > 0.) {
            tot += snum / sdenom;
        }
    }
    return tot;
}

static NPY_INLINE double
bray_curtis_distance_double(const double *u, const double *v, const npy_intp n)
{
    double s1 = 0.0, s2 = 0.0;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        s1 += fabs(u[i] - v[i]);
        s2 += fabs(u[i] + v[i]);
    }
    return s1 / s2;
}

/*
 * Timings with various BLAS implementations (gh-5657) have shown that using
 * OpenBLAS's cblas_ddot here can give some speedup, but only for high-d data
 * and an untuned ATLAS is slower than rolling our own.
 */
static NPY_INLINE double
dot_product(const double *u, const double *v, const npy_intp n)
{
    double s = 0.0;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        s += u[i] * v[i];
    }
    return s;
}

static NPY_INLINE double
mahalanobis_distance(const double *u, const double *v, const double *covinv,
                     double *dimbuf1, double *dimbuf2, const npy_intp n)
{
    npy_intp i;

    for (i = 0; i < n; ++i) {
        dimbuf1[i] = u[i] - v[i];
    }
    /*
     * Note: matrix-vector multiplication (GEMV). Again, OpenBLAS can speed
     * this up for high-d data.
     */
    for (i = 0; i < n; ++i) {
        const double *covrow = covinv + (i * n);
        dimbuf2[i] = dot_product(dimbuf1, covrow, n);
    }
    return sqrt(dot_product(dimbuf1, dimbuf2, n));
}

static NPY_INLINE double
hamming_distance_double(const double *u, const double *v, const npy_intp n, const double *w)
{
    npy_intp i;
    double s = 0;
    double w_sum = 0;

    for (i = 0; i < n; ++i) {
        s += ((double) (u[i] != v[i])) * w[i];
        w_sum += w[i];
    }

    return s / w_sum;
}

static NPY_INLINE double
hamming_distance_char(const char *u, const char *v, const npy_intp n, const double *w)
{
    npy_intp i;
    double s = 0;
    double w_sum = 0;

    for (i = 0; i < n; ++i) {
        s += ((double) (u[i] != v[i])) * w[i];
        w_sum += w[i];
    }

    return s / w_sum;
}

static NPY_INLINE double
yule_distance_char(const char *u, const char *v, const npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, nff = 0, nft = 0, ntf = 0;

    for (i = 0; i < n; ++i) {
        const npy_bool x = (u[i] != 0), y = (v[i] != 0);
        ntt += x & y;
        ntf += x & (!y);
        nft += (!x) & y;
    }
    nff = n - ntt - ntf - nft;
    double half_R = (double)ntf * nft;
    if (half_R == 0.0) {
        return 0.0;
    }
    return (2. * half_R) / ((double)ntt * nff + half_R);
}

static NPY_INLINE double
dice_distance_char(const char *u, const char *v, const npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, ndiff = 0;

    for (i = 0; i < n; ++i) {
        const npy_bool x = (u[i] != 0), y = (v[i] != 0);
        ntt += x & y;
        ndiff += (x != y);
    }
    return ndiff / (2. * ntt + ndiff);
}

static NPY_INLINE double
rogerstanimoto_distance_char(const char *u, const char *v, const npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, ndiff = 0;

    for (i = 0; i < n; ++i) {
        const npy_bool x = (u[i] != 0), y = (v[i] != 0);
        ntt += x & y;
        ndiff += (x != y);
    }
    return (2. * ndiff) / ((double)n + ndiff);
}

static NPY_INLINE double
russellrao_distance_char(const char *u, const char *v, const npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0;

    for (i = 0; i < n; ++i) {
        ntt += (u[i] != 0) & (v[i] != 0);
    }
    return (double)(n - ntt) / n;
}

static NPY_INLINE double
kulsinski_distance_char(const char *u, const char *v, const npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, ndiff = 0;

    for (i = 0; i < n; ++i) {
        const npy_bool x = (u[i] != 0), y = (v[i] != 0);
        ntt += x & y;
        ndiff += (x != y);
    }
    return ((double)ndiff - ntt + n) / ((double)ndiff + n);
}

static NPY_INLINE double
sokalsneath_distance_char(const char *u, const char *v, const npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, ndiff = 0;

    for (i = 0; i < n; ++i) {
        const npy_bool x = (u[i] != 0), y = (v[i] != 0);
        ntt += x & y;
        ndiff += (x != y);
    }
    return (2. * ndiff) / (2. * ndiff + ntt);
}

static NPY_INLINE double
sokalmichener_distance_char(const char *u, const char *v, const npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, ndiff = 0;

    for (i = 0; i < n; ++i) {
        const npy_bool x = (u[i] != 0), y = (v[i] != 0);
        ntt += x & y;
        ndiff += (x != y);
    }
    return (2. * ndiff) / ((double)ndiff + n);
}

static NPY_INLINE double
jaccard_distance_double(const double *u, const double *v, const npy_intp n)
{
    npy_intp denom = 0, num = 0;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        const double x = u[i], y = v[i];
        num += (x != y) & ((x != 0.0) | (y != 0.0));
        denom += (x != 0.0) | (y != 0.0);
    }
    return denom == 0.0 ? 0.0 : (double)num / denom;
}

static NPY_INLINE double
jaccard_distance_char(const char *u, const char *v, const npy_intp n)
{
    npy_intp num = 0, denom = 0;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        const npy_bool x = (u[i] != 0), y = (v[i] != 0);
        num += (x != y);
        denom += x | y;
    }
    return denom == 0.0 ? 0.0 : (double)num / denom;
}


static NPY_INLINE double
jensenshannon_distance_double(const double *p, const double *q, const npy_intp n)
{
    npy_intp i;
    double s = 0.0;
    double p_sum = 0.0;
    double q_sum = 0.0;

    for (i = 0; i < n; ++i) {
        if (p[i] < 0.0 || q[i] < 0.0)
            return HUGE_VAL;
        p_sum += p[i];
        q_sum += q[i];
    }

    if (p_sum == 0.0 || q_sum == 0.0)
        return HUGE_VAL;

    for (i = 0; i < n; ++i) {
        const double p_i = p[i] / p_sum;
        const double q_i = q[i] / q_sum;
        const double m_i = (p_i + q_i) / 2.0;
        if (p_i > 0.0)
            s += p_i * log(p_i / m_i);
        if (q_i > 0.0)
            s += q_i * log(q_i / m_i);
    }

    return sqrt(s / 2.0);
}


static NPY_INLINE double
seuclidean_distance(const double *var, const double *u, const double *v,
                    const npy_intp n)
{
    double s = 0.0;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        const double d = u[i] - v[i];
        s += (d * d) / var[i];
    }
    return sqrt(s);
}

static NPY_INLINE double
city_block_distance_double(const double *u, const double *v, const npy_intp n)
{
    double s = 0.0;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        s += fabs(u[i] - v[i]);
    }
    return s;
}

static NPY_INLINE double
minkowski_distance(const double *u, const double *v, const npy_intp n, const double p)
{
    double s = 0.0;
    npy_intp i;

    for (i = 0; i < n; ++i) {
        const double d = fabs(u[i] - v[i]);
        s += pow(d, p);
    }
    return pow(s, 1.0 / p);
}

static NPY_INLINE double
weighted_minkowski_distance(const double *u, const double *v, const npy_intp n,
                            const double p, const double *w)
{
    npy_intp i = 0;
    double s = 0.0;
    for (i = 0; i < n; ++i) {
        const double d = fabs(u[i] - v[i]);
        s += pow(d, p) * w[i];
    }
    return pow(s, 1.0 / p);
}

#if 0   /* XXX unused */
static void
compute_mean_vector(double *res, const double *X, npy_intp num_rows, const npy_intp n)
{
    npy_intp i, j;
    const double *v;
    for (i = 0; i < n; ++i) {
        res[i] = 0.0;
    }
    for (j = 0; j < num_rows; ++j) {

        v = X + (j * n);
        for (i = 0; i < n; ++i) {
            res[i] += v[i];
        }
    }
    for (i = 0; i < n; ++i) {
        res[i] /= (double)num_rows;
    }
}
#endif

#define DEFINE_CDIST(name, type) \
    static int cdist_ ## name ## _ ## type(const type *XA, const type *XB, \
                                           double *dm,                     \
                                           const npy_intp num_rowsA,       \
                                           const npy_intp num_rowsB,       \
                                           const npy_intp num_cols)        \
    {                                                                      \
        Py_ssize_t i, j;                                                   \
        for (i = 0; i < num_rowsA; ++i) {                                  \
            const type *u = XA + num_cols * i;                             \
            for (j = 0; j < num_rowsB; ++j, ++dm) {                        \
                const type *v = XB + num_cols * j;                         \
                *dm = name ## _distance_ ## type(u, v, num_cols);          \
            }                                                              \
        }                                                                  \
        return 0;\
    }

DEFINE_CDIST(bray_curtis, double)
DEFINE_CDIST(canberra, double)
DEFINE_CDIST(chebyshev, double)
DEFINE_CDIST(city_block, double)
DEFINE_CDIST(euclidean, double)
DEFINE_CDIST(jaccard, double)
DEFINE_CDIST(jensenshannon, double)
DEFINE_CDIST(sqeuclidean, double)

DEFINE_CDIST(dice, char)
DEFINE_CDIST(jaccard, char)
DEFINE_CDIST(kulsinski, char)
DEFINE_CDIST(rogerstanimoto, char)
DEFINE_CDIST(russellrao, char)
DEFINE_CDIST(sokalmichener, char)
DEFINE_CDIST(sokalsneath, char)
DEFINE_CDIST(yule, char)


#define DEFINE_PDIST(name, type) \
    static int pdist_ ## name ## _ ## type(const type *X, double *dm,       \
                                           const npy_intp num_rows,         \
                                           const npy_intp num_cols)         \
    {                                                                       \
        Py_ssize_t i, j;                                                    \
        double *it = dm;                                                    \
        for (i = 0; i < num_rows; ++i) {                                    \
            const type *u = X + num_cols * i;                               \
            for (j = i + 1; j < num_rows; ++j, it++) {                      \
                const type *v = X + num_cols * j;                           \
                *it = name ## _distance_ ## type(u, v, num_cols);           \
            }                                                               \
        }                                                                   \
        return 0; \
    }

DEFINE_PDIST(bray_curtis, double)
DEFINE_PDIST(canberra, double)
DEFINE_PDIST(chebyshev, double)
DEFINE_PDIST(city_block, double)
DEFINE_PDIST(euclidean, double)
DEFINE_PDIST(jaccard, double)
DEFINE_PDIST(jensenshannon, double)
DEFINE_PDIST(sqeuclidean, double)

DEFINE_PDIST(dice, char)
DEFINE_PDIST(jaccard, char)
DEFINE_PDIST(kulsinski, char)
DEFINE_PDIST(rogerstanimoto, char)
DEFINE_PDIST(russellrao, char)
DEFINE_PDIST(sokalmichener, char)
DEFINE_PDIST(sokalsneath, char)
DEFINE_PDIST(yule, char)

static NPY_INLINE int
pdist_mahalanobis(const double *X, double *dm, const npy_intp num_rows,
                  const npy_intp num_cols, const double *covinv)
{
    npy_intp i, j;
    double *dimbuf1 = calloc(2 * num_cols, sizeof(double));
    double *dimbuf2;
    if (!dimbuf1) {
        return -1;
    }

    dimbuf2 = dimbuf1 + num_cols;

    for (i = 0; i < num_rows; ++i) {
        const double *u = X + (num_cols * i);
        for (j = i + 1; j < num_rows; ++j, ++dm) {
            const double *v = X + (num_cols * j);
            *dm = mahalanobis_distance(u, v, covinv, dimbuf1, dimbuf2, num_cols);
        }
    }
    free(dimbuf1);
    return 0;
}

static NPY_INLINE int
pdist_weighted_chebyshev(const double *X, double *dm, npy_intp num_rows,
                         const npy_intp num_cols, const double *w)
{
    npy_intp i, j;

    for (i = 0; i < num_rows; ++i) {
        const double *u = X + (num_cols * i);
        for (j = i + 1; j < num_rows; ++j, ++dm) {
            const double *v = X + (num_cols * j);
            *dm = weighted_chebyshev_distance_double(u, v, num_cols, w);
        }
    }
    return 0;
}

static NPY_INLINE int
pdist_cosine(const double *X, double *dm, const npy_intp num_rows,
             const npy_intp num_cols)
{
    double cosine;
    npy_intp i, j;

    double * norms_buff = calloc(num_rows, sizeof(double));
    if (!norms_buff)
        return -1;

    _row_norms(X, num_rows, num_cols, norms_buff);

    for (i = 0; i < num_rows; ++i) {
        const double *u = X + (num_cols * i);
        for (j = i + 1; j < num_rows; ++j, ++dm) {
            const double *v = X + (num_cols * j);
            cosine = dot_product(u, v, num_cols) / (norms_buff[i] * norms_buff[j]);
            if (fabs(cosine) > 1.) {
                /* Clip to correct rounding error. */
                cosine = npy_copysign(1, cosine);
            }
            *dm = 1. - cosine;
        }
    }
    free(norms_buff);
    return 0;
}

static NPY_INLINE int
pdist_seuclidean(const double *X, const double *var, double *dm,
                 const npy_intp num_rows, const npy_intp num_cols)
{
    npy_intp i, j;

    for (i = 0; i < num_rows; ++i) {
        const double *u = X + (num_cols * i);
        for (j = i + 1; j < num_rows; ++j, ++dm) {
            const double *v = X + (num_cols * j);
            *dm = seuclidean_distance(var, u, v, num_cols);
        }
    }
    return 0;
}

static NPY_INLINE int
pdist_minkowski(const double *X, double *dm, npy_intp num_rows,
                const npy_intp num_cols, const double p)
{
    npy_intp i, j;
    if (p == 1.0) {
        return pdist_city_block_double(X, dm, num_rows, num_cols);
    }
    if (p == 2.0) {
        return pdist_euclidean_double(X, dm, num_rows, num_cols);
    }
    if (sc_isinf(p)) {
        return pdist_chebyshev_double(X, dm, num_rows, num_cols);
    }

    for (i = 0; i < num_rows; ++i) {
        const double *u = X + (num_cols * i);
        for (j = i + 1; j < num_rows; ++j, ++dm) {
            const double *v = X + (num_cols * j);
            *dm = minkowski_distance(u, v, num_cols, p);
        }
    }
    return 0;
}

/* Old weighting type which is inconsistent with other distance metrics.
   Remove in SciPy 1.8 along with wminkowski. */
static NPY_INLINE int
pdist_old_weighted_minkowski(const double *X, double *dm, npy_intp num_rows,
                             const npy_intp num_cols, const double p, const double *w)
{
    npy_intp i, j;

    // Covert from old style weights to new weights
    double * new_weights = malloc(num_cols * sizeof(double));
    if (!new_weights) {
        return 1;
    }
    for (i = 0; i < num_cols; ++i) {
        new_weights[i] = pow(w[i], p);
    }

    for (i = 0; i < num_rows; ++i) {
        const double *u = X + (num_cols * i);
        for (j = i + 1; j < num_rows; ++j, ++dm) {
            const double *v = X + (num_cols * j);
            *dm = weighted_minkowski_distance(u, v, num_cols, p, new_weights);
        }
    }
    free(new_weights);
    return 0;
}

static NPY_INLINE int
pdist_weighted_minkowski(const double *X, double *dm, npy_intp num_rows,
                         const npy_intp num_cols, const double p, const double *w)
{
    npy_intp i, j;

    if (sc_isinf(p)) {
        return pdist_weighted_chebyshev(X, dm, num_rows, num_cols, w);
    }

    for (i = 0; i < num_rows; ++i) {
        const double *u = X + (num_cols * i);
        for (j = i + 1; j < num_rows; ++j, ++dm) {
            const double *v = X + (num_cols * j);
            *dm = weighted_minkowski_distance(u, v, num_cols, p, w);
        }
    }
    return 0;
}

static NPY_INLINE int
pdist_hamming_double(const double *X, double *dm, npy_intp num_rows,
                         const npy_intp num_cols, const double *w)
{
    npy_intp i, j;

    for (i = 0; i < num_rows; ++i) {
        const double *u = X + (num_cols * i);
        for (j = i + 1; j < num_rows; ++j, ++dm) {
            const double *v = X + (num_cols * j);
            *dm = hamming_distance_double(u, v, num_cols, w);
        }
    }
    return 0;
}

static NPY_INLINE int
pdist_hamming_char(const char *X, double *dm, npy_intp num_rows,
                         const npy_intp num_cols, const double *w)
{
    npy_intp i, j;

    for (i = 0; i < num_rows; ++i) {
        const char *u = X + (num_cols * i);
        for (j = i + 1; j < num_rows; ++j, ++dm) {
            const char *v = X + (num_cols * j);
            *dm = hamming_distance_char(u, v, num_cols, w);
        }
    }
    return 0;
}

static NPY_INLINE void
dist_to_squareform_from_vector_generic(char *M, const char *v, const npy_intp n,
                                       npy_intp s)
{
    char *it1 = M + s;
    char *it2;
    npy_intp i, j;

    for (i = 1; i < n; ++i) {
        memcpy(it1, v, (n - i) * s);
        it1 += (n + 1) * s;

        it2 = M + i * (n + 1) * s - s;
        for (j = i; j < n; ++j) {
            memcpy(it2, v, s);
            v += s;
            it2 += n*s;
        }
    }
}

static NPY_INLINE void
dist_to_squareform_from_vector_double(double *M, const double *v, const npy_intp n)
{
    double *it1 = M + 1;
    double *it2;
    npy_intp i, j;

    for (i = 1; i < n; ++i, it1 += n+1) {
        memcpy(it1, v, (n - i) * sizeof(double));

        it2 = M + i * (n + 1) - 1;
        for (j = i; j < n; ++j, ++v, it2 += n) {
            *it2 = *v;
        }
    }
}

static NPY_INLINE void
dist_to_vector_from_squareform(const char *M, char *v, const npy_intp n, npy_intp s)
{
    const char *cit = M + s;
    npy_intp i;

    for (i = 1; i < n; ++i) {
        const npy_intp len = (n - i) * s;
        memcpy(v, cit, len);
        v += len;
        cit += (n + 1) * s;
    }
}

static NPY_INLINE int
cdist_weighted_chebyshev(const double *XA, const double *XB, double *dm,
                         const npy_intp num_rowsA, const npy_intp num_rowsB,
                         const npy_intp num_cols, const double *w)
{
    npy_intp i, j;

    for (i = 0; i < num_rowsA; ++i) {
        const double *u = XA + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            const double *v = XB + (num_cols * j);
            *dm = weighted_chebyshev_distance_double(u, v, num_cols, w);
        }
    }
    return 0;
}


/** cdist */
static NPY_INLINE int
cdist_cosine(const double *XA, const double *XB, double *dm, const npy_intp num_rowsA,
             const npy_intp num_rowsB, const npy_intp num_cols)
{
    double cosine;
    npy_intp i, j;

    double * norms_buffA = calloc(num_rowsA + num_rowsB, sizeof(double));
    double * norms_buffB;
    if (!norms_buffA)
        return -1;

    norms_buffB = norms_buffA + num_rowsA;

    _row_norms(XA, num_rowsA, num_cols, norms_buffA);
    _row_norms(XB, num_rowsB, num_cols, norms_buffB);

    for (i = 0; i < num_rowsA; ++i) {
        const double *u = XA + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            const double *v = XB + (num_cols * j);
            cosine = dot_product(u, v, num_cols) / (norms_buffA[i] * norms_buffB[j]);
            if (fabs(cosine) > 1.) {
                /* Clip to correct rounding error. */
                cosine = npy_copysign(1, cosine);
            }
            *dm = 1. - cosine;
        }
    }
    free(norms_buffA);
    return 0;
}

static NPY_INLINE int
cdist_mahalanobis(const double *XA, const double *XB, double *dm,
                  const npy_intp num_rowsA, const npy_intp num_rowsB,
                  const npy_intp num_cols, const double *covinv)
{
    npy_intp i, j;

    double *dimbuf1 = calloc(2 * num_cols, sizeof(double));
    double *dimbuf2;
    if (!dimbuf1) {
        return -1;
    }
    dimbuf2 = dimbuf1 + num_cols;

    for (i = 0; i < num_rowsA; ++i) {
        const double *u = XA + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            const double *v = XB + (num_cols * j);
            *dm = mahalanobis_distance(u, v, covinv, dimbuf1, dimbuf2, num_cols);
        }
    }
    free(dimbuf1);
    return 0;
}

static NPY_INLINE int
cdist_seuclidean(const double *XA, const double *XB, const double *var,
                 double *dm, const npy_intp num_rowsA, const npy_intp num_rowsB,
                 const npy_intp num_cols)
{
    npy_intp i, j;

    for (i = 0; i < num_rowsA; ++i) {
        const double *u = XA + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            const double *v = XB + (num_cols * j);
            *dm = seuclidean_distance(var, u, v, num_cols);
        }
    }
    return 0;
}

static NPY_INLINE int
cdist_minkowski(const double *XA, const double *XB, double *dm,
                const npy_intp num_rowsA, const npy_intp num_rowsB,
                const npy_intp num_cols, const double p)
{
    npy_intp i, j;
    if (p == 1.0) {
        return cdist_city_block_double(XA, XB, dm, num_rowsA, num_rowsB, num_cols);
    }
    if (p == 2.0) {
        return cdist_euclidean_double(XA, XB, dm, num_rowsA, num_rowsB, num_cols);
    }
    if (sc_isinf(p)) {
        return cdist_chebyshev_double(XA, XB, dm, num_rowsA, num_rowsB, num_cols);
    }

    for (i = 0; i < num_rowsA; ++i) {
        const double *u = XA + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            const double *v = XB + (num_cols * j);
            *dm = minkowski_distance(u, v, num_cols, p);
        }
    }
    return 0;
}

/* Old weighting type which is inconsistent with other distance metrics.
   Remove in SciPy 1.8 along with wminkowski. */
static NPY_INLINE int
cdist_old_weighted_minkowski(const double *XA, const double *XB, double *dm,
                             const npy_intp num_rowsA, const npy_intp num_rowsB,
                             const npy_intp num_cols, const double p,
                             const double *w)
{
    npy_intp i, j;

    // Covert from old style weights to new weights
    double * new_weights = malloc(num_cols * sizeof(double));
    if (!new_weights) {
      return 1;
    }

    for (i = 0; i < num_cols; ++i) {
      new_weights[i] = pow(w[i], p);
    }

    for (i = 0; i < num_rowsA; ++i) {
        const double *u = XA + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            const double *v = XB + (num_cols * j);
            *dm = weighted_minkowski_distance(u, v, num_cols, p, new_weights);
        }
    }
    free(new_weights);
    return 0;
}

static NPY_INLINE int
cdist_weighted_minkowski(const double *XA, const double *XB, double *dm,
                         const npy_intp num_rowsA, const npy_intp num_rowsB,
                         const npy_intp num_cols, const double p,
                         const double *w)
{
    npy_intp i, j;

    if (sc_isinf(p)) {
        return cdist_weighted_chebyshev(XA, XB, dm, num_rowsA, num_rowsB, num_cols, w);
    }

    for (i = 0; i < num_rowsA; ++i) {
        const double *u = XA + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            const double *v = XB + (num_cols * j);
            *dm = weighted_minkowski_distance(u, v, num_cols, p, w);
        }
    }
    return 0;
}

static NPY_INLINE int
cdist_hamming_double(const double *XA, const double *XB, double *dm,
                         const npy_intp num_rowsA, const npy_intp num_rowsB,
                         const npy_intp num_cols,
                         const double *w)
{
    npy_intp i, j;

    for (i = 0; i < num_rowsA; ++i) {
        const double *u = XA + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            const double *v = XB + (num_cols * j);
            *dm = hamming_distance_double(u, v, num_cols, w);
        }
    }
    return 0;
}

static NPY_INLINE int
cdist_hamming_char(const char *XA, const char *XB, double *dm,
                         const npy_intp num_rowsA, const npy_intp num_rowsB,
                         const npy_intp num_cols,
                         const double *w)
{
    npy_intp i, j;

    for (i = 0; i < num_rowsA; ++i) {
        const char *u = XA + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            const char *v = XB + (num_cols * j);
            *dm = hamming_distance_char(u, v, num_cols, w);
        }
    }
    return 0;
}
