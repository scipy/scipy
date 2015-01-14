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

#ifdef __SSE2__
#include <emmintrin.h>

static NPY_INLINE double horizontal_add(__m128d a)
{
    double r;
    __m128d ha = _mm_unpackhi_pd(a, a);
    a = _mm_add_sd(a, ha);
    _mm_store_sd(&r, a);
    return r;
}
#endif

static NPY_INLINE double
sqeuclidean_distance_double(const double *u, const double *v, npy_intp n)
{
    double s = 0.0;
    npy_intp i = 0;
#ifdef __SSE2__
    __m128d sv = _mm_setzero_pd();

    for (; i < n - (n & 1); i+=2) {
        __m128d uv = _mm_loadu_pd(&u[i]);
        __m128d vv = _mm_loadu_pd(&v[i]);
        __m128d dv = _mm_sub_pd(uv, vv);
        dv = _mm_mul_pd(dv, dv);
        sv = _mm_add_pd(sv, dv);
    }
    s = horizontal_add(sv);
#endif

    for (; i < n; i++) {
        double d = u[i] - v[i];
        s += d * d;
    }

    return s;
}

static NPY_INLINE double
euclidean_distance_double(const double *u, const double *v, npy_intp n)
{
    return sqrt(sqeuclidean_distance_double(u, v, n));
}

#if 0   /* XXX unused */
static NPY_INLINE double
ess_distance(const double *u, const double *v, npy_intp n)
{
    double s = 0.0, d;
    npy_intp i;

    for (i = 0; i < n; i++) {
        d = fabs(u[i] - v[i]);
        s += d * d;
    }
    return s;
}
#endif

static NPY_INLINE double
chebyshev_distance_double(const double *u, const double *v, npy_intp n)
{
    double d, maxv = 0.0;
    npy_intp i;

    for (i = 0; i < n; i++) {
        d = fabs(u[i] - v[i]);
        if (d > maxv) {
            maxv = d;
        }
    }
    return maxv;
}

static NPY_INLINE double
canberra_distance_double(const double *u, const double *v, npy_intp n)
{
    double snum = 0.0, sdenom = 0.0, tot = 0.0;
    npy_intp i;

    for (i = 0; i < n; i++) {
        snum = fabs(u[i] - v[i]);
        sdenom = fabs(u[i]) + fabs(v[i]);
        if (sdenom > 0.0) {
            tot += snum / sdenom;
        }
    }
    return tot;
}

static NPY_INLINE double
bray_curtis_distance_double(const double *u, const double *v, npy_intp n)
{
    double s1 = 0.0, s2 = 0.0;
    npy_intp i;

    for (i = 0; i < n; i++) {
        s1 += fabs(u[i] - v[i]);
        s2 += fabs(u[i] + v[i]);
    }
    return s1 / s2;
}

static NPY_INLINE double
mahalanobis_distance(const double *u, const double *v, const double *covinv,
                     double *dimbuf1, double *dimbuf2, npy_intp n)
{
    const double *covrow = covinv;
    double s;
    npy_intp i, j;

    for (i = 0; i < n; i++) {
        dimbuf1[i] = u[i] - v[i];
    }
    for (i = 0; i < n; i++) {
        covrow = covinv + (i * n);
        s = 0.0;
        for (j = 0; j < n; j++) {
            s += dimbuf1[j] * covrow[j];
        }
        dimbuf2[i] = s;
    }
    s = 0.0;
    for (i = 0; i < n; i++) {
        s += dimbuf1[i] * dimbuf2[i];
    }
    return sqrt(s);
}

static NPY_INLINE double
hamming_distance_double(const double *u, const double *v, npy_intp n)
{
    double s = 0.0;
    npy_intp i = 0;
#ifdef __SSE2__
    __m128d sv = _mm_setzero_pd();
    __m128d one = _mm_set1_pd(1.);

    for (; i < n - (n & 1); i+=2) {
        __m128d uv = _mm_loadu_pd(&u[i]);
        __m128d vv = _mm_loadu_pd(&v[i]);
        /* all bit 1 if != else all bit 0 */
        __m128d dv = _mm_cmpneq_pd(uv, vv);
        /* mask 1 to 0 if all bit 0, else keep 1 */
        __m128d oneorzero = _mm_and_pd(dv, one);
        sv = _mm_add_pd(sv, oneorzero);
    }
    s = horizontal_add(sv);
#endif

    for (; i < n; i++) {
        s += (u[i] != v[i]);
    }

    return s / n;
}

static NPY_INLINE double
hamming_distance_char(const char *u, const char *v, npy_intp n)
{
    double s = 0.0;
    npy_intp i;

    for (i = 0; i < n; i++) {
        s += (u[i] != v[i]);
    }
    return s / n;
}

static NPY_INLINE double
yule_bool_distance_char(const char *u, const char *v, npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, nff = 0, nft = 0, ntf = 0;

    for (i = 0; i < n; i++) {
        ntt += (u[i] && v[i]);
        ntf += (u[i] && !v[i]);
        nft += (!u[i] && v[i]);
        nff += (!u[i] && !v[i]);
    }
    return (2.0 * ntf * nft) / ((double)ntt * nff + (double)ntf * nft);  
}

static NPY_INLINE double
matching_distance_char(const char *u, const char *v, npy_intp n)
{
    npy_intp i;
    npy_intp nft = 0, ntf = 0;

    for (i = 0; i < n; i++) {
        ntf += (u[i] && !v[i]);
        nft += (!u[i] && v[i]);
    }
    return (double)(ntf + nft) / (double)(n);
}

static NPY_INLINE double
dice_distance_char(const char *u, const char *v, npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, nft = 0, ntf = 0;

    for (i = 0; i < n; i++) {
        ntt += (u[i] && v[i]);
        ntf += (u[i] && !v[i]);
        nft += (!u[i] && v[i]);
    }
    return (double)(nft + ntf) / (2.0 * ntt + ntf + nft);
}


static NPY_INLINE double
rogerstanimoto_distance_char(const char *u, const char *v, npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, nff = 0, nft = 0, ntf = 0;

    for (i = 0; i < n; i++) {
        ntt += (u[i] && v[i]);
        ntf += (u[i] && !v[i]);
        nft += (!u[i] && v[i]);
        nff += (!u[i] && !v[i]);
    }
    return (2.0 * (ntf + nft)) / ((double)ntt + nff + (2.0 * (ntf + nft)));
}

static NPY_INLINE double
russellrao_distance_char(const char *u, const char *v, npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0;

    for (i = 0; i < n; i++) {
        ntt += (u[i] && v[i]);
    }
    return (double)(n - ntt) / n;
}

static NPY_INLINE double
kulsinski_distance_char(const char *u, const char *v, npy_intp n)
{
    npy_intp _i;
    npy_intp ntt = 0, nft = 0, ntf = 0, nff = 0;

    for (_i = 0; _i < n; _i++) {
        ntt += (u[_i] && v[_i]);
        ntf += (u[_i] && !v[_i]);
        nft += (!u[_i] && v[_i]);
        nff += (!u[_i] && !v[_i]);
    }
    return (double)(ntf + nft - ntt + n) / (double)(ntf + nft + n);
}

static NPY_INLINE double
sokalsneath_distance_char(const char *u, const char *v, npy_intp n)
{
    npy_intp _i;
    npy_intp ntt = 0, nft = 0, ntf = 0;

    for (_i = 0; _i < n; _i++) {
        ntt += (u[_i] && v[_i]);
        ntf += (u[_i] && !v[_i]);
        nft += (!u[_i] && v[_i]);
    }
    return (2.0 * (ntf + nft)) / (2.0 * (ntf + nft) + ntt);
}

static NPY_INLINE double
sokalmichener_distance_char(const char *u, const char *v, npy_intp n)
{
    npy_intp i;
    npy_intp ntt = 0, nft = 0, ntf = 0, nff = 0;

    for (i = 0; i < n; i++) {
        ntt += (u[i] && v[i]);
        nff += (!u[i] && !v[i]);
        ntf += (u[i] && !v[i]);
        nft += (!u[i] && v[i]);
    }
    return (2.0 * (ntf + nft)) / (2.0 * (ntf + nft) + ntt + nff);
}

static NPY_INLINE double
jaccard_distance_double(const double *u, const double *v, npy_intp n)
{
    double denom = 0.0, num = 0.0;
    npy_intp i;

    for (i = 0; i < n; i++) {
        num += (u[i] != v[i]) & ((u[i] != 0.0) | (v[i] != 0.0));
        denom += (u[i] != 0.0) | (v[i] != 0.0);
    }
    return num / denom;
}

static NPY_INLINE double
jaccard_distance_char(const char *u, const char *v, npy_intp n)
{
    double num = 0.0, denom = 0.0;
    npy_intp i;

    for (i = 0; i < n; i++) {
        num += (u[i] != v[i]) & ((u[i] != 0) | (v[i] != 0));
        denom += (u[i] != 0) | (v[i] != 0);
    }
    return num / denom;
}

/* XXX shouldn't we use BLAS for this? */
static NPY_INLINE double
dot_product(const double *u, const double *v, npy_intp n)
{
    double s = 0.0;
    npy_intp i = 0;
#ifdef __SSE2__
    __m128d sv = _mm_setzero_pd();

    for (; i < n - (n & 1); i+=2) {
        __m128d uv = _mm_loadu_pd(&u[i]);
        __m128d vv = _mm_loadu_pd(&v[i]);
        __m128d dv = _mm_mul_pd(uv, vv);
        sv = _mm_add_pd(sv, dv);
    }
    s = horizontal_add(sv);
#endif

    for (; i < n; i++) {
        s += u[i] * v[i];
    }

    return s;
}

static NPY_INLINE double
cosine_distance(const double *u, const double *v, npy_intp n,
                double nu, double nv)
{
    return 1.0 - (dot_product(u, v, n) / (nu * nv));
}

static NPY_INLINE double
seuclidean_distance(const double *var, const double *u, const double *v,
                    npy_intp n)
{
    double s = 0.0;
    npy_intp i = 0;
#ifdef __SSE2__
    __m128d sv = _mm_setzero_pd();

    for (; i < n - (n & 1); i+=2) {
        __m128d uv = _mm_loadu_pd(&u[i]);
        __m128d vv = _mm_loadu_pd(&v[i]);
        __m128d varv = _mm_loadu_pd(&var[i]);
        __m128d dv = _mm_sub_pd(uv, vv);
        dv = _mm_mul_pd(dv, dv);
        dv = _mm_div_pd(dv, varv);
        sv = _mm_add_pd(sv, dv);
    }
    s = horizontal_add(sv);
#endif

    for (; i < n; i++) {
        double d = u[i] - v[i];
        s = s + (d * d) / var[i];
    }

    return sqrt(s);
}

static NPY_INLINE double
city_block_distance_double(const double *u, const double *v, npy_intp n)
{
    double s = 0.0, d;
    npy_intp i;

    for (i = 0; i < n; i++) {
        d = fabs(u[i] - v[i]);
        s = s + d;
    }
    return s;
}

static NPY_INLINE double
minkowski_distance(const double *u, const double *v, npy_intp n, double p)
{
    double s = 0.0, d;
    npy_intp i;

    for (i = 0; i < n; i++) {
        d = fabs(u[i] - v[i]);
        s = s + pow(d, p);
    }
    return pow(s, 1.0 / p);
}

static NPY_INLINE double
weighted_minkowski_distance(const double *u, const double *v, npy_intp n,
                            double p, const double *w)
{
  npy_intp i = 0;
  double s = 0.0, d;
  for (i = 0; i < n; i++) {
    d = fabs(u[i] - v[i]) * w[i];
    s = s + pow(d, p);
  }
  return pow(s, 1.0 / p);
}

#if 0   /* XXX unused */
static void
compute_mean_vector(double *res, const double *X, npy_intp m, npy_intp n)
{
    npy_intp i, j;
    const double *v;
    for (i = 0; i < n; i++) {
        res[i] = 0.0;
    }
    for (j = 0; j < m; j++) {

        v = X + (j * n);
        for (i = 0; i < n; i++) {
            res[i] += v[i];
        }
    }
    for (i = 0; i < n; i++) {
        res[i] /= (double)m;
    }
}
#endif

static NPY_INLINE void
pdist_mahalanobis(const double *X, const double *covinv, double *dimbuf,
                  double *dm, npy_intp m, npy_intp n)
{
    npy_intp i, j;
    const double *u, *v;

    double *dimbuf1 = dimbuf;
    double *dimbuf2 = dimbuf + n;

    for (i = 0; i < m; i++) {
        for (j = i + 1; j < m; j++, dm++) {
            u = X + (n * i);
            v = X + (n * j);
            *dm = mahalanobis_distance(u, v, covinv, dimbuf1, dimbuf2, n);
        }
    }
    dimbuf2 = 0;
}

static NPY_INLINE void
pdist_cosine(const double *X, double *dm, npy_intp m, npy_intp n,
             const double *norms)
{
    const double *u, *v;
    npy_intp i, j;

    for (i = 0; i < m; i++) {
        for (j = i + 1; j < m; j++, dm++) {
            u = X + (n * i);
            v = X + (n * j);
            *dm = cosine_distance(u, v, n, norms[i], norms[j]);
        }
    }
}

static NPY_INLINE void
pdist_seuclidean(const double *X, const double *var, double *dm,
                 npy_intp m, npy_intp n)
{
    const double *u, *v;
    npy_intp i, j;

    for (i = 0; i < m; i++) {
        for (j = i + 1; j < m; j++, dm++) {
            u = X + (n * i);
            v = X + (n * j);
            *dm = seuclidean_distance(var, u, v, n);
        }
    }
}

static NPY_INLINE void
pdist_minkowski(const double *X, double *dm, npy_intp m, npy_intp n, double p)
{
    const double *u, *v;
    npy_intp i, j;

    for (i = 0; i < m; i++) {
        for (j = i + 1; j < m; j++, dm++) {
            u = X + (n * i);
            v = X + (n * j);
            *dm = minkowski_distance(u, v, n, p);
        }
    }
}

static NPY_INLINE void
pdist_weighted_minkowski(const double *X, double *dm, npy_intp m, npy_intp n,
                         double p, const double *w)
{
    const double *u, *v;
    npy_intp i, j;

    for (i = 0; i < m; i++) {
        for (j = i + 1; j < m; j++, dm++) {
            u = X + (n * i);
            v = X + (n * j);
            *dm = weighted_minkowski_distance(u, v, n, p, w);
        }
    }
}

static NPY_INLINE void
dist_to_squareform_from_vector(double *M, const double *v, npy_intp n)
{
    double *it;
    npy_intp i, j;

    for (i = 0; i < n - 1; i++) {
        it = M + (i * n) + i + 1;
        for (j = i + 1; j < n; j++, it++, v++) {
            *it = *v;
        }
    }
}

static NPY_INLINE void
dist_to_vector_from_squareform(const double *M, double *v, npy_intp n)
{
    const double *cit;
    npy_intp i, j;

    for (i = 0; i < n - 1; i++) {
        cit = M + (i * n) + i + 1;
        for (j = i + 1; j < n; j++, v++, cit++) {
            *v = *cit;
        }
    }
}


/** cdist */

static NPY_INLINE void
cdist_mahalanobis(const double *XA, const double *XB,
                  const double *covinv, double *dimbuf,
                  double *dm, npy_intp mA, npy_intp mB, npy_intp n)
{
    npy_intp i, j;
    const double *u, *v;

    double *dimbuf1 = dimbuf;
    double *dimbuf2 = dimbuf + n;

    for (i = 0; i < mA; i++) {
        for (j = 0; j < mB; j++, dm++) {
            u = XA + (n * i);
            v = XB + (n * j);
            *dm = mahalanobis_distance(u, v, covinv, dimbuf1, dimbuf2, n);
        }
    }
    dimbuf2 = 0;
}

static NPY_INLINE void
cdist_cosine(const double *XA, const double *XB, double *dm,
             npy_intp mA, npy_intp mB, npy_intp n,
             const double *normsA, const double *normsB)
{
    const double *u, *v;
    npy_intp i, j;

    for (i = 0; i < mA; i++) {
        for (j = 0; j < mB; j++, dm++) {
            u = XA + (n * i);
            v = XB + (n * j);
            *dm = cosine_distance(u, v, n, normsA[i], normsB[j]);
        }
    }
}

static NPY_INLINE void
cdist_seuclidean(const double *XA, const double *XB, const double *var,
                 double *dm, npy_intp mA, npy_intp mB, npy_intp n)
{
    const double *u, *v;
    npy_intp i, j;

    for (i = 0; i < mA; i++) {
        for (j = 0; j < mB; j++, dm++) {
            u = XA + (n * i);
            v = XB + (n * j);
            *dm = seuclidean_distance(var, u, v, n);
        }
    }
}

static NPY_INLINE void
cdist_minkowski(const double *XA, const double *XB, double *dm,
                npy_intp mA, npy_intp mB, npy_intp n, double p)
{
    const double *u, *v;
    npy_intp i, j;

    for (i = 0; i < mA; i++) {
        for (j = 0; j < mB; j++, dm++) {
            u = XA + (n * i);
            v = XB + (n * j);
            *dm = minkowski_distance(u, v, n, p);
        }
    }
}

static NPY_INLINE void
cdist_weighted_minkowski(const double *XA, const double *XB, double *dm,
                         npy_intp mA, npy_intp mB, npy_intp n, double p,
                         const double *w)
{
    npy_intp i, j;
    const double *u, *v;

    for (i = 0; i < mA; i++) {
        for (j = 0; j < mB; j++, dm++) {
            u = XA + (n * i);
            v = XB + (n * j);
            *dm = weighted_minkowski_distance(u, v, n, p, w);
        }
    }
}

#define DEFINE_CDIST(name, type) \
    static void cdist_ ## name ## _ ## type(const type *XA, const type *XB, \
                                            double *dm,                     \
                                            npy_intp mA, npy_intp mB,       \
                                            npy_intp n)                     \
    {                                                                       \
        Py_ssize_t i, j;                                                    \
        const type *u, *v;                                                  \
        for (i = 0; i < mA; i++) {                                          \
            for (j = 0; j < mB; j++, dm++) {                                \
                u = XA + n * i;                                             \
                v = XB + n * j;                                             \
                *dm = name ## _distance_ ## type(u, v, n);                  \
            }                                                               \
        }                                                                   \
    }

DEFINE_CDIST(bray_curtis, double)
DEFINE_CDIST(canberra, double)
DEFINE_CDIST(chebyshev, double)
DEFINE_CDIST(city_block, double)
DEFINE_CDIST(euclidean, double)
DEFINE_CDIST(hamming, double)
DEFINE_CDIST(jaccard, double)
DEFINE_CDIST(sqeuclidean, double)

DEFINE_CDIST(dice, char)
DEFINE_CDIST(hamming, char)
DEFINE_CDIST(jaccard, char)
DEFINE_CDIST(kulsinski, char)
DEFINE_CDIST(matching, char)
DEFINE_CDIST(rogerstanimoto, char)
DEFINE_CDIST(russellrao, char)
DEFINE_CDIST(sokalmichener, char)
DEFINE_CDIST(sokalsneath, char)
DEFINE_CDIST(yule_bool, char)


#define DEFINE_PDIST(name, type) \
    static void pdist_ ## name ## _ ## type(const type *X, double *dm,      \
                                            npy_intp m, npy_intp n)         \
    {                                                                       \
        Py_ssize_t i, j;                                                    \
        const type *u, *v;                                                  \
        double *it = dm;                                                    \
        for (i = 0; i < m; i++) {                                           \
            for (j = i + 1; j < m; j++, it++) {                             \
                u = X + n * i;                                              \
                v = X + n * j;                                              \
                *it = name ## _distance_ ## type(u, v, n);                  \
            }                                                               \
        }                                                                   \
    }

DEFINE_PDIST(bray_curtis, double)
DEFINE_PDIST(canberra, double)
DEFINE_PDIST(chebyshev, double)
DEFINE_PDIST(city_block, double)
DEFINE_PDIST(euclidean, double)
DEFINE_PDIST(hamming, double)
DEFINE_PDIST(jaccard, double)
DEFINE_PDIST(sqeuclidean, double)

DEFINE_PDIST(dice, char)
DEFINE_PDIST(hamming, char)
DEFINE_PDIST(jaccard, char)
DEFINE_PDIST(kulsinski, char)
DEFINE_PDIST(matching, char)
DEFINE_PDIST(rogerstanimoto, char)
DEFINE_PDIST(russellrao, char)
DEFINE_PDIST(sokalmichener, char)
DEFINE_PDIST(sokalsneath, char)
DEFINE_PDIST(yule_bool, char)
