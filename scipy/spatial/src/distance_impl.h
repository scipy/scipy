/*
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

#include <math.h>
#include <numpy/npy_math.h>

#if !defined(__clang__) && defined(__GNUC__) && defined(__GNUC_MINOR__)
#if __GNUC__ >= 5 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
/* enable auto-vectorizer */
#pragma GCC optimize("tree-vectorize")
/* float associativity required to vectorize reductions */
#pragma GCC optimize("unsafe-math-optimizations")
/* maybe 5% gain, manual unrolling with more accumulators would be better */
#pragma GCC optimize("unroll-loops")
#endif
#endif

static NPY_INLINE double
sqeuclidean_distance_double(const double *u, const double *v, npy_intp n)
{
    double s = 0.0, d;
    npy_intp i;

    for (i = 0; i < n; i++) {
        d = u[i] - v[i];
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
    npy_intp i;

    for (i = 0; i < n; i++) {
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
yule_distance_char(const char *u, const char *v, npy_intp n)
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
    npy_intp i;

    for (i = 0; i < n; i++) {
        s += u[i] * v[i];
    }
    return s;
}

static NPY_INLINE double
seuclidean_distance(const double *var, const double *u, const double *v,
                    npy_intp n)
{
    double s = 0.0, d;
    npy_intp i;

    for (i = 0; i < n; i++) {
        d = u[i] - v[i];
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
pdist_cosine(const double *X, double *dm, npy_intp m, npy_intp n,
             const double *norms)
{
    const double *u, *v;
    double cosine;
    npy_intp i, j;

    for (i = 0; i < m; i++) {
        for (j = i + 1; j < m; j++, dm++) {
            u = X + (n * i);
            v = X + (n * j);
            cosine = dot_product(u, v, n) / (norms[i] * norms[j]);
            if (fabs(cosine) > 1.) {
                /* Clip to correct rounding error. */
                cosine = npy_copysign(1, cosine);
            }
            *dm = 1. - cosine;
        }
    }
}

static NPY_INLINE void
dist_to_squareform_from_vector(double *M, const double *v, npy_intp n)
{
    double *it, *it2;
    npy_intp i, j;

    for (i = 0; i < n - 1; i++) {
        it  = M + (i * n) + i + 1;
        it2 = M + (i+1)*n + i;
        for (j = i + 1; j < n; j++, it++, it2 += n, v++) {
            *it = *v;
            *it2 = *v;
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
