/**
 * distance.c
 *
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
#include <Python.h>
#include <numpy/ndarrayobject.h>

#include <math.h>
#include <stdlib.h>
#include "distance.h"

static NPY_INLINE double sqeuclidean_distance(const double *u, const double *v, int n) {
  int i = 0;
  double s = 0.0, d;
  for (i = 0; i < n; i++) {
    d = u[i] - v[i];
    s = s + d * d;
  }
  return s;
}

static NPY_INLINE double euclidean_distance(const double *u, const double *v, int n) {
    return sqrt(sqeuclidean_distance(u, v, n));
}

#if 0   /* XXX unused */
static NPY_INLINE double ess_distance(const double *u, const double *v, int n) {
  int i = 0;
  double s = 0.0, d;
  for (i = 0; i < n; i++) {
    d = fabs(u[i] - v[i]);
    s = s + d * d;
  }
  return s;
}
#endif

static NPY_INLINE double chebyshev_distance(const double *u, const double *v, int n) {
  int i = 0;
  double d, maxv = 0.0;
  for (i = 0; i < n; i++) {
    d = fabs(u[i] - v[i]);
    if (d > maxv) {
      maxv = d;
    }
  }
  return maxv;
}

static NPY_INLINE double canberra_distance(const double *u, const double *v, int n) {
  int i;
  double snum = 0.0, sdenom = 0.0, tot = 0.0;
  for (i = 0; i < n; i++) {
    snum = fabs(u[i] - v[i]);
    sdenom = fabs(u[i]) + fabs(v[i]);
    if (sdenom > 0.0) {
        tot += snum / sdenom;
    }
  }
  return tot;
}

static NPY_INLINE double bray_curtis_distance(const double *u, const double *v, int n) {
  int i;
  double s1 = 0.0, s2 = 0.0;
  for (i = 0; i < n; i++) {
    s1 += fabs(u[i] - v[i]);
    s2 += fabs(u[i] + v[i]);
  }
  return s1 / s2;
}

static NPY_INLINE double mahalanobis_distance(const double *u, const double *v,
			    const double *covinv, double *dimbuf1,
			    double *dimbuf2, int n) {
  int i, j;
  double s;
  const double *covrow = covinv;
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

double hamming_distance(const double *u, const double *v, int n) {
  int i = 0;
  double s = 0.0;
  for (i = 0; i < n; i++) {
    s = s + (u[i] != v[i]);
  }
  return s / (double)n;
}

static NPY_INLINE double hamming_distance_bool(const char *u, const char *v, int n) {
  int i = 0;
  double s = 0.0;
  for (i = 0; i < n; i++) {
    s = s + (u[i] != v[i]);
  }
  return s / (double)n;
}

static NPY_INLINE double yule_distance_bool(const char *u, const char *v, int n) {
  int i = 0;
  int ntt = 0, nff = 0, nft = 0, ntf = 0;
  for (i = 0; i < n; i++) {
    ntt += (u[i] && v[i]);
    ntf += (u[i] && !v[i]);
    nft += (!u[i] && v[i]);
    nff += (!u[i] && !v[i]);
  }
  return (2.0 * ntf * nft) / (double)(ntt * nff + ntf * nft);  
}

static NPY_INLINE double matching_distance_bool(const char *u, const char *v, int n) {
  int i = 0;
  int nft = 0, ntf = 0;
  for (i = 0; i < n; i++) {
    ntf += (u[i] && !v[i]);
    nft += (!u[i] && v[i]);
  }
  return (double)(ntf + nft) / (double)(n);
}

static NPY_INLINE double dice_distance_bool(const char *u, const char *v, int n) {
  int i = 0;
  int ntt = 0, nft = 0, ntf = 0;
  for (i = 0; i < n; i++) {
    ntt += (u[i] && v[i]);
    ntf += (u[i] && !v[i]);
    nft += (!u[i] && v[i]);
  }
  return (double)(nft + ntf) / (double)(2.0 * ntt + ntf + nft);
}


static NPY_INLINE double rogerstanimoto_distance_bool(const char *u, const char *v, int n) {
  int i = 0;
  int ntt = 0, nff = 0, nft = 0, ntf = 0;
  for (i = 0; i < n; i++) {
    ntt += (u[i] && v[i]);
    ntf += (u[i] && !v[i]);
    nft += (!u[i] && v[i]);
    nff += (!u[i] && !v[i]);
  }
  return (2.0 * (ntf + nft)) / ((double)ntt + nff + (2.0 * (ntf + nft)));
}

static NPY_INLINE double russellrao_distance_bool(const char *u, const char *v, int n) {
  int i = 0;
  /**  int nff = 0, nft = 0, ntf = 0;**/
  int ntt = 0;
  for (i = 0; i < n; i++) {
    /**    nff += (!u[i] && !v[i]);
    ntf += (u[i] && !v[i]);
    nft += (!u[i] && v[i]);**/
    ntt += (u[i] && v[i]);
  }
  /**  return (double)(ntf + nft + nff) / (double)n;**/
  return (double) (n - ntt) / (double) n;
}

static NPY_INLINE double kulsinski_distance_bool(const char *u, const char *v, int n) {
  int _i = 0;
  int ntt = 0, nft = 0, ntf = 0, nff = 0;
  for (_i = 0; _i < n; _i++) {
    ntt += (u[_i] && v[_i]);
    ntf += (u[_i] && !v[_i]);
    nft += (!u[_i] && v[_i]);
    nff += (!u[_i] && !v[_i]);
  }
  return ((double)(ntf + nft - ntt + n)) / ((double)(ntf + nft + n));
}

static NPY_INLINE double sokalsneath_distance_bool(const char *u, const char *v, int n) {
  int _i = 0;
  int ntt = 0, nft = 0, ntf = 0;
  for (_i = 0; _i < n; _i++) {
    ntt += (u[_i] && v[_i]);
    ntf += (u[_i] && !v[_i]);
    nft += (!u[_i] && v[_i]);
  }
  return (2.0 * (ntf + nft))/(2.0 * (ntf + nft) + ntt);
}

static NPY_INLINE double sokalmichener_distance_bool(const char *u, const char *v, int n) {
  int _i = 0;
  int ntt = 0, nft = 0, ntf = 0, nff = 0;
  for (_i = 0; _i < n; _i++) {
    ntt += (u[_i] && v[_i]);
    nff += (!u[_i] && !v[_i]);
    ntf += (u[_i] && !v[_i]);
    nft += (!u[_i] && v[_i]);
  }
  return (2.0 * (ntf + nft))/(2.0 * (ntf + nft) + ntt + nff);
}

static NPY_INLINE double jaccard_distance(const double *u, const double *v, int n) {
  int i = 0;
  double denom = 0.0, num = 0.0;
  for (i = 0; i < n; i++) {
    num += (u[i] != v[i]) && ((u[i] != 0.0) || (v[i] != 0.0));
    denom += (u[i] != 0.0) || (v[i] != 0.0);
  }
  return num / denom;
}

static NPY_INLINE double jaccard_distance_bool(const char *u, const char *v, int n) {
  int i = 0;
  double num = 0.0, denom = 0.0;
  for (i = 0; i < n; i++) {
    num += (u[i] != v[i]) && ((u[i] != 0) || (v[i] != 0));
    denom += (u[i] != 0) || (v[i] != 0);
  }
  return num / denom;
}

static NPY_INLINE double dot_product(const double *u, const double *v, int n) {
  int i;
  double s = 0.0;
  for (i = 0; i < n; i++) {
    s += u[i] * v[i];
  }
  return s;
}

static NPY_INLINE double cosine_distance(const double *u, const double *v, int n,
		       const double nu, const double nv) {
  return 1.0 - (dot_product(u, v, n) / (nu * nv));
}

static NPY_INLINE double seuclidean_distance(const double *var,
			   const double *u, const double *v, int n) {
  int i = 0;
  double s = 0.0, d;
  for (i = 0; i < n; i++) {
    d = u[i] - v[i];
    s = s + (d * d) / var[i];
  }
  return sqrt(s);
}

static NPY_INLINE double city_block_distance(const double *u, const double *v, int n) {
  int i = 0;
  double s = 0.0, d;
  for (i = 0; i < n; i++) {
    d = fabs(u[i] - v[i]);
    s = s + d;
  }
  return s;
}

double minkowski_distance(const double *u, const double *v, int n, double p) {
  int i = 0;
  double s = 0.0, d;
  for (i = 0; i < n; i++) {
    d = fabs(u[i] - v[i]);
    s = s + pow(d, p);
  }
  return pow(s, 1.0 / p);
}

double weighted_minkowski_distance(const double *u, const double *v, int n, double p, const double *w) {
  int i = 0;
  double s = 0.0, d;
  for (i = 0; i < n; i++) {
    d = fabs(u[i] - v[i]) * w[i];
    s = s + pow(d, p);
  }
  return pow(s, 1.0 / p);
}

void compute_mean_vector(double *res, const double *X, int m, int n) {
  int i, j;
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

void pdist_euclidean(const double *X, double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = euclidean_distance(u, v, n);
    }
  }
}

void pdist_sqeuclidean(const double *X, double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = sqeuclidean_distance(u, v, n);
    }
  }
}

void pdist_mahalanobis(const double *X, const double *covinv,
		       double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  double *dimbuf1, *dimbuf2;
  dimbuf1 = (double*)malloc(sizeof(double) * 2 * n);
  dimbuf2 = dimbuf1 + n;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = mahalanobis_distance(u, v, covinv, dimbuf1, dimbuf2, n);
    }
  }
  dimbuf2 = 0;
  free(dimbuf1);
}

void pdist_bray_curtis(const double *X, double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = bray_curtis_distance(u, v, n);
    }
  }
}

void pdist_canberra(const double *X, double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = canberra_distance(u, v, n);
    }
  }
}

void pdist_hamming(const double *X, double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = hamming_distance(u, v, n);
    }
  }
}

void pdist_hamming_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = hamming_distance_bool(u, v, n);
    }
  }
}

void pdist_jaccard(const double *X, double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = jaccard_distance(u, v, n);
    }
  }
}

void pdist_jaccard_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = jaccard_distance_bool(u, v, n);
    }
  }
}


void pdist_chebyshev(const double *X, double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = chebyshev_distance(u, v, n);
    }
  }
}

void pdist_cosine(const double *X, double *dm, int m, int n, const double *norms) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = cosine_distance(u, v, n, norms[i], norms[j]);
    }
  }
}

void pdist_seuclidean(const double *X, const double *var,
		     double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = seuclidean_distance(var, u, v, n);
    }
  }
}

void pdist_city_block(const double *X, double *dm, int m, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = city_block_distance(u, v, n);
    }
  }
}

void pdist_minkowski(const double *X, double *dm, int m, int n, double p) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = minkowski_distance(u, v, n, p);
    }
  }
}

void pdist_weighted_minkowski(const double *X, double *dm, int m, int n, double p, const double *w) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = weighted_minkowski_distance(u, v, n, p, w);
    }
  }
}

void pdist_yule_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = yule_distance_bool(u, v, n);
    }
  }
}

void pdist_matching_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = matching_distance_bool(u, v, n);
    }
  }
}

void pdist_dice_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = dice_distance_bool(u, v, n);
    }
  }
}

void pdist_rogerstanimoto_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = rogerstanimoto_distance_bool(u, v, n);
    }
  }
}

void pdist_russellrao_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = russellrao_distance_bool(u, v, n);
    }
  }
}

void pdist_kulsinski_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = kulsinski_distance_bool(u, v, n);
    }
  }
}

void pdist_sokalsneath_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = sokalsneath_distance_bool(u, v, n);
    }
  }
}

void pdist_sokalmichener_bool(const char *X, double *dm, int m, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++, it++) {
      u = X + (n * i);
      v = X + (n * j);
      *it = sokalmichener_distance_bool(u, v, n);
    }
  }
}

void dist_to_squareform_from_vector(double *M, const double *v, int n) {
  double *it;
  const double *cit;
  int i, j;
  cit = v;
  for (i = 0; i < n - 1; i++) {
    it = M + (i * n) + i + 1;
    for (j = i + 1; j < n; j++, it++, cit++) {
      *it = *cit;
    }
  }
}

void dist_to_vector_from_squareform(const double *M, double *v, int n) {
  double *it;
  const double *cit;
  int i, j;
  it = v;
  for (i = 0; i < n - 1; i++) {
    cit = M + (i * n) + i + 1;
    for (j = i + 1; j < n; j++, it++, cit++) {
      *it = *cit;
    }
  }
}


/** cdist */

void cdist_euclidean(const double *XA,
		     const double *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = euclidean_distance(u, v, n);
    }
  }
}

void cdist_sqeuclidean(const double *XA,
		     const double *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = sqeuclidean_distance(u, v, n);
    }
  }
}

void cdist_mahalanobis(const double *XA,
		       const double *XB,
		       const double *covinv,
		       double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  double *dimbuf1, *dimbuf2;
  dimbuf1 = (double*)malloc(sizeof(double) * 2 * n);
  dimbuf2 = dimbuf1 + n;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = mahalanobis_distance(u, v, covinv, dimbuf1, dimbuf2, n);
    }
  }
  dimbuf2 = 0;
  free(dimbuf1);
}

void cdist_bray_curtis(const double *XA, const double *XB,
		       double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = bray_curtis_distance(u, v, n);
    }
  }
}

void cdist_canberra(const double *XA,
		    const double *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = canberra_distance(u, v, n);
    }
  }
}

void cdist_hamming(const double *XA,
		   const double *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = hamming_distance(u, v, n);
    }
  }
}

void cdist_hamming_bool(const char *XA,
			const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = hamming_distance_bool(u, v, n);
    }
  }
}

void cdist_jaccard(const double *XA,
		   const double *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = jaccard_distance(u, v, n);
    }
  }
}

void cdist_jaccard_bool(const char *XA,
			const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = jaccard_distance_bool(u, v, n);
    }
  }
}


void cdist_chebyshev(const double *XA,
		     const double *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = chebyshev_distance(u, v, n);
    }
  }
}

void cdist_cosine(const double *XA,
		  const double *XB, double *dm, int mA, int mB, int n,
		  const double *normsA, const double *normsB) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = cosine_distance(u, v, n, normsA[i], normsB[j]);
    }
  }
}

void cdist_seuclidean(const double *XA,
		      const double *XB,
		      const double *var,
		      double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = seuclidean_distance(var, u, v, n);
    }
  }
}

void cdist_city_block(const double *XA, const double *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = city_block_distance(u, v, n);
    }
  }
}

void cdist_minkowski(const double *XA, const double *XB, double *dm, int mA, int mB, int n, double p) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = minkowski_distance(u, v, n, p);
    }
  }
}

void cdist_weighted_minkowski(const double *XA, const double *XB, double *dm, int mA, int mB, int n, double p, const double *w) {
  int i, j;
  const double *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = weighted_minkowski_distance(u, v, n, p, w);
    }
  }
}

void cdist_yule_bool(const char *XA, const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = yule_distance_bool(u, v, n);
    }
  }
}

void cdist_matching_bool(const char *XA, const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = matching_distance_bool(u, v, n);
    }
  }
}

void cdist_dice_bool(const char *XA, const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = dice_distance_bool(u, v, n);
    }
  }
}

void cdist_rogerstanimoto_bool(const char *XA, const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = rogerstanimoto_distance_bool(u, v, n);
    }
  }
}

void cdist_russellrao_bool(const char *XA, const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = russellrao_distance_bool(u, v, n);
    }
  }
}

void cdist_kulsinski_bool(const char *XA, const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = kulsinski_distance_bool(u, v, n);
    }
  }
}

void cdist_sokalsneath_bool(const char *XA, const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = sokalsneath_distance_bool(u, v, n);
    }
  }
}

void cdist_sokalmichener_bool(const char *XA, const char *XB, double *dm, int mA, int mB, int n) {
  int i, j;
  const char *u, *v;
  double *it = dm;
  for (i = 0; i < mA; i++) {
    for (j = 0; j < mB; j++, it++) {
      u = XA + (n * i);
      v = XB + (n * j);
      *it = sokalmichener_distance_bool(u, v, n);
    }
  }
}
