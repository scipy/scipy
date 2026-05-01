/* Copyright (C) 2026 Max-Planck-Society
   Author: Martin Reinecke */

/* SPDX-License-Identifier: BSD-3-Clause OR GPL-2.0-or-later */

/*
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.
* Neither the name of the copyright holder nor the names of its contributors may
  be used to endorse or promote products derived from this software without
  specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <array>
#include <cmath>
#include <vector>
#include "ducc0/infra/mav.h"
#include "ducc0/infra/simd.h"

namespace ducc0 {

constexpr int PSWF_ERROR = 42;

/* Class for evaluation of the prolate spheroidal wavefunction
   of order zero (Psi_0^c) inside [-1,1], for arbitrary frequency parameter c.
   Computation is done using a basis of Legendre polynomials.
   This implementation is based on work by Libin Lu for FINUFFT.
   The orignal implementation was done by Vladimir Rokhlin and 
   can be found in src/common/specialfunctions/ of the DMK repo
   https://github.com/flatironinstitute/dmk */
class PSWF0 {
private:
  double c;
  std::vector<double> workdata; // Legendre coefficients
  std::vector<double> f1;   // f1 = (j-1.)/j
  double xv0;

  static void prolcoef(double rlam, int k, double c, double &alpha, double &beta,
                       double &gamma) {
    double kf     = k;
    double alpha0 = kf * (kf - 1.) / ((2. * kf + 1.) * (2. * kf - 1.));
    double beta0  = ((kf + 1.) * (kf + 1.) / (2. * kf + 3.) + kf * kf / (2. * kf - 1.)) /
                   (2. * kf + 1.);
    double gamma0 = (kf + 1.) * (kf + 2.) / ((2. * kf + 1.) * (2. * kf + 3.));

    alpha = -c * c * alpha0;
    beta  = rlam - kf * (kf + 1.) - c * c * beta0;
    gamma = -c * c * gamma0;
  }

  // fills as, bs, cs in [0; (n+2)/2]
  static void prolmatr(std::vector<double> &as, std::vector<double> &bs,
                       std::vector<double> &cs, int n, double c, double rlam) {
    for (int k = 0; 2*k <= n + 2; ++k) {
      prolcoef(rlam, 2*k, c, as[k], bs[k], cs[k]);

      if (k != 0)
        as[k] *= std::sqrt((2*k + .5) / (2*k - 1.5));
      cs[k] *= std::sqrt((2*k + .5) / (2*k + 2.5));
    }
  }

  static void prolql1(int n, std::vector<double> &d, std::vector<double> &e) {
    if (n == 1) return;

    for (int i = 1; i < n; ++i)
      e[i - 1] = e[i];
    e[n - 1] = 0.0;

    for (int l = 0; l < n; ++l) {
      int j = 0;
      while (true) {
        int m;
        for (m = l; m < n - 1; ++m) {
          double tst1 = std::abs(d[m]) + std::abs(d[m + 1]);
          double tst2 = tst1 + std::abs(e[m]);
          if (tst2 == tst1) break;
        }

        if (m == l) break;
        if (j == 30) throw int(PSWF_ERROR);
        ++j;

        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        double r = std::sqrt(g * g + 1.0);
        g        = d[m] - d[l] + e[l] / (g + std::copysign(r, g));
        double s = 1.0;
        double c = 1.0;
        double p = 0.0;

        for (int i = m - 1; i >= l; --i) {
          double f = s * e[i];
          double b = c * e[i];
          r        = std::sqrt(f * f + g * g);
          e[i + 1] = r;
          if (r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;
            break;
          }
          s        = f / r;
          c        = g / r;
          g        = d[i + 1] - p;
          r        = (d[i] - g) * s + 2.0 * c * b;
          p        = s * r;
          d[i + 1] = g + p;
          g        = c * r - b;
        }

        if (r == 0.0) break;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }

      if (l == 0) continue;
      for (int i = l; i > 0; --i) {
        if (d[i] >= d[i - 1]) break;
        std::swap(d[i], d[i - 1]);
      }
    }
  }

  static void prolfact(std::vector<double> &a, const std::vector<double> &b,
                       const std::vector<double> &c, int n, std::vector<double> &u,
                       std::vector<double> &v, std::vector<double> &w) {
    // Eliminate down and up, and scale
    for (int i = 0; i + 1 < n; ++i) {
      double d = c[i + 1] / a[i];
      a[i + 1] -= b[i] * d;
      u[i] = d;
      v[i+1] = b[i] / a[i+1];
      w[i+1] = 1. / a[i+1];
    }
    w[0] = 1./a[0];
  }

  static void prolsolv(const std::vector<double> &u, const std::vector<double> &v,
                       const std::vector<double> &w, int n, std::vector<double> &rhs) {
    // Eliminate down
    for (int i = 0; i + 1 < n; ++i) rhs[i + 1] -= u[i] * rhs[i];

    // Eliminate up and scale
    for (int i = n - 1; i > 0; --i) {
      rhs[i - 1] -= rhs[i] * v[i];
      rhs[i] *= w[i];
    }
    rhs[0] *= w[0];
  }

  static void prolfun0(int n, double c, std::vector<double> &xk, double eps) {
    double delta = 1.0e-8;

    xk.resize(n/2 + 3);
    std::vector<double> as(n/2 + 2), bs(n/2 + 2), cs(n/2 + 2), u(n/2 + 2), v(n/2 + 2), w(n/2 + 2);
    prolmatr(as, bs, cs, n, c, 0.);

    prolql1(n / 2, bs, as);

    std::fill(xk.begin(), xk.end(), 1.0);

    double rlam = -bs[n / 2 - 1] + delta;
    prolmatr(as, bs, cs, n, c, rlam);

    prolfact(bs, cs, as, n / 2, u, v, w);

    constexpr int numit = 4;
    for (int ijk = 0; ijk < numit; ++ijk) {
      prolsolv(u, v, w, n / 2, xk);

      double d = 0;
      for (int j = 0; j < n / 2; ++j) d += xk[j] * xk[j];

      d = std::sqrt(d);
      for (int j = 0; j < n / 2; ++j) {
        xk[j] /= d;
        as[j] = xk[j];
      }
    }

    int imax=0;
    for (int i = 0; i < n / 2; ++i) {
      if (std::abs(xk[i]) > eps) imax = i;
      xk[i] *= std::sqrt(i * 2 + .5);
    }
    xk.resize(imax + 1);
  }

  static void prolps0i(double c, std::vector<double> &work) {
    static const std::array<int, 20> ns = {48,  64,  80,  92,  106, 120, 130,
                                           144, 156, 168, 178, 190, 202, 214,
                                           224, 236, 248, 258, 268, 280};

    int i = static_cast<int>(c / 10);
    int n = (i < int(ns.size())) ? ns[i] : static_cast<int>(c * 3) / 2;

    prolfun0(n, c, work, 1e-16);
  }

  template<typename T, size_t N> array<T,N> eval_raw(array<T, N> x) const {
    array<T,N> pjm1, pjm2, val;
    for (size_t n=0; n<N; ++n) {
      pjm1[n] = 0;
      pjm2[n] = 1;
      val[n] = workdata[0];
    }

    for (size_t i=1; i < workdata.size(); ++i) {
      for (size_t n=0; n<N; ++n) {
        pjm1[n] = (f1[2*i-2]+1.) * x[n] * pjm2[n] - f1[2*i-2] * pjm1[n];
        pjm2[n] = (f1[2*i-1]+1.) * x[n] * pjm1[n] - f1[2*i-1] * pjm2[n];
        val[n] += workdata[i] * pjm2[n];
      }
    }
    return val;
  }
  template<typename T> T eval_raw(T x) const {
    T pjm1 = 0;
    T pjm2 = 1;
    T val = workdata[0];

    for (size_t i=1; i < workdata.size(); ++i) {
      pjm1 = (f1[2*i-2]+1.) * x * pjm2 - f1[2*i-2] * pjm1;
      pjm2 = (f1[2*i-1]+1.) * x * pjm1 - f1[2*i-1] * pjm2;
      val += workdata[i] * pjm2;
    }
    return val;
  }

public:
  PSWF0(double c_) : c(c_) {
    prolps0i(c, workdata);
    f1.resize(2 * workdata.size() - 2);
    for (size_t i = 0; i < f1.size(); ++i)
      f1[i] = i / (i+1.);
    xv0 = 1. / eval_raw(0.);
  }

  double operator()(double x) const {
    if (std::abs(x) > 1) return 0.;
    return eval_raw(x) * xv0;
  }

  void multi_eval(const cmav<double,1> &x, const vmav<double,1> &res) const {
    MR_assert(x.size()==res.size(), "array size mismatch");

    using Tv = native_simd<double>;
    constexpr size_t vlen = Tv::size();
    constexpr size_t nvec = 4;
    size_t i = 0;
    for (; i + nvec*vlen <= x.size(); i += nvec*vlen) {
      array<Tv,nvec> xx;
      for (size_t n=0; n<nvec; ++n)
        for (size_t m = 0; m < vlen; ++m)
          xx[n][m] = x(i + n*vlen + m);
      const auto val = eval_raw(xx);
      for (size_t n=0; n<nvec; ++n)
        for (size_t m = 0; m < vlen; ++m)
          res(i + n*vlen + m) = (std::abs(xx[n][m])>1) ? 0 : val[n][m] * xv0;
    }
    for (; i + vlen <= x.size(); i += vlen) {
      Tv xx;
      for (size_t m = 0; m < vlen; ++m) xx[m] = x(i + m);
      const auto val = eval_raw(xx);
      for (size_t m = 0; m < vlen; ++m)
        res(i + m) = (std::abs(xx[m])>1) ? 0 : val[m] * xv0;
    }

    for (; i < x.size(); ++i) res(i) = operator()(x(i));
  }

};

}
