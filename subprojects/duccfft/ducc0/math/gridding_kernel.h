/*
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* Copyright (C) 2020-2024 Max-Planck-Society
   Author: Martin Reinecke */

#ifndef DUCC0_GRIDDING_KERNEL_H
#define DUCC0_GRIDDING_KERNEL_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <vector>
#include <memory>
#include <cmath>
#include <type_traits>
#include <limits>
#include "ducc0/infra/useful_macros.h"
#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/simd.h"
#include "ducc0/infra/threading.h"
#include "ducc0/math/gl_integrator.h"
#include "ducc0/math/constants.h"

namespace ducc0 {

namespace detail_gridding_kernel {

using namespace std;

vector<double> getCoeffs(size_t W, size_t D, const function<double(double)> &func);

/*! A GriddingKernel is considered to be a symmetric real-valued function
    defined on the interval [-1; 1].
    This range is subdivided into W equal-sized parts. */
class GriddingKernel
  {
  public:
    virtual ~GriddingKernel() {}

    virtual size_t support() const = 0;

    /* Computes the correction function at a given coordinate.
       Useful coordinates lie in the range [0; 0.5]. */
    virtual double corfunc(double x) const = 0;

    /* Computes the correction function values at a coordinates
       [0, dx, 2*dx, ..., (n-1)*dx]  */
    virtual vector<double> corfunc(size_t n, double dx, int nthreads=1) const = 0;

    /* Returns the kernel's value at x */
    virtual double eval(double x) const = 0;
  };

class KernelCorrection
  {
  protected:
    using Tsimd = native_simd<double>;
    static constexpr size_t vlen=Tsimd::size();
    vector<double> x, wgtpsi;
    size_t supp;

  public:
    /* Compute correction factors for gridding kernel
       This implementation follows eqs. (3.8) to (3.10) of Barnett et al. 2018 */
    template<typename T> [[gnu::always_inline]] T corfunc(T v) const
      {
      T tmp=0;
      for (size_t i=0; i<x.size(); ++i)
        tmp += T(wgtpsi[i])*cos(T(x[i])*v);
      return T(1)/tmp;
      }
    /* Compute correction factors for gridding kernel
       This implementation follows eqs. (3.8) to (3.10) of Barnett et al. 2018 */
    vector<double> corfunc(size_t n, double dx, int nthreads=1) const
      {
      vector<double> res(n);
      Tsimd itimesdx;
      for (size_t i=0; i<vlen; ++i) itimesdx[i]=i*dx;
      execStatic(n, nthreads, 0, [&](auto &sched)
        {
        while (auto rng=sched.getNext())
          {
          auto i = rng.lo;
          for (; i+vlen<=rng.hi; i+=vlen)
            {
            auto v = corfunc(itimesdx+i*dx);
            v.copy_to(&res[i],element_aligned_tag());
            }
          for(; i<rng.hi; ++i)
            res[i] = corfunc(i*dx);
          }
        });
      return res;
      }
  };

class GLFullCorrection: public KernelCorrection
  {
  public:
    GLFullCorrection(size_t W, const function<double(double)> &func)
      {
      supp = W;
      // 1.5*W is safe according to FINUFFT. It may even be posible to reduce
      // to 1*W.
      size_t p = size_t(1.5*W)+2;
      GL_Integrator integ(2*p);
      x = integ.coordsSymmetric();
      wgtpsi = integ.weightsSymmetric();
      for (size_t i=0; i<x.size(); ++i)
        {
        wgtpsi[i] *= func(x[i])*supp*0.5;
        x[i] *= pi*supp;
        }
      }
  };

/*! This class implements the \a GriddingKernel interface by approximating the
    provided function with \a W polynomials of degree \a D. */
class PolynomialKernel: public GriddingKernel
  {
  private:
    size_t W, D;
    vector<double> coeff;
    KernelCorrection corr;

  public:
    PolynomialKernel(size_t W_, size_t D_, const function<double(double)> &func,
      const KernelCorrection &corr_)
      : W(W_), D(D_),
        coeff(getCoeffs(W_, D_, func)),
        corr(corr_)
      {}

    virtual size_t support() const { return W; }

    virtual double corfunc(double x) const { return corr.corfunc(x); }

    /* Computes the correction function values at a coordinates
       [0, dx, 2*dx, ..., (n-1)*dx]  */
    virtual vector<double> corfunc(size_t n, double dx, int nthreads=1) const
      { return corr.corfunc(n, dx, nthreads); }

    const vector<double> &Coeff() const { return coeff; }
    size_t degree() const { return D; }

    const KernelCorrection &Corr() const { return corr; }

    double eval(double x) const
      {
      if (abs(x)>=1) return 0.;
      double xrel = W*0.5*(x+1.);
      size_t nth = size_t(xrel);
      nth = min<size_t>(nth, W-1);
      double locx = ((xrel-nth)-0.5)*2; // should be in [-1; 1]
      double res = coeff[nth];
      for (size_t i=1; i<=D; ++i)
        res = res*locx+coeff[i*W+nth];
      return res;
      }
  };

/*! This class is initialized with a \a PolynomialKernel object and provides
    low level methods for extremely fast kernel evaluations. */
template<size_t W, typename Tsimd> class TemplateKernel
  {
  private:
    static constexpr auto D=W+3+(W&1);
    using T = typename Tsimd::value_type;
#ifdef DUCC0_HOMEGROWN_SIMD
    using Tvl = typename Tsimd::Tv;
#else
    using Tvl = Tsimd;
#endif
    static constexpr auto vlen = Tsimd::size();
    static constexpr auto nvec = (W+vlen-1)/vlen;
    static constexpr auto nvec_eval = (nvec+1)/2;

    std::array<Tsimd,(D+1)*nvec_eval> coeff;
    const T *scoeff;
    static constexpr auto sstride = nvec_eval*vlen;

    void transferCoeffs(const vector<double> &input, size_t d_input)
      {
      auto ofs = D-d_input;
      if (ofs>0)
        for (size_t i=0; i<nvec_eval; ++i)
          coeff[i] = 0;
      for (size_t j=0; j<=d_input; ++j)
        {
        for (size_t i=0; i<min(W, nvec_eval*vlen); ++i)
          coeff[(j+ofs)*nvec_eval + i/vlen][i%vlen] = T(input[j*W+i]);
        for (size_t i=W; i<nvec_eval*vlen; ++i)
          coeff[(j+ofs)*nvec_eval + i/vlen][i%vlen] = T(0);
        }
      }

  public:
    TemplateKernel(const PolynomialKernel &krn)
      : scoeff(reinterpret_cast<T *>(&coeff[0]))
      {
      MR_assert(W==krn.support(), "support mismatch");
      MR_assert(D>=krn.degree(), "degree mismatch");
      transferCoeffs(krn.Coeff(), krn.degree());
      }

    constexpr size_t support() const { return W; }

    double eval(double x) const
      {
      if (abs(x)>=1) return 0.;
      double xrel = W*0.5*(x+1.);
      size_t nth = size_t(xrel);
      nth = min<size_t>(nth, W-1);
      double locx = ((xrel-nth)-0.5)*2; // should be in [-1; 1]
      if (nth>=nvec_eval*vlen)
        {
        locx*=-1;
        nth=W-1-nth;
        }
      double res = scoeff[nth];
      for (size_t i=1; i<=D; ++i)
        res = res*locx+scoeff[i*sstride+nth];
      return res;
      }

    [[gnu::always_inline]] void eval2s(T x, T y, T z, size_t nth, Tsimd * DUCC0_RESTRICT res) const
      {
      z = (z-nth)*2+(W-1);
      T x2=x*x, y2=y*y, z2=z*z;
      if (nth>=nvec_eval*vlen)
        {
        z*=-1;
        nth=W-1-nth;
        }
      if constexpr((nvec>1)&&((W%vlen)!=0))
        res[nvec-1] = res[2*nvec-1] = 0;

      T zfac;
      {
      Tvl tvalx = coeff[0], tvaly = coeff[0];
      Tvl tvalx2 = coeff[nvec_eval], tvaly2 = coeff[nvec_eval];
      auto ptrz = scoeff+nth;
      auto tvalz = *ptrz, tvalz2 = ptrz[sstride];
      for (size_t j=2; j<D; j+=2)
        {
        tvalx = tvalx*x2 + Tvl(coeff[j*nvec_eval]);
        tvaly = tvaly*y2 + Tvl(coeff[j*nvec_eval]);
        tvalz = tvalz*z2 + ptrz[j*sstride];
        tvalx2 = tvalx2*x2 + Tvl(coeff[(j+1)*nvec_eval]);
        tvaly2 = tvaly2*y2 + Tvl(coeff[(j+1)*nvec_eval]);
        tvalz2 = tvalz2*z2 + ptrz[(j+1)*sstride];
        }
      zfac = tvalz*z+tvalz2;
      res[0] = (tvalx*x+tvalx2)*zfac;
      res[nvec] = tvaly*y+tvaly2;
      auto tmpx = Tsimd(tvalx2-tvalx*x)*zfac;
      auto tmpy = Tsimd(tvaly2-tvaly*y);
      for (size_t j=0, j2=W-1; (j<vlen)&&(j2>=nvec_eval*vlen); ++j,--j2)
        {
        res[j2/vlen][j2%vlen] = T(tmpx[j]);
        res[nvec+j2/vlen][j2%vlen] = T(tmpy[j]);
        }
      }
      for (size_t i=1; i<nvec_eval; ++i)
        {
        Tvl tvalx = coeff[i], tvaly = coeff[i];
        Tvl tvalx2 = coeff[i+nvec_eval], tvaly2 = coeff[i+nvec_eval];
        for (size_t j=2; j<D; j+=2)
          {
          tvalx = tvalx*x2 + Tvl(coeff[i+j*nvec_eval]);
          tvaly = tvaly*y2 + Tvl(coeff[i+j*nvec_eval]);
          tvalx2 = tvalx2*x2 + Tvl(coeff[i+(j+1)*nvec_eval]);
          tvaly2 = tvaly2*y2 + Tvl(coeff[i+(j+1)*nvec_eval]);
          }
        res[i] = (tvalx*x+tvalx2)*zfac;
        res[nvec+i] = tvaly*y+tvaly2;
        auto tmpx = Tsimd(tvalx2-tvalx*x)*zfac;
        auto tmpy = Tsimd(tvaly2-tvaly*y);
        for (size_t j=0, j2=W-1-i*vlen; (j<vlen)&&(j2>=nvec_eval*vlen); ++j,--j2)
          {
          res[j2/vlen][j2%vlen] = T(tmpx[j]);
          res[nvec+j2/vlen][j2%vlen] = T(tmpy[j]);
          }
        }
      }
    [[gnu::always_inline]] void eval1(T x, Tsimd * DUCC0_RESTRICT res) const
      {
      T x2=x*x;

      if constexpr((nvec>1)&&((W%vlen)!=0))
        res[nvec-1] = 0;
      for (size_t i=0; i<nvec_eval; ++i)
        {
        Tvl tvalx = coeff[i];
        Tvl tvalx2 = coeff[i+nvec_eval];
        for (size_t j=2; j<D; j+=2)
          {
          tvalx = tvalx*x2 + Tvl(coeff[i+j*nvec_eval]);
          tvalx2 = tvalx2*x2 + Tvl(coeff[i+(j+1)*nvec_eval]);
          }
        res[i] = tvalx*x+tvalx2;
        auto tmp = Tsimd(tvalx2-tvalx*x);
        for (size_t j=0, j2=W-1-i*vlen; (j<vlen)&&(j2>=nvec_eval*vlen); ++j,--j2)
          res[j2/vlen][j2%vlen] = T(tmp[j]);
        }
      }
    [[gnu::always_inline]] void eval2(T x, T y, Tsimd * DUCC0_RESTRICT res) const
      {
      T x2=x*x, y2=y*y;

      if constexpr((nvec>1)&&((W%vlen)!=0))
        res[nvec-1] = res[2*nvec-1] = 0;
      for (size_t i=0; i<nvec_eval; ++i)
        {
        Tvl tvalx = coeff[i], tvaly = coeff[i];
        Tvl tvalx2 = coeff[i+nvec_eval], tvaly2 = coeff[i+nvec_eval];
        for (size_t j=2; j<D; j+=2)
          {
          tvalx = tvalx*x2 + Tvl(coeff[i+j*nvec_eval]);
          tvaly = tvaly*y2 + Tvl(coeff[i+j*nvec_eval]);
          tvalx2 = tvalx2*x2 + Tvl(coeff[i+(j+1)*nvec_eval]);
          tvaly2 = tvaly2*y2 + Tvl(coeff[i+(j+1)*nvec_eval]);
          }
        res[i] = tvalx*x+tvalx2;
        res[nvec+i] = tvaly*y+tvaly2;
        auto tmpx = Tsimd(tvalx2-tvalx*x);
        auto tmpy = Tsimd(tvaly2-tvaly*y);
        for (size_t j=0, j2=W-1-i*vlen; (j<vlen)&&(j2>=nvec_eval*vlen); ++j,--j2)
          {
          res[j2/vlen][j2%vlen] = T(tmpx[j]);
          res[nvec+j2/vlen][j2%vlen] = T(tmpy[j]);
          }
        }
      }
    [[gnu::always_inline]] void eval3(T x, T y, T z, Tsimd * DUCC0_RESTRICT res) const
      {
      T x2=x*x, y2=y*y, z2=z*z;

      if constexpr((nvec>1)&&((W%vlen)!=0))
        res[nvec-1] = res[2*nvec-1]  = res[3*nvec-1] = 0;
      for (size_t i=0; i<nvec_eval; ++i)
        {
        Tvl tvalx = coeff[i], tvaly = coeff[i], tvalz = coeff[i];
        Tvl tvalx2 = coeff[i+nvec_eval], tvaly2 = coeff[i+nvec_eval], tvalz2 = coeff[i+nvec_eval];
        for (size_t j=2; j<D; j+=2)
          {
          tvalx = tvalx*x2 + Tvl(coeff[i+j*nvec_eval]);
          tvaly = tvaly*y2 + Tvl(coeff[i+j*nvec_eval]);
          tvalz = tvalz*z2 + Tvl(coeff[i+j*nvec_eval]);
          tvalx2 = tvalx2*x2 + Tvl(coeff[i+(j+1)*nvec_eval]);
          tvaly2 = tvaly2*y2 + Tvl(coeff[i+(j+1)*nvec_eval]);
          tvalz2 = tvalz2*z2 + Tvl(coeff[i+(j+1)*nvec_eval]);
          }
        res[i] = tvalx*x+tvalx2;
        res[nvec+i] = tvaly*y+tvaly2;
        res[2*nvec+i] = tvalz*z+tvalz2;
        auto tmpx = Tsimd(tvalx2-tvalx*x);
        auto tmpy = Tsimd(tvaly2-tvaly*y);
        auto tmpz = Tsimd(tvalz2-tvalz*z);
        for (size_t j=0, j2=W-1-i*vlen; (j<vlen)&&(j2>=nvec_eval*vlen); ++j,--j2)
          {
          res[j2/vlen][j2%vlen] = T(tmpx[j]);
          res[nvec+j2/vlen][j2%vlen] = T(tmpy[j]);
          res[2*nvec+j2/vlen][j2%vlen] = T(tmpz[j]);
          }
        }
      }
  };

struct KernelParams
  {
  size_t W;
  double ofactor, epsilon, beta, e0;
  size_t ndim;
  bool singleprec;
  };

shared_ptr<PolynomialKernel> selectKernel(size_t idx);
const KernelParams &getKernel(size_t idx);

template<typename T> constexpr inline size_t Wmax()
  { return is_same<T,float>::value ? 8 : 16; }

extern const vector<KernelParams> KernelDB;

/*! Returns the 2-parameter ES kernel for the given oversampling factor,
 *  dimensionality, and error that has the smallest support. */
template<typename T> auto selectKernel(double ofactor, size_t ndim, double epsilon)
  {
  constexpr bool singleprec = is_same<T, float>::value;
  size_t Wmin = Wmax<T>();
  size_t idx = KernelDB.size();
  for (size_t i=0; i<KernelDB.size(); ++i)
    {
    if  ((KernelDB[i].ndim==ndim) && (KernelDB[i].singleprec==singleprec)
      && (KernelDB[i].ofactor<=ofactor) && (KernelDB[i].epsilon<=epsilon)
      && (KernelDB[i].W<=Wmin))
      {
      idx = i;
      Wmin = KernelDB[i].W;
      }
    }
  return selectKernel(idx);
  }

template<typename T> auto getAvailableKernels(double epsilon,
  size_t ndim, double ofactor_min=1.1, double ofactor_max=2.6)
  {
  constexpr bool singleprec = is_same<T, float>::value;
  vector<double> ofc(20, ofactor_max);
  vector<size_t> idx(20, KernelDB.size());
  size_t Wlim = Wmax<T>();
  for (size_t i=0; i<KernelDB.size(); ++i)
    {
    auto ofactor = KernelDB[i].ofactor;
    size_t W = KernelDB[i].W;
    if ((KernelDB[i].ndim==ndim) && (KernelDB[i].singleprec==singleprec)
      && (W<=Wlim) && (KernelDB[i].epsilon<=epsilon)
      && (ofactor<=ofc[W]) && (ofactor>=ofactor_min))
      {
      ofc[W] = ofactor;
      idx[W] = i;
      }
    }
  vector<size_t> res;
  for (auto v: idx)
    if (v<KernelDB.size()) res.push_back(v);
  MR_assert(!res.empty(),
    "No appropriate kernel found for the specified combination of parameters\n"
    "(epsilon, sigma_min, sigma_max, ndim, floating point precision).");
  return res;
  }

double bestEpsilon(size_t ndim, bool singleprec,
  double ofactor_min=1.1, double ofactor_max=2.6);

}

using detail_gridding_kernel::GriddingKernel;
using detail_gridding_kernel::getKernel;
using detail_gridding_kernel::selectKernel;
using detail_gridding_kernel::getAvailableKernels;
using detail_gridding_kernel::bestEpsilon;
using detail_gridding_kernel::PolynomialKernel;
using detail_gridding_kernel::TemplateKernel;
using detail_gridding_kernel::KernelParams;

}

#endif
