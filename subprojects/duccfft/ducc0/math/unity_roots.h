/* Copyright (C) 2019-2026 Max-Planck-Society
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

#ifndef DUCC0_UNITY_ROOTS_H
#define DUCC0_UNITY_ROOTS_H

#include <cmath>
#include <memory>
#include <cstddef>
#include <type_traits>
#include <vector>
#include <tuple>

namespace ducc0 {

namespace detail_unity_roots {

using namespace std;

template<typename T, typename Tc> class RootCalc
  {
  private:
    static constexpr const auto pi = 3.141592653589793238462643383279502884197L;
    using Thigh = typename conditional<(sizeof(T)>sizeof(double)), T, double>::type;
    size_t n;
    Thigh ang;

  public:
    RootCalc(size_t n_) : n(n_), ang(Thigh(0.25L*pi/n)) {}
#if 1
    Tc calc(size_t x) const
      {
      x<<=3;
      if (x<4*n) // first half
        {
        if (x<2*n) // first quadrant
          {
          if (x<n) return {cos(Thigh(x)*ang), sin(Thigh(x)*ang)};
          return {sin(Thigh(2*n-x)*ang), cos(Thigh(2*n-x)*ang)};
          }
        else // second quadrant
          {
          x-=2*n;
          if (x<n) return {-sin(Thigh(x)*ang), cos(Thigh(x)*ang)};
          return {-cos(Thigh(2*n-x)*ang), sin(Thigh(2*n-x)*ang)};
          }
        }
      else
        {
        x=8*n-x;
        if (x<2*n) // third quadrant
          {
          if (x<n) return {cos(Thigh(x)*ang), -sin(Thigh(x)*ang)};
          return {sin(Thigh(2*n-x)*ang), -cos(Thigh(2*n-x)*ang)};
          }
        else // fourth quadrant
          {
          x-=2*n;
          if (x<n) return {-sin(Thigh(x)*ang), -cos(Thigh(x)*ang)};
          return {-cos(Thigh(2*n-x)*ang), -sin(Thigh(2*n-x)*ang)};
          }
        }
      }
#else
    // alternative version, similar speed, but maybe a bit more accurate
    Tc calc(size_t x) const
      {
      static constexpr Thigh pi = Thigh(3.141592653589793238462643383279502884197L);
      Thigh n4 = Thigh(n<<2);

      x<<=3;
      if (x<4*n) // first half
        {
        if (x<2*n) // first quadrant
          {
          if (x<n)
            {
            auto ang = (x/n4)*pi;
            return {cos(ang), sin(ang)};
            }
          auto ang = ((2*n-x)/n4)*pi;
          return {sin(ang), cos(ang)};
          }
        else // second quadrant
          {
          x-=2*n;
          if (x<n)
            {
            auto ang = (x/n4)*pi;
            return {-sin(ang), cos(ang)};
            }
          auto ang = ((2*n-x)/n4)*pi;
          return {-cos(ang), sin(ang)};
          }
        }
      else
        {
        x=8*n-x;
        if (x<2*n) // third quadrant
          {
          if (x<n)
            {
            auto ang = (x/n4)*pi;
            return {cos(ang), -sin(ang)};
            }
          auto ang = ((2*n-x)/n4)*pi;
          return {sin(ang), -cos(ang)};
          }
        else // fourth quadrant
          {
          x-=2*n;
          if (x<n)
            {
            auto ang = (x/n4)*pi;
            return {-sin(ang), -cos(ang)};
            }
          auto ang = ((2*n-x)/n4)*pi;
          return {-cos(ang), -sin(ang)};
          }
        }
      }
#endif
  };

template<typename T, typename Tc> class UnityRoots
  {
  private:
    using Thigh = typename conditional<(sizeof(T)>sizeof(double)), T, double>::type;
    struct cmplx_ { Thigh r, i; };
    size_t N, mask, shift;
    vector<cmplx_> v1, v2;

  public:
    UnityRoots(size_t n)
      : N(n)
      {
      RootCalc<T, cmplx_> rc(n);
      size_t nval = (n+2)/2;
      shift = 1;
      while((size_t(1)<<shift)*(size_t(1)<<shift) < nval) ++shift;
      mask = (size_t(1)<<shift)-1;
      v1.resize(mask+1);
      v1[0]={Thigh(1), Thigh(0)};
      for (size_t i=1; i<v1.size(); ++i)
        v1[i]=rc.calc(i);
      v2.resize((nval+mask)/(mask+1));
      v2[0]={Thigh(1), Thigh(0)};
      for (size_t i=1; i<v2.size(); ++i)
        v2[i]=rc.calc(i*(mask+1));
      }

    size_t size() const { return N; }
    size_t footprint() const { return (v1.size()+v2.size())*sizeof(cmplx_); }

    Tc operator[](size_t idx) const
      {
      if (2*idx<=N)
        {
        auto x1=v1[idx&mask], x2=v2[idx>>shift];
        return Tc(T(x1.r*x2.r-x1.i*x2.i), T(x1.r*x2.i+x1.i*x2.r));
        }
      idx = N-idx;
      auto x1=v1[idx&mask], x2=v2[idx>>shift];
      return Tc(T(x1.r*x2.r-x1.i*x2.i), -T(x1.r*x2.i+x1.i*x2.r));
      }
  };

// Rational approximation for floating-point numbers
/* f : number to convert.
 * num, denom: returned parts of the rational.
 * md: max denominator value.  Note that machine floating point number
 *     has a finite resolution (1e-16 ish for 64 bit double), so specifying
 *     a "best match with minimal error" is often wrong, because one can
 *     always just retrieve the significand and return that divided by
 *     2**52, which is in a sense accurate, but generally not very useful:
 *     1.0/7.0 would be "2573485501354569/18014398509481984", for example.
 */
template<typename T> static tuple<int64_t, int64_t> rat_approx0(T f, int64_t md)
  {
  static_assert(is_same<T,double>::value || is_same<T,float>::value,
    "unsupported floating point type");

  if (md <= 1) return make_tuple(int64_t(f), int64_t(1));

  // take care of negative numbers
  bool neg = f<0;
  if (neg) f = -f;

  // move number into interval [1;2]
  int64_t sub = int64_t(floor(f))-1;
  f -= sub;

  // Multiply number by a sufficiently large power of 2 to make it integer
  int64_t n = ((int64_t)1)<<53;
  f = ldexp(f,53);
  int64_t d = int64_t(f);

  /* continued fraction and check denominator each step */
  int64_t h[3] = { 0, 1, 0 }, k[3] = { 1, 0, 0 };
  for (int i = 0; i < 64; i++) {
    int64_t a = n ? d / n : 0;
    if (i && !a) break;

    int64_t x = d;
    d = n;
    n = x % n;

    x = a;
    if (k[1] * a + k[0] >= md) {
      x = (md - k[0]) / k[1];
      if ((x*2 >= a) || (k[1]>=md))
        i = 65;
      else
        break;
    }

    h[2] = x * h[1] + h[0]; h[0] = h[1]; h[1] = h[2];
    k[2] = x * k[1] + k[0]; k[0] = k[1]; k[1] = k[2];
  }
  return make_tuple(neg ? -(h[1]+k[1]*sub) : h[1]+k[1]*sub, k[1]);
}

template<typename T, typename Tc> class MultiExp
  {
  private:
    using Thigh = typename conditional<(sizeof(T)>sizeof(double)), T, double>::type;
    struct cmplx_ { Thigh r, i; };
    size_t N, mask, shift;
    vector<cmplx_> v1, v2;

  public:
    MultiExp(T ang0, size_t n)
      : N(n)
      {
      // Check if ang0 is a clean fraction of 2*pi.
      // If yes, use RootCalc internally, because this will be
      // more accurate.
      static constexpr T pi = T(3.141592653589793238462643383279502884197L);
      T frac = ang0/(2*pi);
      frac -= floor(frac);  // bring frac into [0;1] range
      // limit denominator to avoid spending too much time
      auto[num,den] = rat_approx0(frac, 1000000);
      auto ftden = frac*T(den), tnum = T(num);
      if ((tnum==ftden) || (nextafter(tnum,ftden)==ftden))
        {
        RootCalc<T,cmplx_> rc(den);
        size_t nval = n+2;
        shift = 1;
        while((size_t(1)<<shift)*(size_t(1)<<shift) < nval) ++shift;
        mask = (size_t(1)<<shift)-1;
        v1.resize(mask+1);
        v1[0]={Thigh(1), Thigh(0)};
        for (size_t i=1; i<v1.size(); ++i)
          v1[i] = rc.calc((i*num)%den);
        v2.resize((nval+mask)/(mask+1));
        v2[0]={Thigh(1), Thigh(0)};
        size_t inc = ((mask+1)*num)%den;
        for (size_t i=1; i<v2.size(); ++i)
          v2[i] = rc.calc((i*inc)%den);
        }
      else
        {
        Thigh ang = ang0;
        size_t nval = n+2;
        shift = 1;
        while((size_t(1)<<shift)*(size_t(1)<<shift) < nval) ++shift;
        mask = (size_t(1)<<shift)-1;
        v1.resize(mask+1);
        v1[0]={Thigh(1), Thigh(0)};
        for (size_t i=1; i<v1.size(); ++i)
          v1[i] = {cos(i*ang), sin(i*ang)};
        v2.resize((nval+mask)/(mask+1));
        v2[0]={Thigh(1), Thigh(0)};
        for (size_t i=1; i<v2.size(); ++i)
          v2[i] = {cos((i*(mask+1))*ang), sin((i*(mask+1))*ang)};
        }
      }

    size_t size() const { return N; }

    Tc operator[](size_t idx) const
      {
      auto x1=v1[idx&mask], x2=v2[idx>>shift];
      return Tc(T(x1.r*x2.r-x1.i*x2.i), T(x1.r*x2.i+x1.i*x2.r));
      }
  };

}

using detail_unity_roots::UnityRoots;
using detail_unity_roots::MultiExp;

}

#endif
