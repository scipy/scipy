/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file math_utils.h
 *  Various convenience mathematical functions.
 *
 *  Copyright (C) 2002-2021 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_MATH_UTILS_H
#define DUCC0_MATH_UTILS_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace ducc0 {

namespace detail_math_utils {

using namespace std;

/*! Returns \e true if | \a a-b | <= \a epsilon * | \a b |, else \e false. */
template<typename F> inline bool approx (F a, F b, F epsilon=F(1e-5L))
  { return abs(a-b) <= (epsilon*abs(b)); }

/*! Returns \e true if | \a a-b | <= \a epsilon, else \e false. */
template<typename F> inline bool abs_approx (F a, F b, F epsilon=F(1e-5L))
  { return abs(a-b) <= epsilon; }

/*! Returns the largest integer which is smaller than (or equal to) \a arg. */
template<typename I, typename F> inline I ifloor (F arg)
  { return I(floor(arg)); }

/*! Returns the integer which is nearest to \a arg. */
template<typename I, typename F> inline I nearest (F arg)
  { return ifloor<I>(arg+F(0.5)); }

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
template<typename F> inline F fmodulo (F v1, F v2)
  {
  if (v1>=0)
    return (v1<v2) ? v1 : fmod(v1,v2);
  auto tmp=fmod(v1,v2)+v2;
  return (tmp==v2) ? F(0) : tmp;
  }

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
template<typename I> inline I imodulo (I v1, I v2)
  { I v=v1%v2; return (v>=0) ? v : v+v2; }

/*! Returns -1 if \a signvalue is negative, else +1. */
template<typename T> inline T sign (const T& signvalue)
  { return (signvalue>=0) ? 1 : -1; }

/*! Returns \a val*pow(-1,m) */
template<typename T, typename I> inline T xpow (I m, T val)
  { return (m&1) ? -val : val; }

/*! Returns the integer \a n, which fulfills \a n*n<=arg<(n+1)*(n+1). */
template<typename I> inline uint32_t isqrt (I arg)
  {
  if constexpr (sizeof(I)<=4)
    return uint32_t (sqrt(arg+0.5));
  I res = I(sqrt(double(arg)+0.5));
  if (uint64_t(arg)<(uint64_t(1)<<50)) return uint32_t(res);
  if (res*res>arg)
    --res;
  else if ((res+1)*(res+1)<=arg)
    ++res;
  return uint32_t(res);
  }

/*! Returns the largest integer \a n that fulfills \a 2^n<=arg. */
template<typename I> inline unsigned int ilog2 (I arg)
  {
#ifdef __GNUC__
  if (arg==0) return 0;
  if constexpr (sizeof(I)==sizeof(int))
    return 8*sizeof(int)-1-__builtin_clz(arg);
  if constexpr (sizeof(I)==sizeof(long))
    return 8*sizeof(long)-1-__builtin_clzl(arg);
  if constexpr (sizeof(I)==sizeof(long long))
    return 8*sizeof(long long)-1-__builtin_clzll(arg);
#endif
  unsigned int res=0;
  while (arg > 0xFFFF) { res+=16; arg>>=16; }
  if (arg > 0x00FF) { res|=8; arg>>=8; }
  if (arg > 0x000F) { res|=4; arg>>=4; }
  if (arg > 0x0003) { res|=2; arg>>=2; }
  if (arg > 0x0001) { res|=1; }
  return res;
  }

template<typename I> inline unsigned int ilog2_nonnull (I arg)
  {
#ifdef __GNUC__
  if constexpr (sizeof(I)<=sizeof(int))
    return 8*sizeof(int)-1-__builtin_clz(arg);
  if constexpr (sizeof(I)==sizeof(long))
    return 8*sizeof(long)-1-__builtin_clzl(arg);
  if constexpr (sizeof(I)==sizeof(long long))
    return 8*sizeof(long long)-1-__builtin_clzll(arg);
#endif
  return ilog2 (arg);
  }

template<typename I> inline int trailingZeros(I arg)
  {
  if (arg==0) return sizeof(I)<<3;
#ifdef __GNUC__
  if constexpr (sizeof(I)<=sizeof(int))
    return __builtin_ctz(arg);
  if constexpr (sizeof(I)==sizeof(long))
    return __builtin_ctzl(arg);
  if constexpr (sizeof(I)==sizeof(long long))
    return __builtin_ctzll(arg);
#endif
  int res=0;
  while ((arg&0xFFFF)==0) { res+=16; arg>>=16; }
  if ((arg&0x00FF)==0) { res|=8; arg>>=8; }
  if ((arg&0x000F)==0) { res|=4; arg>>=4; }
  if ((arg&0x0003)==0) { res|=2; arg>>=2; }
  if ((arg&0x0001)==0) { res|=1; }
  return res;
  }

template<typename T1, typename T2>
inline bool multiequal (const T1 &a, const T2 &b)
  { return (a==b); }

template<typename T1, typename T2, typename... Args>
inline bool multiequal (const T1 &a, const T2 &b, Args... args)
  { return (a==b) && multiequal (a, args...); }

template<typename T> class tree_adder
  {
  private:
    vector<T> state;
    size_t n;

  public:
    tree_adder(): state(1,T(0)), n(0) {}

    void add (const T &val)
      {
      state[0]+=val; ++n;
      if (n>(size_t(1)<<(state.size()-1)))
        state.push_back(T(0));
      int shift=0;
      while (((n>>shift)&1)==0)
        {
        state[shift+1]+=state[shift];
        state[shift]=T(0);
        ++shift;
        }
      }
    T result() const
      {
      T sum(0);
      for (size_t i=0; i<state.size(); ++i)
        sum+=state[i];
      return sum;
      }
  };

/*! Returns \a atan2(y,x) if \a x!=0 or \a y!=0; else returns 0. */
inline double safe_atan2 (double y, double x)
  { return ((x==0.) && (y==0.)) ? 0.0 : atan2(y,x); }

}

using detail_math_utils::approx;
using detail_math_utils::abs_approx;
using detail_math_utils::ifloor;
using detail_math_utils::nearest;
using detail_math_utils::fmodulo;
using detail_math_utils::imodulo;
using detail_math_utils::sign;
using detail_math_utils::xpow;
using detail_math_utils::isqrt;
using detail_math_utils::ilog2;
using detail_math_utils::ilog2_nonnull;
using detail_math_utils::trailingZeros;
using detail_math_utils::multiequal;
using detail_math_utils::tree_adder;
using detail_math_utils::safe_atan2;

}

#endif
