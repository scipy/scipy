/*
This file is part of pocketfft.

Copyright (C) 2010-2019 Max-Planck-Society
Copyright (C) 2019 Peter Bell

For the odd-sized DCT-IV transforms:
  Copyright (C) 2003, 2007-14 Matteo Frigo
  Copyright (C) 2003, 2007-14 Massachusetts Institute of Technology

Authors: Martin Reinecke, Peter Bell

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

#ifndef POCKETFFT_HDRONLY_H
#define POCKETFFT_HDRONLY_H

#ifndef __cplusplus
#error This file is C++ and requires a C++ compiler.
#endif

#if !(__cplusplus >= 201103L || _MSVC_LANG+0L >= 201103L)
#error This file requires at least C++11 support.
#endif

#ifndef POCKETFFT_CACHE_SIZE
#define POCKETFFT_CACHE_SIZE 16
#endif

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <memory>
#include <vector>
#include <complex>
#if POCKETFFT_CACHE_SIZE!=0
#include <array>
#include <mutex>
#endif
#ifdef POCKETFFT_OPENMP
#include <omp.h>
#endif


#if defined(__GNUC__)
#define POCKETFFT_NOINLINE __attribute__((noinline))
#define POCKETFFT_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define POCKETFFT_NOINLINE __declspec(noinline)
#define POCKETFFT_RESTRICT __restrict
#else
#define POCKETFFT_NOINLINE
#define POCKETFFT_RESTRICT
#endif

namespace pocketfft {

namespace detail {

using namespace std;

using shape_t = vector<size_t>;
using stride_t = vector<ptrdiff_t>;

constexpr bool FORWARD  = true,
               BACKWARD = false;

// only enable vector support for gcc>=5.0 and clang>=5.0
#ifndef POCKETFFT_NO_VECTORS
#define POCKETFFT_NO_VECTORS
#if defined(__INTEL_COMPILER)
// do nothing. This is necessary because this compiler also sets __GNUC__.
#elif defined(__clang__)
#if __clang__>=5
#undef POCKETFFT_NO_VECTORS
#endif
#elif defined(__GNUC__)
#if __GNUC__>=5
#undef POCKETFFT_NO_VECTORS
#endif
#endif
#endif

template<typename T> struct VLEN { static constexpr size_t val=1; };

#ifndef POCKETFFT_NO_VECTORS
#if (defined(__AVX512F__))
template<> struct VLEN<float> { static constexpr size_t val=16; };
template<> struct VLEN<double> { static constexpr size_t val=8; };
#elif (defined(__AVX__))
template<> struct VLEN<float> { static constexpr size_t val=8; };
template<> struct VLEN<double> { static constexpr size_t val=4; };
#elif (defined(__SSE2__))
template<> struct VLEN<float> { static constexpr size_t val=4; };
template<> struct VLEN<double> { static constexpr size_t val=2; };
#elif (defined(__VSX__))
template<> struct VLEN<float> { static constexpr size_t val=4; };
template<> struct VLEN<double> { static constexpr size_t val=2; };
#else
#define POCKETFFT_NO_VECTORS
#endif
#endif

template<typename T> class arr
  {
  private:
    T *p;
    size_t sz;

#if defined(POCKETFFT_NO_VECTORS)
    static T *ralloc(size_t num)
      {
      if (num==0) return nullptr;
      void *res = malloc(num*sizeof(T));
      if (!res) throw bad_alloc();
      return reinterpret_cast<T *>(res);
      }
    static void dealloc(T *ptr)
      { free(ptr); }
#elif __cplusplus >= 201703L
    static T *ralloc(size_t num)
      {
      if (num==0) return nullptr;
      void *res = aligned_alloc(64,num*sizeof(T));
      if (!res) throw bad_alloc();
      return reinterpret_cast<T *>(res);
      }
    static void dealloc(T *ptr)
      { free(ptr); }
#else // portable emulation
    static T *ralloc(size_t num)
      {
      if (num==0) return nullptr;
      void *ptr = malloc(num*sizeof(T)+64);
      if (!ptr) throw bad_alloc();
      T *res = reinterpret_cast<T *>
        ((reinterpret_cast<size_t>(ptr) & ~(size_t(63))) + 64);
      (reinterpret_cast<void**>(res))[-1] = ptr;
      return res;
      }
    static void dealloc(T *ptr)
      { if (ptr) free((reinterpret_cast<void**>(ptr))[-1]); }
#endif

  public:
    arr() : p(0), sz(0) {}
    arr(size_t n) : p(ralloc(n)), sz(n) {}
    arr(arr &&other)
      : p(other.p), sz(other.sz)
      { other.p=nullptr; other.sz=0; }
    ~arr() { dealloc(p); }

    void resize(size_t n)
      {
      if (n==sz) return;
      dealloc(p);
      p = ralloc(n);
      sz = n;
      }

    T &operator[](size_t idx) { return p[idx]; }
    const T &operator[](size_t idx) const { return p[idx]; }

    T *data() { return p; }
    const T *data() const { return p; }

    size_t size() const { return sz; }
  };

template<typename T> struct cmplx {
  T r, i;
  cmplx() {}
  cmplx(T r_, T i_) : r(r_), i(i_) {}
  void Set(T r_, T i_) { r=r_; i=i_; }
  void Set(T r_) { r=r_; i=T(0); }
  void Split(T &r_, T &i_) const { r_=r; i_=i; }
  void SplitConj(T &r_, T &i_) const { r_=r; i_=-i; }
  cmplx &operator+= (const cmplx &other)
    { r+=other.r; i+=other.i; return *this; }
  template<typename T2>cmplx &operator*= (T2 other)
    { r*=other; i*=other; return *this; }
  template<typename T2>cmplx &operator*= (const cmplx<T2> &other)
    {
    T tmp = r*other.r - i*other.i;
    i = r*other.i + i*other.r;
    r = tmp;
    return *this;
    }
  template<typename T2>cmplx &operator+= (const cmplx<T2> &other)
    { r+=other.r; i+=other.i; return *this; }
  template<typename T2>cmplx &operator-= (const cmplx<T2> &other)
    { r-=other.r; i-=other.i; return *this; }
  template<typename T2> auto operator* (const T2 &other) const
    -> cmplx<decltype(r*other)>
    { return {r*other, i*other}; }
  template<typename T2> auto operator+ (const cmplx<T2> &other) const
    -> cmplx<decltype(r+other.r)>
    { return {r+other.r, i+other.i}; }
  template<typename T2> auto operator- (const cmplx<T2> &other) const
    -> cmplx<decltype(r+other.r)>
    { return {r-other.r, i-other.i}; }
  template<typename T2> auto operator* (const cmplx<T2> &other) const
    -> cmplx<decltype(r+other.r)>
    { return {r*other.r-i*other.i, r*other.i + i*other.r}; }
  template<bool fwd, typename T2> auto special_mul (const cmplx<T2> &other) const
    -> cmplx<decltype(r+other.r)>
    {
    using Tres = cmplx<decltype(r+other.r)>;
    return fwd ? Tres(r*other.r+i*other.i, i*other.r-r*other.i)
               : Tres(r*other.r-i*other.i, r*other.i+i*other.r);
    }
};
template<typename T> inline void PMINPLACE(T &a, T &b)
  { T t = a; a+=b; b=t-b; }
template<typename T> inline void MPINPLACE(T &a, T &b)
  { T t = a; a-=b; b=t+b; }
template<typename T> void PMC(T &a, T &b, const T &c, const T &d)
  { a = c+d; b = c-d; }
template<typename T> cmplx<T> conj(const cmplx<T> &a)
  { return {a.r, -a.i}; }

template<typename T> void ROT90(cmplx<T> &a)
  { auto tmp_=a.r; a.r=-a.i; a.i=tmp_; }
template<bool fwd, typename T> void ROTX90(cmplx<T> &a)
  { auto tmp_= fwd ? -a.r : a.r; a.r = fwd ? a.i : -a.i; a.i=tmp_; }

//
// twiddle factor section
//

template<typename T> class sincos_2pibyn
  {
  private:
    using Thigh = typename conditional<(sizeof(T)>sizeof(double)), T, double>::type;
    arr<T> data;

    void my_sincosm1pi (Thigh a_, Thigh *POCKETFFT_RESTRICT res)
      {
      if (sizeof(Thigh)>sizeof(double)) // don't have the code for long double
        {
        constexpr Thigh pi = Thigh(3.141592653589793238462643383279502884197L);
        auto s = sin(pi*a_);
        res[1] = s;
        res[0] = (s*s)/(-sqrt((1-s)*(1+s))-1);
        return;
        }
      // adapted from https://stackoverflow.com/questions/42792939/
      // CAUTION: this function only works for arguments in the range
      //          [-0.25; 0.25]!
      double a = double(a_);
      double s = a * a;
      /* Approximate cos(pi*x)-1 for x in [-0.25,0.25] */
      double r =     -1.0369917389758117e-4;
      r = fma (r, s,  1.9294935641298806e-3);
      r = fma (r, s, -2.5806887942825395e-2);
      r = fma (r, s,  2.3533063028328211e-1);
      r = fma (r, s, -1.3352627688538006e+0);
      r = fma (r, s,  4.0587121264167623e+0);
      r = fma (r, s, -4.9348022005446790e+0);
      double c = r*s;
      /* Approximate sin(pi*x) for x in [-0.25,0.25] */
      r =             4.6151442520157035e-4;
      r = fma (r, s, -7.3700183130883555e-3);
      r = fma (r, s,  8.2145868949323936e-2);
      r = fma (r, s, -5.9926452893214921e-1);
      r = fma (r, s,  2.5501640398732688e+0);
      r = fma (r, s, -5.1677127800499516e+0);
      s = s * a;
      r = r * s;
      s = fma (a, 3.1415926535897931e+0, r);
      res[0] = c;
      res[1] = s;
      }

    POCKETFFT_NOINLINE void calc_first_octant(size_t den,
      T * POCKETFFT_RESTRICT res)
      {
      size_t n = (den+4)>>3;
      if (n==0) return;
      res[0]=1.; res[1]=0.;
      if (n==1) return;
      size_t l1 = size_t(sqrt(n));
      arr<Thigh> tmp(2*l1);
      for (size_t i=1; i<l1; ++i)
        {
        my_sincosm1pi(Thigh(2*i)/Thigh(den),&tmp[2*i]);
        res[2*i  ] = T(tmp[2*i]+1);
        res[2*i+1] = T(tmp[2*i+1]);
        }
      size_t start=l1;
      while(start<n)
        {
        Thigh cs[2];
        my_sincosm1pi((Thigh(2*start))/Thigh(den),cs);
        res[2*start] = T(cs[0]+1);
        res[2*start+1] = T(cs[1]);
        size_t end = l1;
        if (start+end>n) end = n-start;
        for (size_t i=1; i<end; ++i)
          {
          Thigh csx[2]={tmp[2*i], tmp[2*i+1]};
          res[2*(start+i)] = T(((cs[0]*csx[0] - cs[1]*csx[1] + cs[0]) + csx[0]) + 1);
          res[2*(start+i)+1] = T((cs[0]*csx[1] + cs[1]*csx[0]) + cs[1] + csx[1]);
          }
        start += l1;
        }
      }

    void calc_first_quadrant(size_t n, T * POCKETFFT_RESTRICT res)
      {
      T * POCKETFFT_RESTRICT p = res+n;
      calc_first_octant(n<<1, p);
      size_t ndone=(n+2)>>2;
      size_t i=0, idx1=0, idx2=2*ndone-2;
      for (; i+1<ndone; i+=2, idx1+=2, idx2-=2)
        {
        res[idx1] = p[2*i  ]; res[idx1+1] = p[2*i+1];
        res[idx2] = p[2*i+3]; res[idx2+1] = p[2*i+2];
        }
      if (i!=ndone)
        { res[idx1] = p[2*i]; res[idx1+1] = p[2*i+1]; }
      }

    void calc_first_half(size_t n, T * POCKETFFT_RESTRICT res)
      {
      int ndone=int(n+1)>>1;
      T * p = res+n-1;
      calc_first_octant(n<<2, p);
      int i4=0, in=int(n), i=0;
      for (; i4<=in-i4; ++i, i4+=4) // octant 0
        { res[2*i] = p[2*i4]; res[2*i+1] = p[2*i4+1]; }
      for (; i4-in <= 0; ++i, i4+=4) // octant 1
        { auto xm = in-i4; res[2*i] = p[2*xm+1]; res[2*i+1] = p[2*xm]; }
      for (; i4<=3*in-i4; ++i, i4+=4) // octant 2
        { auto xm = i4-in; res[2*i] = -p[2*xm+1]; res[2*i+1] = p[2*xm]; }
      for (; i<ndone; ++i, i4+=4) // octant 3
        { auto xm = 2*in-i4; res[2*i] = -p[2*xm]; res[2*i+1] = p[2*xm+1]; }
      }

    void fill_first_quadrant(size_t n, T * POCKETFFT_RESTRICT res)
      {
      constexpr T hsqt2 = T(0.707106781186547524400844362104849L);
      size_t quart = n>>2;
      if ((n&7)==0)
        res[quart] = res[quart+1] = hsqt2;
      for (size_t i=2, j=2*quart-2; i<quart; i+=2, j-=2)
        { res[j] = res[i+1]; res[j+1] = res[i]; }
      }

    POCKETFFT_NOINLINE void fill_first_half(size_t n, T * POCKETFFT_RESTRICT res)
      {
      size_t half = n>>1;
      if ((n&3)==0)
        for (size_t i=0; i<half; i+=2)
          { res[i+half] = -res[i+1]; res[i+half+1] = res[i]; }
      else
        for (size_t i=2, j=2*half-2; i<half; i+=2, j-=2)
          { res[j] = -res[i]; res[j+1] = res[i+1]; }
      }

    void fill_second_half(size_t n, T * POCKETFFT_RESTRICT res)
      {
      if ((n&1)==0)
        for (size_t i=0; i<n; ++i)
          res[i+n] = -res[i];
      else
        for (size_t i=2, j=2*n-2; i<n; i+=2, j-=2)
          { res[j] = res[i]; res[j+1] = -res[i+1]; }
      }

    POCKETFFT_NOINLINE void sincos_2pibyn_half(size_t n, T * POCKETFFT_RESTRICT res)
      {
      if ((n&3)==0)
        {
        calc_first_octant(n, res);
        fill_first_quadrant(n, res);
        fill_first_half(n, res);
        }
      else if ((n&1)==0)
        {
        calc_first_quadrant(n, res);
        fill_first_half(n, res);
        }
      else
        calc_first_half(n, res);
      }

  public:
    POCKETFFT_NOINLINE sincos_2pibyn(size_t n, bool half)
      : data(2*n)
      {
      sincos_2pibyn_half(n, data.data());
      if (!half) fill_second_half(n, data.data());
      }

    T operator[](size_t idx) const { return data[idx]; }
    const T *rdata() const { return data; }
    const cmplx<T> *cdata() const
      { return reinterpret_cast<const cmplx<T> *>(data.data()); }
  };

struct util // hack to avoid duplicate symbols
  {
  static POCKETFFT_NOINLINE size_t largest_prime_factor (size_t n)
    {
    size_t res=1;
    while ((n&1)==0)
      { res=2; n>>=1; }
    for (size_t x=3; x*x<=n; x+=2)
      while ((n%x)==0)
        { res=x; n/=x; }
    if (n>1) res=n;
    return res;
    }

  static POCKETFFT_NOINLINE double cost_guess (size_t n)
    {
    constexpr double lfp=1.1; // penalty for non-hardcoded larger factors
    size_t ni=n;
    double result=0.;
    while ((n&1)==0)
      { result+=2; n>>=1; }
    for (size_t x=3; x*x<=n; x+=2)
      while ((n%x)==0)
        {
        result+= (x<=5) ? double(x) : lfp*double(x); // penalize larger prime factors
        n/=x;
        }
    if (n>1) result+=(n<=5) ? double(n) : lfp*double(n);
    return result*double(ni);
    }

  /* returns the smallest composite of 2, 3, 5, 7 and 11 which is >= n */
  static POCKETFFT_NOINLINE size_t good_size_cmplx(size_t n)
    {
    if (n<=12) return n;

    size_t bestfac=2*n;
    for (size_t f11=1; f11<bestfac; f11*=11)
      for (size_t f117=f11; f117<bestfac; f117*=7)
        for (size_t f1175=f117; f1175<bestfac; f1175*=5)
          {
          size_t x=f1175;
          while (x<n) x*=2;
          for (;;)
            {
            if (x<n)
              x*=3;
            else if (x>n)
              {
              if (x<bestfac) bestfac=x;
              if (x&1) break;
              x>>=1;
              }
            else
              return n;
            }
          }
    return bestfac;
    }

  /* returns the smallest composite of 2, 3, 5 which is >= n */
  static POCKETFFT_NOINLINE size_t good_size_real(size_t n)
    {
    if (n<=6) return n;

    size_t bestfac=2*n;
    for (size_t f5=1; f5<bestfac; f5*=5)
      {
      size_t x = f5;
      while (x<n) x *= 2;
      for (;;)
        {
        if (x<n)
          x*=3;
        else if (x>n)
          {
          if (x<bestfac) bestfac=x;
          if (x&1) break;
          x>>=1;
          }
        else
          return n;
        }
      }
    return bestfac;
    }

  static size_t prod(const shape_t &shape)
    {
    size_t res=1;
    for (auto sz: shape)
      res*=sz;
    return res;
    }

  static POCKETFFT_NOINLINE void sanity_check(const shape_t &shape,
    const stride_t &stride_in, const stride_t &stride_out, bool inplace)
    {
    auto ndim = shape.size();
    if (ndim<1) throw runtime_error("ndim must be >= 1");
    if ((stride_in.size()!=ndim) || (stride_out.size()!=ndim))
      throw runtime_error("stride dimension mismatch");
    if (inplace && (stride_in!=stride_out))
      throw runtime_error("stride mismatch");
    }

  static POCKETFFT_NOINLINE void sanity_check(const shape_t &shape,
    const stride_t &stride_in, const stride_t &stride_out, bool inplace,
    const shape_t &axes)
    {
    sanity_check(shape, stride_in, stride_out, inplace);
    auto ndim = shape.size();
    shape_t tmp(ndim,0);
    for (auto ax : axes)
      {
      if (ax>=ndim) throw invalid_argument("bad axis number");
      if (++tmp[ax]>1) throw invalid_argument("axis specified repeatedly");
      }
    }

  static POCKETFFT_NOINLINE void sanity_check(const shape_t &shape,
    const stride_t &stride_in, const stride_t &stride_out, bool inplace,
    size_t axis)
    {
    sanity_check(shape, stride_in, stride_out, inplace);
    if (axis>=shape.size()) throw invalid_argument("bad axis number");
    }

#ifdef POCKETFFT_OPENMP
    static size_t nthreads() { return size_t(omp_get_num_threads()); }
    static size_t thread_num() { return size_t(omp_get_thread_num()); }
    static size_t thread_count (size_t nthreads, const shape_t &shape,
      size_t axis)
      {
      if (nthreads==1) return 1;
      if (prod(shape) < 20*shape[axis]) return 1;
      return (nthreads==0) ? size_t(omp_get_max_threads()) : nthreads;
      }
#else
    static constexpr size_t nthreads() { return 1; }
    static constexpr size_t thread_num() { return 0; }
#endif
  };

//
// complex FFTPACK transforms
//

template<typename T0> class cfftp
  {
  private:
    struct fctdata
      {
      size_t fct;
      cmplx<T0> *tw, *tws;
      };

    size_t length;
    arr<cmplx<T0>> mem;
    vector<fctdata> fact;

    void add_factor(size_t factor)
      { fact.push_back({factor, nullptr, nullptr}); }

template<bool fwd, typename T> void pass2 (size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const cmplx<T0> * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=2;

  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto WA = [wa, ido](size_t x, size_t i)
    { return wa[i-1+x*(ido-1)]; };

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      CH(0,k,0) = CC(0,0,k)+CC(0,1,k);
      CH(0,k,1) = CC(0,0,k)-CC(0,1,k);
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      CH(0,k,0) = CC(0,0,k)+CC(0,1,k);
      CH(0,k,1) = CC(0,0,k)-CC(0,1,k);
      for (size_t i=1; i<ido; ++i)
        {
        CH(i,k,0) = CC(i,0,k)+CC(i,1,k);
        CH(i,k,1) = (CC(i,0,k)-CC(i,1,k)).template special_mul<fwd>(WA(0,i));
        }
      }
  }

#define POCKETFFT_PREP3(idx) \
        T t0 = CC(idx,0,k), t1, t2; \
        PMC (t1,t2,CC(idx,1,k),CC(idx,2,k)); \
        CH(idx,k,0)=t0+t1;
#define POCKETFFT_PARTSTEP3a(u1,u2,twr,twi) \
        { \
        T ca,cb; \
        ca=t0+t1*twr; \
        cb=t2*twi; ROT90(cb); \
        PMC(CH(0,k,u1),CH(0,k,u2),ca,cb) ;\
        }
#define POCKETFFT_PARTSTEP3b(u1,u2,twr,twi) \
        { \
        T ca,cb,da,db; \
        ca=t0+t1*twr; \
        cb=t2*twi; ROT90(cb); \
        PMC(da,db,ca,cb); \
        CH(i,k,u1) = da.template special_mul<fwd>(WA(u1-1,i)); \
        CH(i,k,u2) = db.template special_mul<fwd>(WA(u2-1,i)); \
        }
template<bool fwd, typename T> void pass3 (size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const cmplx<T0> * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=3;
  constexpr T0 tw1r=-0.5,
               tw1i= (fwd ? -1: 1) * T0(0.8660254037844386467637231707529362L);

  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto WA = [wa, ido](size_t x, size_t i)
    { return wa[i-1+x*(ido-1)]; };

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      POCKETFFT_PREP3(0)
      POCKETFFT_PARTSTEP3a(1,2,tw1r,tw1i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      POCKETFFT_PREP3(0)
      POCKETFFT_PARTSTEP3a(1,2,tw1r,tw1i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        POCKETFFT_PREP3(i)
        POCKETFFT_PARTSTEP3b(1,2,tw1r,tw1i)
        }
      }
  }

#undef POCKETFFT_PARTSTEP3b
#undef POCKETFFT_PARTSTEP3a
#undef POCKETFFT_PREP3

template<bool fwd, typename T> void pass4 (size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const cmplx<T0> * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=4;

  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto WA = [wa, ido](size_t x, size_t i)
    { return wa[i-1+x*(ido-1)]; };

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      T t1, t2, t3, t4;
      PMC(t2,t1,CC(0,0,k),CC(0,2,k));
      PMC(t3,t4,CC(0,1,k),CC(0,3,k));
      ROTX90<fwd>(t4);
      PMC(CH(0,k,0),CH(0,k,2),t2,t3);
      PMC(CH(0,k,1),CH(0,k,3),t1,t4);
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      T t1, t2, t3, t4;
      PMC(t2,t1,CC(0,0,k),CC(0,2,k));
      PMC(t3,t4,CC(0,1,k),CC(0,3,k));
      ROTX90<fwd>(t4);
      PMC(CH(0,k,0),CH(0,k,2),t2,t3);
      PMC(CH(0,k,1),CH(0,k,3),t1,t4);
      }
      for (size_t i=1; i<ido; ++i)
        {
        T c2, c3, c4, t1, t2, t3, t4;
        T cc0=CC(i,0,k), cc1=CC(i,1,k),cc2=CC(i,2,k),cc3=CC(i,3,k);
        PMC(t2,t1,cc0,cc2);
        PMC(t3,t4,cc1,cc3);
        ROTX90<fwd>(t4);
        PMC(CH(i,k,0),c3,t2,t3);
        PMC(c2,c4,t1,t4);
        CH(i,k,1) = c2.template special_mul<fwd>(WA(0,i));
        CH(i,k,2) = c3.template special_mul<fwd>(WA(1,i));
        CH(i,k,3) = c4.template special_mul<fwd>(WA(2,i));
        }
      }
  }

#define POCKETFFT_PREP5(idx) \
        T t0 = CC(idx,0,k), t1, t2, t3, t4; \
        PMC (t1,t4,CC(idx,1,k),CC(idx,4,k)); \
        PMC (t2,t3,CC(idx,2,k),CC(idx,3,k)); \
        CH(idx,k,0).r=t0.r+t1.r+t2.r; \
        CH(idx,k,0).i=t0.i+t1.i+t2.i;

#define POCKETFFT_PARTSTEP5a(u1,u2,twar,twbr,twai,twbi) \
        { \
        T ca,cb; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(CH(0,k,u1),CH(0,k,u2),ca,cb); \
        }

#define POCKETFFT_PARTSTEP5b(u1,u2,twar,twbr,twai,twbi) \
        { \
        T ca,cb,da,db; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(da,db,ca,cb); \
        CH(i,k,u1) = da.template special_mul<fwd>(WA(u1-1,i)); \
        CH(i,k,u2) = db.template special_mul<fwd>(WA(u2-1,i)); \
        }
template<bool fwd, typename T> void pass5 (size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const cmplx<T0> * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=5;
  constexpr T0 tw1r= T0(0.3090169943749474241022934171828191L),
               tw1i= (fwd ? -1: 1) * T0(0.9510565162951535721164393333793821L),
               tw2r= T0(-0.8090169943749474241022934171828191L),
               tw2i= (fwd ? -1: 1) * T0(0.5877852522924731291687059546390728L);

  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto WA = [wa, ido](size_t x, size_t i)
    { return wa[i-1+x*(ido-1)]; };

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      POCKETFFT_PREP5(0)
      POCKETFFT_PARTSTEP5a(1,4,tw1r,tw2r,+tw1i,+tw2i)
      POCKETFFT_PARTSTEP5a(2,3,tw2r,tw1r,+tw2i,-tw1i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      POCKETFFT_PREP5(0)
      POCKETFFT_PARTSTEP5a(1,4,tw1r,tw2r,+tw1i,+tw2i)
      POCKETFFT_PARTSTEP5a(2,3,tw2r,tw1r,+tw2i,-tw1i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        POCKETFFT_PREP5(i)
        POCKETFFT_PARTSTEP5b(1,4,tw1r,tw2r,+tw1i,+tw2i)
        POCKETFFT_PARTSTEP5b(2,3,tw2r,tw1r,+tw2i,-tw1i)
        }
      }
  }

#undef POCKETFFT_PARTSTEP5b
#undef POCKETFFT_PARTSTEP5a
#undef POCKETFFT_PREP5

#define POCKETFFT_PREP7(idx) \
        T t1 = CC(idx,0,k), t2, t3, t4, t5, t6, t7; \
        PMC (t2,t7,CC(idx,1,k),CC(idx,6,k)); \
        PMC (t3,t6,CC(idx,2,k),CC(idx,5,k)); \
        PMC (t4,t5,CC(idx,3,k),CC(idx,4,k)); \
        CH(idx,k,0).r=t1.r+t2.r+t3.r+t4.r; \
        CH(idx,k,0).i=t1.i+t2.i+t3.i+t4.i;

#define POCKETFFT_PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,out1,out2) \
        { \
        T ca,cb; \
        ca.r=t1.r+x1*t2.r+x2*t3.r+x3*t4.r; \
        ca.i=t1.i+x1*t2.i+x2*t3.i+x3*t4.i; \
        cb.i=y1*t7.r y2*t6.r y3*t5.r; \
        cb.r=-(y1*t7.i y2*t6.i y3*t5.i); \
        PMC(out1,out2,ca,cb); \
        }
#define POCKETFFT_PARTSTEP7a(u1,u2,x1,x2,x3,y1,y2,y3) \
        POCKETFFT_PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,CH(0,k,u1),CH(0,k,u2))
#define POCKETFFT_PARTSTEP7(u1,u2,x1,x2,x3,y1,y2,y3) \
        { \
        T da,db; \
        POCKETFFT_PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,da,db) \
        CH(i,k,u1) = da.template special_mul<fwd>(WA(u1-1,i)); \
        CH(i,k,u2) = db.template special_mul<fwd>(WA(u2-1,i)); \
        }

template<bool fwd, typename T> void pass7(size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const cmplx<T0> * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=7;
  constexpr T0 tw1r= T0(0.6234898018587335305250048840042398L),
               tw1i= (fwd ? -1 : 1) * T0(0.7818314824680298087084445266740578L),
               tw2r= T0(-0.2225209339563144042889025644967948L),
               tw2i= (fwd ? -1 : 1) * T0(0.9749279121818236070181316829939312L),
               tw3r= T0(-0.9009688679024191262361023195074451L),
               tw3i= (fwd ? -1 : 1) * T0(0.433883739117558120475768332848359L);

  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto WA = [wa, ido](size_t x, size_t i)
    { return wa[i-1+x*(ido-1)]; };

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      POCKETFFT_PREP7(0)
      POCKETFFT_PARTSTEP7a(1,6,tw1r,tw2r,tw3r,+tw1i,+tw2i,+tw3i)
      POCKETFFT_PARTSTEP7a(2,5,tw2r,tw3r,tw1r,+tw2i,-tw3i,-tw1i)
      POCKETFFT_PARTSTEP7a(3,4,tw3r,tw1r,tw2r,+tw3i,-tw1i,+tw2i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      POCKETFFT_PREP7(0)
      POCKETFFT_PARTSTEP7a(1,6,tw1r,tw2r,tw3r,+tw1i,+tw2i,+tw3i)
      POCKETFFT_PARTSTEP7a(2,5,tw2r,tw3r,tw1r,+tw2i,-tw3i,-tw1i)
      POCKETFFT_PARTSTEP7a(3,4,tw3r,tw1r,tw2r,+tw3i,-tw1i,+tw2i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        POCKETFFT_PREP7(i)
        POCKETFFT_PARTSTEP7(1,6,tw1r,tw2r,tw3r,+tw1i,+tw2i,+tw3i)
        POCKETFFT_PARTSTEP7(2,5,tw2r,tw3r,tw1r,+tw2i,-tw3i,-tw1i)
        POCKETFFT_PARTSTEP7(3,4,tw3r,tw1r,tw2r,+tw3i,-tw1i,+tw2i)
        }
      }
  }

#undef POCKETFFT_PARTSTEP7
#undef POCKETFFT_PARTSTEP7a0
#undef POCKETFFT_PARTSTEP7a
#undef POCKETFFT_PREP7

template <bool fwd, typename T> void ROTX45(T &a)
  {
  constexpr T0 hsqt2=T0(0.707106781186547524400844362104849L);
  if (fwd)
    { auto tmp_=a.r; a.r=hsqt2*(a.r+a.i); a.i=hsqt2*(a.i-tmp_); }
  else
    { auto tmp_=a.r; a.r=hsqt2*(a.r-a.i); a.i=hsqt2*(a.i+tmp_); }
  }
template <bool fwd, typename T> void ROTX135(T &a)
  {
  constexpr T0 hsqt2=T0(0.707106781186547524400844362104849L);
  if (fwd)
    { auto tmp_=a.r; a.r=hsqt2*(a.i-a.r); a.i=hsqt2*(-tmp_-a.i); }
  else
    { auto tmp_=a.r; a.r=hsqt2*(-a.r-a.i); a.i=hsqt2*(tmp_-a.i); }
  }

template<bool fwd, typename T> void pass8 (size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const cmplx<T0> * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=8;

  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto WA = [wa, ido](size_t x, size_t i)
    { return wa[i-1+x*(ido-1)]; };

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      T a0, a1, a2, a3, a4, a5, a6, a7;
      PMC(a1,a5,CC(0,1,k),CC(0,5,k));
      PMC(a3,a7,CC(0,3,k),CC(0,7,k));
      PMINPLACE(a1,a3);
      ROTX90<fwd>(a3);

      ROTX90<fwd>(a7);
      PMINPLACE(a5,a7);
      ROTX45<fwd>(a5);
      ROTX135<fwd>(a7);

      PMC(a0,a4,CC(0,0,k),CC(0,4,k));
      PMC(a2,a6,CC(0,2,k),CC(0,6,k));
      PMC(CH(0,k,0),CH(0,k,4),a0+a2,a1);
      PMC(CH(0,k,2),CH(0,k,6),a0-a2,a3);
      ROTX90<fwd>(a6);
      PMC(CH(0,k,1),CH(0,k,5),a4+a6,a5);
      PMC(CH(0,k,3),CH(0,k,7),a4-a6,a7);
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      T a0, a1, a2, a3, a4, a5, a6, a7;
      PMC(a1,a5,CC(0,1,k),CC(0,5,k));
      PMC(a3,a7,CC(0,3,k),CC(0,7,k));
      PMINPLACE(a1,a3);
      ROTX90<fwd>(a3);

      ROTX90<fwd>(a7);
      PMINPLACE(a5,a7);
      ROTX45<fwd>(a5);
      ROTX135<fwd>(a7);

      PMC(a0,a4,CC(0,0,k),CC(0,4,k));
      PMC(a2,a6,CC(0,2,k),CC(0,6,k));
      PMC(CH(0,k,0),CH(0,k,4),a0+a2,a1);
      PMC(CH(0,k,2),CH(0,k,6),a0-a2,a3);
      ROTX90<fwd>(a6);
      PMC(CH(0,k,1),CH(0,k,5),a4+a6,a5);
      PMC(CH(0,k,3),CH(0,k,7),a4-a6,a7);
      }
      for (size_t i=1; i<ido; ++i)
        {
        T a0, a1, a2, a3, a4, a5, a6, a7;
        PMC(a1,a5,CC(i,1,k),CC(i,5,k));
        PMC(a3,a7,CC(i,3,k),CC(i,7,k));
        ROTX90<fwd>(a7);
        PMINPLACE(a1,a3);
        ROTX90<fwd>(a3);
        PMINPLACE(a5,a7);
        ROTX45<fwd>(a5);
        ROTX135<fwd>(a7);
        PMC(a0,a4,CC(i,0,k),CC(i,4,k));
        PMC(a2,a6,CC(i,2,k),CC(i,6,k));
        PMINPLACE(a0,a2);
        CH(i,k,0) = a0+a1;
        CH(i,k,4) = (a0-a1).template special_mul<fwd>(WA(3,i));
        CH(i,k,2) = (a2+a3).template special_mul<fwd>(WA(1,i));
        CH(i,k,6) = (a2-a3).template special_mul<fwd>(WA(5,i));
        ROTX90<fwd>(a6);
        PMINPLACE(a4,a6);
        CH(i,k,1) = (a4+a5).template special_mul<fwd>(WA(0,i));
        CH(i,k,5) = (a4-a5).template special_mul<fwd>(WA(4,i));
        CH(i,k,3) = (a6+a7).template special_mul<fwd>(WA(2,i));
        CH(i,k,7) = (a6-a7).template special_mul<fwd>(WA(6,i));
        }
      }
   }


#define POCKETFFT_PREP11(idx) \
        T t1 = CC(idx,0,k), t2, t3, t4, t5, t6, t7, t8, t9, t10, t11; \
        PMC (t2,t11,CC(idx,1,k),CC(idx,10,k)); \
        PMC (t3,t10,CC(idx,2,k),CC(idx, 9,k)); \
        PMC (t4,t9 ,CC(idx,3,k),CC(idx, 8,k)); \
        PMC (t5,t8 ,CC(idx,4,k),CC(idx, 7,k)); \
        PMC (t6,t7 ,CC(idx,5,k),CC(idx, 6,k)); \
        CH(idx,k,0).r=t1.r+t2.r+t3.r+t4.r+t5.r+t6.r; \
        CH(idx,k,0).i=t1.i+t2.i+t3.i+t4.i+t5.i+t6.i;

#define POCKETFFT_PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,out1,out2) \
        { \
        T ca = t1 + t2*x1 + t3*x2 + t4*x3 + t5*x4 +t6*x5, \
          cb; \
        cb.i=y1*t11.r y2*t10.r y3*t9.r y4*t8.r y5*t7.r; \
        cb.r=-(y1*t11.i y2*t10.i y3*t9.i y4*t8.i y5*t7.i ); \
        PMC(out1,out2,ca,cb); \
        }
#define POCKETFFT_PARTSTEP11a(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5) \
        POCKETFFT_PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,CH(0,k,u1),CH(0,k,u2))
#define POCKETFFT_PARTSTEP11(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5) \
        { \
        T da,db; \
        POCKETFFT_PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,da,db) \
        CH(i,k,u1) = da.template special_mul<fwd>(WA(u1-1,i)); \
        CH(i,k,u2) = db.template special_mul<fwd>(WA(u2-1,i)); \
        }

template<bool fwd, typename T> void pass11 (size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const cmplx<T0> * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=11;
  constexpr T0 tw1r= T0(0.8412535328311811688618116489193677L),
               tw1i= (fwd ? -1 : 1) * T0(0.5406408174555975821076359543186917L),
               tw2r= T0(0.4154150130018864255292741492296232L),
               tw2i= (fwd ? -1 : 1) * T0(0.9096319953545183714117153830790285L),
               tw3r= T0(-0.1423148382732851404437926686163697L),
               tw3i= (fwd ? -1 : 1) * T0(0.9898214418809327323760920377767188L),
               tw4r= T0(-0.6548607339452850640569250724662936L),
               tw4i= (fwd ? -1 : 1) * T0(0.7557495743542582837740358439723444L),
               tw5r= T0(-0.9594929736144973898903680570663277L),
               tw5i= (fwd ? -1 : 1) * T0(0.2817325568414296977114179153466169L);

  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto WA = [wa, ido](size_t x, size_t i)
    { return wa[i-1+x*(ido-1)]; };

  if (ido==1)
    for (size_t k=0; k<l1; ++k)
      {
      POCKETFFT_PREP11(0)
      POCKETFFT_PARTSTEP11a(1,10,tw1r,tw2r,tw3r,tw4r,tw5r,+tw1i,+tw2i,+tw3i,+tw4i,+tw5i)
      POCKETFFT_PARTSTEP11a(2, 9,tw2r,tw4r,tw5r,tw3r,tw1r,+tw2i,+tw4i,-tw5i,-tw3i,-tw1i)
      POCKETFFT_PARTSTEP11a(3, 8,tw3r,tw5r,tw2r,tw1r,tw4r,+tw3i,-tw5i,-tw2i,+tw1i,+tw4i)
      POCKETFFT_PARTSTEP11a(4, 7,tw4r,tw3r,tw1r,tw5r,tw2r,+tw4i,-tw3i,+tw1i,+tw5i,-tw2i)
      POCKETFFT_PARTSTEP11a(5, 6,tw5r,tw1r,tw4r,tw2r,tw3r,+tw5i,-tw1i,+tw4i,-tw2i,+tw3i)
      }
  else
    for (size_t k=0; k<l1; ++k)
      {
      {
      POCKETFFT_PREP11(0)
      POCKETFFT_PARTSTEP11a(1,10,tw1r,tw2r,tw3r,tw4r,tw5r,+tw1i,+tw2i,+tw3i,+tw4i,+tw5i)
      POCKETFFT_PARTSTEP11a(2, 9,tw2r,tw4r,tw5r,tw3r,tw1r,+tw2i,+tw4i,-tw5i,-tw3i,-tw1i)
      POCKETFFT_PARTSTEP11a(3, 8,tw3r,tw5r,tw2r,tw1r,tw4r,+tw3i,-tw5i,-tw2i,+tw1i,+tw4i)
      POCKETFFT_PARTSTEP11a(4, 7,tw4r,tw3r,tw1r,tw5r,tw2r,+tw4i,-tw3i,+tw1i,+tw5i,-tw2i)
      POCKETFFT_PARTSTEP11a(5, 6,tw5r,tw1r,tw4r,tw2r,tw3r,+tw5i,-tw1i,+tw4i,-tw2i,+tw3i)
      }
      for (size_t i=1; i<ido; ++i)
        {
        POCKETFFT_PREP11(i)
        POCKETFFT_PARTSTEP11(1,10,tw1r,tw2r,tw3r,tw4r,tw5r,+tw1i,+tw2i,+tw3i,+tw4i,+tw5i)
        POCKETFFT_PARTSTEP11(2, 9,tw2r,tw4r,tw5r,tw3r,tw1r,+tw2i,+tw4i,-tw5i,-tw3i,-tw1i)
        POCKETFFT_PARTSTEP11(3, 8,tw3r,tw5r,tw2r,tw1r,tw4r,+tw3i,-tw5i,-tw2i,+tw1i,+tw4i)
        POCKETFFT_PARTSTEP11(4, 7,tw4r,tw3r,tw1r,tw5r,tw2r,+tw4i,-tw3i,+tw1i,+tw5i,-tw2i)
        POCKETFFT_PARTSTEP11(5, 6,tw5r,tw1r,tw4r,tw2r,tw3r,+tw5i,-tw1i,+tw4i,-tw2i,+tw3i)
        }
      }
  }

#undef PARTSTEP11
#undef PARTSTEP11a0
#undef PARTSTEP11a
#undef POCKETFFT_PREP11

template<bool fwd, typename T> void passg (size_t ido, size_t ip,
  size_t l1, T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const cmplx<T0> * POCKETFFT_RESTRICT wa,
  const cmplx<T0> * POCKETFFT_RESTRICT csarr)
  {
  const size_t cdim=ip;
  size_t ipph = (ip+1)/2;
  size_t idl1 = ido*l1;

  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto CX = [cc, ido, l1](size_t a, size_t b, size_t c) -> T&
    { return cc[a+ido*(b+l1*c)]; };
  auto CX2 = [cc, idl1](size_t a, size_t b) -> T&
    { return cc[a+idl1*b]; };
  auto CH2 = [ch, idl1](size_t a, size_t b) -> const T&
    { return ch[a+idl1*b]; };

  arr<cmplx<T0>> wal(ip);
  wal[0] = cmplx<T0>(1., 0.);
  for (size_t i=1; i<ip; ++i)
    wal[i]=cmplx<T0>(csarr[i].r,fwd ? -csarr[i].i : csarr[i].i);

  for (size_t k=0; k<l1; ++k)
    for (size_t i=0; i<ido; ++i)
      CH(i,k,0) = CC(i,0,k);
  for (size_t j=1, jc=ip-1; j<ipph; ++j, --jc)
    for (size_t k=0; k<l1; ++k)
      for (size_t i=0; i<ido; ++i)
        PMC(CH(i,k,j),CH(i,k,jc),CC(i,j,k),CC(i,jc,k));
  for (size_t k=0; k<l1; ++k)
    for (size_t i=0; i<ido; ++i)
      {
      T tmp = CH(i,k,0);
      for (size_t j=1; j<ipph; ++j)
        tmp+=CH(i,k,j);
      CX(i,k,0) = tmp;
      }
  for (size_t l=1, lc=ip-1; l<ipph; ++l, --lc)
    {
    // j=0
    for (size_t ik=0; ik<idl1; ++ik)
      {
      CX2(ik,l).r = CH2(ik,0).r+wal[l].r*CH2(ik,1).r+wal[2*l].r*CH2(ik,2).r;
      CX2(ik,l).i = CH2(ik,0).i+wal[l].r*CH2(ik,1).i+wal[2*l].r*CH2(ik,2).i;
      CX2(ik,lc).r=-wal[l].i*CH2(ik,ip-1).i-wal[2*l].i*CH2(ik,ip-2).i;
      CX2(ik,lc).i=wal[l].i*CH2(ik,ip-1).r+wal[2*l].i*CH2(ik,ip-2).r;
      }

    size_t iwal=2*l;
    size_t j=3, jc=ip-3;
    for (; j<ipph-1; j+=2, jc-=2)
      {
      iwal+=l; if (iwal>ip) iwal-=ip;
      cmplx<T0> xwal=wal[iwal];
      iwal+=l; if (iwal>ip) iwal-=ip;
      cmplx<T0> xwal2=wal[iwal];
      for (size_t ik=0; ik<idl1; ++ik)
        {
        CX2(ik,l).r += CH2(ik,j).r*xwal.r+CH2(ik,j+1).r*xwal2.r;
        CX2(ik,l).i += CH2(ik,j).i*xwal.r+CH2(ik,j+1).i*xwal2.r;
        CX2(ik,lc).r -= CH2(ik,jc).i*xwal.i+CH2(ik,jc-1).i*xwal2.i;
        CX2(ik,lc).i += CH2(ik,jc).r*xwal.i+CH2(ik,jc-1).r*xwal2.i;
        }
      }
    for (; j<ipph; ++j, --jc)
      {
      iwal+=l; if (iwal>ip) iwal-=ip;
      cmplx<T0> xwal=wal[iwal];
      for (size_t ik=0; ik<idl1; ++ik)
        {
        CX2(ik,l).r += CH2(ik,j).r*xwal.r;
        CX2(ik,l).i += CH2(ik,j).i*xwal.r;
        CX2(ik,lc).r -= CH2(ik,jc).i*xwal.i;
        CX2(ik,lc).i += CH2(ik,jc).r*xwal.i;
        }
      }
    }

  // shuffling and twiddling
  if (ido==1)
    for (size_t j=1, jc=ip-1; j<ipph; ++j, --jc)
      for (size_t ik=0; ik<idl1; ++ik)
        {
        T t1=CX2(ik,j), t2=CX2(ik,jc);
        PMC(CX2(ik,j),CX2(ik,jc),t1,t2);
        }
  else
    {
    for (size_t j=1, jc=ip-1; j<ipph; ++j,--jc)
      for (size_t k=0; k<l1; ++k)
        {
        T t1=CX(0,k,j), t2=CX(0,k,jc);
        PMC(CX(0,k,j),CX(0,k,jc),t1,t2);
        for (size_t i=1; i<ido; ++i)
          {
          T x1, x2;
          PMC(x1,x2,CX(i,k,j),CX(i,k,jc));
          size_t idij=(j-1)*(ido-1)+i-1;
          CX(i,k,j) = x1.template special_mul<fwd>(wa[idij]);
          idij=(jc-1)*(ido-1)+i-1;
          CX(i,k,jc) = x2.template special_mul<fwd>(wa[idij]);
          }
        }
    }
  }

template<bool fwd, typename T> void pass_all(T c[], T0 fct)
  {
  if (length==1) { c[0]*=fct; return; }
  size_t l1=1;
  arr<T> ch(length);
  T *p1=c, *p2=ch.data();

  for(size_t k1=0; k1<fact.size(); k1++)
    {
    size_t ip=fact[k1].fct;
    size_t l2=ip*l1;
    size_t ido = length/l2;
    if     (ip==4)
      pass4<fwd> (ido, l1, p1, p2, fact[k1].tw);
    else if(ip==8)
      pass8<fwd>(ido, l1, p1, p2, fact[k1].tw);
    else if(ip==2)
      pass2<fwd>(ido, l1, p1, p2, fact[k1].tw);
    else if(ip==3)
      pass3<fwd> (ido, l1, p1, p2, fact[k1].tw);
    else if(ip==5)
      pass5<fwd> (ido, l1, p1, p2, fact[k1].tw);
    else if(ip==7)
      pass7<fwd> (ido, l1, p1, p2, fact[k1].tw);
    else if(ip==11)
      pass11<fwd> (ido, l1, p1, p2, fact[k1].tw);
    else
      {
      passg<fwd>(ido, ip, l1, p1, p2, fact[k1].tw, fact[k1].tws);
      swap(p1,p2);
      }
    swap(p1,p2);
    l1=l2;
    }
  if (p1!=c)
    {
    if (fct!=1.)
      for (size_t i=0; i<length; ++i)
        c[i] = ch[i]*fct;
    else
      memcpy (c,p1,length*sizeof(T));
    }
  else
    if (fct!=1.)
      for (size_t i=0; i<length; ++i)
        c[i] *= fct;
  }

  public:
    template<typename T> void forward(T c[], T0 fct)
      { pass_all<true>(c, fct); }

    template<typename T> void backward(T c[], T0 fct)
      { pass_all<false>(c, fct); }

  private:
    POCKETFFT_NOINLINE void factorize()
      {
      size_t len=length;
      while ((len&7)==0)
        { add_factor(8); len>>=3; }
      while ((len&3)==0)
        { add_factor(4); len>>=2; }
      if ((len&1)==0)
        {
        len>>=1;
        // factor 2 should be at the front of the factor list
        add_factor(2);
        swap(fact[0].fct, fact.back().fct);
        }
      for (size_t divisor=3; divisor*divisor<=len; divisor+=2)
        while ((len%divisor)==0)
          {
          add_factor(divisor);
          len/=divisor;
          }
      if (len>1) add_factor(len);
      }

    size_t twsize() const
      {
      size_t twsize=0, l1=1;
      for (size_t k=0; k<fact.size(); ++k)
        {
        size_t ip=fact[k].fct, ido= length/(l1*ip);
        twsize+=(ip-1)*(ido-1);
        if (ip>11)
          twsize+=ip;
        l1*=ip;
        }
      return twsize;
      }

    void comp_twiddle()
      {
      sincos_2pibyn<T0> twid(length, false);
      auto twiddle = twid.cdata();
      size_t l1=1;
      size_t memofs=0;
      for (size_t k=0; k<fact.size(); ++k)
        {
        size_t ip=fact[k].fct, ido=length/(l1*ip);
        fact[k].tw=mem.data()+memofs;
        memofs+=(ip-1)*(ido-1);
        for (size_t j=1; j<ip; ++j)
          for (size_t i=1; i<ido; ++i)
            fact[k].tw[(j-1)*(ido-1)+i-1] = twiddle[j*l1*i];
        if (ip>11)
          {
          fact[k].tws=mem.data()+memofs;
          memofs+=ip;
          for (size_t j=0; j<ip; ++j)
            fact[k].tws[j] = twiddle[j*l1*ido];
          }
        l1*=ip;
        }
      }

  public:
    POCKETFFT_NOINLINE cfftp(size_t length_)
      : length(length_)
      {
      if (length==0) throw runtime_error("zero-length FFT requested");
      if (length==1) return;
      factorize();
      mem.resize(twsize());
      comp_twiddle();
      }
  };

//
// real-valued FFTPACK transforms
//

template<typename T0> class rfftp
  {
  private:
    struct fctdata
      {
      size_t fct;
      T0 *tw, *tws;
      };

    size_t length;
    arr<T0> mem;
    vector<fctdata> fact;

    void add_factor(size_t factor)
      { fact.push_back({factor, nullptr, nullptr}); }

template<typename T> inline void PM(T &a, T &b, T c, T d)
  { a=c+d; b=c-d; }

/* (a+ib) = conj(c+id) * (e+if) */
template<typename T1, typename T2, typename T3> inline void MULPM
  (T1 &a, T1 &b, T2 c, T2 d, T3 e, T3 f)
  {  a=c*e+d*f; b=c*f-d*e; }

template<typename T> void radf2 (size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=2;

  auto WA = [wa,ido](size_t x, size_t i) { return wa[i+x*(ido-1)]; };
  auto CC = [cc,ido,l1](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+l1*c)]; };
  auto CH = [ch,ido,cdim](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+cdim*c)]; };

  for (size_t k=0; k<l1; k++)
    PM (CH(0,0,k),CH(ido-1,1,k),CC(0,k,0),CC(0,k,1));
  if ((ido&1)==0)
    for (size_t k=0; k<l1; k++)
      {
      CH(    0,1,k) = -CC(ido-1,k,1);
      CH(ido-1,0,k) =  CC(ido-1,k,0);
      }
  if (ido<=2) return;
  for (size_t k=0; k<l1; k++)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      T tr2, ti2;
      MULPM (tr2,ti2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1));
      PM (CH(i-1,0,k),CH(ic-1,1,k),CC(i-1,k,0),tr2);
      PM (CH(i  ,0,k),CH(ic  ,1,k),ti2,CC(i  ,k,0));
      }
  }

// a2=a+b; b2=i*(b-a);
#define POCKETFFT_REARRANGE(rx, ix, ry, iy) \
  {\
  auto t1=rx+ry, t2=ry-rx, t3=ix+iy, t4=ix-iy; \
  rx=t1; ix=t3; ry=t4; iy=t2; \
  }

template<typename T> void radf3(size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=3;
  constexpr T0 taur=-0.5, taui=T0(0.8660254037844386467637231707529362L);

  auto WA = [wa,ido](size_t x, size_t i) { return wa[i+x*(ido-1)]; };
  auto CC = [cc,ido,l1](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+l1*c)]; };
  auto CH = [ch,ido,cdim](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+cdim*c)]; };

  for (size_t k=0; k<l1; k++)
    {
    T cr2=CC(0,k,1)+CC(0,k,2);
    CH(0,0,k) = CC(0,k,0)+cr2;
    CH(0,2,k) = taui*(CC(0,k,2)-CC(0,k,1));
    CH(ido-1,1,k) = CC(0,k,0)+taur*cr2;
    }
  if (ido==1) return;
  for (size_t k=0; k<l1; k++)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      T di2, di3, dr2, dr3;
      MULPM (dr2,di2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1)); // d2=conj(WA0)*CC1
      MULPM (dr3,di3,WA(1,i-2),WA(1,i-1),CC(i-1,k,2),CC(i,k,2)); // d3=conj(WA1)*CC2
      POCKETFFT_REARRANGE(dr2, di2, dr3, di3);
      CH(i-1,0,k) = CC(i-1,k,0)+dr2; // c add
      CH(i  ,0,k) = CC(i  ,k,0)+di2;
      T tr2 = CC(i-1,k,0)+taur*dr2; // c add
      T ti2 = CC(i  ,k,0)+taur*di2;
      T tr3 = taui*dr3;  // t3 = taui*i*(d3-d2)?
      T ti3 = taui*di3;
      PM(CH(i-1,2,k),CH(ic-1,1,k),tr2,tr3); // PM(i) = t2+t3
      PM(CH(i  ,2,k),CH(ic  ,1,k),ti3,ti2); // PM(ic) = conj(t2-t3)
      }
  }

template<typename T> void radf4(size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=4;
  constexpr T0 hsqt2=T0(0.707106781186547524400844362104849L);

  auto WA = [wa,ido](size_t x, size_t i) { return wa[i+x*(ido-1)]; };
  auto CC = [cc,ido,l1](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+l1*c)]; };
  auto CH = [ch,ido,cdim](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+cdim*c)]; };

  for (size_t k=0; k<l1; k++)
    {
    T tr1,tr2;
    PM (tr1,CH(0,2,k),CC(0,k,3),CC(0,k,1));
    PM (tr2,CH(ido-1,1,k),CC(0,k,0),CC(0,k,2));
    PM (CH(0,0,k),CH(ido-1,3,k),tr2,tr1);
    }
  if ((ido&1)==0)
    for (size_t k=0; k<l1; k++)
      {
      T ti1=-hsqt2*(CC(ido-1,k,1)+CC(ido-1,k,3));
      T tr1= hsqt2*(CC(ido-1,k,1)-CC(ido-1,k,3));
      PM (CH(ido-1,0,k),CH(ido-1,2,k),CC(ido-1,k,0),tr1);
      PM (CH(    0,3,k),CH(    0,1,k),ti1,CC(ido-1,k,2));
      }
  if (ido<=2) return;
  for (size_t k=0; k<l1; k++)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      T ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
      MULPM(cr2,ci2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1));
      MULPM(cr3,ci3,WA(1,i-2),WA(1,i-1),CC(i-1,k,2),CC(i,k,2));
      MULPM(cr4,ci4,WA(2,i-2),WA(2,i-1),CC(i-1,k,3),CC(i,k,3));
      PM(tr1,tr4,cr4,cr2);
      PM(ti1,ti4,ci2,ci4);
      PM(tr2,tr3,CC(i-1,k,0),cr3);
      PM(ti2,ti3,CC(i  ,k,0),ci3);
      PM(CH(i-1,0,k),CH(ic-1,3,k),tr2,tr1);
      PM(CH(i  ,0,k),CH(ic  ,3,k),ti1,ti2);
      PM(CH(i-1,2,k),CH(ic-1,1,k),tr3,ti4);
      PM(CH(i  ,2,k),CH(ic  ,1,k),tr4,ti3);
      }
  }

template<typename T> void radf5(size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=5;
  constexpr T0 tr11= T0(0.3090169943749474241022934171828191L),
               ti11= T0(0.9510565162951535721164393333793821L),
               tr12= T0(-0.8090169943749474241022934171828191L),
               ti12= T0(0.5877852522924731291687059546390728L);

  auto WA = [wa,ido](size_t x, size_t i) { return wa[i+x*(ido-1)]; };
  auto CC = [cc,ido,l1](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+l1*c)]; };
  auto CH = [ch,ido,cdim](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+cdim*c)]; };

  for (size_t k=0; k<l1; k++)
    {
    T cr2, cr3, ci4, ci5;
    PM (cr2,ci5,CC(0,k,4),CC(0,k,1));
    PM (cr3,ci4,CC(0,k,3),CC(0,k,2));
    CH(0,0,k)=CC(0,k,0)+cr2+cr3;
    CH(ido-1,1,k)=CC(0,k,0)+tr11*cr2+tr12*cr3;
    CH(0,2,k)=ti11*ci5+ti12*ci4;
    CH(ido-1,3,k)=CC(0,k,0)+tr12*cr2+tr11*cr3;
    CH(0,4,k)=ti12*ci5-ti11*ci4;
    }
  if (ido==1) return;
  for (size_t k=0; k<l1;++k)
    for (size_t i=2, ic=ido-2; i<ido; i+=2, ic-=2)
      {
      T di2, di3, di4, di5, dr2, dr3, dr4, dr5;
      MULPM (dr2,di2,WA(0,i-2),WA(0,i-1),CC(i-1,k,1),CC(i,k,1));
      MULPM (dr3,di3,WA(1,i-2),WA(1,i-1),CC(i-1,k,2),CC(i,k,2));
      MULPM (dr4,di4,WA(2,i-2),WA(2,i-1),CC(i-1,k,3),CC(i,k,3));
      MULPM (dr5,di5,WA(3,i-2),WA(3,i-1),CC(i-1,k,4),CC(i,k,4));
      POCKETFFT_REARRANGE(dr2, di2, dr5, di5);
      POCKETFFT_REARRANGE(dr3, di3, dr4, di4);
      CH(i-1,0,k)=CC(i-1,k,0)+dr2+dr3;
      CH(i  ,0,k)=CC(i  ,k,0)+di2+di3;
      T tr2=CC(i-1,k,0)+tr11*dr2+tr12*dr3;
      T ti2=CC(i  ,k,0)+tr11*di2+tr12*di3;
      T tr3=CC(i-1,k,0)+tr12*dr2+tr11*dr3;
      T ti3=CC(i  ,k,0)+tr12*di2+tr11*di3;
      T tr5 = ti11*dr5 + ti12*dr4;
      T ti5 = ti11*di5 + ti12*di4;
      T tr4 = ti12*dr5 - ti11*dr4;
      T ti4 = ti12*di5 - ti11*di4;
      PM(CH(i-1,2,k),CH(ic-1,1,k),tr2,tr5);
      PM(CH(i  ,2,k),CH(ic  ,1,k),ti5,ti2);
      PM(CH(i-1,4,k),CH(ic-1,3,k),tr3,tr4);
      PM(CH(i  ,4,k),CH(ic  ,3,k),ti4,ti3);
      }
  }

#undef POCKETFFT_REARRANGE

template<typename T> void radfg(size_t ido, size_t ip, size_t l1,
  T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa, const T0 * POCKETFFT_RESTRICT csarr)
  {
  const size_t cdim=ip;
  size_t ipph=(ip+1)/2;
  size_t idl1 = ido*l1;

  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> const T&
    { return ch[a+ido*(b+l1*c)]; };
  auto C1 = [cc,ido,l1] (size_t a, size_t b, size_t c) -> T&
    { return cc[a+ido*(b+l1*c)]; };
  auto C2 = [cc,idl1] (size_t a, size_t b) -> T&
    { return cc[a+idl1*b]; };
  auto CH2 = [ch,idl1] (size_t a, size_t b) -> T&
    { return ch[a+idl1*b]; };

  if (ido>1)
    {
    for (size_t j=1, jc=ip-1; j<ipph; ++j,--jc)              // 114
      {
      size_t is=(j-1)*(ido-1),
             is2=(jc-1)*(ido-1);
      for (size_t k=0; k<l1; ++k)                            // 113
        {
        size_t idij=is;
        size_t idij2=is2;
        for (size_t i=1; i<=ido-2; i+=2)                      // 112
          {
          T t1=C1(i,k,j ), t2=C1(i+1,k,j ),
                 t3=C1(i,k,jc), t4=C1(i+1,k,jc);
          T x1=wa[idij]*t1 + wa[idij+1]*t2,
                 x2=wa[idij]*t2 - wa[idij+1]*t1,
                 x3=wa[idij2]*t3 + wa[idij2+1]*t4,
                 x4=wa[idij2]*t4 - wa[idij2+1]*t3;
          C1(i  ,k,j ) = x1+x3;
          C1(i  ,k,jc) = x2-x4;
          C1(i+1,k,j ) = x2+x4;
          C1(i+1,k,jc) = x3-x1;
          idij+=2;
          idij2+=2;
          }
        }
      }
    }

  for (size_t j=1, jc=ip-1; j<ipph; ++j,--jc)                // 123
    for (size_t k=0; k<l1; ++k)                              // 122
      {
      T t1=C1(0,k,j), t2=C1(0,k,jc);
      C1(0,k,j ) = t1+t2;
      C1(0,k,jc) = t2-t1;
      }

//everything in C
//memset(ch,0,ip*l1*ido*sizeof(double));

  for (size_t l=1,lc=ip-1; l<ipph; ++l,--lc)                 // 127
    {
    for (size_t ik=0; ik<idl1; ++ik)                         // 124
      {
      CH2(ik,l ) = C2(ik,0)+csarr[2*l]*C2(ik,1)+csarr[4*l]*C2(ik,2);
      CH2(ik,lc) = csarr[2*l+1]*C2(ik,ip-1)+csarr[4*l+1]*C2(ik,ip-2);
      }
    size_t iang = 2*l;
    size_t j=3, jc=ip-3;
    for (; j<ipph-3; j+=4,jc-=4)              // 126
      {
      iang+=l; if (iang>=ip) iang-=ip;
      T0 ar1=csarr[2*iang], ai1=csarr[2*iang+1];
      iang+=l; if (iang>=ip) iang-=ip;
      T0 ar2=csarr[2*iang], ai2=csarr[2*iang+1];
      iang+=l; if (iang>=ip) iang-=ip;
      T0 ar3=csarr[2*iang], ai3=csarr[2*iang+1];
      iang+=l; if (iang>=ip) iang-=ip;
      T0 ar4=csarr[2*iang], ai4=csarr[2*iang+1];
      for (size_t ik=0; ik<idl1; ++ik)                       // 125
        {
        CH2(ik,l ) += ar1*C2(ik,j )+ar2*C2(ik,j +1)
                     +ar3*C2(ik,j +2)+ar4*C2(ik,j +3);
        CH2(ik,lc) += ai1*C2(ik,jc)+ai2*C2(ik,jc-1)
                     +ai3*C2(ik,jc-2)+ai4*C2(ik,jc-3);
        }
      }
    for (; j<ipph-1; j+=2,jc-=2)              // 126
      {
      iang+=l; if (iang>=ip) iang-=ip;
      T0 ar1=csarr[2*iang], ai1=csarr[2*iang+1];
      iang+=l; if (iang>=ip) iang-=ip;
      T0 ar2=csarr[2*iang], ai2=csarr[2*iang+1];
      for (size_t ik=0; ik<idl1; ++ik)                       // 125
        {
        CH2(ik,l ) += ar1*C2(ik,j )+ar2*C2(ik,j +1);
        CH2(ik,lc) += ai1*C2(ik,jc)+ai2*C2(ik,jc-1);
        }
      }
    for (; j<ipph; ++j,--jc)              // 126
      {
      iang+=l; if (iang>=ip) iang-=ip;
      T0 ar=csarr[2*iang], ai=csarr[2*iang+1];
      for (size_t ik=0; ik<idl1; ++ik)                       // 125
        {
        CH2(ik,l ) += ar*C2(ik,j );
        CH2(ik,lc) += ai*C2(ik,jc);
        }
      }
    }
  for (size_t ik=0; ik<idl1; ++ik)                         // 101
    CH2(ik,0) = C2(ik,0);
  for (size_t j=1; j<ipph; ++j)                              // 129
    for (size_t ik=0; ik<idl1; ++ik)                         // 128
      CH2(ik,0) += C2(ik,j);

// everything in CH at this point!
//memset(cc,0,ip*l1*ido*sizeof(double));

  for (size_t k=0; k<l1; ++k)                                // 131
    for (size_t i=0; i<ido; ++i)                             // 130
      CC(i,0,k) = CH(i,k,0);

  for (size_t j=1, jc=ip-1; j<ipph; ++j,--jc)                // 137
    {
    size_t j2=2*j-1;
    for (size_t k=0; k<l1; ++k)                              // 136
      {
      CC(ido-1,j2,k) = CH(0,k,j);
      CC(0,j2+1,k) = CH(0,k,jc);
      }
    }

  if (ido==1) return;

  for (size_t j=1, jc=ip-1; j<ipph; ++j,--jc)                // 140
    {
    size_t j2=2*j-1;
    for(size_t k=0; k<l1; ++k)                               // 139
      for(size_t i=1, ic=ido-i-2; i<=ido-2; i+=2, ic-=2)      // 138
        {
        CC(i   ,j2+1,k) = CH(i  ,k,j )+CH(i  ,k,jc);
        CC(ic  ,j2  ,k) = CH(i  ,k,j )-CH(i  ,k,jc);
        CC(i+1 ,j2+1,k) = CH(i+1,k,j )+CH(i+1,k,jc);
        CC(ic+1,j2  ,k) = CH(i+1,k,jc)-CH(i+1,k,j );
        }
    }
  }

template<typename T> void radb2(size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=2;

  auto WA = [wa,ido](size_t x, size_t i) { return wa[i+x*(ido-1)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };

  for (size_t k=0; k<l1; k++)
    PM (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(ido-1,1,k));
  if ((ido&1)==0)
    for (size_t k=0; k<l1; k++)
      {
      CH(ido-1,k,0) = 2*CC(ido-1,0,k);
      CH(ido-1,k,1) =-2*CC(0    ,1,k);
      }
  if (ido<=2) return;
  for (size_t k=0; k<l1;++k)
    for (size_t i=2; i<ido; i+=2)
      {
      size_t ic=ido-i;
      T ti2, tr2;
      PM (CH(i-1,k,0),tr2,CC(i-1,0,k),CC(ic-1,1,k));
      PM (ti2,CH(i  ,k,0),CC(i  ,0,k),CC(ic  ,1,k));
      MULPM (CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),ti2,tr2);
      }
  }

template<typename T> void radb3(size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=3;
  constexpr T0 taur=-0.5, taui=T0(0.8660254037844386467637231707529362L);

  auto WA = [wa,ido](size_t x, size_t i) { return wa[i+x*(ido-1)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };

  for (size_t k=0; k<l1; k++)
    {
    T tr2=2*CC(ido-1,1,k);
    T cr2=CC(0,0,k)+taur*tr2;
    CH(0,k,0)=CC(0,0,k)+tr2;
    T ci3=2*taui*CC(0,2,k);
    PM (CH(0,k,2),CH(0,k,1),cr2,ci3);
    }
  if (ido==1) return;
  for (size_t k=0; k<l1; k++)
    for (size_t i=2, ic=ido-2; i<ido; i+=2, ic-=2)
      {
      T tr2=CC(i-1,2,k)+CC(ic-1,1,k); // t2=CC(I) + conj(CC(ic))
      T ti2=CC(i  ,2,k)-CC(ic  ,1,k);
      T cr2=CC(i-1,0,k)+taur*tr2;     // c2=CC +taur*t2
      T ci2=CC(i  ,0,k)+taur*ti2;
      CH(i-1,k,0)=CC(i-1,0,k)+tr2;         // CH=CC+t2
      CH(i  ,k,0)=CC(i  ,0,k)+ti2;
      T cr3=taui*(CC(i-1,2,k)-CC(ic-1,1,k));// c3=taui*(CC(i)-conj(CC(ic)))
      T ci3=taui*(CC(i  ,2,k)+CC(ic  ,1,k));
      T di2, di3, dr2, dr3;
      PM(dr3,dr2,cr2,ci3); // d2= (cr2-ci3, ci2+cr3) = c2+i*c3
      PM(di2,di3,ci2,cr3); // d3= (cr2+ci3, ci2-cr3) = c2-i*c3
      MULPM(CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),di2,dr2); // ch = WA*d2
      MULPM(CH(i,k,2),CH(i-1,k,2),WA(1,i-2),WA(1,i-1),di3,dr3);
      }
  }

template<typename T> void radb4(size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=4;
  constexpr T0 sqrt2=T0(1.414213562373095048801688724209698L);

  auto WA = [wa,ido](size_t x, size_t i) { return wa[i+x*(ido-1)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };

  for (size_t k=0; k<l1; k++)
    {
    T tr1, tr2;
    PM (tr2,tr1,CC(0,0,k),CC(ido-1,3,k));
    T tr3=2*CC(ido-1,1,k);
    T tr4=2*CC(0,2,k);
    PM (CH(0,k,0),CH(0,k,2),tr2,tr3);
    PM (CH(0,k,3),CH(0,k,1),tr1,tr4);
    }
  if ((ido&1)==0)
    for (size_t k=0; k<l1; k++)
      {
      T tr1,tr2,ti1,ti2;
      PM (ti1,ti2,CC(0    ,3,k),CC(0    ,1,k));
      PM (tr2,tr1,CC(ido-1,0,k),CC(ido-1,2,k));
      CH(ido-1,k,0)=tr2+tr2;
      CH(ido-1,k,1)=sqrt2*(tr1-ti1);
      CH(ido-1,k,2)=ti2+ti2;
      CH(ido-1,k,3)=-sqrt2*(tr1+ti1);
      }
  if (ido<=2) return;
  for (size_t k=0; k<l1;++k)
    for (size_t i=2; i<ido; i+=2)
      {
      T ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
      size_t ic=ido-i;
      PM (tr2,tr1,CC(i-1,0,k),CC(ic-1,3,k));
      PM (ti1,ti2,CC(i  ,0,k),CC(ic  ,3,k));
      PM (tr4,ti3,CC(i  ,2,k),CC(ic  ,1,k));
      PM (tr3,ti4,CC(i-1,2,k),CC(ic-1,1,k));
      PM (CH(i-1,k,0),cr3,tr2,tr3);
      PM (CH(i  ,k,0),ci3,ti2,ti3);
      PM (cr4,cr2,tr1,tr4);
      PM (ci2,ci4,ti1,ti4);
      MULPM (CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),ci2,cr2);
      MULPM (CH(i,k,2),CH(i-1,k,2),WA(1,i-2),WA(1,i-1),ci3,cr3);
      MULPM (CH(i,k,3),CH(i-1,k,3),WA(2,i-2),WA(2,i-1),ci4,cr4);
      }
  }

template<typename T> void radb5(size_t ido, size_t l1,
  const T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa)
  {
  constexpr size_t cdim=5;
  constexpr T0 tr11= T0(0.3090169943749474241022934171828191L),
               ti11= T0(0.9510565162951535721164393333793821L),
               tr12= T0(-0.8090169943749474241022934171828191L),
               ti12= T0(0.5877852522924731291687059546390728L);

  auto WA = [wa,ido](size_t x, size_t i) { return wa[i+x*(ido-1)]; };
  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };

  for (size_t k=0; k<l1; k++)
    {
    T ti5=CC(0,2,k)+CC(0,2,k);
    T ti4=CC(0,4,k)+CC(0,4,k);
    T tr2=CC(ido-1,1,k)+CC(ido-1,1,k);
    T tr3=CC(ido-1,3,k)+CC(ido-1,3,k);
    CH(0,k,0)=CC(0,0,k)+tr2+tr3;
    T cr2=CC(0,0,k)+tr11*tr2+tr12*tr3;
    T cr3=CC(0,0,k)+tr12*tr2+tr11*tr3;
    T ci4, ci5;
    MULPM(ci5,ci4,ti5,ti4,ti11,ti12);
    PM(CH(0,k,4),CH(0,k,1),cr2,ci5);
    PM(CH(0,k,3),CH(0,k,2),cr3,ci4);
    }
  if (ido==1) return;
  for (size_t k=0; k<l1;++k)
    for (size_t i=2, ic=ido-2; i<ido; i+=2, ic-=2)
      {
      T tr2, tr3, tr4, tr5, ti2, ti3, ti4, ti5;
      PM(tr2,tr5,CC(i-1,2,k),CC(ic-1,1,k));
      PM(ti5,ti2,CC(i  ,2,k),CC(ic  ,1,k));
      PM(tr3,tr4,CC(i-1,4,k),CC(ic-1,3,k));
      PM(ti4,ti3,CC(i  ,4,k),CC(ic  ,3,k));
      CH(i-1,k,0)=CC(i-1,0,k)+tr2+tr3;
      CH(i  ,k,0)=CC(i  ,0,k)+ti2+ti3;
      T cr2=CC(i-1,0,k)+tr11*tr2+tr12*tr3;
      T ci2=CC(i  ,0,k)+tr11*ti2+tr12*ti3;
      T cr3=CC(i-1,0,k)+tr12*tr2+tr11*tr3;
      T ci3=CC(i  ,0,k)+tr12*ti2+tr11*ti3;
      T ci4, ci5, cr5, cr4;
      MULPM(cr5,cr4,tr5,tr4,ti11,ti12);
      MULPM(ci5,ci4,ti5,ti4,ti11,ti12);
      T dr2, dr3, dr4, dr5, di2, di3, di4, di5;
      PM(dr4,dr3,cr3,ci4);
      PM(di3,di4,ci3,cr4);
      PM(dr5,dr2,cr2,ci5);
      PM(di2,di5,ci2,cr5);
      MULPM(CH(i,k,1),CH(i-1,k,1),WA(0,i-2),WA(0,i-1),di2,dr2);
      MULPM(CH(i,k,2),CH(i-1,k,2),WA(1,i-2),WA(1,i-1),di3,dr3);
      MULPM(CH(i,k,3),CH(i-1,k,3),WA(2,i-2),WA(2,i-1),di4,dr4);
      MULPM(CH(i,k,4),CH(i-1,k,4),WA(3,i-2),WA(3,i-1),di5,dr5);
      }
  }

template<typename T> void radbg(size_t ido, size_t ip, size_t l1,
  T * POCKETFFT_RESTRICT cc, T * POCKETFFT_RESTRICT ch,
  const T0 * POCKETFFT_RESTRICT wa, const T0 * POCKETFFT_RESTRICT csarr)
  {
  const size_t cdim=ip;
  size_t ipph=(ip+1)/ 2;
  size_t idl1 = ido*l1;

  auto CC = [cc,ido,cdim](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+cdim*c)]; };
  auto CH = [ch,ido,l1](size_t a, size_t b, size_t c) -> T&
    { return ch[a+ido*(b+l1*c)]; };
  auto C1 = [cc,ido,l1](size_t a, size_t b, size_t c) -> const T&
    { return cc[a+ido*(b+l1*c)]; };
  auto C2 = [cc,idl1](size_t a, size_t b) -> T&
    { return cc[a+idl1*b]; };
  auto CH2 = [ch,idl1](size_t a, size_t b) -> T&
    { return ch[a+idl1*b]; };

  for (size_t k=0; k<l1; ++k)        // 102
    for (size_t i=0; i<ido; ++i)     // 101
      CH(i,k,0) = CC(i,0,k);
  for (size_t j=1, jc=ip-1; j<ipph; ++j, --jc)   // 108
    {
    size_t j2=2*j-1;
    for (size_t k=0; k<l1; ++k)
      {
      CH(0,k,j ) = 2*CC(ido-1,j2,k);
      CH(0,k,jc) = 2*CC(0,j2+1,k);
      }
    }

  if (ido!=1)
    {
    for (size_t j=1, jc=ip-1; j<ipph; ++j,--jc)   // 111
      {
      size_t j2=2*j-1;
      for (size_t k=0; k<l1; ++k)
        for (size_t i=1, ic=ido-i-2; i<=ido-2; i+=2, ic-=2)      // 109
          {
          CH(i  ,k,j ) = CC(i  ,j2+1,k)+CC(ic  ,j2,k);
          CH(i  ,k,jc) = CC(i  ,j2+1,k)-CC(ic  ,j2,k);
          CH(i+1,k,j ) = CC(i+1,j2+1,k)-CC(ic+1,j2,k);
          CH(i+1,k,jc) = CC(i+1,j2+1,k)+CC(ic+1,j2,k);
          }
      }
    }
  for (size_t l=1,lc=ip-1; l<ipph; ++l,--lc)
    {
    for (size_t ik=0; ik<idl1; ++ik)
      {
      C2(ik,l ) = CH2(ik,0)+csarr[2*l]*CH2(ik,1)+csarr[4*l]*CH2(ik,2);
      C2(ik,lc) = csarr[2*l+1]*CH2(ik,ip-1)+csarr[4*l+1]*CH2(ik,ip-2);
      }
    size_t iang=2*l;
    size_t j=3,jc=ip-3;
    for(; j<ipph-3; j+=4,jc-=4)
      {
      iang+=l; if(iang>ip) iang-=ip;
      T0 ar1=csarr[2*iang], ai1=csarr[2*iang+1];
      iang+=l; if(iang>ip) iang-=ip;
      T0 ar2=csarr[2*iang], ai2=csarr[2*iang+1];
      iang+=l; if(iang>ip) iang-=ip;
      T0 ar3=csarr[2*iang], ai3=csarr[2*iang+1];
      iang+=l; if(iang>ip) iang-=ip;
      T0 ar4=csarr[2*iang], ai4=csarr[2*iang+1];
      for (size_t ik=0; ik<idl1; ++ik)
        {
        C2(ik,l ) += ar1*CH2(ik,j )+ar2*CH2(ik,j +1)
                    +ar3*CH2(ik,j +2)+ar4*CH2(ik,j +3);
        C2(ik,lc) += ai1*CH2(ik,jc)+ai2*CH2(ik,jc-1)
                    +ai3*CH2(ik,jc-2)+ai4*CH2(ik,jc-3);
        }
      }
    for(; j<ipph-1; j+=2,jc-=2)
      {
      iang+=l; if(iang>ip) iang-=ip;
      T0 ar1=csarr[2*iang], ai1=csarr[2*iang+1];
      iang+=l; if(iang>ip) iang-=ip;
      T0 ar2=csarr[2*iang], ai2=csarr[2*iang+1];
      for (size_t ik=0; ik<idl1; ++ik)
        {
        C2(ik,l ) += ar1*CH2(ik,j )+ar2*CH2(ik,j +1);
        C2(ik,lc) += ai1*CH2(ik,jc)+ai2*CH2(ik,jc-1);
        }
      }
    for(; j<ipph; ++j,--jc)
      {
      iang+=l; if(iang>ip) iang-=ip;
      T0 war=csarr[2*iang], wai=csarr[2*iang+1];
      for (size_t ik=0; ik<idl1; ++ik)
        {
        C2(ik,l ) += war*CH2(ik,j );
        C2(ik,lc) += wai*CH2(ik,jc);
        }
      }
    }
  for (size_t j=1; j<ipph; ++j)
    for (size_t ik=0; ik<idl1; ++ik)
      CH2(ik,0) += CH2(ik,j);
  for (size_t j=1, jc=ip-1; j<ipph; ++j,--jc)   // 124
    for (size_t k=0; k<l1; ++k)
      {
      CH(0,k,j ) = C1(0,k,j)-C1(0,k,jc);
      CH(0,k,jc) = C1(0,k,j)+C1(0,k,jc);
      }

  if (ido==1) return;

  for (size_t j=1, jc=ip-1; j<ipph; ++j, --jc)  // 127
    for (size_t k=0; k<l1; ++k)
      for (size_t i=1; i<=ido-2; i+=2)
        {
        CH(i  ,k,j ) = C1(i  ,k,j)-C1(i+1,k,jc);
        CH(i  ,k,jc) = C1(i  ,k,j)+C1(i+1,k,jc);
        CH(i+1,k,j ) = C1(i+1,k,j)+C1(i  ,k,jc);
        CH(i+1,k,jc) = C1(i+1,k,j)-C1(i  ,k,jc);
        }

// All in CH

  for (size_t j=1; j<ip; ++j)
    {
    size_t is = (j-1)*(ido-1);
    for (size_t k=0; k<l1; ++k)
      {
      size_t idij = is;
      for (size_t i=1; i<=ido-2; i+=2)
        {
        T t1=CH(i,k,j), t2=CH(i+1,k,j);
        CH(i  ,k,j) = wa[idij]*t1-wa[idij+1]*t2;
        CH(i+1,k,j) = wa[idij]*t2+wa[idij+1]*t1;
        idij+=2;
        }
      }
    }
  }

    template<typename T> void copy_and_norm(T *c, T *p1, size_t n, T0 fct)
      {
      if (p1!=c)
        {
        if (fct!=1.)
          for (size_t i=0; i<n; ++i)
            c[i] = fct*p1[i];
        else
          memcpy (c,p1,n*sizeof(T));
        }
      else
        if (fct!=1.)
          for (size_t i=0; i<n; ++i)
            c[i] *= fct;
      }

  public:
    template<typename T> void forward(T c[], T0 fct)
      {
      if (length==1) { c[0]*=fct; return; }
      size_t n=length;
      size_t l1=n, nf=fact.size();
      arr<T> ch(n);
      T *p1=c, *p2=ch.data();

      for(size_t k1=0; k1<nf;++k1)
        {
        size_t k=nf-k1-1;
        size_t ip=fact[k].fct;
        size_t ido=n / l1;
        l1 /= ip;
        if(ip==4)
          radf4(ido, l1, p1, p2, fact[k].tw);
        else if(ip==2)
          radf2(ido, l1, p1, p2, fact[k].tw);
        else if(ip==3)
          radf3(ido, l1, p1, p2, fact[k].tw);
        else if(ip==5)
          radf5(ido, l1, p1, p2, fact[k].tw);
        else
          { radfg(ido, ip, l1, p1, p2, fact[k].tw, fact[k].tws); swap (p1,p2); }
        swap (p1,p2);
        }
      copy_and_norm(c,p1,n,fct);
      }

    template<typename T> void backward(T c[], T0 fct)
      {
      if (length==1) { c[0]*=fct; return; }
      size_t n=length;
      size_t l1=1, nf=fact.size();
      arr<T> ch(n);
      T *p1=c, *p2=ch.data();

      for(size_t k=0; k<nf; k++)
        {
        size_t ip = fact[k].fct,
               ido= n/(ip*l1);
        if(ip==4)
          radb4(ido, l1, p1, p2, fact[k].tw);
        else if(ip==2)
          radb2(ido, l1, p1, p2, fact[k].tw);
        else if(ip==3)
          radb3(ido, l1, p1, p2, fact[k].tw);
        else if(ip==5)
          radb5(ido, l1, p1, p2, fact[k].tw);
        else
          radbg(ido, ip, l1, p1, p2, fact[k].tw, fact[k].tws);
        swap (p1,p2);
        l1*=ip;
        }
      copy_and_norm(c,p1,n,fct);
      }

  private:
    void factorize()
      {
      size_t len=length;
      while ((len%4)==0)
        { add_factor(4); len>>=2; }
      if ((len%2)==0)
        {
        len>>=1;
        // factor 2 should be at the front of the factor list
        add_factor(2);
        swap(fact[0].fct, fact.back().fct);
        }
      for (size_t divisor=3; divisor*divisor<=len; divisor+=2)
        while ((len%divisor)==0)
          {
          add_factor(divisor);
          len/=divisor;
          }
      if (len>1) add_factor(len);
      }

    size_t twsize() const
      {
      size_t twsz=0, l1=1;
      for (size_t k=0; k<fact.size(); ++k)
        {
        size_t ip=fact[k].fct, ido=length/(l1*ip);
        twsz+=(ip-1)*(ido-1);
        if (ip>5) twsz+=2*ip;
        l1*=ip;
        }
      return twsz;
      }

    void comp_twiddle()
      {
      sincos_2pibyn<T0> twid(length, true);
      size_t l1=1;
      T0 *ptr=mem.data();
      for (size_t k=0; k<fact.size(); ++k)
        {
        size_t ip=fact[k].fct, ido=length/(l1*ip);
        if (k<fact.size()-1) // last factor doesn't need twiddles
          {
          fact[k].tw=ptr; ptr+=(ip-1)*(ido-1);
          for (size_t j=1; j<ip; ++j)
            for (size_t i=1; i<=(ido-1)/2; ++i)
              {
              fact[k].tw[(j-1)*(ido-1)+2*i-2] = twid[2*j*l1*i];
              fact[k].tw[(j-1)*(ido-1)+2*i-1] = twid[2*j*l1*i+1];
              }
          }
        if (ip>5) // special factors required by *g functions
          {
          fact[k].tws=ptr; ptr+=2*ip;
          fact[k].tws[0] = 1.;
          fact[k].tws[1] = 0.;
          for (size_t i=2, ic=2*ip-2; i<=ic; i+=2, ic-=2)
            {
            fact[k].tws[i  ] = twid[i*(length/ip)];
            fact[k].tws[i+1] = twid[i*(length/ip)+1];
            fact[k].tws[ic]   = twid[i*(length/ip)];
            fact[k].tws[ic+1] = -twid[i*(length/ip)+1];
            }
          }
        l1*=ip;
        }
      }

  public:
    POCKETFFT_NOINLINE rfftp(size_t length_)
      : length(length_)
      {
      if (length==0) throw runtime_error("zero-length FFT requested");
      if (length==1) return;
      factorize();
      mem.resize(twsize());
      comp_twiddle();
      }
};

//
// complex Bluestein transforms
//

template<typename T0> class fftblue
  {
  private:
    size_t n, n2;
    cfftp<T0> plan;
    arr<cmplx<T0>> mem;
    cmplx<T0> *bk, *bkf;

    template<bool fwd, typename T> void fft(cmplx<T> c[], T0 fct)
      {
      arr<cmplx<T>> akf(n2);

      /* initialize a_k and FFT it */
      for (size_t m=0; m<n; ++m)
        akf[m] = c[m].template special_mul<fwd>(bk[m]);
      auto zero = akf[0]*T0(0);
      for (size_t m=n; m<n2; ++m)
        akf[m]=zero;

      plan.forward (akf.data(),1.);

      /* do the convolution */
      for (size_t m=0; m<n2; ++m)
        akf[m] = akf[m].template special_mul<!fwd>(bkf[m]);

      /* inverse FFT */
      plan.backward (akf.data(),1.);

      /* multiply by b_k */
      for (size_t m=0; m<n; ++m)
        c[m] = akf[m].template special_mul<fwd>(bk[m])*fct;
      }

  public:
    POCKETFFT_NOINLINE fftblue(size_t length)
      : n(length), n2(util::good_size_cmplx(n*2-1)), plan(n2), mem(n+n2),
        bk(mem.data()), bkf(mem.data()+n)
      {
      /* initialize b_k */
      sincos_2pibyn<T0> tmp_(2*n, false);
      auto tmp = tmp_.cdata();
      bk[0].Set(1, 0);

      size_t coeff=0;
      for (size_t m=1; m<n; ++m)
        {
        coeff+=2*m-1;
        if (coeff>=2*n) coeff-=2*n;
        bk[m] = tmp[coeff];
        }

      /* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
      T0 xn2 = T0(1)/T0(n2);
      bkf[0] = bk[0]*xn2;
      for (size_t m=1; m<n; ++m)
        bkf[m] = bkf[n2-m] = bk[m]*xn2;
      for (size_t m=n;m<=(n2-n);++m)
        bkf[m].Set(0.,0.);
      plan.forward(bkf,1.);
      }

    template<typename T> void backward(cmplx<T> c[], T0 fct)
      { fft<false>(c,fct); }

    template<typename T> void forward(cmplx<T> c[], T0 fct)
      { fft<true>(c,fct); }

    template<typename T> void backward_r(T c[], T0 fct)
      {
      arr<cmplx<T>> tmp(n);
      tmp[0].Set(c[0],c[0]*0);
      memcpy (reinterpret_cast<void *>(tmp.data()+1),
              reinterpret_cast<void *>(c+1), (n-1)*sizeof(T));
      if ((n&1)==0) tmp[n/2].i=T0(0)*c[0];
      for (size_t m=1; 2*m<n; ++m)
        tmp[n-m].Set(tmp[m].r, -tmp[m].i);
      fft<false>(tmp.data(),fct);
      for (size_t m=0; m<n; ++m)
        c[m] = tmp[m].r;
      }

    template<typename T> void forward_r(T c[], T0 fct)
      {
      arr<cmplx<T>> tmp(n);
      auto zero = T0(0)*c[0];
      for (size_t m=0; m<n; ++m)
        tmp[m].Set(c[m], zero);
      fft<true>(tmp.data(),fct);
      c[0] = tmp[0].r;
      memcpy (c+1, tmp.data()+1, (n-1)*sizeof(T));
      }
    };

//
// flexible (FFTPACK/Bluestein) complex 1D transform
//

template<typename T0> class pocketfft_c
  {
  private:
    unique_ptr<cfftp<T0>> packplan;
    unique_ptr<fftblue<T0>> blueplan;
    size_t len;

  public:
    POCKETFFT_NOINLINE pocketfft_c(size_t length)
      : len(length)
      {
      if (length==0) throw runtime_error("zero-length FFT requested");
      size_t tmp = (length<50) ? 0 : util::largest_prime_factor(length);
      if (tmp*tmp <= length)
        {
        packplan=unique_ptr<cfftp<T0>>(new cfftp<T0>(length));
        return;
        }
      double comp1 = util::cost_guess(length);
      double comp2 = 2*util::cost_guess(util::good_size_cmplx(2*length-1));
      comp2*=1.5; /* fudge factor that appears to give good overall performance */
      if (comp2<comp1) // use Bluestein
        blueplan=unique_ptr<fftblue<T0>>(new fftblue<T0>(length));
      else
        packplan=unique_ptr<cfftp<T0>>(new cfftp<T0>(length));
      }

    template<typename T> POCKETFFT_NOINLINE void backward(cmplx<T> c[], T0 fct) const
      { packplan ? packplan->backward(c,fct) : blueplan->backward(c,fct); }

    template<typename T> POCKETFFT_NOINLINE void forward(cmplx<T> c[], T0 fct) const
      { packplan ? packplan->forward(c,fct) : blueplan->forward(c,fct); }

    size_t length() const { return len; }
  };

//
// flexible (FFTPACK/Bluestein) real-valued 1D transform
//

template<typename T0> class pocketfft_r
  {
  private:
    unique_ptr<rfftp<T0>> packplan;
    unique_ptr<fftblue<T0>> blueplan;
    size_t len;

  public:
    POCKETFFT_NOINLINE pocketfft_r(size_t length)
      : len(length)
      {
      if (length==0) throw runtime_error("zero-length FFT requested");
      size_t tmp = (length<50) ? 0 : util::largest_prime_factor(length);
      if (tmp*tmp <= length)
        {
        packplan=unique_ptr<rfftp<T0>>(new rfftp<T0>(length));
        return;
        }
      double comp1 = 0.5*util::cost_guess(length);
      double comp2 = 2*util::cost_guess(util::good_size_cmplx(2*length-1));
      comp2*=1.5; /* fudge factor that appears to give good overall performance */
      if (comp2<comp1) // use Bluestein
        blueplan=unique_ptr<fftblue<T0>>(new fftblue<T0>(length));
      else
        packplan=unique_ptr<rfftp<T0>>(new rfftp<T0>(length));
      }

    template<typename T> POCKETFFT_NOINLINE void backward(T c[], T0 fct) const
      {
      packplan ? packplan->backward(c,fct)
               : blueplan->backward_r(c,fct);
      }

    template<typename T> POCKETFFT_NOINLINE void forward(T c[], T0 fct) const
      {
      packplan ? packplan->forward(c,fct)
               : blueplan->forward_r(c,fct);
      }

    size_t length() const { return len; }
  };


//
// sine/cosine transforms
//

template<typename T0> class T_dct1
  {
  private:
    pocketfft_r<T0> fftplan;

  public:
    POCKETFFT_NOINLINE T_dct1(size_t length)
      : fftplan(2*(length-1)) {}

    template<typename T> POCKETFFT_NOINLINE void exec(T c[], T0 fct, bool ortho,
      int /*type*/, bool /*cosine*/) const
      {
      constexpr T0 sqrt2=T0(1.414213562373095048801688724209698L);
      size_t N=fftplan.length(), n=N/2+1;
      if (ortho)
        { c[0]*=sqrt2; c[n-1]*=sqrt2; }
      arr<T> tmp(N);
      tmp[0] = c[0];
      for (size_t i=1; i<n; ++i)
        tmp[i] = tmp[N-i] = c[i];
      fftplan.forward(tmp.data(), fct);
      c[0] = tmp[0];
      for (size_t i=1; i<n; ++i)
        c[i] = tmp[2*i-1];
      if (ortho)
        { c[0]/=sqrt2; c[n-1]/=sqrt2; }
      }

    size_t length() const { return fftplan.length()/2+1; }
  };

template<typename T0> class T_dst1
  {
  private:
    pocketfft_r<T0> fftplan;

  public:
    POCKETFFT_NOINLINE T_dst1(size_t length)
      : fftplan(2*(length+1)) {}

    template<typename T> POCKETFFT_NOINLINE void exec(T c[], T0 fct,
      bool /*ortho*/, int /*type*/, bool /*cosine*/) const
      {
      size_t N=fftplan.length(), n=N/2-1;
      arr<T> tmp(N);
      tmp[0] = tmp[n+1] = c[0]*0;
      for (size_t i=0; i<n; ++i)
        { tmp[i+1]=c[i]; tmp[N-1-i]=-c[i]; }
      fftplan.forward(tmp.data(), fct);
      for (size_t i=0; i<n; ++i)
        c[i] = -tmp[2*i+2];
      }

    size_t length() const { return fftplan.length()/2-1; }
  };

template<typename T0> class T_dcst23
  {
  private:
    pocketfft_r<T0> fftplan;
    vector<T0> twiddle;

  public:
    POCKETFFT_NOINLINE T_dcst23(size_t length)
      : fftplan(length), twiddle(length)
      {
      constexpr T0 pi = T0(3.141592653589793238462643383279502884197L);
      for (size_t i=0; i<length; ++i)
        twiddle[i] = T0(cos(0.5*pi*T0(i+1)/T0(length)));
      }

    template<typename T> POCKETFFT_NOINLINE void exec(T c[], T0 fct, bool ortho,
      int type, bool cosine) const
      {
      constexpr T0 sqrt2=T0(1.414213562373095048801688724209698L);
      size_t N=length();
      size_t NS2 = (N+1)/2;
      if (type==2)
        {
        if (!cosine)
          for (size_t k=1; k<N; k+=2)
            c[k] = -c[k];
        c[0] *= 2;
        if ((N&1)==0) c[N-1]*=2;
        for (size_t k=1; k<N-1; k+=2)
          MPINPLACE(c[k+1], c[k]);
        fftplan.backward(c, fct);
        for (size_t k=1, kc=N-1; k<NS2; ++k, --kc)
          {
          T t1 = twiddle[k-1]*c[kc]+twiddle[kc-1]*c[k];
          T t2 = twiddle[k-1]*c[k]-twiddle[kc-1]*c[kc];
          c[k] = T0(0.5)*(t1+t2); c[kc]=T0(0.5)*(t1-t2);
          }
        if ((N&1)==0)
          c[NS2] *= twiddle[NS2-1];
        if (!cosine)
          for (size_t k=0, kc=N-1; k<kc; ++k, --kc)
            swap(c[k], c[kc]);
        if (ortho) c[0]/=sqrt2;
        }
      else
        {
        if (ortho) c[0]*=sqrt2;
        if (!cosine)
          for (size_t k=0, kc=N-1; k<NS2; ++k, --kc)
            swap(c[k], c[kc]);
        for (size_t k=1, kc=N-1; k<NS2; ++k, --kc)
          {
          T t1=c[k]+c[kc], t2=c[k]-c[kc];
          c[k] = twiddle[k-1]*t2+twiddle[kc-1]*t1;
          c[kc]= twiddle[k-1]*t1-twiddle[kc-1]*t2;
          }
        if ((N&1)==0)
          c[NS2] *= 2*twiddle[NS2-1];
        fftplan.forward(c, fct);
        for (size_t k=1; k<N-1; k+=2)
          MPINPLACE(c[k], c[k+1]);
        if (!cosine)
          for (size_t k=1; k<N; k+=2)
            c[k] = -c[k];
        }
      }

    size_t length() const { return fftplan.length(); }
  };

template<typename T0> class T_dcst4
  {
  // even length algorithm from
  // https://www.appletonaudio.com/blog/2013/derivation-of-fast-dct-4-algorithm-based-on-dft/
  private:
    size_t N;
    unique_ptr<pocketfft_c<T0>> fft;
    unique_ptr<pocketfft_r<T0>> rfft;
    arr<cmplx<T0>> C2;

  public:
    POCKETFFT_NOINLINE T_dcst4(size_t length)
      : N(length),
        fft((N&1) ? nullptr : new pocketfft_c<T0>(N/2)),
        rfft((N&1)? new pocketfft_r<T0>(N) : nullptr),
        C2((N&1) ? 0 : N/2)
      {
      constexpr T0 pi = T0(3.141592653589793238462643383279502884197L);
      if ((N&1)==0)
        for (size_t i=0; i<N/2; ++i)
          {
          T0 ang = -pi/T0(N)*(T0(i)+T0(0.125));
          C2[i].Set(cos(ang), sin(ang));
          }
      }

    template<typename T> POCKETFFT_NOINLINE void exec(T c[], T0 fct,
      bool /*ortho*/, int /*type*/, bool cosine) const
      {
      constexpr T0 sqrt2=T0(1.414213562373095048801688724209698L);
      size_t n2 = N/2;
      if (!cosine)
        for (size_t k=0, kc=N-1; k<n2; ++k, --kc)
          swap(c[k], c[kc]);
      if (N&1)
        {
        // The following code is derived from the FFTW3 function apply_re11()
        // and is released under the 3-clause BSD license with friendly
        // permission of Matteo Frigo.

        auto SGN_SET = [](T x, size_t i) {return (i%2) ? -x : x;};
        arr<T> y(N);
        {
        size_t i=0, m=n2;
        for (; m<N; ++i, m+=4)
          y[i] = c[m];
        for (; m<2*N; ++i, m+=4)
          y[i] = -c[2*N-m-1];
        for (; m<3*N; ++i, m+=4)
          y[i] = -c[m-2*N];
        for (; m<4*N; ++i, m+=4)
          y[i] = c[4*N-m-1];
        m -= 4*N;
        for (; i<N; ++i, m+=4)
          y[i] = c[m];
        }
        rfft->forward(y.data(), fct);
        {
        size_t i=0;
        for (; i+i+1<n2; ++i)
          {
          size_t k = i+i+1;
          T c1=y[2*k-1], s1=y[2*k], c2=y[2*k+1], s2=y[2*k+2];
          c[i] = sqrt2 * (SGN_SET(c1, (i+1)/2) + SGN_SET(s1, i/2));
          c[N-(i+1)] = sqrt2 * (SGN_SET(c1, (N-i)/2) - SGN_SET(s1, (N-(i+1))/2));
          c[n2-(i+1)] = sqrt2 * (SGN_SET(c2, (n2-i)/2) - SGN_SET(s2, (n2-(i+1))/2));
          c[n2+(i+1)] = sqrt2 * (SGN_SET(c2, (n2+i+2)/2) + SGN_SET(s2, (n2+(i+1))/2));
          }
        if (i+i+1 == n2)
          {
          T cx=y[2*n2-1], sx=y[2*n2];
          c[i] = sqrt2 * (SGN_SET(cx, (i+1)/2) + SGN_SET(sx, i/2));
          c[N-(i+1)] = sqrt2 * (SGN_SET(cx, (i+2)/2) + SGN_SET(sx, (i+1)/2));
          }
        }
        c[n2] = sqrt2 * SGN_SET(y[0], (n2+1)/2);

        // FFTW-derived code ends here
        }
      else
        {
        // even length algorithm from
        // https://www.appletonaudio.com/blog/2013/derivation-of-fast-dct-4-algorithm-based-on-dft/
        arr<cmplx<T>> y(n2);
        for(size_t i=0; i<n2; ++i)
          {
          y[i].Set(c[2*i],c[N-1-2*i]);
          y[i] *= C2[i];
          }
        fft->forward(y.data(), fct);
        for(size_t i=0, ic=n2-1; i<n2; ++i, --ic)
          {
          c[2*i  ] =  2*(y[i ].r*C2[i ].r-y[i ].i*C2[i ].i);
          c[2*i+1] = -2*(y[ic].i*C2[ic].r+y[ic].r*C2[ic].i);
          }
        }
      if (!cosine)
        for (size_t k=1; k<N; k+=2)
          c[k] = -c[k];
      }

    size_t length() const { return N; }
  };


//
// multi-D infrastructure
//

template<typename T> shared_ptr<T> get_plan(size_t length)
  {
#if POCKETFFT_CACHE_SIZE==0
  return make_shared<T>(length);
#else
  constexpr size_t nmax=POCKETFFT_CACHE_SIZE;
  static array<shared_ptr<T>, nmax> cache;
  static array<size_t, nmax> last_access{{0}};
  static size_t access_counter = 0;
  static mutex mut;

  auto find_in_cache = [&]() -> shared_ptr<T>
    {
    for (size_t i=0; i<nmax; ++i)
      if (cache[i] && (cache[i]->length()==length))
        {
        // no need to update if this is already the most recent entry
        if (last_access[i]!=access_counter)
          {
          last_access[i] = ++access_counter;
          // Guard against overflow
          if (access_counter == 0)
            last_access.fill(0);
          }
        return cache[i];
        }

    return nullptr;
    };

  {
  lock_guard<mutex> lock(mut);
  auto p = find_in_cache();
  if (p) return p;
  }
  auto plan = make_shared<T>(length);
  {
  lock_guard<mutex> lock(mut);
  auto p = find_in_cache();
  if (p) return p;

  size_t lru = 0;
  for (size_t i=1; i<nmax; ++i)
    if (last_access[i] < last_access[lru])
      lru = i;

  cache[lru] = plan;
  last_access[lru] = ++access_counter;
  }
  return plan;
#endif
  }

class arr_info
  {
  protected:
    shape_t shp;
    stride_t str;

  public:
    arr_info(const shape_t &shape_, const stride_t &stride_)
      : shp(shape_), str(stride_) {}
    size_t ndim() const { return shp.size(); }
    size_t size() const { return util::prod(shp); }
    const shape_t &shape() const { return shp; }
    size_t shape(size_t i) const { return shp[i]; }
    const stride_t &stride() const { return str; }
    const ptrdiff_t &stride(size_t i) const { return str[i]; }
  };

template<typename T> class cndarr: public arr_info
  {
  protected:
    const char *d;

  public:
    cndarr(const void *data_, const shape_t &shape_, const stride_t &stride_)
      : arr_info(shape_, stride_),
        d(reinterpret_cast<const char *>(data_)) {}
    const T &operator[](ptrdiff_t ofs) const
      { return *reinterpret_cast<const T *>(d+ofs); }
  };

template<typename T> class ndarr: public cndarr<T>
  {
  public:
    ndarr(void *data_, const shape_t &shape_, const stride_t &stride_)
      : cndarr<T>::cndarr(const_cast<const void *>(data_), shape_, stride_)
      {}
    T &operator[](ptrdiff_t ofs)
      { return *reinterpret_cast<T *>(const_cast<char *>(cndarr<T>::d+ofs)); }
  };

template<size_t N> class multi_iter
  {
  private:
    shape_t pos;
    const arr_info &iarr, &oarr;
    ptrdiff_t p_ii, p_i[N], str_i, p_oi, p_o[N], str_o;
    size_t idim, rem;

    void advance_i()
      {
      for (int i_=int(pos.size())-1; i_>=0; --i_)
        {
        auto i = size_t(i_);
        if (i==idim) continue;
        p_ii += iarr.stride(i);
        p_oi += oarr.stride(i);
        if (++pos[i] < iarr.shape(i))
          return;
        pos[i] = 0;
        p_ii -= ptrdiff_t(iarr.shape(i))*iarr.stride(i);
        p_oi -= ptrdiff_t(oarr.shape(i))*oarr.stride(i);
        }
      }

  public:
    multi_iter(const arr_info &iarr_, const arr_info &oarr_, size_t idim_)
      : pos(iarr_.ndim(), 0), iarr(iarr_), oarr(oarr_), p_ii(0),
        str_i(iarr.stride(idim_)), p_oi(0), str_o(oarr.stride(idim_)),
        idim(idim_), rem(iarr.size()/iarr.shape(idim))
      {
      auto nshares = util::nthreads();
      if (nshares==1) return;
      if (nshares==0) throw runtime_error("can't run with zero threads");
      auto myshare = util::thread_num();
      if (myshare>=nshares) throw runtime_error("impossible share requested");
      size_t nbase = rem/nshares;
      size_t additional = rem%nshares;
      size_t lo = myshare*nbase + ((myshare<additional) ? myshare : additional);
      size_t hi = lo+nbase+(myshare<additional);
      size_t todo = hi-lo;

      size_t chunk = rem;
      for (size_t i=0; i<pos.size(); ++i)
        {
        if (i==idim) continue;
        chunk /= iarr.shape(i);
        size_t n_advance = lo/chunk;
        pos[i] += n_advance;
        p_ii += ptrdiff_t(n_advance)*iarr.stride(i);
        p_oi += ptrdiff_t(n_advance)*oarr.stride(i);
        lo -= n_advance*chunk;
        }
      rem = todo;
      }
    void advance(size_t n)
      {
      if (rem<n) throw runtime_error("underrun");
      for (size_t i=0; i<n; ++i)
        {
        p_i[i] = p_ii;
        p_o[i] = p_oi;
        advance_i();
        }
      rem -= n;
      }
    ptrdiff_t iofs(size_t i) const { return p_i[0] + ptrdiff_t(i)*str_i; }
    ptrdiff_t iofs(size_t j, size_t i) const { return p_i[j] + ptrdiff_t(i)*str_i; }
    ptrdiff_t oofs(size_t i) const { return p_o[0] + ptrdiff_t(i)*str_o; }
    ptrdiff_t oofs(size_t j, size_t i) const { return p_o[j] + ptrdiff_t(i)*str_o; }
    size_t length_in() const { return iarr.shape(idim); }
    size_t length_out() const { return oarr.shape(idim); }
    ptrdiff_t stride_in() const { return str_i; }
    ptrdiff_t stride_out() const { return str_o; }
    size_t remaining() const { return rem; }
  };

class simple_iter
  {
  private:
    shape_t pos;
    const arr_info &arr;
    ptrdiff_t p;
    size_t rem;

  public:
    simple_iter(const arr_info &arr_)
      : pos(arr_.ndim(), 0), arr(arr_), p(0), rem(arr_.size()) {}
    void advance()
      {
      --rem;
      for (int i_=int(pos.size())-1; i_>=0; --i_)
        {
        auto i = size_t(i_);
        p += arr.stride(i);
        if (++pos[i] < arr.shape(i))
          return;
        pos[i] = 0;
        p -= ptrdiff_t(arr.shape(i))*arr.stride(i);
        }
      }
    ptrdiff_t ofs() const { return p; }
    size_t remaining() const { return rem; }
  };

class rev_iter
  {
  private:
    shape_t pos;
    const arr_info &arr;
    vector<char> rev_axis;
    vector<char> rev_jump;
    size_t last_axis, last_size;
    shape_t shp;
    ptrdiff_t p, rp;
    size_t rem;

  public:
    rev_iter(const arr_info &arr_, const shape_t &axes)
      : pos(arr_.ndim(), 0), arr(arr_), rev_axis(arr_.ndim(), 0),
        rev_jump(arr_.ndim(), 1), p(0), rp(0)
      {
      for (auto ax: axes)
        rev_axis[ax]=1;
      last_axis = axes.back();
      last_size = arr.shape(last_axis)/2 + 1;
      shp = arr.shape();
      shp[last_axis] = last_size;
      rem=1;
      for (auto i: shp)
        rem *= i;
      }
    void advance()
      {
      --rem;
      for (int i_=int(pos.size())-1; i_>=0; --i_)
        {
        auto i = size_t(i_);
        p += arr.stride(i);
        if (!rev_axis[i])
          rp += arr.stride(i);
        else
          {
          rp -= arr.stride(i);
          if (rev_jump[i])
            {
            rp += ptrdiff_t(arr.shape(i))*arr.stride(i);
            rev_jump[i] = 0;
            }
          }
        if (++pos[i] < shp[i])
          return;
        pos[i] = 0;
        p -= ptrdiff_t(shp[i])*arr.stride(i);
        if (rev_axis[i])
          {
          rp -= ptrdiff_t(arr.shape(i)-shp[i])*arr.stride(i);
          rev_jump[i] = 1;
          }
        else
          rp -= ptrdiff_t(shp[i])*arr.stride(i);
        }
      }
    ptrdiff_t ofs() const { return p; }
    ptrdiff_t rev_ofs() const { return rp; }
    size_t remaining() const { return rem; }
  };

#ifndef POCKETFFT_NO_VECTORS
template<typename T> struct VTYPE {};
template<> struct VTYPE<float>
  {
  using type = float __attribute__ ((vector_size (VLEN<float>::val*sizeof(float))));
  };
template<> struct VTYPE<double>
  {
  using type = double __attribute__ ((vector_size (VLEN<double>::val*sizeof(double))));
  };
template<> struct VTYPE<long double>
  {
  using type = long double __attribute__ ((vector_size (VLEN<long double>::val*sizeof(long double))));
  };
#endif

template<typename T> arr<char> alloc_tmp(const shape_t &shape,
  size_t axsize, size_t elemsize)
  {
  auto othersize = util::prod(shape)/axsize;
  auto tmpsize = axsize*((othersize>=VLEN<T>::val) ? VLEN<T>::val : 1);
  return arr<char>(tmpsize*elemsize);
  }
template<typename T> arr<char> alloc_tmp(const shape_t &shape,
  const shape_t &axes, size_t elemsize)
  {
  size_t fullsize=util::prod(shape);
  size_t tmpsize=0;
  for (size_t i=0; i<axes.size(); ++i)
    {
    auto axsize = shape[axes[i]];
    auto othersize = fullsize/axsize;
    auto sz = axsize*((othersize>=VLEN<T>::val) ? VLEN<T>::val : 1);
    if (sz>tmpsize) tmpsize=sz;
    }
  return arr<char>(tmpsize*elemsize);
  }

#ifdef POCKETFFT_OPENMP
#define POCKETFFT_NTHREADS nthreads
#else
#define POCKETFFT_NTHREADS
#endif

template<typename T> POCKETFFT_NOINLINE void general_c(
  const cndarr<cmplx<T>> &in, ndarr<cmplx<T>> &out,
  const shape_t &axes, bool forward, T fct, size_t POCKETFFT_NTHREADS)
  {
  shared_ptr<pocketfft_c<T>> plan;

  for (size_t iax=0; iax<axes.size(); ++iax)
    {
    constexpr auto vlen = VLEN<T>::val;
    size_t len=in.shape(axes[iax]);
    if ((!plan) || (len!=plan->length()))
      plan = get_plan<pocketfft_c<T>>(len);

#ifdef POCKETFFT_OPENMP
#pragma omp parallel num_threads(util::thread_count(nthreads, in.shape(), axes[iax]))
#endif
{
    auto storage = alloc_tmp<T>(in.shape(), len, sizeof(cmplx<T>));
    const auto &tin(iax==0? in : out);
    multi_iter<vlen> it(tin, out, axes[iax]);
#ifndef POCKETFFT_NO_VECTORS
    if (vlen>1)
      while (it.remaining()>=vlen)
        {
        using vtype = typename VTYPE<T>::type;
        it.advance(vlen);
        auto tdatav = reinterpret_cast<cmplx<vtype> *>(storage.data());
        for (size_t i=0; i<len; ++i)
          for (size_t j=0; j<vlen; ++j)
            {
            tdatav[i].r[j] = tin[it.iofs(j,i)].r;
            tdatav[i].i[j] = tin[it.iofs(j,i)].i;
            }
        forward ? plan->forward (tdatav, fct) : plan->backward(tdatav, fct);
        for (size_t i=0; i<len; ++i)
          for (size_t j=0; j<vlen; ++j)
            out[it.oofs(j,i)].Set(tdatav[i].r[j],tdatav[i].i[j]);
        }
#endif
    while (it.remaining()>0)
      {
      it.advance(1);
      auto buf = it.stride_out() == sizeof(cmplx<T>) ?
        &out[it.oofs(0)] : reinterpret_cast<cmplx<T> *>(storage.data());

      if (buf != &tin[it.iofs(0)])
        for (size_t i=0; i<len; ++i)
          buf[i] = tin[it.iofs(i)];
      forward ? plan->forward (buf, fct) : plan->backward(buf, fct);
      if (buf != &out[it.oofs(0)])
        for (size_t i=0; i<len; ++i)
          out[it.oofs(i)] = buf[i];
      }
} // end of parallel region
    fct = T(1); // factor has been applied, use 1 for remaining axes
    }
  }

template<typename T> POCKETFFT_NOINLINE void general_hartley(
  const cndarr<T> &in, ndarr<T> &out, const shape_t &axes, T fct,
  size_t POCKETFFT_NTHREADS)
  {
  shared_ptr<pocketfft_r<T>> plan;

  for (size_t iax=0; iax<axes.size(); ++iax)
    {
    constexpr auto vlen = VLEN<T>::val;
    size_t len=in.shape(axes[iax]);
    if ((!plan) || (len!=plan->length()))
      plan = get_plan<pocketfft_r<T>>(len);

#ifdef POCKETFFT_OPENMP
#pragma omp parallel num_threads(util::thread_count(nthreads, in.shape(), axes[iax]))
#endif
{
    auto storage = alloc_tmp<T>(in.shape(), len, sizeof(T));
    const auto &tin(iax==0 ? in : out);
    multi_iter<vlen> it(tin, out, axes[iax]);
#ifndef POCKETFFT_NO_VECTORS
    if (vlen>1)
      while (it.remaining()>=vlen)
        {
        using vtype = typename VTYPE<T>::type;
        it.advance(vlen);
        auto tdatav = reinterpret_cast<vtype *>(storage.data());
        for (size_t i=0; i<len; ++i)
          for (size_t j=0; j<vlen; ++j)
            tdatav[i][j] = tin[it.iofs(j,i)];
        plan->forward(tdatav, fct);
        for (size_t j=0; j<vlen; ++j)
          out[it.oofs(j,0)] = tdatav[0][j];
        size_t i=1, i1=1, i2=len-1;
        for (i=1; i<len-1; i+=2, ++i1, --i2)
          for (size_t j=0; j<vlen; ++j)
            {
            out[it.oofs(j,i1)] = tdatav[i][j]+tdatav[i+1][j];
            out[it.oofs(j,i2)] = tdatav[i][j]-tdatav[i+1][j];
            }
        if (i<len)
          for (size_t j=0; j<vlen; ++j)
            out[it.oofs(j,i1)] = tdatav[i][j];
        }
#endif
    while (it.remaining()>0)
      {
      it.advance(1);
      auto tdata = reinterpret_cast<T *>(storage.data());
      for (size_t i=0; i<len; ++i)
        tdata[i] = tin[it.iofs(i)];
      plan->forward(tdata, fct);
      // Hartley order
      out[it.oofs(0)] = tdata[0];
      size_t i=1, i1=1, i2=len-1;
      for (i=1; i<len-1; i+=2, ++i1, --i2)
        {
        out[it.oofs(i1)] = tdata[i]+tdata[i+1];
        out[it.oofs(i2)] = tdata[i]-tdata[i+1];
        }
      if (i<len)
        out[it.oofs(i1)] = tdata[i];
      }
} // end of parallel region
    fct = T(1); // factor has been applied, use 1 for remaining axes
    }
  }

template<typename Trafo, typename T> POCKETFFT_NOINLINE void general_dcst(
  const cndarr<T> &in, ndarr<T> &out, const shape_t &axes,
  T fct, bool ortho, int type, bool cosine, size_t POCKETFFT_NTHREADS)
  {
  shared_ptr<Trafo> plan;

  for (size_t iax=0; iax<axes.size(); ++iax)
    {
    constexpr auto vlen = VLEN<T>::val;
    size_t len=in.shape(axes[iax]);
    if ((!plan) || (len!=plan->length()))
      plan = get_plan<Trafo>(len);

#ifdef POCKETFFT_OPENMP
#pragma omp parallel num_threads(util::thread_count(nthreads, in.shape(), axes[iax]))
#endif
{
    auto storage = alloc_tmp<T>(in.shape(), len, sizeof(T));
    const auto &tin(iax==0 ? in : out);
    multi_iter<vlen> it(tin, out, axes[iax]);
#ifndef POCKETFFT_NO_VECTORS
    if (vlen>1)
      while (it.remaining()>=vlen)
        {
        using vtype = typename VTYPE<T>::type;
        it.advance(vlen);
        auto tdatav = reinterpret_cast<vtype *>(storage.data());
        for (size_t i=0; i<len; ++i)
          for (size_t j=0; j<vlen; ++j)
            tdatav[i][j] = tin[it.iofs(j,i)];
        plan->exec(tdatav, fct, ortho, type, cosine);
        for (size_t i=0; i<len; ++i)
          for (size_t j=0; j<vlen; ++j)
            out[it.oofs(j,i)] = tdatav[i][j];
        }
#endif
    while (it.remaining()>0)
      {
      it.advance(1);
      auto buf = it.stride_out() == sizeof(T) ? &out[it.oofs(0)]
        : reinterpret_cast<T *>(storage.data());

      if (buf != &tin[it.iofs(0)])
        for (size_t i=0; i<len; ++i)
          buf[i] = tin[it.iofs(i)];
      plan->exec(buf, fct, ortho, type, cosine);
      if (buf != &out[it.oofs(0)])
        for (size_t i=0; i<len; ++i)
          out[it.oofs(i)] = buf[i];
      }
} // end of parallel region
    fct = T(1); // factor has been applied, use 1 for remaining axes
    }
  }

template<typename T> POCKETFFT_NOINLINE void general_r2c(
  const cndarr<T> &in, ndarr<cmplx<T>> &out, size_t axis, bool forward, T fct,
  size_t POCKETFFT_NTHREADS)
  {
  auto plan = get_plan<pocketfft_r<T>>(in.shape(axis));
  constexpr auto vlen = VLEN<T>::val;
  size_t len=in.shape(axis);
#ifdef POCKETFFT_OPENMP
#pragma omp parallel num_threads(util::thread_count(nthreads, in.shape(), axis))
#endif
{
  auto storage = alloc_tmp<T>(in.shape(), len, sizeof(T));
  multi_iter<vlen> it(in, out, axis);
#ifndef POCKETFFT_NO_VECTORS
  if (vlen>1)
    while (it.remaining()>=vlen)
      {
      using vtype = typename VTYPE<T>::type;
      it.advance(vlen);
      auto tdatav = reinterpret_cast<vtype *>(storage.data());
      for (size_t i=0; i<len; ++i)
        for (size_t j=0; j<vlen; ++j)
          tdatav[i][j] = in[it.iofs(j,i)];
      plan->forward(tdatav, fct);
      for (size_t j=0; j<vlen; ++j)
        out[it.oofs(j,0)].Set(tdatav[0][j]);
      size_t i=1, ii=1;
      if (forward)
        for (; i<len-1; i+=2, ++ii)
          for (size_t j=0; j<vlen; ++j)
            out[it.oofs(j,ii)].Set(tdatav[i][j], tdatav[i+1][j]);
      else
        for (; i<len-1; i+=2, ++ii)
          for (size_t j=0; j<vlen; ++j)
            out[it.oofs(j,ii)].Set(tdatav[i][j], -tdatav[i+1][j]);
      if (i<len)
        for (size_t j=0; j<vlen; ++j)
          out[it.oofs(j,ii)].Set(tdatav[i][j]);
      }
#endif
  while (it.remaining()>0)
    {
    it.advance(1);
    auto tdata = reinterpret_cast<T *>(storage.data());
    for (size_t i=0; i<len; ++i)
      tdata[i] = in[it.iofs(i)];
    plan->forward(tdata, fct);
    out[it.oofs(0)].Set(tdata[0]);
    size_t i=1, ii=1;
    if (forward)
      for (; i<len-1; i+=2, ++ii)
        out[it.oofs(ii)].Set(tdata[i], tdata[i+1]);
    else
      for (; i<len-1; i+=2, ++ii)
        out[it.oofs(ii)].Set(tdata[i], -tdata[i+1]);
    if (i<len)
      out[it.oofs(ii)].Set(tdata[i]);
    }
} // end of parallel region
  }
template<typename T> POCKETFFT_NOINLINE void general_c2r(
  const cndarr<cmplx<T>> &in, ndarr<T> &out, size_t axis, bool forward, T fct,
  size_t POCKETFFT_NTHREADS)
  {
  auto plan = get_plan<pocketfft_r<T>>(out.shape(axis));
  constexpr auto vlen = VLEN<T>::val;
  size_t len=out.shape(axis);
#ifdef POCKETFFT_OPENMP
#pragma omp parallel num_threads(util::thread_count(nthreads, in.shape(), axis))
#endif
{
  auto storage = alloc_tmp<T>(out.shape(), len, sizeof(T));
  multi_iter<vlen> it(in, out, axis);
#ifndef POCKETFFT_NO_VECTORS
  if (vlen>1)
    while (it.remaining()>=vlen)
      {
      using vtype = typename VTYPE<T>::type;
      it.advance(vlen);
      auto tdatav = reinterpret_cast<vtype *>(storage.data());
      for (size_t j=0; j<vlen; ++j)
        tdatav[0][j]=in[it.iofs(j,0)].r;
      {
      size_t i=1, ii=1;
      if (forward)
        for (; i<len-1; i+=2, ++ii)
          for (size_t j=0; j<vlen; ++j)
            in[it.iofs(j,ii)].SplitConj(tdatav[i][j], tdatav[i+1][j]);
      else
        for (; i<len-1; i+=2, ++ii)
          for (size_t j=0; j<vlen; ++j)
            in[it.iofs(j,ii)].Split(tdatav[i][j], tdatav[i+1][j]);
      if (i<len)
        for (size_t j=0; j<vlen; ++j)
          tdatav[i][j] = in[it.iofs(j,ii)].r;
      }
      plan->backward(tdatav, fct);
      for (size_t i=0; i<len; ++i)
        for (size_t j=0; j<vlen; ++j)
          out[it.oofs(j,i)] = tdatav[i][j];
      }
#endif
  while (it.remaining()>0)
    {
    it.advance(1);
    auto tdata = reinterpret_cast<T *>(storage.data());
    tdata[0]=in[it.iofs(0)].r;
    {
    size_t i=1, ii=1;
    if (forward)
      for (; i<len-1; i+=2, ++ii)
        in[it.iofs(ii)].SplitConj(tdata[i], tdata[i+1]);
    else
      for (; i<len-1; i+=2, ++ii)
        in[it.iofs(ii)].Split(tdata[i], tdata[i+1]);
    if (i<len)
      tdata[i] = in[it.iofs(ii)].r;
    }
    plan->backward(tdata, fct);
    for (size_t i=0; i<len; ++i)
      out[it.oofs(i)] = tdata[i];
    }
} // end of parallel region
  }

template<typename T> POCKETFFT_NOINLINE void general_r(
  const cndarr<T> &in, ndarr<T> &out, const shape_t &axes, bool r2c,
  bool forward, T fct, size_t POCKETFFT_NTHREADS)
  {
  shared_ptr<pocketfft_r<T>> plan;

  for (size_t iax=0; iax<axes.size(); ++iax)
    {
    constexpr auto vlen = VLEN<T>::val;
    size_t len=in.shape(axes[iax]);
    if ((!plan) || (len!=plan->length()))
      plan = get_plan<pocketfft_r<T>>(len);

#ifdef POCKETFFT_OPENMP
#pragma omp parallel num_threads(util::thread_count(nthreads, in.shape(), axes[iax]))
#endif
{
    auto storage = alloc_tmp<T>(in.shape(), len, sizeof(T));
    const auto &tin(iax==0 ? in : out);
    multi_iter<vlen> it(tin, out, axes[iax]);
#ifndef POCKETFFT_NO_VECTORS
    if (vlen>1)
      while (it.remaining()>=vlen)
        {
        using vtype = typename VTYPE<T>::type;
        it.advance(vlen);
        auto tdatav = reinterpret_cast<vtype *>(storage.data());
        for (size_t i=0; i<len; ++i)
          for (size_t j=0; j<vlen; ++j)
            tdatav[i][j] = tin[it.iofs(j,i)];
        if ((!r2c) && forward)
          for (size_t i=2; i<len; i+=2)
            tdatav[i] = -tdatav[i];
        forward ? plan->forward (tdatav, fct)
                : plan->backward(tdatav, fct);
        if (r2c && (!forward))
          for (size_t i=2; i<len; i+=2)
            tdatav[i] = -tdatav[i];
        for (size_t i=0; i<len; ++i)
          for (size_t j=0; j<vlen; ++j)
            out[it.oofs(j,i)] = tdatav[i][j];
        }
#endif
    while (it.remaining()>0)
      {
      it.advance(1);
      auto buf = it.stride_out() == sizeof(T) ?
        &out[it.oofs(0)] : reinterpret_cast<T *>(storage.data());

      if (buf != &tin[it.iofs(0)])
        for (size_t i=0; i<len; ++i)
          buf[i] = tin[it.iofs(i)];
      if ((!r2c) && forward)
        for (size_t i=2; i<len; i+=2)
          buf[i] = -buf[i];
      forward ? plan->forward(buf, fct) : plan->backward(buf, fct);
      if (r2c && (!forward))
        for (size_t i=2; i<len; i+=2)
          buf[i] = -buf[i];
      if (buf != &out[it.oofs(0)])
        for (size_t i=0; i<len; ++i)
          out[it.oofs(i)] = buf[i];
      }
} // end of parallel region
    fct = T(1); // factor has been applied, use 1 for remaining axes
    }
  }

#undef POCKETFFT_NTHREADS

template<typename T> void c2c(const shape_t &shape, const stride_t &stride_in,
  const stride_t &stride_out, const shape_t &axes, bool forward,
  const complex<T> *data_in, complex<T> *data_out, T fct,
  size_t nthreads=1)
  {
  if (util::prod(shape)==0) return;
  util::sanity_check(shape, stride_in, stride_out, data_in==data_out, axes);
  cndarr<cmplx<T>> ain(data_in, shape, stride_in);
  ndarr<cmplx<T>> aout(data_out, shape, stride_out);
  general_c(ain, aout, axes, forward, fct, nthreads);
  }

template<typename T> void dct(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  int type, const T *data_in, T *data_out, T fct, bool ortho, size_t nthreads=1)
  {
  if ((type<1) || (type>4)) throw invalid_argument("invalid DCT type");
  if (util::prod(shape)==0) return;
  util::sanity_check(shape, stride_in, stride_out, data_in==data_out, axes);
  cndarr<T> ain(data_in, shape, stride_in);
  ndarr<T> aout(data_out, shape, stride_out);
  if (type==1)
    general_dcst<T_dct1<T>>(ain, aout, axes, fct, ortho, type, true, nthreads);
  else if (type==4)
    general_dcst<T_dcst4<T>>(ain, aout, axes, fct, ortho, type, true, nthreads);
  else
    general_dcst<T_dcst23<T>>(ain, aout, axes, fct, ortho, type, true, nthreads);
  }

template<typename T> void dst(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  int type, const T *data_in, T *data_out, T fct, bool ortho, size_t nthreads=1)
  {
  if ((type<1) || (type>4)) throw invalid_argument("invalid DST type");
  if (util::prod(shape)==0) return;
  util::sanity_check(shape, stride_in, stride_out, data_in==data_out, axes);
  cndarr<T> ain(data_in, shape, stride_in);
  ndarr<T> aout(data_out, shape, stride_out);
  if (type==1)
    general_dcst<T_dst1<T>>(ain, aout, axes, fct, ortho, type, false, nthreads);
  else if (type==4)
    general_dcst<T_dcst4<T>>(ain, aout, axes, fct, ortho, type, false, nthreads);
  else
    general_dcst<T_dcst23<T>>(ain, aout, axes, fct, ortho, type, false, nthreads);
  }

template<typename T> void r2c(const shape_t &shape_in,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  bool forward, const T *data_in, complex<T> *data_out, T fct,
  size_t nthreads=1)
  {
  if (util::prod(shape_in)==0) return;
  util::sanity_check(shape_in, stride_in, stride_out, false, axis);
  cndarr<T> ain(data_in, shape_in, stride_in);
  shape_t shape_out(shape_in);
  shape_out[axis] = shape_in[axis]/2 + 1;
  ndarr<cmplx<T>> aout(data_out, shape_out, stride_out);
  general_r2c(ain, aout, axis, forward, fct, nthreads);
  }

template<typename T> void r2c(const shape_t &shape_in,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  bool forward, const T *data_in, complex<T> *data_out, T fct,
  size_t nthreads=1)
  {
  if (util::prod(shape_in)==0) return;
  util::sanity_check(shape_in, stride_in, stride_out, false, axes);
  r2c(shape_in, stride_in, stride_out, axes.back(), forward, data_in, data_out,
    fct, nthreads);
  if (axes.size()==1) return;

  shape_t shape_out(shape_in);
  shape_out[axes.back()] = shape_in[axes.back()]/2 + 1;
  auto newaxes = shape_t{axes.begin(), --axes.end()};
  c2c(shape_out, stride_out, stride_out, newaxes, forward, data_out, data_out,
    T(1), nthreads);
  }

template<typename T> void c2r(const shape_t &shape_out,
  const stride_t &stride_in, const stride_t &stride_out, size_t axis,
  bool forward, const complex<T> *data_in, T *data_out, T fct,
  size_t nthreads=1)
  {
  if (util::prod(shape_out)==0) return;
  util::sanity_check(shape_out, stride_in, stride_out, false, axis);
  shape_t shape_in(shape_out);
  shape_in[axis] = shape_out[axis]/2 + 1;
  cndarr<cmplx<T>> ain(data_in, shape_in, stride_in);
  ndarr<T> aout(data_out, shape_out, stride_out);
  general_c2r(ain, aout, axis, forward, fct, nthreads);
  }

template<typename T> void c2r(const shape_t &shape_out,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  bool forward, const complex<T> *data_in, T *data_out, T fct,
  size_t nthreads=1)
  {
  if (util::prod(shape_out)==0) return;
  if (axes.size()==1)
    return c2r(shape_out, stride_in, stride_out, axes[0], forward,
      data_in, data_out, fct, nthreads);
  util::sanity_check(shape_out, stride_in, stride_out, false, axes);
  auto shape_in = shape_out;
  shape_in[axes.back()] = shape_out[axes.back()]/2 + 1;
  auto nval = util::prod(shape_in);
  stride_t stride_inter(shape_in.size());
  stride_inter.back() = sizeof(cmplx<T>);
  for (int i=int(shape_in.size())-2; i>=0; --i)
    stride_inter[size_t(i)] =
      stride_inter[size_t(i+1)]*ptrdiff_t(shape_in[size_t(i+1)]);
  arr<complex<T>> tmp(nval);
  auto newaxes = shape_t({axes.begin(), --axes.end()});
  c2c(shape_in, stride_in, stride_inter, newaxes, forward, data_in, tmp.data(),
    T(1), nthreads);
  c2r(shape_out, stride_inter, stride_out, axes.back(), forward,
    tmp.data(), data_out, fct, nthreads);
  }

template<typename T> void r2r_fftpack(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  bool real2hermitian, bool forward, const T *data_in, T *data_out, T fct,
  size_t nthreads=1)
  {
  if (util::prod(shape)==0) return;
  util::sanity_check(shape, stride_in, stride_out, data_in==data_out, axes);
  cndarr<T> ain(data_in, shape, stride_in);
  ndarr<T> aout(data_out, shape, stride_out);
  general_r(ain, aout, axes, real2hermitian, forward, fct, nthreads);
  }

template<typename T> void r2r_separable_hartley(const shape_t &shape,
  const stride_t &stride_in, const stride_t &stride_out, const shape_t &axes,
  const T *data_in, T *data_out, T fct, size_t nthreads=1)
  {
  if (util::prod(shape)==0) return;
  util::sanity_check(shape, stride_in, stride_out, data_in==data_out, axes);
  cndarr<T> ain(data_in, shape, stride_in);
  ndarr<T> aout(data_out, shape, stride_out);
  general_hartley(ain, aout, axes, fct, nthreads);
  }

} // namespace detail

using detail::FORWARD;
using detail::BACKWARD;
using detail::shape_t;
using detail::stride_t;
using detail::c2c;
using detail::c2r;
using detail::r2c;
using detail::r2r_fftpack;
using detail::r2r_separable_hartley;
using detail::dct;
using detail::dst;

} // namespace pocketfft

#undef POCKETFFT_NOINLINE
#undef POCKETFFT_RESTRICT

#endif // POCKETFFT_HDRONLY_H
