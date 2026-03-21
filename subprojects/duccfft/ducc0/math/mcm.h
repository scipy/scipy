/*
 *  This file is part of DUCC.
 *
 *  DUCC is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  DUCC is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DUCC; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  DUCC is being developed at the Max-Planck-Institut fuer Astrophysik
 */

/*
 *  Copyright (C) 2023-2026 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef DUCC0_MCM_H
#define DUCC0_MCM_H

#include "ducc0/infra/simd.h"
#include "ducc0/infra/mav.h"
#include "ducc0/math/constants.h"
#include "ducc0/math/wigner3j.h"

namespace ducc0 {

using namespace std;

static bool neededForToeplitz(int l1, int l2, int l_exact, int l_toeplitz, int dl_band)
  {
  if (l_exact<0) return true;  // we want everything
  if (l1<l2) swap(l1,l2);
  if (l2<=l_exact) return true;
  if (l2==l_toeplitz) return true;
  if ((l2<l_toeplitz) && (l1-l2<=dl_band)) return true;
  if (l1==l2) return true;
  return false;
  }

template<typename T> static void toeplitz_fill(const vmav<T,2> &mat,
  const cmav<T,1> &diag0, int l_exact, int l_toeplitz, int dl_band, size_t nthreads)
  {
  int s0=mat.shape(0), s1=mat.shape(1);
  int lmax=max(s0-1,s1-1);
  int lsmall=min(s0-1,s1-1);
  if ((l_exact<0) || (l_exact>=lmax))  // nothing to do
    return;
  MR_assert(l_toeplitz<=lsmall, "l_toeplitz too large");
  MR_assert(l_exact<=lsmall, "l_exact too large");
  vector<double> diag(lmax+1);
  for (int l=0; l<=lmax; ++l)
    diag[l] = sqrt(double(diag0(l)));
  vector<double> row_toep(lmax+1-l_toeplitz);
  vector<double> row_exact(lmax+1-l_exact);
  if (s1>=s0)
    {
    for (int l=l_toeplitz; l<=lmax; ++l)
      row_toep[l-l_toeplitz] = mat(l_toeplitz,l)/(diag[l_toeplitz]*diag[l]);
    for (int l=l_exact; l<=lmax; ++l)
      row_exact[l-l_exact] = mat(l_exact,l)/(diag[l_exact]*diag[l]);
    }
  else
    {
    for (int l=l_toeplitz; l<=lmax; ++l)
      row_toep[l-l_toeplitz] = mat(l,l_toeplitz)/(diag[l_toeplitz]*diag[l]);
    for (int l=l_exact; l<=lmax; ++l)
      row_exact[l-l_exact] = mat(l,l_exact)/(diag[l_exact]*diag[l]);
    }

  execDynamic(lsmall-l_exact, nthreads, 8, [&](Scheduler &sched)
    {
    while (auto rng=sched.getNext())
      for (size_t xl1=rng.lo; xl1<rng.hi; ++xl1)
        {
        int l1 = int(xl1)+l_exact+1;
        for (int l2=l1+1; l2<=lmax; ++l2)
          {
          if ((l1<l_toeplitz)&&((l2-l1)<=dl_band))  // already computed exactly
            continue;
          if (l1==l_toeplitz)  // already computed exactly
            continue;
          T val=0;
          if ((l2-l1)<=(lmax-l_toeplitz))  // use column at l_toeplitz
            val = T(row_toep[l2-l1]*diag[l1]*diag[l2]);
          else  // use column at l_exact
            val = T(row_exact[l2-l1]*diag[l1]*diag[l2]);
          if (l2<s1) mat(l1,l2) = val;
          if (l2<s0) mat(l2,l1) = val;
          }
        }
    });
  }

template<typename Tsimd, size_t nspec> inline array<Tsimd,nspec> sum_wig00
  (int el1, int el2, int lmax_spec, const Tsimd * DUCC0_RESTRICT res, const vmav<double,2> &spec2)
  {
  array<Tsimd,nspec> val;
  for (size_t ispec=0; ispec<nspec; ++ispec)
    val[ispec]=0;
  int el3min = el2-el1;
  int max_i = min(el1+el2, int(lmax_spec)) - el3min;
  for (int i=0, i2=0; i<=max_i; i+=2, ++i2)
    {
    int el3 = el3min+i;
    for (size_t ispec=0; ispec<nspec; ++ispec)
      val[ispec] += res[i2]*Tsimd(&spec2(ispec,el3), element_aligned_tag());
    }
  return val;
  }
template<typename Tsimd, typename Tspec, typename Tval> inline void sum_wig00
  (int el1, int el2, int lmax_spec, size_t nspec, const Tsimd * DUCC0_RESTRICT res, const Tspec &spec2, Tval &val)
  {
  for (size_t ispec=0; ispec<nspec; ++ispec)
    val[ispec]=0;
  int el3min = el2-el1;
  int max_i = min(el1+el2, int(lmax_spec)) - el3min;
  for (int i=0, i2=0; i<=max_i; i+=2, ++i2)
    {
    int el3 = el3min+i;
    for (size_t ispec=0; ispec<nspec; ++ispec)
      val[ispec] += res[i2]*Tsimd(&spec2(ispec,el3), element_aligned_tag());
    }
  }

template<typename Tsimd, typename Tout, size_t nspec> inline void store_mat
  (int el1, int el2, size_t s1_, size_t s2_, const vmav<Tout,3> &mat, const vmav<Tout,2> &diag, const array<Tsimd, nspec> &val)
  {
  int s1=int(s1_), s2=int(s2_);
  constexpr auto vlen=int(Tsimd::size());
  for (size_t ispec=0; ispec<nspec; ++ispec)
    {
    for (int k=0; k<vlen; ++k)
      {
      if ((el1<s1) && (el2+k<s2)) mat(ispec, el1, el2+k) = Tout(val[ispec][k]);
      if ((el2+k<s1) && (el1<s2)) mat(ispec, el2+k, el1) = Tout(val[ispec][k]);
      }
    if (el1==el2) diag(ispec, el1) = Tout(val[ispec][0]);
    }
  }
template<typename Tsimd, typename Tout, typename Tval> inline void store_mat
  (int el1, int el2, size_t s1_, size_t s2_, size_t nspec, const vmav<Tout,3> &mat, const vmav<Tout,2> &diag, const Tval &val)
  {
  int s1=int(s1_), s2=int(s2_);
  constexpr auto vlen=int(Tsimd::size());
  for (size_t ispec=0; ispec<nspec; ++ispec)
    {
    for (int k=0; k<vlen; ++k)
      {
      if ((el1<s1) && (el2+k<s2)) mat(ispec, el1, el2+k) = Tout(val[ispec][k]);
      if ((el2+k<s1) && (el1<s2)) mat(ispec, el2+k, el1) = Tout(val[ispec][k]);
      }
    if (el1==el2) diag(ispec, el1) = Tout(val[ispec][0]);
    }
  }
template<typename Tsimd, typename Tout> inline void zero_mat
  (int el1, int el2, size_t s1_, size_t s2_, const vmav<Tout,3> &mat, const vmav<Tout,2> &diag)
  {
  int s1=int(s1_), s2=int(s2_);
  size_t nmat=mat.shape(0);
  constexpr auto vlen=int(Tsimd::size());
  for (size_t imat=0; imat<nmat; ++imat)
    {
    for (int k=0; k<vlen; ++k)
      {
      if ((el1<s1) && (el2+k<s2)) mat(imat, el1, el2+k) = Tout(0);
      if ((el2+k<s1) && (el1<s2)) mat(imat, el2+k, el1) = Tout(0);
      }
    if (el1==el2) diag(imat, el1) = Tout(0);
    }
  }

template<typename Tout> void coupling_matrix_spin0_square(const cmav<double,2> &spec,
  const vmav<Tout,3> &mat, int l_exact, int l_toeplitz, int dl_band, size_t nthreads)
  {
  size_t nspec=spec.shape(0);
  MR_assert(spec.shape(1)>=1, "spec.shape[1] is too small.");
  auto lmax_spec = spec.shape(1)-1;
  MR_assert(mat.shape(0)==nspec, "number of spectra and matrices mismatch");
  MR_assert(mat.size()>0, "matrix must not be zero-sized");
  size_t s1=mat.shape(1), s2=mat.shape(2);
  size_t lmax = max(s1-1, s2-1);
  size_t lsmall = min(s1-1, s2-1);
  using Tsimd = native_simd<double>;
  constexpr size_t vlen = Tsimd::size();
  auto lmax_spec_used = min(2*lmax, lmax_spec);
  auto spec2(vmav<double,2>::build_noncritical({nspec, lmax_spec_used+1+vlen-1}, PAGE_IN(nthreads)));
  for (size_t l=0; l<=lmax_spec_used; ++l)
    for (size_t i=0; i<nspec; ++i)
      spec2(i,l) = spec(i,l)/ducc0::fourpi*(2.*l+1.);
  for (size_t l=lmax_spec_used+1; l<spec2.shape(1); ++l)
    for (size_t i=0; i<nspec; ++i)
      spec2(i,l) = 0.;
  auto diag(vmav<Tout,2>::build_noncritical({nspec, lmax+1}));
  execDynamic(lmax+1, nthreads, 1, [&](ducc0::Scheduler &sched)
    {
    vmav<Tsimd,1> resfullv({lmax+1});
    vmav<Tsimd,1> val_({nspec});
    Tsimd * DUCC0_RESTRICT val = val_.data();
    Tsimd lofs;
    for (size_t k=0; k<vlen; ++k)
      lofs[k]=double(k);
    while (auto rng=sched.getNext()) for(int el1=int(rng.lo); el1<int(rng.hi); ++el1)
      {
      for (int el2=el1; el2<=((el1<=int(lsmall))?int(lmax):el1); el2+=vlen)
        {
        bool necessary=false;
        for (size_t i=0; i<vlen; ++i)
          if (neededForToeplitz(el1, el2+i, l_exact, l_toeplitz, dl_band))
            necessary=true;
        if (!necessary) continue;

        int el3min = el2-el1;
        if (el3min<=int(lmax_spec))
          {
          wigner3j_00_vec_squared_compact(Tsimd(el1), Tsimd(el2)+lofs,
            subarray<1>(resfullv, {{size_t(0), size_t(el1+1)}}));
          const Tsimd * DUCC0_RESTRICT res = resfullv.data();

          if (nspec==1)
            store_mat(el1, el2, s1, s2, mat, diag, sum_wig00<Tsimd,1>(el1, el2, lmax_spec, res, spec2));
          else if (nspec==2)
            store_mat(el1, el2, s1, s2, mat, diag, sum_wig00<Tsimd,2>(el1, el2, lmax_spec, res, spec2));
          else if (nspec==3)
            store_mat(el1, el2, s1, s2, mat, diag, sum_wig00<Tsimd,3>(el1, el2, lmax_spec, res, spec2));
          else if (nspec==4)
            store_mat(el1, el2, s1, s2, mat, diag, sum_wig00<Tsimd,4>(el1, el2, lmax_spec, res, spec2));
          else if (nspec<=50)
            {
            array<Tsimd,50> val;
            sum_wig00<Tsimd>(el1, el2, lmax_spec, nspec, res, spec2, val); 
            store_mat<Tsimd>(el1, el2, s1, s2, nspec, mat, diag, val);
            }
          else
            {
            sum_wig00<Tsimd>(el1, el2, lmax_spec, nspec, res, spec2, val); 
            store_mat<Tsimd>(el1, el2, s1, s2, nspec, mat, diag, val);
            }
          }
        else
          zero_mat<Tsimd>(el1, el2, s1, s2, mat, diag);
        }
      }
    });
  if (l_exact>=0)
    for (size_t imat=0; imat<mat.shape(0); ++imat)
      toeplitz_fill(subarray<2>(mat,{{imat},{},{}}), subarray<1>(diag,{{imat},{}}), l_exact, l_toeplitz, dl_band, nthreads);
  }

template<size_t opmask, typename Tsimd, size_t nspec> inline array<array<Tsimd,4>,nspec> sum_wig02
  (int el1, int el2, size_t lmax_spec, const Tsimd * DUCC0_RESTRICT wp0,
   const Tsimd * DUCC0_RESTRICT wp1, const vmav<double,2> &spec2)
  {
  int el3min = el2-el1;
  int el3max = el2+el1;
  int maxidx = min(el3max, int(lmax_spec));
  array<array<Tsimd,4>,nspec> val;
  for (size_t ispec=0; ispec<nspec; ++ispec)
    for (size_t j=0; j<4; ++j)
      val[ispec][j]=0;
  for (int el3=el3min; el3<=maxidx; el3+=2)
    {
    const Tsimd w0=wp0[el3], w1=wp1[el3];
    const Tsimd w00=w0*w0, w01=w0*w1, w11=w1*w1;
    const Tsimd w11p1=wp1[el3+1]*wp1[el3+1];
    for (size_t ispec=0; ispec<nspec; ++ispec)
      {
      Tsimd sp = Tsimd(&spec2(ispec,el3), element_aligned_tag());
      if constexpr (opmask&1)
        val[ispec][0] += w00*sp;
      if constexpr (opmask&2)
        val[ispec][1] += w01*sp;
      if constexpr (opmask&4)
        val[ispec][2] += w11*sp;
      if constexpr (opmask&8)
        val[ispec][3] += w11p1*Tsimd(&spec2(ispec,el3+1), element_aligned_tag());
      }
    }
  return val;
  }
template<size_t opmask, typename Tsimd, typename Tval> inline void sum_wig02
  (int el1, int el2, size_t lmax_spec, size_t nspec, const Tsimd * DUCC0_RESTRICT wp0,
   const Tsimd * DUCC0_RESTRICT wp1, const vmav<double,2> &spec2, Tval &val)
  {
  int el3min = el2-el1;
  int el3max = el2+el1;
  int maxidx = min(el3max, int(lmax_spec));
  for (size_t ispec=0; ispec<nspec; ++ispec)
    for (size_t j=0; j<4; ++j)
      val[ispec][j]=0;
  for (int el3=el3min; el3<=maxidx; el3+=2)
    {
    const Tsimd w0=wp0[el3], w1=wp1[el3];
    const Tsimd w00=w0*w0, w01=w0*w1, w11=w1*w1;
    const Tsimd w11p1=wp1[el3+1]*wp1[el3+1];
    for (size_t ispec=0; ispec<nspec; ++ispec)
      {
      Tsimd sp = Tsimd(&spec2(ispec,el3), element_aligned_tag());
      if constexpr (opmask&1)
        val[ispec][0] += w00*sp;
      if constexpr (opmask&2)
        val[ispec][1] += w01*sp;
      if constexpr (opmask&4)
        val[ispec][2] += w11*sp;
      if constexpr (opmask&8)
        val[ispec][3] += w11p1*Tsimd(&spec2(ispec,el3+1), element_aligned_tag());
      }
    }
  }
template<typename Tsimd, typename Tout, typename Tval> inline void store_mat02
  (int el1, int el2, size_t s1_, size_t s2_, const vector<int> &optype,
   const vmav<Tout,3> &mat, const vmav<Tout,2> &diag, const Tval &val)
  {
  int s1=int(s1_), s2=int(s2_);
  constexpr auto vlen=int(Tsimd::size());
  for (size_t ispec=0, imat=0; ispec<optype.size(); ++ispec)
    {
    auto op = optype[ispec];
    for (int k=0; k<vlen; ++k)
      {
      if ((el2+k<s1) && (el1<s2))
        {
        if (op==4)
          {
          mat(imat, el2+k, el1) = Tout(val[ispec][2][k]);
          mat(imat+1, el2+k, el1) = Tout(val[ispec][3][k]);
          }
        else
          mat(imat, el2+k, el1) = Tout(val[ispec][op][k]);
        }
      if ((el1<s1) && (el2+k<s2))
        {
        if (op==4)
          {
          mat(imat, el1, el2+k) = Tout(val[ispec][2][k]);
          mat(imat+1, el1, el2+k) = Tout(val[ispec][3][k]);
          }
        else
          mat(imat, el1, el2+k) = Tout(val[ispec][op][k]);
        }
      }
    if (el1==el2)
      {
      if (op==4)
        {
        diag(imat, el1) = Tout(val[ispec][2][0]);
        diag(imat+1, el1) = Tout(val[ispec][3][0]);
        }
      else
        diag(imat, el1) = Tout(val[ispec][op][0]);
      }
    imat += (op==4) ? 2 : 1;
    }
  }

template<size_t opmask, typename Tout> void coupling_matrix_rect(
  const cmav<double,2> &spec, const vmav<Tout,3> &mat,
  const vector<int> &optype, int l_exact, int l_toeplitz, int dl_band,
  size_t nthreads)
  {
  if constexpr ((opmask&14)==0)
    return coupling_matrix_spin0_square(spec, mat,
      l_exact, l_toeplitz, dl_band, nthreads);

  size_t nspec=spec.shape(0);
  MR_assert(optype.size()==nspec, "incorrect optype size");
  MR_assert(spec.shape(1)>=1, "lmax_spec is too small.");
  size_t nmat=0;
  for (auto op: optype)
    {
    MR_assert((op>=0) && (op<5), "bad optype entry");
    nmat += (op<4) ? 1 : 2;
    }
  MR_assert(mat.shape(0)==nmat, "incorrect number of output matrices");
  MR_assert(mat.size()>0, "matrix must not be zero-sized");
  size_t s1=mat.shape(1), s2=mat.shape(2);
  size_t lmax = max(s1-1, s2-1);
  size_t lsmall = min(s1-1, s2-1);
  auto lmax_spec = spec.shape(1)-1;
  using Tsimd = native_simd<double>;
  constexpr size_t vlen = Tsimd::size();
  auto lmax_spec_used = min(2*lmax, lmax_spec);
  auto spec2(vmav<double,2>::build_noncritical
    ({nspec, lmax_spec_used+1+vlen-1+1}, PAGE_IN(nthreads)));
  for (size_t l=0; l<=lmax_spec_used; ++l)
    for (size_t i=0; i<nspec; ++i)
      spec2(i,l) = spec(i,l)/ducc0::fourpi*(2.*l+1.);
  for (size_t l=lmax_spec_used+1; l<spec2.shape(1); ++l)
    for (size_t i=0; i<nspec; ++i)
      spec2(i,l) = 0.;
  auto diag(vmav<Tout,2>::build_noncritical({nmat, lmax+1}));
  execDynamic(lmax+1, nthreads, 1, [&](ducc0::Scheduler &sched)
    {
// res arrays are one larger to make loops simpler below
    vmav<Tsimd,2> wig({2, 2*lmax+1+1});
    vmav<array<Tsimd,4>,1> val_({nspec});
    array<Tsimd,4> * DUCC0_RESTRICT val = val_.data();
    Tsimd lofs;
    for (size_t k=0; k<vlen; ++k)
      lofs[k]=double(k);
    while (auto rng=sched.getNext()) for(int el1=int(rng.lo); el1<int(rng.hi); ++el1)
      {
      for (int el2=el1; el2<=((el1<=int(lsmall))?int(lmax):el1); el2+=vlen)
        {
        bool necessary=false;
        for (size_t i=0; i<vlen; ++i)
          if (neededForToeplitz(el1, el2+i, l_exact, l_toeplitz, dl_band))
            necessary=true;
        if (!necessary) continue;

        int el3min = el2-el1;
        int el3max = el2+el1;
        if (el3min<=int(lmax_spec))
          {
          auto tmp=subarray<2>(wig,{{},{size_t(el3min), size_t(el3max+2)}});
          if constexpr(opmask&3)
            flexible_wigner3j_vec(Tsimd(el1), Tsimd(el2)+lofs, 0, 0,
              Tsimd(el3min)+lofs, subarray<1>(tmp, {{0}, {}}));
          if constexpr(opmask&14)
            flexible_wigner3j_vec(Tsimd(el1), Tsimd(el2)+lofs, -2, 2,
              Tsimd(el3min)+lofs, subarray<1>(tmp, {{1}, {}}));
          const Tsimd * DUCC0_RESTRICT wp0 = &wig(0,0);
          const Tsimd * DUCC0_RESTRICT wp1 = &wig(1,0);

          if (nspec==1)
            {
            auto val = sum_wig02<opmask, Tsimd, 1> (el1, el2, lmax_spec, wp0, wp1, spec2);
            store_mat02<Tsimd, Tout> (el1, el2, s1, s2, optype, mat, diag, val);
            }
          else if (nspec==2)
            {
            auto val = sum_wig02<opmask, Tsimd, 2> (el1, el2, lmax_spec, wp0, wp1, spec2);
            store_mat02<Tsimd, Tout> (el1, el2, s1, s2, optype, mat, diag, val);
            }
          else if (nspec==3)
            {
            auto val = sum_wig02<opmask, Tsimd, 3> (el1, el2, lmax_spec, wp0, wp1, spec2);
            store_mat02<Tsimd, Tout> (el1, el2, s1, s2, optype, mat, diag, val);
            }
          else if (nspec==4)
            {
            auto val = sum_wig02<opmask, Tsimd, 4> (el1, el2, lmax_spec, wp0, wp1, spec2);
            store_mat02<Tsimd, Tout> (el1, el2, s1, s2, optype, mat, diag, val);
            }
          else if (nspec<=50)
            {
            array<array<Tsimd,4>,50> val;
            sum_wig02<opmask, Tsimd> (el1, el2, lmax_spec, nspec, wp0, wp1, spec2, val);
            store_mat02<Tsimd, Tout> (el1, el2, s1, s2, optype, mat, diag, val);
            }
          else
            {
            sum_wig02<opmask, Tsimd> (el1, el2, lmax_spec, nspec, wp0, wp1, spec2, val);
            store_mat02<Tsimd, Tout> (el1, el2, s1, s2, optype, mat, diag, val);
            }
          }
        else
          zero_mat<Tsimd, Tout> (el1, el2, s1, s2, mat, diag);
        }
      }
    });
  if (l_exact>=0)
    for (size_t imat=0; imat<mat.shape(0); ++imat)
      toeplitz_fill(subarray<2>(mat,{{imat},{},{}}), subarray<1>(diag,{{imat},{}}), l_exact, l_toeplitz, dl_band, nthreads);
  }

#if 0
template<typename Tout> void coupling_matrix_spin0and2_pure(const cmav<double,3> &spec,
  size_t lmax, const vmav<Tout,4> &mat, size_t nthreads)
  {
  using Tsimd = native_simd<double>;
  constexpr size_t vlen=Tsimd::size();
  constexpr size_t ncomp_spec=4;
  constexpr size_t ncomp_mat=4;
  size_t nspec=spec.shape(0);
  MR_assert(spec.shape(1)==ncomp_spec, "spec.shape[1] must be 4.");
  MR_assert(spec.shape(2)>=1, "lmax_spec is too small.");
  MR_assert(mat.shape(0)==nspec, "number of spectra and matrices mismatch");
  MR_assert(mat.shape(1)==ncomp_mat, "bad number of matrix components");
  MR_assert(mat.shape(2)==lmax+1, "bad number of matrix entries");
  MR_assert(mat.shape(3)==lmax+1, "bad number of matrix entries");
  auto lmax_spec = spec.shape(2)-1;
  auto lmax_spec_used = min(2*lmax, lmax_spec);
  auto spec2(vmav<double,3>::build_noncritical
    ({nspec, ncomp_spec, lmax_spec_used+1+vlen-1+1}, PAGE_IN(nthreads)));
  for (size_t l=0; l<=lmax_spec_used; ++l)
    for (size_t j=0; j<ncomp_spec; ++j)
      for (size_t i=0; i<nspec; ++i)
        spec2(i,j,l) = spec(i,j,l)/ducc0::fourpi*(2.*l+1.);
  for (size_t l=lmax_spec_used+1; l<spec2.shape(2); ++l)
    for (size_t j=0; j<ncomp_spec; ++j)
      for (size_t i=0; i<nspec; ++i)
        spec2(i,j,l) = 0.;
  vector<double> nom1(2*lmax+1+vlen-1+1), nom2(2*lmax+1+vlen-1+1);
  for (size_t el3=0; el3<nom1.size(); ++el3)
    {
    nom1[el3] = 2.*sqrt((el3+1.)*el3);
    nom2[el3] = sqrt((el3+2.)*(el3+1.)*el3*(el3-1.));
    }
  execDynamic(lmax+1, nthreads, 1, [&](ducc0::Scheduler &sched)
    {
    // res arrays are one larger to make loops simpler below
    vmav<Tsimd,2> wig({6, 2*lmax+1+1});
    constexpr size_t nvcomp = 7;
    vmav<array<Tsimd,nvcomp>,1> val_({nspec});
    array<Tsimd,nvcomp> * DUCC0_RESTRICT val = val_.data();
    Tsimd lofs;
    for (size_t k=0; k<vlen; ++k)
      lofs[k]=double(k);
    while (auto rng=sched.getNext()) for(int el1=int(rng.lo); el1<int(rng.hi); ++el1)
      {
      for (int xel2=el1; xel2<=int(lmax); xel2+=vlen)
        {
        Tsimd el2=Tsimd(xel2)+lofs;
        int el3min = abs(xel2-el1);
        int el3max = el1+xel2;
        Tsimd xdenom1 = blend(el2>Tsimd(1.), sqrt(Tsimd(1.) / ((el2-1.)*(el2+2.))), Tsimd(0.)),
              xdenom2 = blend(el2>Tsimd(1.), sqrt(Tsimd(1.) / ((el2+2.)*(el2+1.)*el2*(el2-1.))), Tsimd(0.));
        double xxdenom1 = (el1>1) ? sqrt(1. / ((el1-1.)*(el1+2.))) : 0,
               xxdenom2 = (el1>1) ? sqrt(1. / ((el1+2.)*(el1+1.)*el1*(el1-1.))): 0;
        if (el3min<=int(lmax_spec))
          {
          {
          auto tmp = subarray<2>(wig, {{}, {size_t(el3min), size_t(el3max+2)}});
          constexpr array<int,6> m1 {{0, -2, -2, -2, -2, -2}};
          constexpr array<int,6> m2 {{0,  2,  1,  0,  1,  0}};
          Tsimd tel1(el1);
          array<Tsimd,6> xl1 {{tel1, tel1, tel1, tel1,  el2,  el2}};
          array<Tsimd,6> xl2 {{ el2,  el2,  el2,  el2, tel1, tel1}};
          for (size_t ii=0; ii<6; ++ii)
            flexible_wigner3j_vec(xl1[ii], xl2[ii], m1[ii], m2[ii],
              Tsimd(el3min)+lofs, subarray<1>(tmp, {{ii}, {}}));
          }

          for (size_t ispec=0; ispec<nspec; ++ispec)
            for (size_t j=0; j<nvcomp; ++j)
              val[ispec][j]=0;
          int maxidx = min(el3max, int(lmax_spec));
          for (int el3=el3min; el3<=maxidx; el3+=2)
            {
            Tsimd fac_b = Tsimd(&nom1[el3],element_aligned_tag())*xdenom1,
                  fac_c = Tsimd(&nom2[el3],element_aligned_tag())*xdenom2,
                  xfac_b = Tsimd(&nom1[el3],element_aligned_tag())*xxdenom1,
                  xfac_c = Tsimd(&nom2[el3],element_aligned_tag())*xxdenom2;
//                  fac_b2 = Tsimd(&nom1[el3+1],element_aligned_tag())*xdenom1,
//                  fac_c2 = Tsimd(&nom2[el3+1],element_aligned_tag())*xdenom2,
//                  xfac_b2 = Tsimd(&nom1[el3+1],element_aligned_tag())*xxdenom1,
//                  xfac_c2 = Tsimd(&nom2[el3+1],element_aligned_tag())*xxdenom2;
            for (size_t ispec=0; ispec<nspec; ++ispec)
              {
              const Tsimd s0(&spec2(ispec,0,el3), element_aligned_tag()),
                          s1(&spec2(ispec,1,el3), element_aligned_tag()),
                          s2(&spec2(ispec,2,el3), element_aligned_tag()),
                          s3(&spec2(ispec,3,el3), element_aligned_tag());
              val[ispec][0] += wig(0,el3)*wig(0,el3)*s0;
              auto combin = wig(1,el3) + fac_b*wig(2,el3) + fac_c*wig(3,el3);
              val[ispec][1] += wig(0,el3)*combin*s1;
              val[ispec][2] += wig(0,el3)*combin*s2;
              val[ispec][3] += combin*combin*Tsimd(&spec2(ispec,3,el3), element_aligned_tag());
              auto xcombin = wig(1,el3) + xfac_b*wig(4,el3) + xfac_c*wig(5,el3);
              val[ispec][4] += wig(0,el3)*xcombin*s1;
              val[ispec][5] += wig(0,el3)*xcombin*s2;
              val[ispec][6] += xcombin*xcombin*s3;
//              auto combin2 = wig(1,el3+1) + fac_b2*wig(2,el3+1) + fac_c2*wig(3,el3+1);
//              val[ispec][7] += combin2*combin2*Tsimd(&spec2(ispec,3,el3+1), element_aligned_tag());
//              auto  xcombin2 = wig(1,el3+1) + xfac_b2*wig(4,el3+1) + xfac_c2*wig(5,el3+1);
//              val[ispec][8] += xcombin2*xcombin2*Tsimd(&spec2(ispec,3,el3+1), element_aligned_tag());
              }
            }
          for (size_t ispec=0; ispec<nspec; ++ispec)
            for (size_t k=0; k<vlen; ++k)
              if (el2[k]<=lmax)
                {
                mat(ispec, 0, xel2+k, el1) = Tout((2*el1+1.)*val[ispec][0][k]);
                mat(ispec, 1, xel2+k, el1) = Tout((2*el1+1.)*val[ispec][1][k]);
                mat(ispec, 2, xel2+k, el1) = Tout((2*el1+1.)*val[ispec][2][k]);
                mat(ispec, 3, xel2+k, el1) = Tout((2*el1+1.)*val[ispec][3][k]);
                mat(ispec, 0, el1, xel2+k) = Tout((2*el2[k]+1.)*val[ispec][0][k]);
                mat(ispec, 1, el1, xel2+k) = Tout((2*el2[k]+1.)*val[ispec][4][k]);
                mat(ispec, 2, el1, xel2+k) = Tout((2*el2[k]+1.)*val[ispec][5][k]);
                mat(ispec, 3, el1, xel2+k) = Tout((2*el2[k]+1.)*val[ispec][6][k]);
//                mat(ispec, 4, xel2+k, el1) = (2*el1+1.)*val[ispec][4][k];
//                mat(ispec, 4, el1, xel2+k) = (2*el2[k]+1.)*val[ispec][8][k];
                }
          }
        else
          for (size_t ispec=0; ispec<nspec; ++ispec)
            for (size_t j=0; j<ncomp_mat; ++j)
              for (size_t k=0; k<vlen; ++k)
                if (el2[k]<=lmax)
                  mat(ispec, j, xel2+k, el1) = mat(ispec, j, el1, xel2+k) = 0.;
        }
      }
    });
  }

#endif

}

#endif
