/*! \file ducc0/infra/mav.cc
 *  Classes for dealing with multidimensional arrays
 *
 *  \copyright Copyright (C) 2019-2023 Max-Planck-Society
 *  \author Martin Reinecke
 *  */

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

#include <cmath>
#include <tuple>
#include <algorithm>
#include "ducc0/infra/mav.h"

namespace ducc0 {

namespace detail_mav {

using namespace std;

DUCC0_NOINLINE void opt_shp_str(fmav_info::shape_t &shp, vector<fmav_info::stride_t> &str)
  {
  if (shp.size()>1)
    {
    // sort dimensions in order of descending stride, as far as possible
    vector<size_t> strcrit(shp.size(),~size_t(0));
    for (const auto &curstr: str)
      for (size_t i=0; i<curstr.size(); ++i)
        strcrit[i] = min(strcrit[i],size_t(abs(curstr[i])));

    for (size_t lastdim=shp.size(); lastdim>1; --lastdim)
      {
      auto dim = size_t(min_element(strcrit.begin(),strcrit.begin()+lastdim)
                        -strcrit.begin());
      if ((dim+1!=lastdim) && (strcrit[dim]<strcrit[lastdim-1]))
        {
        swap(strcrit[dim], strcrit[lastdim-1]);
        swap(shp[dim], shp[lastdim-1]);
        for (auto &curstr: str)
          swap(curstr[dim], curstr[lastdim-1]);
        }
      }
    // try merging dimensions
    size_t ndim = shp.size();
    if (ndim>1)
      for (size_t d0=ndim-2; d0+1>0; --d0)
        {
        bool can_merge = true;
        for (const auto &curstr: str)
          can_merge &= curstr[d0] == ptrdiff_t(shp[d0+1])*curstr[d0+1];
        if (can_merge)
          {
          for (auto &curstr: str)
            curstr.erase(curstr.begin()+d0);
          shp[d0+1] *= shp[d0];
          shp.erase(shp.begin()+d0);
          }
        }
    }
  }

DUCC0_NOINLINE tuple<fmav_info::shape_t, vector<fmav_info::stride_t>>
  multiprep(const vector<fmav_info> &info)
  {
  auto narr = info.size();
  MR_assert(narr>=1, "need at least one array");
  for (size_t i=1; i<narr; ++i)
    MR_assert(info[i].shape()==info[0].shape(), "shape mismatch");
  fmav_info::shape_t shp;
  vector<fmav_info::stride_t> str(narr);
  for (size_t i=0; i<info[0].ndim(); ++i)
    if (info[0].shape(i)!=1) // remove axes of length 1
      {
      shp.push_back(info[0].shape(i));
      for (size_t j=0; j<narr; ++j)
        str[j].push_back(info[j].stride(i));
      }
  opt_shp_str(shp, str);
  return make_tuple(shp, str);
  }

DUCC0_NOINLINE tuple<fmav_info::shape_t, vector<fmav_info::stride_t>, size_t, size_t>
  multiprep(const vector<fmav_info> &info, const vector<size_t> &tsizes)
  {
  auto narr = info.size();
  MR_assert(narr>=1, "need at least one array");
  MR_assert(tsizes.size()==narr, "tsizes has wrong length");
  for (size_t i=1; i<narr; ++i)
    MR_assert(info[i].shape()==info[0].shape(), "shape mismatch");
  fmav_info::shape_t shp;
  vector<fmav_info::stride_t> str(narr);
  for (size_t i=0; i<info[0].ndim(); ++i)
    if (info[0].shape(i)!=1) // remove axes of length 1
      {
      shp.push_back(info[0].shape(i));
      for (size_t j=0; j<narr; ++j)
        str[j].push_back(info[j].stride(i));
      }
  opt_shp_str(shp, str);

  bool transpose=false;
  size_t ndim=shp.size();
  if (ndim<2) return make_tuple(shp, str, 0, 0);

  for (size_t j=0; j<narr; ++j)
    if (abs(str[j][ndim-2])<abs(str[j][ndim-1]))
      transpose=true;
  if (!transpose) return make_tuple(shp, str, 0, 0);
  bool crit0=false, crit1=false;
  for (size_t j=0; j<narr; ++j)
    {
    auto str0 = abs(str[j][ndim-2]);
    auto str1 = abs(str[j][ndim-1]);
    if (((tsizes[j]*str0)%4096)==0) crit0=true;
    if (((tsizes[j]*str1)%4096)==0) crit1=true;
    }
  if (crit0 && crit1) return make_tuple(shp, str, 4, 4);
  size_t load=0;
  for (size_t i=0; i<narr; ++i)
    {
    auto tmp = tsizes[i]*min<size_t>(abs(str[i][ndim-2]), abs(str[i][ndim-1]));
    load += min<size_t>(tmp, 64);
    }
  if (crit0) return make_tuple(shp, str, 4, max<size_t>(8, 16384/(4*load)));
  if (crit1) return make_tuple(shp, str, max<size_t>(8, 16384/(4*load)), 4);
  size_t blk = max<size_t>(8, size_t(sqrt(16384./load)));
  return make_tuple(shp, str, blk, blk);
  }

DUCC0_NOINLINE tuple<fmav_info::shape_t, vector<fmav_info::stride_t>>
  multiprep_noopt(const vector<fmav_info> &info)
  {
  auto narr = info.size();
  MR_assert(narr>=1, "need at least one array");
  for (size_t i=1; i<narr; ++i)
    MR_assert(info[i].shape()==info[0].shape(), "shape mismatch");
  fmav_info::shape_t shp;
  vector<fmav_info::stride_t> str(narr);
  for (size_t i=0; i<info[0].ndim(); ++i)
    {
    shp.push_back(info[0].shape(i));
    for (size_t j=0; j<narr; ++j)
      str[j].push_back(info[j].stride(i));
    }
  return make_tuple(shp, str);
  }

}}
