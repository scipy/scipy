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

/** \file ducc0/infra/bucket_sort.h
 *
 *  \copyright Copyright (C) 2020-2021 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_BUCKET_SORT_H
#define DUCC0_BUCKET_SORT_H

#include <cstdint>
#include <algorithm>
#include <array>
#include <cstddef>
#include <type_traits>
#include <vector>
#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/threading.h"
#include "ducc0/infra/aligned_array.h"
#include "ducc0/math/math_utils.h"

namespace ducc0 {

namespace detail_bucket_sort {

using namespace std;

template<typename RAidx, typename Tkey, typename Tidx> void subsort
  (RAidx idx, quick_array<Tkey> &keys, size_t keybits, size_t lo,
   size_t hi, vector<Tidx> &numbers, quick_array<Tidx> &idxbak,
   quick_array<Tkey> &keybak)
  {
  auto nval = hi-lo;
  if (nval<=1) return;
  size_t keyshift = (keybits<=8) ? 0 : keybits-8;
  size_t nkeys = min<size_t>(size_t(1)<<keybits, 256);
  size_t keymask = nkeys-1;
  if (keybak.size()<nval) keybak.resize(nval);
  if (idxbak.size()<nval) idxbak.resize(nval);
  if (numbers.size()<nkeys) numbers.resize(nkeys);
  for (size_t i=0; i<nkeys; ++i) numbers[i]=0;
  for (size_t i=0; i<nval; ++i)
    {
    keybak[i] = keys[i+lo];
    idxbak[i] = idx[i+lo];
    ++numbers[(keys[i+lo]>>keyshift)&keymask];
    }
  Tidx ofs=0;
  for (auto &x: numbers)
    {
    auto tmp = x;
    x = ofs;
    ofs += tmp;
    }
  for (size_t i=0; i<nval; ++i)
    {
    auto loc = (keybak[i]>>keyshift)&keymask;
    keys[lo+numbers[loc]] = keybak[i];
    idx[lo+numbers[loc]] = idxbak[i];
    ++numbers[loc];
    }
  if (keyshift==0) return;
  keybits -= 8;
  vector<Tidx> newnumbers;
  for (size_t i=0; i<nkeys; ++i)
    subsort(idx, keys, keybits, lo + ((i==0) ? 0 : numbers[i-1]),
      lo+numbers[i], newnumbers, idxbak, keybak);
  }

/// Efficient bucket sort to determine an ascending order for the given set
/// of keys.
/** \param keys
 *         Random access iterator to a set of keys to be processed in ascending
 *         order.
 *  \param res
 *         Random access iterator to a writable integer-like sequence which will
 *         hold the key indices in the desired order.
 *  \param nval
 *         The number of keys to be sorted
 *  \param max_key
 *         An upper limit for the maximum key value. The closer this is to the
 *         actual maximum key, the better the algorithm will perform.
 *  \param nthreads
 *         The number of threads to use for the sorting
 */
template<typename RAidx, typename RAkey> void bucket_sort
  (RAkey keys, RAidx res, size_t nval, size_t max_key, size_t nthreads)
  {
  nthreads = min(nthreads, thread_pool_size());
  using Tidx = typename remove_reference<decltype(*res)>::type;
  using Tkey = typename remove_reference<decltype(*keys)>::type;
  // align members with cache lines
  struct alignas(64) vbuf { vector<Tidx> v; };
  vector<vbuf> numbers(nthreads);
  auto keybits = ilog2(max_key)+1;
  size_t keyshift = (keybits<=8) ? 0 : keybits-8;
  size_t nkeys = min<size_t>(size_t(1)<<keybits, 256);
  execParallel(nval, nthreads, [&](size_t tid, size_t lo, size_t hi)
    {
    auto &mybuf(numbers[tid].v);
    mybuf.resize(nkeys,0);
    for (size_t i=lo; i<hi; ++i)
      {
      MR_assert(keys[i]<=max_key, "key too large");
      ++mybuf[(keys[i]>>keyshift)];
      }
    });
  size_t ofs=0;
  for (size_t i=0; i<numbers[0].v.size(); ++i)
    for (size_t t=0; t<nthreads; ++t)
      {
      auto tmp=numbers[t].v[i];
      numbers[t].v[i]=ofs;
      ofs+=tmp;
      }
  quick_array<Tkey> keys2(nval);
  execParallel(nval, nthreads, [&](size_t tid, size_t lo, size_t hi)
    {
    auto &mybuf(numbers[tid].v);
    for (size_t i=lo; i<hi; ++i)
      {
      auto loc = (keys[i]>>keyshift);
      res[mybuf[loc]] = i;
      keys2[mybuf[loc]] = keys[i];
      ++mybuf[loc];
      }
    });
  if (keyshift==0) return;
  keybits -= 8;
  execDynamic(nkeys, nthreads, 1, [&](Scheduler &sched)
    {
    vector<Tidx> newnumbers;
    quick_array<Tkey> keybak;
    quick_array<Tidx> idxbak;
    while (auto rng=sched.getNext())
      for(auto i=rng.lo; i<rng.hi; ++i)
        subsort(res, keys2, keybits, (i==0) ? 0 : numbers[nthreads-1].v[i-1],
          numbers[nthreads-1].v[i], newnumbers, idxbak, keybak);
    });
  }

/// Efficient bucket sort to determine an ascending order for the given set
/// of keys.
/** \param keys
 *         array of keys to be processed in ascending order.
 *         NOTE: this will be changed unpredictably on output!
 *  \param idx
 *         output array to hold the computed indices. Will be resized if needed.
 *  \param max_key
 *         An upper limit for the maximum key value. The closer this is to the
 *         actual maximum key, the better the algorithm will perform.
 *  \param nthreads
 *         The number of threads to use for the sorting
 */
template<typename Tidx, typename Tkey> void bucket_sort2
  (quick_array<Tkey> &keys, quick_array<Tidx> &idx, size_t max_key,
   size_t nthreads)
  {
  auto nval = keys.size();
  idx.realloc(nval);
  nthreads = min(nthreads, thread_pool_size());
  auto sizelimit = max<Tidx>(1, nval/nthreads);
  // align members with cache lines
  struct alignas(64) vbuf { vector<Tidx> v; };
  vector<vbuf> numbers(nthreads);
  auto keybits = ilog2(max_key)+1;
  bool last_level = keybits<=8;
  size_t keyshift = last_level ? 0 : keybits-8;
  size_t nkeys = min<size_t>(size_t(1)<<keybits, 256);
  execParallel(nval, nthreads, [&](size_t tid, size_t lo, size_t hi)
    {
    auto &mybuf(numbers[tid].v);
    mybuf.resize(nkeys,0);
    for (size_t i=lo; i<hi; ++i)
      {
      MR_assert(keys[i]<=max_key, "key too large");
      ++mybuf[(keys[i]>>keyshift)];
      }
    });
  size_t ofs=0;
  for (size_t i=0; i<numbers[0].v.size(); ++i)
    for (size_t t=0; t<nthreads; ++t)
      {
      auto tmp=numbers[t].v[i];
      numbers[t].v[i]=ofs;
      ofs+=tmp;
      }
  if (last_level)
    {
    execParallel(nval, nthreads, [&](size_t tid, size_t lo, size_t hi)
      {
      auto &mybuf(numbers[tid].v);
      for (size_t i=lo; i<hi; ++i)
        {
        auto loc = keys[i];
        idx[mybuf[loc]] = i;
        ++mybuf[loc];
        }
      });
    return;
    }
  quick_array<Tkey> keys2(nval);
  quick_array<Tidx> idx2(nval);
  execParallel(nval, nthreads, [&](size_t tid, size_t lo, size_t hi)
    {
    auto &mybuf(numbers[tid].v);
    for (size_t i=lo; i<hi; ++i)
      {
      auto loc = (keys[i]>>keyshift);
      idx2[mybuf[loc]] = i;
      keys2[mybuf[loc]] = keys[i];
      ++mybuf[loc];
      }
    });
  keybits -= 8;

  struct Workitem
    {
    size_t lo, hi;
    size_t keybits;
    bool src_in_2;
    };

  vector<Workitem> items;
  items.reserve(nkeys);
  for (size_t i=0; i<nkeys; ++i)
    {
    auto lo = (i==0) ? Tidx(0) : numbers[nthreads-1].v[i-1];
    auto hi = numbers[nthreads-1].v[i];
    if (hi-lo>1)
      items.emplace_back(Workitem{lo, hi, keybits, true});
    else if (hi-lo==1)
      idx[lo] = idx2[lo];
    }

  function<void(const Workitem &, function<void(const Workitem &)>&)> do_work = [&keys,&keys2,&idx,&idx2](const Workitem &item, function<void(const Workitem &)> &ifunc)
    {
    const auto &keys_in(item.src_in_2 ? keys2 : keys);
    const auto &idx_in(item.src_in_2 ? idx2 : idx);
    auto &keys_out(item.src_in_2 ? keys : keys2);
    auto &idx_out(item.src_in_2 ? idx : idx2);

    auto lo = item.lo;
    auto hi = item.hi;
    auto keybits = item.keybits;
    auto src_in_2 = item.src_in_2;

    auto nval = hi-lo;
    if (nval<=1)
      {
      if (src_in_2 && (nval==1))
        idx[lo] = idx2[lo];
      return;
      }

    bool last_level = keybits<=8;
    size_t keyshift = last_level ? 0 : keybits-8;
    size_t nkeys = min<size_t>(size_t(1)<<keybits, 256);
    size_t keymask = nkeys-1;
    array<Tidx, 256> numbers;
    for (size_t i=0; i<nkeys; ++i) numbers[i]=0;
    bool all_keys_equal=true;
    bool only_one_bin=true;
    for (size_t i=0; i<nval; ++i)
      {
      all_keys_equal =  all_keys_equal && (keys_in[i+lo]==keys_in[lo]);
      only_one_bin = only_one_bin && ((keys_in[i+lo]>>keyshift)==(keys_in[lo]>>keyshift));
      ++numbers[(keys_in[i+lo]>>keyshift)&keymask];
      }
    if (all_keys_equal)  // this subrange is completely sorted
      {
      if (item.src_in_2)
        for (size_t i=lo; i<hi; ++i)
          idx[i] = idx2[i];
      return;
      }

    if (only_one_bin)  // no reordering needed in this step
      {
      ifunc({item.lo, item.hi, item.keybits-8, item.src_in_2});
      return;
      }

    Tidx ofs=0;
    for (size_t i=0; i<nkeys; ++i)
      {
      auto tmp = numbers[i];
      numbers[i] = ofs;
      ofs += tmp;
      }

    if (last_level)  // don't need to do all the work
      {
      for (size_t i=0; i<nval; ++i)
        {
        auto loc = (keys_in[lo+i]>>keyshift)&keymask;
        idx_out[lo+numbers[loc]] = idx_in[lo+i];
        ++numbers[loc];
        }
      if (!item.src_in_2)
        for (size_t i=lo; i<hi; ++i)
          idx[i] = idx2[i];
      return;
      }

    for (size_t i=0; i<nval; ++i)
      {
      auto loc = (keys_in[lo+i]>>keyshift)&keymask;
      keys_out[lo+numbers[loc]] = keys_in[lo+i];
      idx_out[lo+numbers[loc]] = idx_in[lo+i];
      ++numbers[loc];
      }

    for (size_t i=0; i<nkeys; ++i)
      {
      auto llo = (i==0) ? Tidx(0) : numbers[i-1];
      auto lhi = numbers[i];
      auto lsz = lhi-llo;
      if (lsz<=1)
        {
        if ((!item.src_in_2) && (lsz==1))
          idx[item.lo+llo] = idx2[item.lo+llo];
        }
      else
        ifunc(Workitem{item.lo+llo, item.lo+lhi, item.keybits-8, !item.src_in_2});
      }
    };

  function<void(const Workitem &)> subprocess = [&subprocess, &do_work](const Workitem &item)
    {
    do_work(item, subprocess);
    };

  execWorklist(nthreads, items, [sizelimit, &subprocess, &do_work](const Workitem &item, auto insert)
    {
    auto ifunc = (item.hi-item.lo>sizelimit) ? insert : subprocess;
    do_work(item, ifunc);
    });
  }

}

using detail_bucket_sort::bucket_sort;
using detail_bucket_sort::bucket_sort2;

}

#endif
