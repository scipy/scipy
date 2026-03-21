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

/*! \file rangeset.h
 *  Class for storing sets of ranges of integer numbers
 *
 *  Copyright (C) 2011-2021 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_RANGESET_H
#define DUCC0_RANGESET_H

#include <algorithm>
#include <vector>
#include <cstddef>
#include <iostream>
#include "ducc0/infra/error_handling.h"
#include "ducc0/math/math_utils.h"

namespace ducc0 {

/// Class for storing sets of ranges of integer numbers
template<typename T> class rangeset
  {
  private:
    typedef std::vector<T> rtype;
    typedef typename rtype::iterator iterator;
    typedef typename rtype::const_iterator c_iterator;
    rtype r;

    /// Returns the index of the last number in \c r which is <= \a val.
    /// If \a val is smaller than all numbers in \c r, returns -1.
    ptrdiff_t iiv (const T &val) const
      { return ptrdiff_t(std::upper_bound(r.begin(),r.end(),val)-r.begin())-1; }

    void addRemove (T a, T b, ptrdiff_t v)
      {
      ptrdiff_t pos1=iiv(a), pos2=iiv(b);
      if ((pos1>=0) && (r[size_t(pos1)]==a)) --pos1;
      // first to delete is at pos1+1; last is at pos2
      bool insert_a = (pos1&1)==v;
      bool insert_b = (pos2&1)==v;
      ptrdiff_t rmstart=pos1+1+(insert_a ? 1 : 0);
      ptrdiff_t rmend  =pos2-(insert_b ? 1 : 0);

      MR_assert((rmend-rmstart)&1,"cannot happen");

      if (insert_a && insert_b && (pos1+1>pos2)) // insert
        {
        r.insert(r.begin()+pos1+1,2,a);
        r[size_t(pos1+2)]=b;
        }
      else
        {
        if (insert_a) r[size_t(pos1+1)]=a;
        if (insert_b) r[size_t(pos2)]=b;
        r.erase(r.begin()+rmstart,r.begin()+rmend+1);
        }
      }

    /*! Estimate a good strategy for set operations involving two rangesets. */
    static int strategy (size_t sza, size_t szb)
      {
      const double fct1 = 1.;
      const double fct2 = 1.;
      size_t slo = sza<szb ? sza : szb,
            shi = sza<szb ? szb : sza;
      double cost1 = fct1 * (sza+szb);
      double cost2 = fct2 * slo * std::max(1u,ilog2(shi));
      return (cost1<=cost2) ? 1 : (slo==sza) ? 2 : 3;
      }

    void generalUnion1 (const rangeset &a, const rangeset &b,
      bool flip_a, bool flip_b)
      {
      bool state_a=flip_a, state_b=flip_b, state_res=state_a||state_b;
      size_t ia=0, ea=a.r.size(), ib=0, eb=b.r.size();
      bool runa = ia!=ea, runb = ib!=eb;
      while(runa||runb)
        {
        T va = runa ? a.r[ia] : T(0),
          vb = runb ? b.r[ib] : T(0);
        bool adv_a = runa && (!runb || (va<=vb)),
             adv_b = runb && (!runa || (vb<=va));
        if (adv_a) { state_a=!state_a; ++ia; runa = ia!=ea; }
        if (adv_b) { state_b=!state_b; ++ib; runb = ib!=eb; }
        if ((state_a||state_b)!=state_res)
          { r.push_back(adv_a ? va : vb); state_res = !state_res; }
        }
      }
    void generalUnion2 (const rangeset &a, const rangeset &b,
      bool flip_a, bool flip_b)
      {
      ptrdiff_t iva = flip_a ? 0 : -1;
      ptrdiff_t asz=ptrdiff_t(a.r.size()), bsz=ptrdiff_t(b.r.size());
      while (iva<asz)
        {
        ptrdiff_t ivb = (iva==-1) ? -1 : b.iiv(a.r[iva]);
        bool state_b = flip_b^((ivb&1)==0);
        if ((iva>-1) && (!state_b)) r.push_back(a.r[iva]);
        while((ivb<bsz-1)&&((iva==asz-1)||(b.r[ivb+1]<a.r[iva+1])))
          { ++ivb; state_b=!state_b; r.push_back(b.r[ivb]); }
        if ((iva<asz-1)&&(!state_b)) r.push_back(a.r[iva+1]);
        iva+=2;
        }
      }
    void generalUnion (const rangeset &a, const rangeset &b,
      bool flip_a, bool flip_b)
      {
      MR_assert((this!=&a)&&(this!=&b), "cannot overwrite the rangeset");
      if (a.r.empty())
        {
        if (flip_a) clear(); else *this=b;
        return;
        }
      if (b.r.empty())
        {
        if (flip_b) clear(); else *this=a;
        return;
        }

      clear();
      int strat = strategy (a.nranges(), b.nranges());
      (strat==1) ? generalUnion1(a,b,flip_a,flip_b) :
        ((strat==2) ? generalUnion2(a,b,flip_a,flip_b)
                    : generalUnion2(b,a,flip_b,flip_a));
      }
    void generalXor (const rangeset &a, const rangeset &b)
      {
      clear();
      bool state_a=false, state_b=false, state_res=state_a||state_b;
      size_t ia=0, ea=a.r.size(), ib=0, eb=b.r.size();
      bool runa = ia!=ea, runb = ib!=eb;
      while(runa||runb)
        {
        T va = runa ? a.r[ia] : T(0),
          vb = runb ? b.r[ib] : T(0);
        bool adv_a = runa && (!runb || (va<=vb)),
             adv_b = runb && (!runa || (vb<=va));
        if (adv_a) { state_a=!state_a; ++ia; runa = ia!=ea; }
        if (adv_b) { state_b=!state_b; ++ib; runb = ib!=eb; }
        if ((state_a^state_b)!=state_res)
          { r.push_back(adv_a ? va : vb); state_res = !state_res; }
        }
      }

    static bool generalAllOrNothing1 (const rangeset &a, const rangeset &b,
      bool flip_a, bool flip_b)
      {
      bool state_a=flip_a, state_b=flip_b, state_res=state_a||state_b;
      size_t ia=0, ea=a.r.size(), ib=0, eb=b.r.size();
      bool runa = ia!=ea, runb = ib!=eb;
      while(runa||runb)
        {
        T va = runa ? a.r[ia] : T(0),
          vb = runb ? b.r[ib] : T(0);
        bool adv_a = runa && (!runb || (va<=vb)),
             adv_b = runb && (!runa || (vb<=va));
        if (adv_a) { state_a=!state_a; ++ia; runa = ia!=ea; }
        if (adv_b) { state_b=!state_b; ++ib; runb = ib!=eb; }
        if ((state_a||state_b)!=state_res)
          return false;
        }
      return true;
      }
    static bool generalAllOrNothing2 (const rangeset &a, const rangeset &b,
      bool flip_a, bool flip_b)
      {
      ptrdiff_t iva = flip_a ? 0 : -1;
      ptrdiff_t asz=ptrdiff_t(a.r.size()), bsz=ptrdiff_t(b.r.size());
      while (iva<asz)
        {
        if (iva==-1) // implies that flip_a==false
          { if ((!flip_b)||(b.r[0]<a.r[0])) return false; }
        else if (iva==asz-1) // implies that flip_a==false
          { if ((!flip_b)||(b.r[bsz-1]>a.r[asz-1])) return false; }
        else
          {
          ptrdiff_t ivb=b.iiv(a.r[iva]);
          if ((ivb!=bsz-1)&&(b.r[ivb+1]<a.r[iva+1])) return false;
          if (flip_b==((ivb&1)==0)) return false;
          }
        iva+=2;
        }
      return true;
      }
    static bool generalAllOrNothing (const rangeset &a, const rangeset &b,
      bool flip_a, bool flip_b)
      {
      if (a.r.empty())
        return flip_a ? true : b.r.empty();
      if (b.r.empty())
        return flip_b ? true : a.r.empty();
      int strat = strategy (a.nranges(), b.nranges());
      return (strat==1) ? generalAllOrNothing1(a,b,flip_a,flip_b) :
               ((strat==2) ? generalAllOrNothing2(a,b,flip_a,flip_b)
                           : generalAllOrNothing2(b,a,flip_b,flip_a));
      }

  public:
    /// Removes all entries.
    void clear() { r.clear(); }
    /// Reserves space for \a n ranges.
    void reserve(size_t n) { r.reserve(2*n); }
    /// Returns the current number of ranges.
    size_t nranges() const { return r.size()>>1; }
    /// Returns the current number of ranges.
    size_t size() const { return nranges(); }
    /// Returns `true` iff there are no ranges in the rangeset.
    bool empty() const { return r.empty(); }
    /// Returns the current vector of ranges.
    const rtype &data() const { return r; }
    void checkConsistency() const
      {
      MR_assert((r.size()&1)==0,"invalid number of entries");
      for (size_t i=1; i<r.size(); ++i)
        MR_assert(r[i]>r[i-1],"inconsistent entries");
      }
    void setData (const rtype &inp)
      {
      r=inp;
      checkConsistency();
      }

    /// Returns the first value of range \a i.
    const T &ivbegin (size_t i) const { return r[2*i]; }
    /// Returns the one-past-last value of range \a i.
    const T &ivend (size_t i) const { return r[2*i+1]; }
    /// Returns the length of range \a i.
    T ivlen (size_t i) const { return r[2*i+1]-r[2*i]; }

    /// Appends \a [v1;v2[ to the rangeset. \a v1 must be larger
    /// than the minimum of the last range in the rangeset.
    void append(const T &v1, const T &v2)
      {
      if (v2<=v1) return;
      if ((!r.empty()) && (v1<=r.back()))
        {
        MR_assert (v1>=r[r.size()-2],"bad append operation");
        if (v2>r.back()) r.back()=v2;
        }
      else
        { r.push_back(v1); r.push_back(v2); }
      }
    /// Appends \a [v;v+1[ to the rangeset. \a v must be larger
    /// than the minimum of the last range in the rangeset.
    void append(const T &v)
      { append(v,v+1); }

    /// Appends \a other to the rangeset. All values in \a other must be larger
    /// than the minimum of the last range in the rangeset.
    void append (const rangeset &other)
      {
      for (size_t j=0; j<other.nranges(); ++j)
        append(other.ivbegin(j),other.ivend(j));
      }

    /// After this operation, the rangeset contains the union of itself
    /// with \a [v1;v2[.
    void add(const T &v1, const T &v2)
      {
      if (v2<=v1) return;
      if (r.empty() || (v1>=r[r.size()-2])) append(v1,v2);
      addRemove(v1,v2,1);
      }
    /// After this operation, the rangeset contains the union of itself
    /// with \a [v;v+1[.
    void add(const T &v) { add(v,v+1); }

    /// Removes all values within \a [v1;v2[ from the rangeset.
    void remove(const T &v1, const T &v2)
      {
      if (v2<=v1) return;
      if (r.empty()) return;
      if ((v2<=r[0])||(v1>=r.back())) return;
      if ((v1<=r[0]) && (v2>=r.back())) { r.clear(); return; }
      addRemove(v1,v2,0);
      }
    /// Removes the value \a v from the rangeset.
    void remove(const T &v) { remove(v,v+1); }

    /// Removes all values not within \a [a;b[ from the rangeset.
    void intersect (const T &a, const T &b)
      {
      if (r.empty()) return; // nothing to remove
      if ((b<=r[0]) || (a>=r.back())) { r.clear(); return; } // no overlap
      if ((a<=r[0]) && (b>=r.back())) return; // full rangeset in interval

      ptrdiff_t pos2=iiv(b);
      if ((pos2>=0) && (r[size_t(pos2)]==b)) --pos2;
      bool insert_b = (pos2&1)==0;
      r.erase(r.begin()+pos2+1,r.end());
      if (insert_b) r.push_back(b);

      ptrdiff_t pos1=iiv(a);
      bool insert_a = (pos1&1)==0;
      if (insert_a) r[size_t(pos1--)]=a;
      if (pos1>=0)
        r.erase(r.begin(),r.begin()+pos1+1);
      }

    /// Returns the total number of elements in the rangeset.
    size_t nval() const
      {
      size_t result(0);
      for (size_t i=0; i<r.size(); i+=2)
        result+=size_t(r[i+1]-r[i]);
      return result;
      }

    /// After this operation, \a res contains all elements of the rangeset
    /// in ascending order.
    void toVector (std::vector<T> &res) const
      {
      res.clear();
      res.reserve(nval());
      for (size_t i=0; i<r.size(); i+=2)
        for (T m(r[i]); m<r[i+1]; ++m)
          res.push_back(m);
      }

    /// Returns a vector containing all elements of the rangeset in ascending
    /// order.
    std::vector<T> toVector() const
      {
      std::vector<T> res;
      toVector(res);
      return res;
      }

    /// Returns the union of this rangeset and \a other.
    rangeset op_or (const rangeset &other) const
      {
      rangeset res;
      res.generalUnion (*this,other,false,false);
      return res;
      }
    /// Returns the intersection of this rangeset and \a other.
    rangeset op_and (const rangeset &other) const
      {
      rangeset res;
      res.generalUnion (*this,other,true,true);
      return res;
      }
    /// Returns the part of this rangeset which is not in \a other.
    rangeset op_andnot (const rangeset &other) const
      {
      rangeset res;
      res.generalUnion (*this,other,true,false);
      return res;
      }
    /// Returns the parts of this rangeset and \a other, which are not in
    /// both rangesets.
    rangeset op_xor (const rangeset &other) const
      {
      rangeset res;
      res.generalXor (*this,other);
      return res;
      }

    /// Returns the index of the interval containing \a v; if no such interval
    /// exists, -1 is returned.
    ptrdiff_t findInterval (const T &v) const
      {
      ptrdiff_t res = iiv(v);
      return (res&1) ? -1 : res>>1;
      }

    /// Returns \a true iff the rangeset is identical to \a other.
    bool operator== (const rangeset &other) const
      { return r==other.r; }

    /// Returns \a true iff the rangeset contains all values in the range
    /// \a [a;b[.
    bool contains (T a,T b) const
      {
      ptrdiff_t res=iiv(a);
      if (res&1) return false;
      return (b<=r[res+1]);
      }
    /// Returns \a true iff the rangeset contains the value \a v.
    bool contains (T v) const
      { return !(iiv(v)&1); }
    /// Returns \a true iff the rangeset contains all values stored in \a other.
    bool contains (const rangeset &other) const
      { return generalAllOrNothing(*this,other,false,true); }
    /// Returns true if any of the numbers [a;b[ are contained in the set.
    bool overlaps (T a,T b) const
      {
      ptrdiff_t res=iiv(a);
      if ((res&1)==0) return true;
      if (res==ptrdiff_t(r.size())-1) return false; // beyond the end of the set
      return (r[res+1]<b);
      }
    /// Returns true if there is overlap between the set and \a other.
    bool overlaps (const rangeset &other) const
      { return !generalAllOrNothing(*this,other,true,true); }
  };

template<typename T> inline std::ostream &operator<< (std::ostream &os,
  const rangeset<T> &rs)
  {
  os << "{ ";
  for (size_t i=0; i<rs.nranges(); ++i)
    os << "["<<rs.ivbegin(i)<<";"<<rs.ivend(i)<<"[ ";
  return os << "}";
  }

}

#endif
