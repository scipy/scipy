/*
 *  This file is part of the MR utility library.
 *
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

/*
 *  Copyright (C) 2009-2020 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_COMMUNICATION_H
#define DUCC0_COMMUNICATION_H

//#define DUCC0_USE_MPI

#include <cstddef>
#include <typeindex>
#include <vector>
#ifdef DUCC0_USE_MPI
#include <mpi.h>
#endif

#include "ducc0/infra/types.h"
#include "ducc0/infra/mav.h"

namespace ducc0 {

namespace detail_communication {

using namespace std;

class Communication
  {
  public:
    static void init();
    static bool initialized();
    static void finalize();
    static void abort();
  };

class Communicator
  {
  public:
    enum redOp { Sum, Min, Max, Prod };
#ifdef DUCC0_USE_MPI
    using CommType = MPI_Comm;
#else
    using CommType = struct{};
#endif

  private:
    CommType comm_;
    int rank_, num_ranks_;

    void sendrecvRawVoid (const void *sendbuf, size_t sendcnt,
      size_t dest, void *recvbuf, size_t recvcnt, size_t src, type_index type) const;
    void sendrecv_replaceRawVoid (void *data, type_index type, size_t num,
      size_t dest, size_t src) const;
    void allreduceRawVoid (const void *in, void *out, type_index type, size_t num,
      redOp op) const;
    void allgatherRawVoid (const void *in, void *out, type_index type, size_t num)
      const;
    void allgathervRawVoid (const void *in, int numin, void *out,
      const int *numout, const int *disout, type_index type) const;
    /*! NB: \a num refers to the <i>total</i> number of items in the arrays;
        the individual message size is \a num/num_ranks(). */
    void all2allRawVoid (const void *in, void *out, type_index type, size_t num) const;
    void all2allvRawVoid (const void *in, const int *numin, const int *disin,
      void *out, const int *numout, const int *disout, type_index type) const;
    void bcastRawVoid (void *data, type_index type, size_t num, int root) const;
    void redistributeRawVoid (const fmav_info &iin, const void *in,
                              const fmav_info &iout, void *out,
                              size_t axin, size_t axout, type_index type) const;

  public:
    Communicator();
    Communicator(CommType comm);
    ~Communicator();
    Communicator(const Communicator &other) = default;

    int num_ranks() const { return num_ranks_; }
    int rank() const { return rank_; }
    bool master() const { return rank_==0; }
    CommType comm() const { return comm_; }

    void barrier() const;

    Communicator split(size_t subgroup) const;

    template<typename T> void sendrecvRaw (const T *sendbuf, size_t sendcnt,
      size_t dest, T *recvbuf, size_t recvcnt, size_t src) const
      {
      sendrecvRawVoid(sendbuf, sendcnt, dest, recvbuf, recvcnt, src, tidx<T>());
      }
    template<typename T> void sendrecv_replaceRaw (T *data, size_t num,
      size_t dest, size_t src) const
      { sendrecv_replaceRawVoid(data, tidx<T>(), num, dest, src); }
    template<typename T> void allreduceRaw (const T *in, T *out, size_t num,
      redOp op) const
      { allreduceRawVoid (in, out, tidx<T>(), num, op); }
    template<typename T> void allgatherRaw (const T *in, T *out, size_t num)
      const
      { allgatherRawVoid (in, out, tidx<T>(), num); }
    template<typename T> void allgathervRaw (const T *in, int numin, T *out,
      const int *numout, const int *disout) const
      { allgathervRawVoid (in, numin, out, numout, disout, tidx<T>()); }
    template<typename T> vector<T> allgatherVec (const T &in) const
      {
      vector<T> res(num_ranks_);
      allgatherRaw(&in, res.data(), 1);
      return res;
      }

    template<typename T> T allreduce(const T &in, redOp op) const
      {
      T out;
      allreduceRaw (&in, &out, 1, op);
      return out;
      }
    template<typename T> std::vector<T> allreduceVec
      (const std::vector<T> &in, redOp op) const
      {
      std::vector<T> out(in.size());
      allreduceRaw (in.data(), out.data(), in.size(), op);
      return out;
      }
    template<typename T> void sendrecvVec(const vector<T> &sendbuf, size_t dest,
      vector<T> &recvbuf, size_t src) const
      {
      sendrecvRaw(sendbuf.data(), sendbuf.size(), dest, recvbuf.data(), recvbuf.size(), src);
      }
    /*! NB: \a num refers to the <i>total</i> number of items in the arrays;
        the individual message size is \a num/num_ranks(). */
    template<typename T> void all2allRaw (const T *in, T *out, size_t num) const
      { all2allRawVoid (in, out, tidx<T>(), num); }

    template<typename T> void all2allvRaw (const T *in, const int *numin,
      const int *disin, T *out, const int *numout, const int *disout) const
      { all2allvRawVoid (in,numin,disin,out,numout,disout,tidx<T>()); }

    template<typename T> void bcastRaw (T *data, size_t num, int root=0) const
      { bcastRawVoid (data, tidx<T>(), num, root); }

    template<typename T> void redistribute (const cfmav<T> &in,
      const vfmav<T> &out, size_t axin, size_t axout) const
      {
      if (num_ranks()==1)
        mav_apply([](const T &i, T &o){ o=i; }, 1, in, out);
      else
        redistributeRawVoid(in, in.data(), out, out.data(), axin, axout, tidx<T>());
      }
  };

}

using detail_communication::Communication;
using detail_communication::Communicator;

}

#endif
