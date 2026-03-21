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

/*
 *  Copyright (C) 2020 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <cstdlib>
#include <cstring>
#include <numeric>
#include <unordered_map>
#include "ducc0/infra/communication.h"
#include "ducc0/infra/error_handling.h"

namespace ducc0 {

namespace detail_communication {

using namespace std;

void assert_unequal (const void *a, const void *b)
  { MR_assert (a!=b, "input and output buffers must not be identical"); }

#ifdef DUCC0_USE_MPI

class Typemap: public TypeMapper<MPI_Datatype>
  {
  public:
    Typemap()
      {
      add<double>(MPI_DOUBLE);
      add<float>(MPI_FLOAT);
      add<int>(MPI_INT);
      add<long>(MPI_LONG);
      add<char>(MPI_CHAR);
      add<unsigned char>(MPI_BYTE);
      // etc.
      }
   };

Typemap typemap;

MPI_Datatype ndt2mpi (type_index type)
  { return typemap[type]; }

MPI_Op op2mop (Communicator::redOp op)
  {
  switch (op)
    {
    case Communicator::Min : return MPI_MIN;
    case Communicator::Max : return MPI_MAX;
    case Communicator::Sum : return MPI_SUM;
    case Communicator::Prod: return MPI_PROD;
    default: MR_fail ("unsupported reduction operation");
    }
  }

//static
void Communication::init()
  {
  MPI_Init(0,0);
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
  }

//static
bool Communication::initialized()
  {
  int flag=0;
  MPI_Initialized(&flag);
  return flag;
  }

//static
void Communication::finalize()
  { MPI_Finalize(); }

//static
void Communication::abort()
  {
  if (initialized())
    MPI_Abort(MPI_COMM_WORLD, 1);
  else
    exit(1);
  }

Communicator::Communicator(CommType comm)
  : comm_(comm)
  {
  MPI_Comm_size(comm_, &num_ranks_);
  MPI_Comm_rank(comm_, &rank_);
  }

Communicator::Communicator()
  : Communicator(MPI_COMM_WORLD) {}

Communicator::~Communicator()
  {
  if (comm_!=MPI_COMM_WORLD)
     MPI_Comm_free(&comm_);
  }

void Communicator::barrier() const
  { MPI_Barrier(comm_); }

Communicator Communicator::split(size_t color) const
  {
  MPI_Comm comm;
  MPI_Comm_split (comm_, color, rank_, &comm);
  return Communicator(comm);
  }

void Communicator::sendrecvRawVoid (const void *sendbuf, size_t sendcnt,
  size_t dest, void *recvbuf, size_t recvcnt, size_t src, type_index type) const
  {
  if ((sendcnt>0)&&(recvcnt>0)) assert_unequal(sendbuf,recvbuf);

  MPI_Datatype dtype = ndt2mpi(type);
  MPI_Sendrecv (const_cast<void *>(sendbuf),sendcnt,dtype,dest,0,
    recvbuf,recvcnt,dtype,src,0,comm_,MPI_STATUS_IGNORE);
  }
void Communicator::sendrecv_replaceRawVoid (void *data, type_index type, size_t num,
  size_t dest, size_t src) const
  {
  MPI_Sendrecv_replace (data,num,ndt2mpi(type),dest,0,src,0,comm_,
    MPI_STATUS_IGNORE);
  }

void Communicator::allreduceRawVoid (const void *in, void *out, type_index type,
  size_t num, redOp op) const
  {
  void *in2 = (in==out) ? MPI_IN_PLACE : const_cast<void *>(in);
  MPI_Allreduce (in2,out,num,ndt2mpi(type),op2mop(op),comm_);
  }
void Communicator::allgatherRawVoid (const void *in, void *out, type_index type,
  size_t num) const
  {
  if (num>0) assert_unequal(in,out);
  MPI_Datatype tp = ndt2mpi(type);
  MPI_Allgather (const_cast<void *>(in),num,tp,out,num,tp,comm_);
  }
void Communicator::allgathervRawVoid (const void *in, int numin, void *out,
  const int *numout, const int *disout, type_index type) const
  {
  if (numin>0) assert_unequal(in,out);
  MR_assert(numin==numout[rank_],"inconsistent arguments");
  MPI_Datatype tp = ndt2mpi(type);
  MPI_Allgatherv (const_cast<void *>(in),numin,tp,out,const_cast<int *>(numout),
    const_cast<int *>(disout),tp,comm_);
  }
void Communicator::all2allRawVoid (const void *in, void *out, type_index type,
  size_t num) const
  {
  void *in2 = (in==out) ? MPI_IN_PLACE : const_cast<void *>(in);
  MR_assert (num%num_ranks_==0,
    "array size is not divisible by number of ranks");
  MPI_Datatype tp = ndt2mpi(type);
  MPI_Alltoall (in2,num/num_ranks_,tp,out,num/num_ranks_,tp,comm_);
  }
void Communicator::all2allvRawVoid (const void *in, const int *numin,
  const int *disin, void *out, const int *numout, const int *disout, type_index type)
  const
  {
  long commsz=disin[num_ranks_-1]+numin[num_ranks_-1]
             +disout[num_ranks_-1]+numout[num_ranks_-1];
  if (commsz>0) assert_unequal(in,out);
  MPI_Datatype tp = ndt2mpi(type);
  MPI_Alltoallv (const_cast<void *>(in), const_cast<int *>(numin),
    const_cast<int *>(disin), tp, out, const_cast<int *>(numout),
    const_cast<int *>(disout), tp, comm_);
  }

void Communicator::bcastRawVoid (void *data, type_index type, size_t num, int root) const
  { MPI_Bcast (data,num,ndt2mpi(type),root,comm_); }

MPI_Datatype fmav2mpidt(const fmav_info &info, type_index type)
  {
  size_t ndim = info.ndim();
  auto origtype = ndt2mpi(type);
  if (info.size()==0)
    {
    MPI_Datatype res;
    MPI_Type_contiguous(0, origtype, &res);
    MPI_Type_commit(&res);
    return res;
    }
  auto res=origtype;
  bool free_type=false;
  for (size_t i=0; i<ndim; ++i)
    {
    MPI_Datatype tmptype;
    MPI_Type_vector(info.shape(i), 1, info.stride(i), res, &tmptype);
    MPI_Type_commit(&tmptype);
    if (free_type)
      MPI_Type_free(&res);
    free_type = true;
    MPI_Type_create_resized(tmptype, 0, typesize(type), &res);
    MPI_Type_commit(&res);
    MPI_Type_free(&tmptype);
    }
  return res;
  }

void Communicator::redistributeRawVoid(const fmav_info &iin, const void *in,
  const fmav_info &iout, void *out, size_t axin, size_t axout, type_index type) const
  {
  auto ndim = iin.ndim();
  auto nranks = size_t(num_ranks());
  MR_assert(ndim==iout.ndim(), "array dimensions must be equal");
  MR_assert(axin<ndim, "invalid axin");
  MR_assert(axout<ndim, "invalid axout");
  MR_assert(axin!=axout, "axin and axout muyst be different");
  for (size_t i=0; i<ndim; ++i)
    if ((i!=axin) && (i!=axout))
      MR_assert(iin.shape(i)==iout.shape(i), "shape mismatch");
  auto s_in = allgatherVec(int(iin.shape(axin)));
  MR_assert(int(iout.shape(axin))==reduce(s_in.begin(), s_in.end()), "inconsistency");
  auto s_out = allgatherVec(int(iout.shape(axout)));
  MR_assert(int(iin.shape(axout))==reduce(s_out.begin(), s_out.end()), "inconsistency");

  vector<MPI_Datatype> v_in(nranks), v_out(nranks);
  for (size_t i=0; i<nranks; ++i)
    {
    auto tmp = iin.shape();
    tmp[axout] = s_out[i];
    v_in[i] = fmav2mpidt(fmav_info(tmp, iin.stride()), type);
    tmp = iout.shape();
    tmp[axin] = s_in[i];
    v_out[i] = fmav2mpidt(fmav_info(tmp, iout.stride()), type);
    }

  vector<int> disp_in(nranks), disp_out(nranks);
  for (size_t i=0; i<nranks; ++i)
    {
    disp_in[i] = (i==0) ? 0 : disp_in[i-1]+s_out[i-1]*iin.stride(axout)*typesize(type);
    disp_out[i] = (i==0) ? 0 : disp_out[i-1]+s_in[i-1]*iout.stride(axin)*typesize(type);
    }
  vector<int> num(nranks, 1);

  MPI_Alltoallw(in, num.data(), disp_in.data(), v_in.data(),
                out, num.data(), disp_out.data(), v_out.data(), comm_);

  for (auto &t: v_in) MPI_Type_free(&t);
  for (auto &t: v_out) MPI_Type_free(&t);
  }

#else

//static
void Communication::init() {}

//static
bool Communication::initialized()
  { return true; }

//static
void Communication::finalize() {}

//static
void Communication::abort()
  { exit(1); }

Communicator::Communicator()
  : rank_(0), num_ranks_(1) {}

Communicator::~Communicator() {}

void Communicator::barrier() const {}

Communicator Communicator::split(size_t /*color*/) const
  { return *this; }

void Communicator::sendrecvRawVoid (const void *sendbuf, size_t sendcnt,
  size_t dest, void *recvbuf, size_t recvcnt, size_t src, type_index type) const
  {
  MR_assert ((dest==0) && (src==0), "inconsistent call");
  MR_assert (sendcnt==recvcnt, "inconsistent call");
  if (sendcnt>0) assert_unequal(sendbuf,recvbuf);
  memcpy (recvbuf, sendbuf, sendcnt*typesize(type));
  }
void Communicator::sendrecv_replaceRawVoid (void *, type_index, size_t, size_t dest,
  size_t src) const
  { MR_assert ((dest==0) && (src==0), "inconsistent call"); }

void Communicator::allreduceRawVoid (const void *in, void *out, type_index type,
  size_t num, redOp /*op*/) const
  {
  if (in==out) return;
  memcpy (out, in, num*typesize(type));
  }
void Communicator::allgatherRawVoid (const void *in, void *out, type_index type,
  size_t num) const
  { if (num>0) assert_unequal(in,out); memcpy (out, in, num*typesize(type)); }
void Communicator::all2allRawVoid (const void *in, void *out, type_index type,
  size_t num) const
  {
  if (in==out) return;
  memcpy (out, in, num*typesize(type));
  }
void Communicator::allgathervRawVoid (const void *in, int numin, void *out,
  const int *numout, const int *disout, type_index type) const
  {
  if (numin>0) assert_unequal(in,out);
  MR_assert(numin==numout[0],"inconsistent call");
  memcpy (reinterpret_cast<char *>(out)+disout[0]*typesize(type), in,
    numin*typesize(type));
  }
void Communicator::all2allvRawVoid (const void *in, const int *numin,
  const int *disin, void *out, const int *numout, const int *disout, type_index type)
  const
  {
  if (numin[0]>0) assert_unequal(in,out);
  MR_assert (numin[0]==numout[0],"message size mismatch");
  const char *in2 = static_cast<const char *>(in);
  char *out2 = static_cast<char *>(out);
  size_t st=typesize(type);
  memcpy (out2+disout[0]*st,in2+disin[0]*st,numin[0]*st);
  }

void Communicator::bcastRawVoid (void *, type_index, size_t, int) const
  {}

void Communicator::redistributeRawVoid(const fmav_info &, const void *,
  const fmav_info &, void *, size_t, size_t, type_index) const
  { MR_fail("must not get here"); }

#endif

}}
