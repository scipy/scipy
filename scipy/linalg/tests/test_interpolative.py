#*******************************************************************************
#   Copyright (C) 2013 Kenneth L. Ho
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#
#   Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer. Redistributions in binary
#   form must reproduce the above copyright notice, this list of conditions and
#   the following disclaimer in the documentation and/or other materials
#   provided with the distribution.
#
#   None of the names of the copyright holders may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#   POSSIBILITY OF SUCH DAMAGE.
#*******************************************************************************

# FIXME: this needs to be integrated with scipy's test suite

import sys
import scipy.linalg.interpolative as interpolative
import numpy as np
import time

if __name__ == '__main__':
  """
  Test ID routines on a Hilbert matrix.
  """
  # construct Hilbert matrix
  m = n = 1000
  A = np.empty((m, n), order='F')
  for j in range(n):
    for i in range(m):
      A[i,j] = 1. / (i + j + 1)

  # define multiplication functions
  def matvec (x): return np.dot(A,   x)
  def matveca(x): return np.dot(A.T, x)

  # set relative precision
  eps = 1e-12

  # find true rank
  S = np.linalg.svd(A, compute_uv=False)
  try:    rank = np.nonzero(S < 1e-12)[0][0]
  except: rank = n

  # print input summary
  print "Hilbert matrix dimension:        %8i" % n
  print "Working precision:               %8.2e" % eps
  print "Rank to working precision:       %8i" % rank
  print "-----------------------------------------"

  # set print format
  fmt = "%8.2e (s) / %5s"

  # test ID routines
  print "ID routines"
  print "-----------------------------------------"

  # fixed precision
  print "Calling iddp_id  ...",
  t0 = time.clock()
  k, idx, proj = interpolative.compute_id(A, eps, rand=False)
  t = time.clock() - t0
  B = interpolative.reconstruct_skeleton(A, k, idx)
  C = interpolative.reconstruct_matrix_from_id(B, idx, proj)
  print fmt % (t, np.allclose(A, C, eps))

  print "Calling iddp_aid ...",
  t0 = time.clock()
  k, idx, proj = interpolative.compute_id(A, eps)
  t = time.clock() - t0
  B = interpolative.reconstruct_skeleton(A, k, idx)
  C = interpolative.reconstruct_matrix_from_id(B, idx, proj)
  print fmt % (t, np.allclose(A, C, eps))

  print "Calling iddp_rid ...",
  t0 = time.clock()
  k, idx, proj = interpolative.compute_id(m, n, matveca, eps)
  t = time.clock() - t0
  B = interpolative.reconstruct_skeleton(A, k, idx)
  C = interpolative.reconstruct_matrix_from_id(B, idx, proj)
  print fmt % (t, np.allclose(A, C, eps))

  # fixed rank
  k = rank

  print "Calling iddr_id  ...",
  t0 = time.clock()
  idx, proj = interpolative.compute_id(A, k, rand=False)
  t = time.clock() - t0
  B = interpolative.reconstruct_skeleton(A, k, idx)
  C = interpolative.reconstruct_matrix_from_id(B, idx, proj)
  print fmt % (t, np.allclose(A, C, eps))

  print "Calling iddr_aid ...",
  t0 = time.clock()
  idx, proj = interpolative.compute_id(A, k)
  t = time.clock() - t0
  B = interpolative.reconstruct_skeleton(A, k, idx)
  C = interpolative.reconstruct_matrix_from_id(B, idx, proj)
  print fmt % (t, np.allclose(A, C, eps))

  print "Calling iddr_rid ...",
  t0 = time.clock()
  idx, proj = interpolative.compute_id(m, n, matveca, k)
  t = time.clock() - t0
  B = interpolative.reconstruct_skeleton(A, k, idx)
  C = interpolative.reconstruct_matrix_from_id(B, idx, proj)
  print fmt % (t, np.allclose(A, C, eps))

  # test SVD routines
  print "-----------------------------------------"
  print "SVD routines"
  print "-----------------------------------------"

  # fixed precision
  print "Calling iddp_svd ...",
  t0 = time.clock()
  U, S, V = interpolative.svd(A, eps, rand=False)
  t = time.clock() - t0
  B = np.dot(U, np.dot(np.diag(S), V.T))
  print fmt % (t, np.allclose(A, B, eps))

  print "Calling iddp_asvd...",
  t0 = time.clock()
  U, S, V = interpolative.svd(A, eps)
  t = time.clock() - t0
  B = np.dot(U, np.dot(np.diag(S), V.T))
  print fmt % (t, np.allclose(A, B, eps))

  print "Calling iddp_rsvd...",
  t0 = time.clock()
  U, S, V = interpolative.svd(m, n, matvec, matveca, eps)
  t = time.clock() - t0
  B = np.dot(U, np.dot(np.diag(S), V.T))
  print fmt % (t, np.allclose(A, B, eps))

  # fixed rank
  k = rank

  print "Calling iddr_svd ...",
  t0 = time.clock()
  U, S, V = interpolative.svd(A, k, rand=False)
  t = time.clock() - t0
  B = np.dot(U, np.dot(np.diag(S), V.T))
  print fmt % (t, np.allclose(A, B, eps))

  print "Calling iddr_asvd...",
  t0 = time.clock()
  U, S, V = interpolative.svd(A, k)
  t = time.clock() - t0
  B = np.dot(U, np.dot(np.diag(S), V.T))
  print fmt % (t, np.allclose(A, B, eps))

  print "Calling iddr_rsvd...",
  t0 = time.clock()
  U, S, V = interpolative.svd(m, n, matvec, matveca, k)
  t = time.clock() - t0
  B = np.dot(U, np.dot(np.diag(S), V.T))
  print fmt % (t, np.allclose(A, B, eps))
