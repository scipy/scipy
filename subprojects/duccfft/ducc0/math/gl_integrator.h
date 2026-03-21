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

/** \file ducc0/math/gl_integrator.h
 *  Functionality for Gauss-Legendre quadrature
 *
 *  \copyright Copyright (C) 2019-2023 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_GL_INTEGRATOR_H
#define DUCC0_GL_INTEGRATOR_H

#include <cmath>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>
#include "ducc0/math/constants.h"
#include "ducc0/infra/error_handling.h"

namespace ducc0 {

namespace detail_gl_integrator {

using namespace std;

/// Class for computing Gauss-Legendre abscissas, weights and integrals
class GL_Integrator
  {
  private:
    size_t n_;
    vector<double> x, w, th;

  public:
    /// Creates an integrator for \a n abscissas
    /** \note The \a nthreads parameter is obsolescent and ignored. */
    GL_Integrator(size_t n, size_t /*nthreads*/=1);

    /// Returns the approximated integral of \a func in [-1;1]
    template<typename Func> auto integrate(Func f) -> decltype(f(0.))
      {
      using T = decltype(f(0.));
      T res=0;
      size_t istart=0;
      if (n_&1)
        {
        res = f(x[0])*w[0];
        istart=1;
        }
      for (size_t i=istart; i<x.size(); ++i)
        res += (f(x[i])+f(-x[i]))*w[i];
      return res;
      }

    /// Returns the approximated integral of \a func in [-1;1], where \a func
    /// is symmetric with respect to x=0.
    template<typename Func> auto integrateSymmetric(Func f) -> decltype(f(0.))
      {
      using T = decltype(f(0.));
      T res=f(x[0])*w[0];
      if (n_&1) res *= 0.5;
      for (size_t i=1; i<x.size(); ++i)
        res += f(x[i])*w[i];
      return res*2;
      }

    /// Returns the Gauss-Legendre abscissas.
    vector<double> coords() const
      {
      vector<double> res(n_);
      for (size_t i=0; i<x.size(); ++i)
        {
        res[i]=-x[x.size()-1-i];
        res[n_-1-i] = x[x.size()-1-i];
        }
      return res;
      }
    /// Returns the non-negative Gauss-Legendre abscissas.
    const vector<double> &coordsSymmetric() const
      { return x; }

    /// Returns the Gauss-Legendre weights.
    vector<double> weights() const
      {
      vector<double> res(n_);
      for (size_t i=0; i<w.size(); ++i)
        res[i]=res[n_-1-i]=w[w.size()-1-i];
      return res;
      }
    /// Returns the Gauss-Legendre weights for the non-negative abscissas,
    /// with an additional factor of 2 for positive abscissas.
    vector<double> weightsSymmetric() const
      {
      auto res = w;
      if (n_&1) res[0]*=0.5;
      for (auto &v:res) v*=2;
      return res;
      }

    vector<double> thetas() const
      {
      vector<double> res(n_);
      for (size_t i=0; i<th.size(); ++i)
        {
        res[i]= pi-th[th.size()-1-i];
        res[n_-1-i] = th[th.size()-1-i];
        }
      return res;
      }
  };

}

using detail_gl_integrator::GL_Integrator;

}

#endif
