/*
 *  This file is part of ducc0.
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

/** \file ducc0/math/wigner3j.h
 *  Computation of Wigner-3j symbols
 *  Algorithm implemented according to Schulten & Gordon:
 *  J. Math. Phys. 16, p. 10 (1975)
 *
 *  Copyright (C) 2009-2023 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_WIGNER3J_H
#define DUCC0_WIGNER3J_H

#include <vector>
#include "ducc0/infra/mav.h"

namespace ducc0 {

namespace detail_wigner3j {

using namespace std;

/**
  Compute the Wigner 3j symbols for parameters \a l2, \a l3, \a m2, \a m3
  following the algorithm in the SLATEC algorithm DRC3JJ.
  The results are returned in \a res, which is resized appropriately.
 */
void wigner3j (double l2, double l3, double m2, double m3, vector<double> &res);
void wigner3j (double l2, double l3, double m2, double m3, const vmav<double,1> &res);

void wigner3j_00_squared_compact (double l2, double l3, const vmav<double,1> &res);

template<typename Tsimd> void wigner3j_00_vec_squared_compact (Tsimd l2, Tsimd l3, const vmav<Tsimd,1> &res);

void flexible_wigner3j (double l2, double l3, double m2, double m3, double l1min, const vmav<double,1> &res);
template<typename Tsimd> void flexible_wigner3j_vec
  (Tsimd l2, Tsimd l3, double m2, double m3, Tsimd l1min, const vmav<Tsimd,1> &res);

/**
  Compute the Wigner 3j symbols for parameters \a l2, \a l3, \a m2, \a m3
  following the algorithm in the SLATEC algorithm DRC3JJ.
  The results are returned in \a res, which is resized appropriately.
  The l1 value of the first entry in \a res is returned in \a l1min.
 */
void wigner3j_int (int l2, int l3, int m2, int m3, int &l1min, vector<double> &res);
void wigner3j_int (int l2, int l3, int m2, int m3, int &l1min, const vmav<double,1> &res);
int wigner3j_ncoef_int(int l2, int l3, int m2, int m3);
}

using detail_wigner3j::wigner3j;
using detail_wigner3j::wigner3j_int;
using detail_wigner3j::wigner3j_ncoef_int;

using detail_wigner3j::wigner3j_00_vec_squared_compact;

using detail_wigner3j::flexible_wigner3j;
using detail_wigner3j::flexible_wigner3j_vec;
}

#endif
