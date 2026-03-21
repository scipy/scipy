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

/*! \file geom_utils.h
 *  Geometric utility functions.
 *
 *  Copyright (C) 2003-2020 Max-Planck-Society
 *  \author Martin Reinecke
 *  \author Reinhard Hell
 */

#ifndef DUCC0_GEOM_UTILS_H
#define DUCC0_GEOM_UTILS_H

#include <vector>
#include <cmath>
#include "ducc0/math/math_utils.h"
#include "ducc0/math/vec3.h"

namespace ducc0 {

namespace detail_geom_utils {

/*! Returns the orientation when looking from point \a loc on the unit
    sphere in the direction \a dir. \a loc must be normalized. The result
    ranges from -pi to pi, is 0 for North and pi/2 for West, i.e. the angle
    is given in mathematically positive sense.

    If \a loc is the North or South pole, the returned angle is
    \a atan2(dir.y,dir.x). */
inline double orientation (const vec3 &loc, const vec3 &dir)
  {
// FIXME: here is still optimization potential
  if (loc.x==0 && loc.y==0)
    return (loc.z>0) ? safe_atan2(dir.y,-dir.x) : safe_atan2(dir.y,dir.x);
  vec3 east (-loc.y, loc.x, 0);
  vec3 north = crossprod(loc,east);
  return safe_atan2(-dotprod(dir,east),dotprod(dir,north));
  }

/*! Returns the angle between \a v1 and \a v2 in radians. */
inline double v_angle (const vec3 &v1, const vec3 &v2)
  {
  using namespace std;
  return atan2 (crossprod(v1,v2).Length(), dotprod(v1,v2));
  }

/*! Returns the cosine of the angle between the two points on the sphere defined
    by (\a z1, \a phi1) and (\a z2, \a phi2), respectively. \a z is the cosine
    of the colatitude, and \a phi is the longitude. */
inline double cosdist_zphi (double z1, double phi1, double z2, double phi2)
  {
  using namespace std;
  return z1*z2+cos(phi1-phi2)*std::sqrt((1.-z1*z1)*(1.-z2*z2));
  }

/*! Finds the smallest enclosing cone for a point set on the sphere according to
    Barequet & Elber: Information Processing Letters 93(2005), p.83.
    All points are expected to be passed as unit vectors.
    The enclosing cone must have an opening angle <pi/2. */
void find_enclosing_circle (const std::vector<vec3> &point, vec3 &center,
  double &cosrad);

}

using detail_geom_utils::orientation;
using detail_geom_utils::v_angle;
using detail_geom_utils::cosdist_zphi;
using detail_geom_utils::find_enclosing_circle;

}

#endif
