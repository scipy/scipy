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

/*! \file ducc0/math/quaternion.h
 *  Quaternion class for rotation transforms in 3D space
 *
 *  \copyright Copyright (C) 2005-2020 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_QUATERNION_H
#define DUCC0_QUATERNION_H

#include <cmath>
#include "ducc0/math/vec3.h"

namespace ducc0 {

namespace detail_quaternion {

using namespace std;

/*! \defgroup quaterniongroup Quaternions */
/*! \{ */

/*! Quaternion class for rotation transforms in 3D space */
template<typename T> class quaternion_t
  {
  public:
    T x, y, z, w;

    /// Creates an uninitialized quaternion.
    quaternion_t() {}
    /// Creates a quaternion with the components X, Y, Z, and W.
    quaternion_t(T X, T Y, T Z, T W)
      : x(X), y(Y), z(Z), w(W) {}
    /// Creates a quaternion which is equivalent to a rotation characterized by
    /// \a axis and \a angle.
    quaternion_t(const vec3_t<T> &axis, T angle)
      {
      angle*=T(0.5);
      T sa=sin(angle);
      x = sa*axis.x;
      y = sa*axis.y;
      z = sa*axis.z;
      w = cos(angle);
      }

    /// Returns the norm of \a *this.
    T norm() const
      { return x*x+y*y+z*z+w*w; }

    /// Returns \a *this, divided by its norm.
    quaternion_t normalized() const
      {
      T fct = sqrt(T(1)/norm());
      return quaternion_t(x*fct, y*fct, z*fct, w*fct);
      }

    /// Returns \a -(*this)
    quaternion_t operator-() const
      { return quaternion_t(-x,-y,-z,-w); }

    /// Flips the signs of all components.
    void flip()
      { x=-x; y=-y; z=-z; w=-w; }

    /// Returns the conjugate of \a *this.
    quaternion_t conj() const
      { return quaternion_t(-x,-y,-z, w); }

    /// Returns the inverse of \a *this.
    quaternion_t inverse() const
      {
      T fct = T(1)/norm();
      return quaternion_t(-x*fct, -y*fct, -z*fct, w*fct);
      }

    /// Replaces \a *this by \a (*this)*b.
    quaternion_t &operator*= (const quaternion_t &b)
      {
      quaternion_t a=*this;
      x = a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y;
      y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x;
      z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w;
      w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z;
      return *this;
      }

    /// Returns \a (*this)*b.
    quaternion_t operator* (const quaternion_t &b) const
      {
      return quaternion_t (w*b.x + x*b.w + y*b.z - z*b.y,
                           w*b.y - x*b.z + y*b.w + z*b.x,
                           w*b.z + x*b.y - y*b.x + z*b.w,
                           w*b.w - x*b.x - y*b.y - z*b.z);
      }

    /// Returns the dot product of \a (*this) and \a other.
    T dot(const quaternion_t &other) const
      { return x*other.x + y*other.y + z*other.z + w*other.w; }

    /// Returns a tuple of a rotation axis and an angle which are equivalent
    /// to \a (*this).
    auto toAxisAngle() const
      {
      T norm = x*x + y*y + z*z;
      if (norm==T(0))
        return make_tuple(vec3_t<T>(0,0,1), T(0));
      norm = sqrt(norm);
      T inorm = T(1)/norm;
      return make_tuple(vec3_t<T>(x*inorm,y*inorm,z*inorm), 2*atan2(norm, w));
      }
  };

/*! \} */

}

using detail_quaternion::quaternion_t;

}

#endif
