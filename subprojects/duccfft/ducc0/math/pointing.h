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

/*! \file pointing.h
 *  Class representing a direction in 3D space
 *
 *  Copyright (C) 2003-2021 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef DUCC0_POINTING_H
#define DUCC0_POINTING_H

#include <cmath>
#include <iostream>
#include "ducc0/math/vec3.h"

namespace ducc0 {

/*! \defgroup pointinggroup Pointings */
/*! \{ */

/// Class representing a direction in 3D space or a location on the
/// unit sphere. All angles in radians.
class pointing
  {
  public:
    /// Colatitude of the pointing (i.e. the North pole is at \a theta=0).
    double theta;
    /// Longitude of the pointing.
    double phi;

    /// Default constructor. \a theta and \a phi are not initialized.
    pointing() {}
    /// Creates a pointing with \a Theta and \a Phi.
    pointing (double Theta, double Phi) : theta(Theta), phi(Phi) {}

    /// Creates a pointing from the vector \a inp. \a inp need not be
    /// normalized.
    explicit pointing (const vec3 &inp)
      { from_vec3(inp); }
    /// Returns a normalized vector pointing in the same direction.
    explicit operator vec3() const
      { return to_vec3(); }
    /// Returns a normalized vector pointing in the same direction.
    vec3 to_vec3() const;
    /// Converts \a inp to \a ptg. \a inp need not be normalized.
    void from_vec3 (const vec3 &inp);
    /// Changes the angles so that \a 0<=theta<=pi.
    void normalize_theta();
    /// Changes the angles so that \a 0<=theta<=pi and \a 0<=phi<2*pi.
    void normalize();
  };

/// Writes \a p to \a os.
/// \relates pointing
std::ostream &operator<< (std::ostream &os, const pointing &p);

/*! \} */

}

#endif
