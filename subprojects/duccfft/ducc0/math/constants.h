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

/*! \file constants.h
 *  Mathematical constants.
 */

#ifndef DUCC0_CONSTANTS_H
#define DUCC0_CONSTANTS_H

namespace ducc0 {

/*! \defgroup mathconstgroup Mathematical constants */
/*! \{ */

/// pi
constexpr double pi=3.141592653589793238462643383279502884197;
/// 2*pi
constexpr double twopi=6.283185307179586476925286766559005768394;
/// 1./(2*pi)
constexpr double inv_twopi=1.0/twopi;
/// 4*pi
constexpr double fourpi=12.56637061435917295385057353311801153679;
/// pi/2
constexpr double halfpi=1.570796326794896619231321691639751442099;
/// 2/pi
constexpr double inv_halfpi=0.6366197723675813430755350534900574;
/// 1/sqrt(4*pi)
constexpr double inv_sqrt4pi = 0.2820947917738781434740397257803862929220;

/// natural logarithm of 2
constexpr double ln2 = 0.6931471805599453094172321214581766;
/// 1/ln2
constexpr double inv_ln2 = 1.4426950408889634073599246810018921;
/// natural logarithm of 10
constexpr double ln10 = 2.3025850929940456840179914546843642;

/// 1./3
constexpr double onethird=1.0/3.0;
/// 2./3
constexpr double twothird=2.0/3.0;
/// 4./3
constexpr double fourthird=4.0/3.0;

/// degree -> radian conversion factor
constexpr double degr2rad=pi/180.0;
/// arcmin -> radian conversion factor
constexpr double arcmin2rad=degr2rad/60;
/// radian -> degree conversion factor
constexpr double rad2degr=180.0/pi;

//! Ratio between FWHM and sigma of a Gaussian (\f$\sqrt{8\ln2}\f$).
constexpr double sigma2fwhm=2.3548200450309493; // sqrt(8*log(2.))
//! Ratio between sigma and FWHM of a Gaussian (\f$1/\sqrt{8\ln2}\f$).
constexpr double fwhm2sigma=1/sigma2fwhm;

/*! \} */

}

#endif
