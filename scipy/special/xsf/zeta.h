#pragma once

#include "cephes/zeta.h"
#include "cephes/zetac.h"

namespace xsf {

XSF_HOST_DEVICE inline double riemann_zeta(double x) { return cephes::riemann_zeta(x); }

XSF_HOST_DEVICE inline float riemann_zeta(float x) { return riemann_zeta(static_cast<double>(x)); }

XSF_HOST_DEVICE inline double zeta(double x, double q) { return cephes::zeta(x, q); }

XSF_HOST_DEVICE inline float zeta(float x, float q) { return zeta(static_cast<double>(x), static_cast<double>(q)); }

XSF_HOST_DEVICE inline double zetac(double x) { return cephes::zetac(x); }

XSF_HOST_DEVICE inline float zetac(float x) { return zetac(static_cast<double>(x)); }

} // namespace xsf
