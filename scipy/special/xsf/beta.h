#pragma once

#include "cephes/beta.h"

namespace xsf {

XSF_HOST_DEVICE inline double beta(double a, double b) { return cephes::beta(a, b); }

XSF_HOST_DEVICE inline float beta(float a, float b) { return beta(static_cast<double>(a), static_cast<double>(b)); }

XSF_HOST_DEVICE inline double betaln(double a, double b) { return cephes::lbeta(a, b); }

XSF_HOST_DEVICE inline float betaln(float a, float b) { return betaln(static_cast<double>(a), static_cast<double>(b)); }

} // namespace xsf
