#pragma once

#include "cephes/ellie.h"
#include "cephes/ellik.h"
#include "cephes/ellpe.h"
#include "cephes/ellpj.h"
#include "cephes/ellpk.h"
#include "config.h"

namespace xsf {

inline double ellipe(double m) { return cephes::ellpe(m); }

inline float ellipe(float m) { return ellipe(static_cast<double>(m)); }

inline double ellipeinc(double phi, double m) { return cephes::ellie(phi, m); }

inline float ellipeinc(float phi, float m) { return ellipeinc(static_cast<double>(phi), static_cast<double>(m)); }

inline void ellipj(double u, double m, double &sn, double &cn, double &dn, double &ph) {
    cephes::ellpj(u, m, &sn, &cn, &dn, &ph);
}

inline void ellipj(float u, float m, float &sn, float &cn, float &dn, float &ph) {
    double sn_double;
    double cn_double;
    double dn_double;
    double ph_double;
    ellipj(static_cast<double>(u), static_cast<double>(m), sn_double, cn_double, dn_double, ph_double);

    sn = sn_double;
    cn = cn_double;
    dn = dn_double;
    ph = ph_double;
}

inline double ellipkinc(double phi, double m) { return cephes::ellik(phi, m); }

inline float ellipkinc(float phi, float m) { return ellipkinc(static_cast<double>(phi), static_cast<double>(m)); }

XSF_HOST_DEVICE inline double ellipk(double m) { return cephes::ellpk(1.0 - m); }

XSF_HOST_DEVICE inline float ellipk(float m) { return ellipk(static_cast<double>(m)); }

inline double ellipkm1(double p) { return cephes::ellpk(p); }

inline float ellipkm1(float p) { return ellipkm1(static_cast<double>(p)); }

} // namespace xsf
