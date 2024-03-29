#pragma once

#include "specfun.h"

namespace special {

namespace detail {

    template <typename T>
    void klvna(T x, T *ber, T *bei, T *ger, T *gei, T *der, T *dei, T *her, T *hei) {

        // ======================================================
        // Purpose: Compute Kelvin functions ber x, bei x, ker x
        //          and kei x, and their derivatives  ( x > 0 )
        // Input :  x   --- Argument of Kelvin functions
        // Output:  BER --- ber x
        //          BEI --- bei x
        //          GER --- ker x
        //          GEI --- kei x
        //          DER --- ber'x
        //          DEI --- bei'x
        //          HER --- ker'x
        //          HEI --- kei'x
        // ================================================

        int k, km, m;
        T gs, r, x2, x4, pp1, pn1, qp1, qn1, r1, pp0, pn0, qp0, qn0, r0, fac, xt, cs, ss, xd, xe1, xe2, xc1, xc2, cp0,
            cn0, sp0, sn0, rc, rs;
        const T pi = 3.141592653589793;
        const T el = 0.5772156649015329;
        const T eps = 1.0e-15;

        if (x == 0.0) {
            *ber = 1.0;
            *bei = 0.0;
            *ger = 1.0e+300;
            *gei = -0.25 * pi;
            *der = 0.0;
            *dei = 0.0;
            *her = -1.0e+300;
            *hei = 0.0;
            return;
        }

        x2 = 0.25 * x * x;
        x4 = x2 * x2;

        if (fabs(x) < 10.0) {
            *ber = 1.0;
            r = 1.0;
            for (m = 1; m <= 60; m++) {
                r = -0.25 * r / (m * m) / pow((2.0 * m - 1.0), 2) * x4;
                *ber += r;
                if (fabs(r) < fabs(*ber) * eps) {
                    break;
                }
            }

            *bei = x2;
            r = x2;
            for (m = 1; m <= 60; m++) {
                r = -0.25 * r / (m * m) / pow((2.0 * m + 1.0), 2) * x4;
                *bei += r;
                if (fabs(r) < fabs(*bei) * eps) {
                    break;
                }
            }

            *ger = -(log(x / 2.0) + el) * (*ber) + 0.25 * pi * (*bei);

            r = 1.0;
            gs = 0.0;
            for (m = 1; m <= 60; m++) {
                r = -0.25 * r / (m * m) / pow((2.0 * m - 1.0), 2) * x4;
                gs = gs + 1.0 / (2.0 * m - 1.0) + 1.0 / (2.0 * m);
                *ger += r * gs;
                if (fabs(r * gs) < fabs(*ger) * eps) {
                    break;
                }
            }

            *gei = x2 - (log(x / 2.0) + el) * (*bei) - 0.25 * pi * (*ber);

            r = x2;
            gs = 1.0;
            for (m = 1; m <= 60; m++) {
                r = -0.25 * r / (m * m) / pow((2.0 * m + 1.0), 2) * x4;
                gs = gs + 1.0 / (2.0 * m) + 1.0 / (2.0 * m + 1.0);
                *gei += r * gs;
                if (fabs(r * gs) < fabs(*gei) * eps) {
                    break;
                }
            }

            *der = -0.25 * x * x2;
            r = *der;
            for (m = 1; m <= 60; m++) {
                r = -0.25 * r / m / (m + 1.0) / pow((2.0 * m + 1.0), 2) * x4;
                *der += r;
                if (fabs(r) < fabs(*der) * eps) {
                    break;
                }
            }

            *dei = 0.5 * x;
            r = *dei;
            for (m = 1; m <= 60; m++) {
                r = -0.25 * r / (m * m) / (2.0 * m - 1.0) / (2.0 * m + 1.0) * x4;
                *dei += r;
                if (fabs(r) < fabs(*dei) * eps) {
                    break;
                }
            }

            r = -0.25 * x * x2;
            gs = 1.5;
            *her = 1.5 * r - (*ber) / x - (log(x / 2.0) + el) * (*der) + 0.25 * pi * (*dei);
            for (m = 1; m <= 60; m++) {
                r = -0.25 * r / m / (m + 1.0) / pow((2.0 * m + 1.0), 2) * x4;
                gs = gs + 1.0 / (2.0 * m + 1.0) + 1.0 / (2 * m + 2.0);
                *her += r * gs;
                if (fabs(r * gs) < fabs(*her) * eps) {
                    break;
                }
            }

            r = 0.5 * x;
            gs = 1.0;
            *hei = 0.5 * x - (*bei) / x - (log(x / 2.0) + el) * (*dei) - 0.25 * pi * (*der);
            for (m = 1; m <= 60; m++) {
                r = -0.25 * r / (m * m) / (2 * m - 1.0) / (2 * m + 1.0) * x4;
                gs = gs + 1.0 / (2.0 * m) + 1.0 / (2 * m + 1.0);
                *hei += r * gs;
                if (fabs(r * gs) < fabs(*hei) * eps) {
                    return;
                }
            }
        } else {
            pp0 = 1.0;
            pn0 = 1.0;
            qp0 = 0.0;
            qn0 = 0.0;
            r0 = 1.0;
            km = 18;
            if (fabs(x) >= 40.0)
                km = 10;
            fac = 1.0;
            for (k = 1; k <= km; k++) {
                fac = -fac;
                xt = 0.25 * k * pi - trunc(0.125 * k) * 2.0 * pi;
                cs = cos(xt);
                ss = sin(xt);
                r0 = 0.125 * r0 * pow((2.0 * k - 1.0), 2) / k / x;
                rc = r0 * cs;
                rs = r0 * ss;
                pp0 += rc;
                pn0 += fac * rc;
                qp0 += rs;
                qn0 += fac * rs;
            }

            xd = x / sqrt(2.0);
            xe1 = exp(xd);
            xe2 = exp(-xd);
            xc1 = 1.0 / sqrt(2.0 * pi * x);
            xc2 = sqrt(0.5 * pi / x);
            cp0 = cos(xd + 0.125 * pi);
            cn0 = cos(xd - 0.125 * pi);
            sp0 = sin(xd + 0.125 * pi);
            sn0 = sin(xd - 0.125 * pi);

            *ger = xc2 * xe2 * (pn0 * cp0 - qn0 * sp0);
            *gei = xc2 * xe2 * (-pn0 * sp0 - qn0 * cp0);
            *ber = xc1 * xe1 * (pp0 * cn0 + qp0 * sn0) - (*gei) / pi;
            *bei = xc1 * xe1 * (pp0 * sn0 - qp0 * cn0) + (*ger) / pi;

            pp1 = 1.0;
            pn1 = 1.0;
            qp1 = 0.0;
            qn1 = 0.0;
            r1 = 1.0;
            fac = 1.0;
            for (int k = 1; k <= km; k++) {
                fac = -fac;
                xt = 0.25 * k * pi - (int) (0.125 * k) * 2.0 * pi;
                cs = cos(xt);
                ss = sin(xt);
                r1 = 0.125 * r1 * (4.0 - pow(2.0 * k - 1.0, 2)) / (k * x);
                rc = r1 * cs;
                rs = r1 * ss;
                pp1 += fac * rc;
                pn1 += rc;
                qp1 += fac * rs;
                qn1 += rs;
            }
            *her = xc2 * xe2 * (-pn1 * cn0 + qn1 * sn0);
            *hei = xc2 * xe2 * (pn1 * sn0 + qn1 * cn0);
            *der = xc1 * xe1 * (pp1 * cp0 + qp1 * sp0) - (*hei) / pi;
            *dei = xc1 * xe1 * (pp1 * sp0 - qp1 * cp0) + (*her) / pi;
        }
        return;
    }

} // namespace detail

template <typename T>
T ber(T x) {
    std::complex<T> Be;
    T ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        x = -x;
    }

    detail::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Be.real(ber);
    Be.imag(bei);
    SPECFUN_ZCONVINF("ber", Be);
    return Be.real();
}

template <typename T>
T bei(T x) {
    std::complex<T> Be;
    T ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        x = -x;
    }

    detail::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Be.real(ber);
    Be.imag(bei);
    SPECFUN_ZCONVINF("bei", Be);
    return Be.imag();
}

template <typename T>
T ker(T x) {
    std::complex<T> Ke;
    T ber, bei, ger, gei, der, dei, her, hei;
    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    detail::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Ke.real(ger);
    Ke.imag(gei);
    SPECFUN_ZCONVINF("ker", Ke);
    return Ke.real();
}

template <typename T>
T kei(T x) {
    std::complex<T> Ke;
    T ber, bei, ger, gei, der, dei, her, hei;
    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    detail::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Ke.real(ger);
    Ke.imag(gei);
    SPECFUN_ZCONVINF("kei", Ke);
    return Ke.imag();
}

template <typename T>
T berp(T x) {
    std::complex<T> Bep;
    T ber, bei, ger, gei, der, dei, her, hei;
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }

    detail::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Bep.real(der);
    Bep.imag(dei);
    SPECFUN_ZCONVINF("berp", Bep);
    if (flag) {
        return -Bep.real();
    }
    return Bep.real();
}

template <typename T>
T beip(T x) {
    std::complex<T> Bep;
    T ber, bei, ger, gei, der, dei, her, hei;
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    detail::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Bep.real(der);
    Bep.imag(dei);
    SPECFUN_ZCONVINF("beip", Bep);
    if (flag) {
        return -Bep.imag();
    }
    return Bep.imag();
}

template <typename T>
T kerp(T x) {
    std::complex<T> Kep;
    T ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    detail::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Kep.real(her);
    Kep.imag(hei);
    SPECFUN_ZCONVINF("kerp", Kep);
    return Kep.real();
}

template <typename T>
T keip(T x) {
    std::complex<T> Kep;
    T ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    detail::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Kep.real(her);
    Kep.imag(hei);
    SPECFUN_ZCONVINF("keip", Kep);
    return Kep.imag();
}

template <typename T>
void kelvin(T x, std::complex<T> *Be, std::complex<T> *Ke, std::complex<T> *Bep, std::complex<T> *Kep) {
    int flag = 0;
    T ber, bei, ger, gei, der, dei, her, hei;
    if (x < 0) {
        x = -x;
        flag = 1;
    }

    detail::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Be->real(ber);
    Be->imag(bei);
    Ke->real(ger);
    Ke->imag(gei);
    Bep->real(der);
    Bep->imag(dei);
    Kep->real(her);
    Kep->imag(hei);

    SPECFUN_ZCONVINF("klvna", *Be);
    SPECFUN_ZCONVINF("klvna", *Ke);
    SPECFUN_ZCONVINF("klvna", *Bep);
    SPECFUN_ZCONVINF("klvna", *Kep);
    if (flag) {
        Bep->real(-Bep->real());
        Bep->imag(-Bep->imag());
        Ke->real(std::numeric_limits<T>::quiet_NaN());
        Ke->imag(std::numeric_limits<T>::quiet_NaN());
        Kep->real(std::numeric_limits<T>::quiet_NaN());
        Kep->imag(std::numeric_limits<T>::quiet_NaN());
    }
}

inline void klvnzo(int nt, int kd, double *zo) {

    // ====================================================
    // Purpose: Compute the zeros of Kelvin functions
    // Input :  NT  --- Total number of zeros
    //          KD  --- Function code
    //          KD=1 to 8 for ber x, bei x, ker x, kei x,
    //                    ber'x, bei'x, ker'x and kei'x,
    //                    respectively.
    // Output:  ZO(M) --- the M-th zero of Kelvin function
    //                    for code KD
    // Routine called:
    //          KLVNA for computing Kelvin functions and
    //          their derivatives
    // ====================================================

    double ber, bei, ger, gei, der, dei, her, hei;
    double rt0[9] = {0.0, 2.84891, 5.02622, 1.71854, 3.91467, 6.03871, 3.77268, 2.66584, 4.93181};
    double rt = rt0[kd];

    for (int m = 1; m <= nt; m++) {
        while (1) {
            detail::klvna(rt, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
            if (kd == 1) {
                rt -= ber / der;
            } else if (kd == 2) {
                rt -= bei / dei;
            } else if (kd == 3) {
                rt -= ger / her;
            } else if (kd == 4) {
                rt -= gei / hei;
            } else if (kd == 5) {
                rt -= der / (-bei - der / rt);
            } else if (kd == 6) {
                rt -= dei / (ber - dei / rt);
            } else if (kd == 7) {
                rt -= her / (-gei - her / rt);
            } else {
                rt -= hei / (ger - hei / rt);
            }

            if (fabs(rt - rt0[kd]) <= 5e-10) {
                break;
            } else {
                rt0[kd] = rt;
            }
        }
        zo[m - 1] = rt;
        rt += 4.44;
    }
}

} // namespace special
