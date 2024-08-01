/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/* Copyright Benjamin Sobotta 2012
 *
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)
 */

/*
 * Reference:
 * Mike Patefield, David Tandy
 * FAST AND ACCURATE CALCULATION OF OWEN'S T-FUNCTION
 * Journal of Statistical Software, 5 (5), 1-25
 */
#pragma once

#include "../config.h"

#include "ndtr.h"
#include "unity.h"

namespace xsf {
namespace cephes {

    namespace detail {

        constexpr int owens_t_SELECT_METHOD[] = {
            0,  0,  1,  12, 12, 12, 12, 12, 12, 12, 12, 15, 15, 15, 8,  0,  1,  1,  2,  2,  4,  4,  13, 13,
            14, 14, 15, 15, 15, 8,  1,  1,  2,  2,  2,  4,  4,  14, 14, 14, 14, 15, 15, 15, 9,  1,  1,  2,
            4,  4,  4,  4,  6,  6,  15, 15, 15, 15, 15, 9,  1,  2,  2,  4,  4,  5,  5,  7,  7,  16, 16, 16,
            11, 11, 10, 1,  2,  4,  4,  4,  5,  5,  7,  7,  16, 16, 16, 11, 11, 11, 1,  2,  3,  3,  5,  5,
            7,  7,  16, 16, 16, 16, 16, 11, 11, 1,  2,  3,  3,  5,  5,  17, 17, 17, 17, 16, 16, 16, 11, 11};

        constexpr double owens_t_HRANGE[] = {0.02, 0.06, 0.09, 0.125, 0.26, 0.4, 0.6,
                                             1.6,  1.7,  2.33, 2.4,   3.36, 3.4, 4.8};

        constexpr double owens_t_ARANGE[] = {0.025, 0.09, 0.15, 0.36, 0.5, 0.9, 0.99999};

        constexpr double owens_t_ORD[] = {2, 3, 4, 5, 7, 10, 12, 18, 10, 20, 30, 0, 4, 7, 8, 20, 0, 0};

        constexpr int owens_t_METHODS[] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6};

        constexpr double owens_t_C[] = {
            1.0,
            -1.0,
            1.0,
            -0.9999999999999998,
            0.9999999999999839,
            -0.9999999999993063,
            0.9999999999797337,
            -0.9999999995749584,
            0.9999999933226235,
            -0.9999999188923242,
            0.9999992195143483,
            -0.9999939351372067,
            0.9999613559769055,
            -0.9997955636651394,
            0.9990927896296171,
            -0.9965938374119182,
            0.9891001713838613,
            -0.9700785580406933,
            0.9291143868326319,
            -0.8542058695956156,
            0.737965260330301,
            -0.585234698828374,
            0.4159977761456763,
            -0.25882108752419436,
            0.13755358251638927,
            -0.060795276632595575,
            0.021633768329987153,
            -0.005934056934551867,
            0.0011743414818332946,
            -0.0001489155613350369,
            9.072354320794358e-06,
        };

        constexpr double owens_t_PTS[] = {
            0.35082039676451715489E-02, 0.31279042338030753740E-01, 0.85266826283219451090E-01,
            0.16245071730812277011E+00, 0.25851196049125434828E+00, 0.36807553840697533536E+00,
            0.48501092905604697475E+00, 0.60277514152618576821E+00, 0.71477884217753226516E+00,
            0.81475510988760098605E+00, 0.89711029755948965867E+00, 0.95723808085944261843E+00,
            0.99178832974629703586E+00};

        constexpr double owens_t_WTS[] = {
            0.18831438115323502887E-01, 0.18567086243977649478E-01, 0.18042093461223385584E-01,
            0.17263829606398753364E-01, 0.16243219975989856730E-01, 0.14994592034116704829E-01,
            0.13535474469662088392E-01, 0.11886351605820165233E-01, 0.10070377242777431897E-01,
            0.81130545742299586629E-02, 0.60419009528470238773E-02, 0.38862217010742057883E-02,
            0.16793031084546090448E-02};

        XSF_HOST_DEVICE inline int get_method(double h, double a) {
            int ihint, iaint, i;

            ihint = 14;
            iaint = 7;

            for (i = 0; i < 14; i++) {
                if (h <= owens_t_HRANGE[i]) {
                    ihint = i;
                    break;
                }
            }

            for (i = 0; i < 7; i++) {
                if (a <= owens_t_ARANGE[i]) {
                    iaint = i;
                    break;
                }
            }
            return owens_t_SELECT_METHOD[iaint * 15 + ihint];
        }

        XSF_HOST_DEVICE inline double owens_t_norm1(double x) { return xsf::cephes::erf(x / std::sqrt(2)) / 2; }

        XSF_HOST_DEVICE inline double owens_t_norm2(double x) {
            return xsf::cephes::erfc(x / std::sqrt(2)) / 2;
        }

        XSF_HOST_DEVICE inline double owensT1(double h, double a, double m) {
            int j = 1;
            int jj = 1;

            double hs = -0.5 * h * h;
            double dhs = std::exp(hs);
            double as = a * a;
            double aj = a / (2 * M_PI);
            double dj = xsf::cephes::expm1(hs);
            double gj = hs * dhs;

            double val = std::atan(a) / (2 * M_PI);

            while (1) {
                val += dj * aj / jj;

                if (m <= j) {
                    break;
                }
                j++;
                jj += 2;
                aj *= as;
                dj = gj - dj;
                gj *= hs / j;
            }

            return val;
        }

        XSF_HOST_DEVICE inline double owensT2(double h, double a, double ah, double m) {
            int i = 1;
            int maxi = 2 * m + 1;
            double hs = h * h;
            double as = -a * a;
            double y = 1.0 / hs;
            double val = 0.0;
            double vi = a * std::exp(-0.5 * ah * ah) / std::sqrt(2 * M_PI);
            double z = (xsf::cephes::ndtr(ah) - 0.5) / h;

            while (1) {
                val += z;
                if (maxi <= i) {
                    break;
                }
                z = y * (vi - i * z);
                vi *= as;
                i += 2;
            }
            val *= std::exp(-0.5 * hs) / std::sqrt(2 * M_PI);

            return val;
        }

        XSF_HOST_DEVICE inline double owensT3(double h, double a, double ah) {
            double aa, hh, y, vi, zi, result;
            int i;

            aa = a * a;
            hh = h * h;
            y = 1 / hh;

            vi = a * std::exp(-ah * ah / 2) / std::sqrt(2 * M_PI);
            zi = owens_t_norm1(ah) / h;
            result = 0;

            for (i = 0; i <= 30; i++) {
                result += zi * owens_t_C[i];
                zi = y * ((2 * i + 1) * zi - vi);
                vi *= aa;
            }

            result *= std::exp(-hh / 2) / std::sqrt(2 * M_PI);

            return result;
        }

        XSF_HOST_DEVICE inline double owensT4(double h, double a, double m) {
            double maxi, hh, naa, ai, yi, result;
            int i;

            maxi = 2 * m + 1;
            hh = h * h;
            naa = -a * a;

            i = 1;
            ai = a * std::exp(-hh * (1 - naa) / 2) / (2 * M_PI);
            yi = 1;
            result = 0;

            while (1) {
                result += ai * yi;

                if (maxi <= i) {
                    break;
                }

                i += 2;
                yi = (1 - hh * yi) / i;
                ai *= naa;
            }

            return result;
        }

        XSF_HOST_DEVICE inline double owensT5(double h, double a) {
            double result, r, aa, nhh;
            int i;

            result = 0;
            r = 0;
            aa = a * a;
            nhh = -0.5 * h * h;

            for (i = 1; i < 14; i++) {
                r = 1 + aa * owens_t_PTS[i - 1];
                result += owens_t_WTS[i - 1] * std::exp(nhh * r) / r;
            }

            result *= a;

            return result;
        }

        XSF_HOST_DEVICE inline double owensT6(double h, double a) {
            double normh, y, r, result;

            normh = owens_t_norm2(h);
            y = 1 - a;
            r = std::atan2(y, (1 + a));
            result = normh * (1 - normh) / 2;

            if (r != 0) {
                result -= r * std::exp(-y * h * h / (2 * r)) / (2 * M_PI);
            }

            return result;
        }

        XSF_HOST_DEVICE inline double owens_t_dispatch(double h, double a, double ah) {
            int index, meth_code;
            double m, result;

            if (h == 0) {
                return std::atan(a) / (2 * M_PI);
            }
            if (a == 0) {
                return 0;
            }
            if (a == 1) {
                return owens_t_norm2(-h) * owens_t_norm2(h) / 2;
            }

            index = get_method(h, a);
            m = owens_t_ORD[index];
            meth_code = owens_t_METHODS[index];

            switch (meth_code) {
            case 1:
                result = owensT1(h, a, m);
                break;
            case 2:
                result = owensT2(h, a, ah, m);
                break;
            case 3:
                result = owensT3(h, a, ah);
                break;
            case 4:
                result = owensT4(h, a, m);
                break;
            case 5:
                result = owensT5(h, a);
                break;
            case 6:
                result = owensT6(h, a);
                break;
            default:
                result = std::numeric_limits<double>::quiet_NaN();
            }

            return result;
        }

    } // namespace detail

    XSF_HOST_DEVICE inline double owens_t(double h, double a) {
        double result, fabs_a, fabs_ah, normh, normah;

        if (std::isnan(h) || std::isnan(a)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        /* exploit that T(-h,a) == T(h,a) */
        h = std::abs(h);

        /*
         * Use equation (2) in the paper to remap the arguments such that
         * h >= 0 and 0 <= a <= 1 for the call of the actual computation
         * routine.
         */
        fabs_a = std::abs(a);
        fabs_ah = fabs_a * h;

        if (fabs_a == std::numeric_limits<double>::infinity()) {
            /* See page 13 in the paper */
            result = 0.5 * detail::owens_t_norm2(h);
        } else if (h == std::numeric_limits<double>::infinity()) {
            result = 0;
        } else if (fabs_a <= 1) {
            result = detail::owens_t_dispatch(h, fabs_a, fabs_ah);
        } else {
            if (fabs_ah <= 0.67) {
                normh = detail::owens_t_norm1(h);
                normah = detail::owens_t_norm1(fabs_ah);
                result = 0.25 - normh * normah - detail::owens_t_dispatch(fabs_ah, (1 / fabs_a), h);
            } else {
                normh = detail::owens_t_norm2(h);
                normah = detail::owens_t_norm2(fabs_ah);
                result = (normh + normah) / 2 - normh * normah - detail::owens_t_dispatch(fabs_ah, (1 / fabs_a), h);
            }
        }

        if (a < 0) {
            /* exploit that T(h,-a) == -T(h,a) */
            return -result;
        }

        return result;
    }

} // namespace cephes
} // namespace xsf
