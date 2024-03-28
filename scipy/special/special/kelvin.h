#pragma once

#include "specfun.h"

namespace special {

template <typename T>
T ber(T x) {
    std::complex<T> Be;
    T ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        x = -x;
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
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

    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
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
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
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
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
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
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
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
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
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
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
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
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
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

    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
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

} // namespace special
