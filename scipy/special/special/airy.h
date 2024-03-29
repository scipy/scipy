#pragma once

#include "specfun.h"

namespace special {

template <typename T>
void itairy(T x, T *apt, T *bpt, T *ant, T *bnt) {
    T tmp;
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::itairy(x, apt, bpt, ant, bnt);
    if (flag) { /* negative limit -- switch signs and roles */
        tmp = *apt;
        *apt = -*ant;
        *ant = -tmp;
        tmp = *bpt;
        *bpt = -*bnt;
        *bnt = -tmp;
    }
}

} // namespace special
