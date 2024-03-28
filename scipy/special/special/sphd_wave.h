#pragma once

#include "specfun.h"

namespace special {

template <typename T>
void oblate_aswfa_nocv(T m, T n, T c, T x, T *s1, T *s1d) {
    *s1 = oblate_aswfa_nocv(m, n, c, x, s1d);
}

} // namespace special
