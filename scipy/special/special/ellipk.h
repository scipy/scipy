#pragma once

#include "config.h"
#include "cephes/ellpk.h"

namespace special {

    SPECFUN_HOST_DEVICE inline double ellipk(double m) {
        return cephes::ellpk(1.0 - m);
    }
}
