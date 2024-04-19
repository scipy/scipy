#pragma once

#include "amos/amos.h"
#include "error.h"

namespace special {

inline sf_error_t ierr_to_sferr(int nz, int ierr) {
    /* Return sf_error equivalents for amos ierr values */

    if (nz != 0) {
        return SF_ERROR_UNDERFLOW;
    }

    switch (ierr) {
    case 1:
        return SF_ERROR_DOMAIN;
    case 2:
        return SF_ERROR_OVERFLOW;
    case 3:
        return SF_ERROR_LOSS;
    case 4:
        return SF_ERROR_NO_RESULT;
    case 5: /* Algorithm termination condition not met */
        return SF_ERROR_NO_RESULT;
    }

    return SF_ERROR_OK;
}
} // namespace special
