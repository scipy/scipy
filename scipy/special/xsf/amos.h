#pragma once

#include "amos/amos.h"
#include "error.h"

namespace xsf {

//
// Return sf_error equivalents for AMOS ierr values.
// 'ierr' refers to the last parameter of the AMOS functions
// airy(), besh(), besi(), besj(), besk(), besy(), and biry().
//
inline sf_error_t ierr_to_sferr(int nz, int ierr) {
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
    case 6: /* Memory allocation failed */
        return SF_ERROR_MEMORY;
    }

    return SF_ERROR_OK;
}
} // namespace xsf
