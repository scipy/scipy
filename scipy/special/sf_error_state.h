#pragma once

#include "scipy_dll.h"
#include "xsf/error.h"


#ifdef __cplusplus
extern "C" {
#endif

    typedef enum {
        SF_ERROR_IGNORE = 0,  /* Ignore errors */
        SF_ERROR_WARN,        /* Warn on errors */
        SF_ERROR_RAISE        /* Raise on errors */
    } sf_action_t;
    
    SCIPY_DLL void scipy_sf_error_set_action(sf_error_t code, sf_action_t action);

    SCIPY_DLL sf_action_t scipy_sf_error_get_action(sf_error_t code);

#ifdef __cplusplus
}
#endif
