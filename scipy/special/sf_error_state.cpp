#include <stdlib.h>
#include "sf_error_state.h"

#ifdef __MINGW32__
#include <mutex>
std::mutex err_mutex;

#define THREAD_LOCAL
#else
#define THREAD_LOCAL thread_local
#endif



/* If this isn't volatile clang tries to optimize it away */
static THREAD_LOCAL sf_action_t sf_error_actions[] = {
    SF_ERROR_IGNORE, /* SF_ERROR_OK */
    SF_ERROR_IGNORE, /* SF_ERROR_SINGULAR */
    SF_ERROR_IGNORE, /* SF_ERROR_UNDERFLOW */
    SF_ERROR_IGNORE, /* SF_ERROR_OVERFLOW */
    SF_ERROR_IGNORE, /* SF_ERROR_SLOW */
    SF_ERROR_IGNORE, /* SF_ERROR_LOSS */
    SF_ERROR_IGNORE, /* SF_ERROR_NO_RESULT */
    SF_ERROR_IGNORE, /* SF_ERROR_DOMAIN */
    SF_ERROR_IGNORE, /* SF_ERROR_ARG */
    SF_ERROR_IGNORE, /* SF_ERROR_OTHER */
    SF_ERROR_RAISE,  /* SF_ERROR_MEMORY */
    SF_ERROR_IGNORE  /* SF_ERROR__LAST */
};


SCIPY_DLL void scipy_sf_error_set_action(sf_error_t code, sf_action_t action)
{
    #ifdef __MINGW32__
    std::lock_guard<std::mutex> guard(err_mutex);
    #endif
    sf_error_actions[(int)code] = action;
}


SCIPY_DLL sf_action_t scipy_sf_error_get_action(sf_error_t code)
{
    #ifdef __MINGW32__
    std::lock_guard<std::mutex> guard(err_mutex);
    #endif
    return sf_error_actions[(int)code];
}
