#pragma once

#include <xsf/error.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    SF_ERROR_IGNORE = 0,  /* Ignore errors */
    SF_ERROR_WARN,        /* Warn on errors */
    SF_ERROR_RAISE        /* Raise on errors */
} sf_action_t;

extern const char *sf_error_messages[];

void sf_error_set_action(sf_error_t code, sf_action_t action);

sf_action_t sf_error_get_action(sf_error_t code);

void sf_error(const char *func_name, sf_error_t code, const char *fmt, ...);

void sf_error_check_fpe(const char *func_name);

#ifdef __cplusplus
}
#endif
