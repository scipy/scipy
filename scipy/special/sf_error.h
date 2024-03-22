#ifndef SF_ERROR_H_
#define SF_ERROR_H_

#include "special/error.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef enum {
    SF_ERROR_IGNORE = 0,  /* Ignore errors */
    SF_ERROR_WARN,        /* Warn on errors */
    SF_ERROR_RAISE        /* Raise on errors */
} sf_action_t;

typedef void (*sf_callback_t)(sf_action_t, const char *);

typedef int (*sf_callback_fpe_t)();

extern const char *sf_error_messages[];
void sf_error(const char *func_name, sf_error_t code, const char *fmt, ...);
void sf_error_check_fpe(const char *func_name);
void sf_error_set_action(sf_error_t code, sf_action_t action);
sf_action_t sf_error_get_action(sf_error_t code);

void sf_error_set_callback(sf_callback_t callback);
void sf_error_set_callback_fpe(sf_callback_fpe_t callback_fpe);

#ifdef __cplusplus
}
#endif

#endif /* SF_ERROR_H_ */
