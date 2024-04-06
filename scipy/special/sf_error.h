#pragma once

#include "sf_error_state.h"

#ifdef __cplusplus
extern "C" {
#endif

extern const char *sf_error_messages[];
void sf_error(const char *func_name, sf_error_t code, const char *fmt, ...);
void sf_error_check_fpe(const char *func_name);


#ifdef __cplusplus
}
#endif
