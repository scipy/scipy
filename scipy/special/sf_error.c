#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "sf_error.h"

const char *sf_error_messages[] = {"no error",
                                   "singularity",
                                   "underflow",
                                   "overflow",
                                   "too slow convergence",
                                   "loss of precision",
                                   "no result obtained",
                                   "domain error",
                                   "invalid input argument",
                                   "other error",
                                   NULL};

/* If this isn't volatile clang tries to optimize it away */
static volatile sf_action_t sf_error_actions[] = {
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
    SF_ERROR_IGNORE  /* SF_ERROR__LAST */
};

static sf_callback_t sf_error_callback = NULL;

static sf_callback_fpe_t sf_error_callback_fpe = NULL;

void sf_error_set_action(sf_error_t code, sf_action_t action) { sf_error_actions[(int) code] = action; }

sf_action_t sf_error_get_action(sf_error_t code) { return sf_error_actions[(int) code]; }

void sf_error_set_callback(sf_callback_t callback) { sf_error_callback = callback; }

void sf_error_set_callback_fpe(sf_callback_fpe_t callback_fpe) { sf_error_callback_fpe = callback_fpe; }

void sf_error_v(const char *func_name, sf_error_t code, const char *fmt, va_list ap) {
    char msg[2048], info[1024];
    sf_action_t action;

    if ((int) code < 0 || (int) code >= 10) {
        code = SF_ERROR_OTHER;
    }
    action = sf_error_get_action(code);
    if (action == SF_ERROR_IGNORE) {
        return;
    }

    if (func_name == NULL) {
        func_name = "?";
    }

    if (fmt != NULL && fmt[0] != '\0') {
        vsnprintf(info, 1024, fmt, ap);
        snprintf(msg, 2048, "scipy.special/%s: (%s) %s", func_name, sf_error_messages[(int) code], info);
    } else {
        snprintf(msg, 2048, "scipy.special/%s: %s", func_name, sf_error_messages[(int) code]);
    }

    if (sf_error_callback != NULL) {
        sf_error_callback(action, msg);
    }
}

void sf_error(const char *func_name, sf_error_t code, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    sf_error_v(func_name, code, fmt, ap);
    va_end(ap);
}

#define UFUNC_FPE_DIVIDEBYZERO 1
#define UFUNC_FPE_OVERFLOW 2
#define UFUNC_FPE_UNDERFLOW 4
#define UFUNC_FPE_INVALID 8

void sf_error_check_fpe(const char *func_name) {
    if (sf_error_callback_fpe != NULL) {
        int status = sf_error_callback_fpe();
        if (status & UFUNC_FPE_DIVIDEBYZERO) {
            sf_error(func_name, SF_ERROR_SINGULAR, "floating point division by zero");
        }
        if (status & UFUNC_FPE_UNDERFLOW) {
            sf_error(func_name, SF_ERROR_UNDERFLOW, "floating point underflow");
        }
        if (status & UFUNC_FPE_OVERFLOW) {
            sf_error(func_name, SF_ERROR_OVERFLOW, "floating point overflow");
        }
        if (status & UFUNC_FPE_INVALID) {
            sf_error(func_name, SF_ERROR_DOMAIN, "floating point invalid value");
        }
    }
}
