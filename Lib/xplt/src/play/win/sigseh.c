/*
 * sigseh.c -- $Id$
 * signal handing using Microsoft structured exception handling
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#include "playw.h"
#include "pstdlib.h"

static int w_catch(int code);
static int w_quitting = 0;
static int w_catch_count = 0;

int
w_protect(int (*run)(void))
{
  int result = 0;
  do {
    __try {
      if (p_signalling) w_caught(), w_catch_count = 0;
      result = run();
      w_quitting = 1;
    } __except (EXCEPTION_EXECUTE_HANDLER) {
      p_signalling = w_catch(GetExceptionCode());
      if (w_catch_count++ > 6) w_quitting = 1;
    }
  } while (!w_quitting);
  return result;
}

void (*w_abort_hook)(void) = 0;

void
p_abort(void)
{
  if (!p_signalling) p_signalling = PSIG_SOFT;
  if (w_abort_hook) w_abort_hook();
  RaiseException(0, 0, 0,0);  /* exception code 0 not used by Windows */
}

void
w_siginit(void)
{
}

static int
w_catch(int code)
{
  if (!code && p_signalling)                      code = p_signalling;
  else if (!code && !p_signalling)                code = PSIG_SOFT;
  else if (code==EXCEPTION_FLT_DIVIDE_BY_ZERO ||
           code==EXCEPTION_FLT_OVERFLOW ||
           code==EXCEPTION_FLT_INVALID_OPERATION ||
           code==EXCEPTION_INT_DIVIDE_BY_ZERO ||
           code==EXCEPTION_INT_OVERFLOW ||
           code==EXCEPTION_FLT_UNDERFLOW ||
           code==EXCEPTION_FLT_INEXACT_RESULT ||
           code==EXCEPTION_FLT_DENORMAL_OPERAND ||
           code==EXCEPTION_FLT_STACK_CHECK)       code = PSIG_FPE;
  else if (code==EXCEPTION_ACCESS_VIOLATION ||
           code==EXCEPTION_IN_PAGE_ERROR ||
           code==EXCEPTION_ARRAY_BOUNDS_EXCEEDED) code = PSIG_SEGV;
  else if (code==EXCEPTION_ILLEGAL_INSTRUCTION ||
           code==EXCEPTION_PRIV_INSTRUCTION ||
           code==EXCEPTION_BREAKPOINT ||
           code==EXCEPTION_SINGLE_STEP)           code = PSIG_ILL;
  else if (code==EXCEPTION_DATATYPE_MISALIGNMENT) code = PSIG_BUS;
  else if (code==EXCEPTION_NONCONTINUABLE_EXCEPTION ||
           code==EXCEPTION_INVALID_DISPOSITION ||
           code==EXCEPTION_STACK_OVERFLOW)        code = PSIG_OTHER;
  else                                            code = PSIG_OTHER;

  if (code==PSIG_FPE) w_fpu_setup();

  return code;
}
