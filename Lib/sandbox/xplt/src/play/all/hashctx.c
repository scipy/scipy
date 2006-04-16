/*
 * hashctx.c -- $Id$
 * generic pointer<->context association functions
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "phash.h"

static p_hashtab *ctx_table = 0;

void
p_setctx(void *ptr, void *context)
{
  if (!ctx_table) ctx_table = p_halloc(64);
  p_hinsert(ctx_table, P_PHASH(ptr), context);
}

void *
p_getctx(void *ptr)
{
  return ctx_table? p_hfind(ctx_table, P_PHASH(ptr)) : 0;
}
