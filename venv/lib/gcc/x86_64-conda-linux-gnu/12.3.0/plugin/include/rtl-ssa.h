// On-the-side RTL SSA representation                               -*- C++ -*-
// Copyright (C) 2020-2022 Free Software Foundation, Inc.
//
// This file is part of GCC.
//
// GCC is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 3, or (at your option) any later
// version.
//
// GCC is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// along with GCC; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.

#ifndef GCC_RTL_SSA_H
#define GCC_RTL_SSA_H 1

// This is an aggregation header file.  This means it should contain only
// other include files.

#if 0
// Files that use this one should first have:
#define INCLUDE_ALGORITHM
#define INCLUDE_FUNCTIONAL
#include "config.h"
#include "system.h"
#include "coretypes.h"
#include "backend.h"
#include "rtl.h"
#include "df.h"
#endif

// Needed by splay-tree-utils.h and directly by rtl-ssa.
#include "pretty-print.h"

// Needed directly by recog.h.
#include "insn-config.h"

// Needed directly by rtl-ssa.
#include "splay-tree-utils.h"
#include "recog.h"
#include "regs.h"
#include "function-abi.h"
#include "obstack-utils.h"
#include "mux-utils.h"
#include "rtlanal.h"

// Provides the global crtl->ssa.
#include "memmodel.h"
#include "tm_p.h"
#include "emit-rtl.h"

// The rtl-ssa files themselves.
#include "rtl-ssa/accesses.h"
#include "rtl-ssa/insns.h"
#include "rtl-ssa/blocks.h"
#include "rtl-ssa/changes.h"
#include "rtl-ssa/functions.h"
#include "rtl-ssa/is-a.inl"
#include "rtl-ssa/access-utils.h"
#include "rtl-ssa/insn-utils.h"
#include "rtl-ssa/movement.h"
#include "rtl-ssa/change-utils.h"
#include "rtl-ssa/member-fns.inl"

#endif
