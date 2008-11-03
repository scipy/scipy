/***************************************************************************
 * blitz/tau.h       Integration with the Tau profiling package.
 *                   See http://www.acl.lanl.gov/tau/
 *
 * $Id: tau.h 1414 2005-11-01 22:04:59Z cookedm $
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This code was relicensed under the modified BSD license for use in SciPy
 * by Todd Veldhuizen (see LICENSE.txt in the weave directory).
 *
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ***************************************************************************/

#ifndef BZ_TAU_H
#define BZ_TAU_H

#ifdef BZ_TAU_PROFILING
 #define TAU_BLITZ  TAU_USER1
 #include <Profile/Profiler.h>

#else
 #define TYPE_STRING(profileString, str)
 #define PROFILED_BLOCK(name, type)
 #define TAU_TYPE_STRING(profileString, str)
 #define TAU_PROFILE(name, type, group)
 #define TAU_PROFILE_TIMER(var, name, type, group)
 #define TAU_PROFILE_START(var)
 #define TAU_PROFILE_STOP(var)
 #define TAU_PROFILE_STMT(stmt)
 #define TAU_PROFILE_EXIT(msg)
 #define TAU_PROFILE_INIT(argc, argv)
 #define TAU_PROFILE_SET_NODE(node)
 #define CT(obj)
#endif // ! BZ_TAU_PROFILING

#endif // BZ_TAU_H
