/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsInt.h
 * @brief The definition for the integer type to use
 */

#ifndef UTIL_HIGHS_INT_H_
#define UTIL_HIGHS_INT_H_

#include <stdint.h>

#ifdef __cplusplus
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#endif
#include <inttypes.h>

#include "HConfig.h"

#ifdef HIGHSINT64
typedef int64_t HighsInt;
typedef uint64_t HighsUInt;
#define HIGHSINT_FORMAT PRId64
#else
typedef int HighsInt;
typedef unsigned int HighsUInt;
#define HIGHSINT_FORMAT "d"
#endif

#endif
