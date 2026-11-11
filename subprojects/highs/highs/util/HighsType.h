/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsType.h
 * @brief The definition for basic types to use
 */

#ifndef UTIL_HIGHS_TYPE_H_
#define UTIL_HIGHS_TYPE_H_

#include "util/HighsInt.h"

// vector<bool> is not thread-safe, so HiGHS uses vector<HighsBool>

typedef uint8_t HighsBool;

#endif
