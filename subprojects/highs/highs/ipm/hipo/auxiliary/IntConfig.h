#ifndef HIPO_INT_CONFIG_H
#define HIPO_INT_CONFIG_H

#include "lp_data/HConst.h"
#include "util/HighsType.h"

namespace hipo {

// Generic integer type from HiGHS
typedef HighsInt Int;

// Integer type for factorisation
typedef int64_t Int64;

// The matrix (AS or NE) is formed using Int, so it must have fewer than
// kHighsIInf nonzero entries.
//
// The factorisation uses a combination of Int and Int64. Int64 is used only
// where needed, i.e. when dealing with nonzeros of the factor or when checking
// overflows.
// BLAS is 32-bit, so the vectors used by BLAS must be addressable with 32-bit
// integers.
//

}  // namespace hipo

#endif