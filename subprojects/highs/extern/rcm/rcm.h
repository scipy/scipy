#ifndef REVERSE_CUTHILL_MCKEE_H
#define REVERSE_CUTHILL_MCKEE_H

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

// Function to compute Reverse Cuthill-McKee ordering.
// Taken from sparsepak:
// https://people.sc.fsu.edu/~jburkardt/f77_src/sparsepak/sparsepak.html
// available under MIT license.
// Changes:
// - Type int substituted with HighsInt
// - Added return codes for errors.
// - Input matrix made const.
// - Permutation returned with zero-based numbering.
// - Input matrix provided with zero-base numbering.
//

HighsInt Highs_genrcm(HighsInt node_num, HighsInt adj_num,
                      const HighsInt adj_row[], const HighsInt adj[],
                      HighsInt perm[]);

#endif