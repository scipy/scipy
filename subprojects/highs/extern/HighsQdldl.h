/*
 * Interface to the QDLDL-derived kernels bundled with HiGHS, which are
 * compiled only when HIPO is enabled.
 *
 * This header lives outside extern/qdldl so that it is available even when
 * extern/qdldl is absent (e.g. SciPy trims it from its sdist), which keeps
 * HighsExtrasApi layout-compatible with an extras library built with HIPO.
 * It is also included by extern/qdldl/qdldl.h, so is the only definition of
 * the QDLDL interface used by HiGHS.
 */
#ifndef HIGHS_QDLDL_INTERFACE_H
#define HIGHS_QDLDL_INTERFACE_H

#include <cstdint>

#include "HConfig.h"

#ifdef HIGHSINT64
typedef int64_t HighsInt;
#else
typedef int HighsInt;
#endif

// Elimination tree (parent) and column counts of L, including the diagonal,
// of the n-by-n upper triangular matrix (ptr, rows). work has size n.
void Highs_qdldl_etree(HighsInt n, const HighsInt* ptr, const HighsInt* rows,
                       HighsInt* parent, int64_t* colcount, HighsInt* work);

// Computes row `row` of L from column `row` of the upper triangular matrix
// (ptr, rows, val), given rows 0..row-1 of L. Entries are appended to columns
// of L at next[col]; the diagonal of L at ptrL[col] holds D^-1. mark (all zero
// on entry and return), stack, revpattern and x (all zero on entry and return)
// have size n. Adds the number of entries of L scanned to *work. Returns the
// pivot d of the row; the caller stores 1/d and the diagonal index.
double Highs_qdldl_factor_row(HighsInt row, const HighsInt* ptr,
                              const HighsInt* rows, const double* val,
                              const HighsInt* parent, const int64_t* ptrL,
                              HighsInt* rowsL, double* valL, int64_t* next,
                              uint8_t* mark, HighsInt* stack,
                              HighsInt* revpattern, double* x, int64_t* work);

#endif /* HIGHS_QDLDL_INTERFACE_H */
