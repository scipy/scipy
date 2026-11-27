/*
 * Adapted from QDLDL (https://github.com/osqp/qdldl).
 *
 * Copyright 2018, Paul Goulart, Bartolomeo Stellato, Goran Banjac, The OSQP
 * developers
 *
 * Licensed under the Apache License, Version 2.0. See qdldl.h for the notice
 * and the list of changes.
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#include "qdldl/qdldl.h"

#include <cassert>

void Highs_qdldl_etree(HighsInt n, const HighsInt* ptr, const HighsInt* rows,
                       HighsInt* parent, int64_t* colcount, HighsInt* work) {
  for (HighsInt i = 0; i < n; ++i) {
    parent[i] = -1;
    colcount[i] = 1;
    work[i] = 0;
  }

  for (HighsInt col = 0; col < n; ++col) {
    work[col] = col;
    for (HighsInt el = ptr[col]; el < ptr[col + 1]; ++el) {
      HighsInt row = rows[el];
      while (work[row] != col) {
        if (parent[row] < 0) parent[row] = col;
        ++colcount[row];
        work[row] = col;
        row = parent[row];
      }
    }
  }
}

double Highs_qdldl_factor_row(HighsInt row, const HighsInt* ptr,
                              const HighsInt* rows, const double* val,
                              const HighsInt* parent, const int64_t* ptrL,
                              HighsInt* rowsL, double* valL, int64_t* next,
                              uint8_t* mark, HighsInt* stack,
                              HighsInt* revpattern, double* x, int64_t* work) {
  // Sparsity pattern of the row of L, by traversing the row subtree.
  HighsInt rowcount = 0;
  for (HighsInt el = ptr[row]; el < ptr[row + 1]; ++el) {
    HighsInt i = rows[el];
    HighsInt top = 0;

    assert(i <= row);

    while (i >= 0 && i < row && !mark[i]) {
      mark[i] = 1;
      stack[top++] = i;
      i = parent[i];
    }

    while (top > 0) revpattern[rowcount++] = stack[--top];
  }

  // Values of the row of L: l = Dinv * L^-1 * b, computed in x.
  x[row] = 0.0;
  for (HighsInt el = ptr[row]; el < ptr[row + 1]; ++el) x[rows[el]] = val[el];
  double d = x[row];
  x[row] = 0.0;

  for (HighsInt i = 0; i < rowcount; ++i) {
    const HighsInt col = revpattern[rowcount - 1 - i];

    double value = x[col];
    x[col] = 0.0;
    assert(rowsL[ptrL[col]] == col);
    for (int64_t el = ptrL[col] + 1; el < next[col]; ++el)
      x[rowsL[el]] -= valL[el] * value;
    *work += next[col] - ptrL[col];

    const double Dinv = valL[ptrL[col]];
    rowsL[next[col]] = row;
    valL[next[col]] = value * Dinv;

    d -= value * valL[next[col]];
    next[col]++;

    mark[col] = 0;
  }

  return d;
}
