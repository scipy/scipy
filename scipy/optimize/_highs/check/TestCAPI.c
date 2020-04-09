#include "interfaces/highs_c_api.h"

#include <stdio.h>
#include <stdlib.h>
// Force asserts to be checked always.
#undef NDEBUG
#include <assert.h>

void minimal_api() {
  int numcol = 2;
  int numrow = 2;
  int nnz = 4;
  int i;

  double cc[2] = {1.0, -2.0};
  double cl[2] = {0.0, 0.0};
  double cu[2] = {10.0, 10.0};
  double rl[2] = {0.0, 0.0};
  double ru[2] = {2.0, 1.0};
  int astart[3] = {0, 2, 4};
  int aindex[4] = {0, 1, 0, 1};
  double avalue[4] = {1.0, 2.0, 1.0, 3.0};

  double* cv = (double*)malloc(sizeof(double) * numcol);
  double* cd = (double*)malloc(sizeof(double) * numcol);
  double* rv = (double*)malloc(sizeof(double) * numrow);
  double* rd = (double*)malloc(sizeof(double) * numrow);

  int* cbs = (int*)malloc(sizeof(int) * numcol);
  int* rbs = (int*)malloc(sizeof(int) * numrow);

  int modelstatus;

  int status = Highs_call(numcol, numrow, nnz, cc, cl, cu, rl, ru, astart, aindex, avalue, cv,
            cd, rv, rd, cbs, rbs, &modelstatus);
  assert(status == 0);

  for (i = 0; i < numcol; i++) {
    printf("x%d = %lf\n", i, cv[i]);
  }

  free(cv);
  free(cd);
  free(rv);
  free(rd);
  free(cbs);
  free(rbs);
}

void full_api() {
  void* highs;

  highs = Highs_create();

  double cc[2] = {1.0, -2.0};
  double cl[2] = {0.0, 0.0};
  double cu[2] = {10.0, 10.0};
  double rl[2] = {0.0, 0.0};
  double ru[2] = {2.0, 1.0};
  int astart[3] = {0, 2, 4};
  int aindex[4] = {0, 1, 0, 1};
  double avalue[4] = {1.0, 2.0, 1.0, 3.0};

  assert( Highs_addCols(highs, 2, cc, cl, cu, 0, NULL, NULL, NULL) );
  assert( Highs_addRows(highs, 2, rl, ru,  4, astart, aindex, avalue) );

  Highs_run(highs);
  Highs_destroy(highs);
}

void options() {
  void* highs = Highs_create();

  int simplex_scale_strategy;
  Highs_setHighsIntOptionValue(highs, "simplex_scale_strategy", 0);
  Highs_getHighsIntOptionValue(highs, "simplex_scale_strategy", &simplex_scale_strategy);
  assert( simplex_scale_strategy == 0 );

  double primal_feasibility_tolerance;
  Highs_setHighsDoubleOptionValue(highs, "primal_feasibility_tolerance", 2.0);
  Highs_getHighsDoubleOptionValue(highs, "primal_feasibility_tolerance", &primal_feasibility_tolerance);
  assert( primal_feasibility_tolerance == 2.0 );

  Highs_destroy(highs);
}

int main() {
  minimal_api();
  full_api();
  options();
  return 0;
}
