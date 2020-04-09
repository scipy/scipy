#include "interfaces/highs_c_api.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// gcc call_highs_from_c.c -o highstest -I ../build/install_folder/include/ -L ../build/install_folder/lib/ -lhighs

void minimal_api() {
  // Form and solve the LP
  // Min    f  = 2x_0 + 3x_1
  // s.t.                x_1 <= 6
  //       10 <=  x_0 + 2x_1 <= 14
  //        8 <= 2x_0 +  x_1
  // 0 <= x_0 <= 3; 1 <= x_1

  int numcol = 2;
  int numrow = 3;
  int nnz = 5;
  int i;

  // Define the column costs, lower bounds and upper bounds
  double cc[2] = {2.0, 3.0};
  double cl[2] = {0.0, 1.0};
  double cu[2] = {3.0, 1.0e30};
  // Define the row lower bounds and upper bounds
  double rl[3] = {-1.0e30, 10.0, 8.0};
  double ru[3] = {6.0, 14.0, 1.0e30};
  // Define the constraint matrix column-wise
  int astart[3] = {0, 2, 5};
  int aindex[5] = {1, 2, 0, 1, 2};
  double avalue[5] = {1.0, 2.0, 1.0, 2.0, 1.0};

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

  printf("Run status = %d; Model status = %d\n", status, modelstatus);

  if (modelstatus == 9) {
    // Report the column primal and dual values, and basis status
    for (i = 0; i < numcol; i++) {
      printf("Col%d = %lf; dual = %lf; status = %d; \n", i, cv[i], cd[i], cbs[i]);
    }
    // Report the row primal and dual values, and basis status
    for (i = 0; i < numrow; i++) {
      printf("Row%d = %lf; dual = %lf; status = %d; \n", i, rv[i], rd[i], rbs[i]);
    }
  }

  free(cv);
  free(cd);
  free(rv);
  free(rd);
  free(cbs);
  free(rbs);
}

void full_api() {
  // Form and solve the LP
  // Min    f  = 2x_0 + 3x_1
  // s.t.                x_1 <= 6
  //       10 <=  x_0 + 2x_1 <= 14
  //        8 <= 2x_0 +  x_1
  // 0 <= x_0 <= 3; 1 <= x_1

  void* highs;

  highs = Highs_create();

  int numcol = 2;
  int numrow = 3;
  int nnz = 5;
  int i;

  // Define the column costs, lower bounds and upper bounds
  double cc[2] = {2.0, 3.0};
  double cl[2] = {0.0, 1.0};
  double cu[2] = {3.0, 1.0e30};
  // Define the row lower bounds and upper bounds
  double rl[3] = {-1.0e30, 10.0, 8.0};
  double ru[3] = {6.0, 14.0, 1.0e30};
  // Define the constraint matrix row-wise, as it is added to the LP
  // with the rows
  int arstart[4] = {0, 1, 3, 5};
  int arindex[5] = {1, 0, 1, 0, 1};
  double arvalue[5] = {1.0, 1.0, 2.0, 2.0, 1.0};

  double* cv = (double*)malloc(sizeof(double) * numcol);
  double* cd = (double*)malloc(sizeof(double) * numcol);
  double* rv = (double*)malloc(sizeof(double) * numrow);
  double* rd = (double*)malloc(sizeof(double) * numrow);

  int* cbs = (int*)malloc(sizeof(int) * numcol);
  int* rbs = (int*)malloc(sizeof(int) * numrow);

  // Add two columns to the empty LP
  assert( Highs_addCols(highs, 2, cc, cl, cu, 0, NULL, NULL, NULL) );
  // Add three rows to the 2-column LP
  assert( Highs_addRows(highs, 3, rl, ru,  5, arstart, arindex, arvalue) );

  int* sense;
  Highs_getObjectiveSense(highs, sense);
  printf("LP problem has objective sense = %d\n", *sense);

  *sense *= -1;
  Highs_changeObjectiveSense(highs, *sense);

  *sense *= -1;
  Highs_changeObjectiveSense(highs, *sense);

  Highs_getObjectiveSense(highs, sense);
  printf("LP problem has old objective sense = %d\n", *sense);

  int simplex_scale_strategy;
  Highs_getHighsIntOptionValue(highs, "simplex_scale_strategy", &simplex_scale_strategy);
  printf("simplex_scale_strategy = %d: setting it to 3\n", simplex_scale_strategy);
  simplex_scale_strategy = 3;
  Highs_setHighsIntOptionValue(highs, "simplex_scale_strategy", simplex_scale_strategy);

  double primal_feasibility_tolerance;
  Highs_getHighsDoubleOptionValue(highs, "primal_feasibility_tolerance", &primal_feasibility_tolerance);
  printf("primal_feasibility_tolerance = %g: setting it to 1e-6\n", primal_feasibility_tolerance);
  primal_feasibility_tolerance = 1e-6;
  Highs_setHighsDoubleOptionValue(highs, "primal_feasibility_tolerance", primal_feasibility_tolerance);

  int status = Highs_run(highs);
  // Get the model status
  const int scaled_model = 0;
  int modelstatus = Highs_getModelStatus(highs, scaled_model);

  printf("Run status = %d; Model status = %d\n", status, modelstatus);

  double objective_function_value;
  Highs_getHighsDoubleInfoValue(highs, "objective_function_value", &objective_function_value);
  int simplex_iteration_count = 0;
  Highs_getHighsIntInfoValue(highs, "simplex_iteration_count", &simplex_iteration_count);

  printf("Objective value = %g; Iteration count = %d\n", objective_function_value, simplex_iteration_count);
  if (modelstatus == 9) {
    // Get the primal and dual solution 
    Highs_getSolution(highs, cv, cd, rv, rd);
    // Get the basis
    Highs_getBasis(highs, cbs, rbs);
    // Report the column primal and dual values, and basis status
    for (i = 0; i < numcol; i++) {
      printf("Col%d = %lf; dual = %lf; status = %d; \n", i, cv[i], cd[i], cbs[i]);
    }
    // Report the row primal and dual values, and basis status
    for (i = 0; i < numrow; i++) {
      printf("Row%d = %lf; dual = %lf; status = %d; \n", i, rv[i], rd[i], rbs[i]);
    }
  }

  free(cv);
  free(cd);
  free(rv);
  free(rd);
  free(cbs);
  free(rbs);

  Highs_destroy(highs);
}

int main() { 
  minimal_api();
  full_api();
  return 0; 
}
