#include "interfaces/highs_c_api.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// gcc call_highs_from_c.c -o highstest -I ../build/install_folder/include/ -L ../build/install_folder/lib/ -lhighs

void minimal_api() {
  // This illustrates the use of Highs_call, the simple C interface to
  // HiGHS. It's designed to solve the general LP problem
  //
  // Min c^Tx subject to L <= Ax <= U; l <= x <= u
  //
  // where A is a matrix with m rows and n columns
  //
  // The scalar n is numcol
  // The scalar m is numrow
  //
  // The vector c is colcost
  // The vector l is collower
  // The vector u is colupper
  // The vector L is rowlower
  // The vector U is rowupper
  //
  // The matrix A is represented in packed column-wise form: only its
  // nonzeros are stored
  // 
  // * The number of nonzeros in A is nnz
  // 
  // * The row indices of the nonnzeros in A are stored column-by-column
  // in aindex
  // 
  // * The values of the nonnzeros in A are stored column-by-column in
  // avalue
  // 
  // * The position in aindex/avalue of the index/value of the first
  // nonzero in each column is stored in astart
  // 
  // Note that astart[0] must be zero
  //
  // After a successful call to Highs_call, the primal and dual
  // solution, and the simplex basis are returned as follows
  //
  // The vector x is colvalue
  // The vector Ax is rowvalue
  // The vector of dual values for the variables x is coldual
  // The vector of dual values for the variables Ax is rowdual
  // The basic/nonbasic status of the variables x is colbasisstatus
  // The basic/nonbasic status of the variables Ax is rowbasisstatus
  //
  // The status of the solution obtained is modelstatus
  //
  // To solve maximization problems, the values in c must be negated
  //
  // The use of Highs_call is illustrated for the LP
  //
  // Min    f  = 2x_0 + 3x_1
  // s.t.                x_1 <= 6
  //       10 <=  x_0 + 2x_1 <= 14
  //        8 <= 2x_0 +  x_1
  // 0 <= x_0 <= 3; 1 <= x_1

  int numcol = 2;
  int numrow = 3;
  int nnz = 5;

  // Define the column costs, lower bounds and upper bounds
  double colcost[numcol] = {2.0, 3.0};
  double collower[numcol] = {0.0, 1.0};
  double colupper[numcol] = {3.0, 1.0e30};
  // Define the row lower bounds and upper bounds
  double rowlower[numrow] = {-1.0e30, 10.0, 8.0};
  double rowupper[numrow] = {6.0, 14.0, 1.0e30};
  // Define the constraint matrix column-wise
  int astart[numcol] = {0, 2};
  int aindex[nnz] = {1, 2, 0, 1, 2};
  double avalue[nnz] = {1.0, 2.0, 1.0, 2.0, 1.0};

  double* colvalue = (double*)malloc(sizeof(double) * numcol);
  double* coldual = (double*)malloc(sizeof(double) * numcol);
  double* rowvalue = (double*)malloc(sizeof(double) * numrow);
  double* rowdual = (double*)malloc(sizeof(double) * numrow);

  int* colbasisstatus = (int*)malloc(sizeof(int) * numcol);
  int* rowbasisstatus = (int*)malloc(sizeof(int) * numrow);

  int modelstatus; 

  int status = Highs_call(numcol, numrow, nnz,
			  colcost, collower, colupper,
			  rowlower, rowupper,
			  astart, aindex, avalue,
			  colvalue, coldual, rowvalue, rowdual,
			  colbasisstatus, rowbasisstatus,
			  &modelstatus);
            
  assert(status == 0);

  printf("Run status = %d; Model status = %d\n", status, modelstatus);

  int i;
  if (modelstatus == 9) {
    // Report the column primal and dual values, and basis status
    for (i = 0; i < numcol; i++) {
      printf("Col%d = %lf; dual = %lf; status = %d; \n", i, colvalue[i], coldual[i], colbasisstatus[i]);
    }
    // Report the row primal and dual values, and basis status
    for (i = 0; i < numrow; i++) {
      printf("Row%d = %lf; dual = %lf; status = %d; \n", i, rowvalue[i], rowdual[i], rowbasisstatus[i]);
    }
  }

  free(colvalue);
  free(coldual);
  free(rowvalue);
  free(rowdual);
  free(colbasisstatus);
  free(rowbasisstatus);
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
  double colcost[numcol] = {2.0, 3.0};
  double collower[numcol] = {0.0, 1.0};
  double colupper[numcol] = {3.0, 1.0e30};
  // Define the row lower bounds and upper bounds
  double rowlower[numrow] = {-1.0e30, 10.0, 8.0};
  double rowupper[numrow] = {6.0, 14.0, 1.0e30};
  // Define the constraint matrix row-wise, as it is added to the LP
  // with the rows
  int arstart[numrow] = {0, 1, 3};
  int arindex[nnz] = {1, 0, 1, 0, 1};
  double arvalue[nnz] = {1.0, 1.0, 2.0, 2.0, 1.0};

  double* colvalue = (double*)malloc(sizeof(double) * numcol);
  double* coldual = (double*)malloc(sizeof(double) * numcol);
  double* rowvalue = (double*)malloc(sizeof(double) * numrow);
  double* rowdual = (double*)malloc(sizeof(double) * numrow);

  int* colbasisstatus = (int*)malloc(sizeof(int) * numcol);
  int* rowbasisstatus = (int*)malloc(sizeof(int) * numrow);

  // Add two columns to the empty LP
  assert( Highs_addCols(highs, numcol, colcost, collower, colupper, 0, NULL, NULL, NULL) );
  // Add three rows to the 2-column LP
  assert( Highs_addRows(highs, numrow, rowlower, rowupper, nnz, arstart, arindex, arvalue) );

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
    Highs_getSolution(highs, colvalue, coldual, rowvalue, rowdual);
    // Get the basis
    Highs_getBasis(highs, colbasisstatus, rowbasisstatus);
    // Report the column primal and dual values, and basis status
    for (i = 0; i < numcol; i++) {
      printf("Col%d = %lf; dual = %lf; status = %d; \n", i, colvalue[i], coldual[i], colbasisstatus[i]);
    }
    // Report the row primal and dual values, and basis status
    for (i = 0; i < numrow; i++) {
      printf("Row%d = %lf; dual = %lf; status = %d; \n", i, rowvalue[i], rowdual[i], rowbasisstatus[i]);
    }
  }

  free(colvalue);
  free(coldual);
  free(rowvalue);
  free(rowdual);
  free(colbasisstatus);
  free(rowbasisstatus);

  Highs_destroy(highs);
}

int main() { 
  minimal_api();
  full_api();
  return 0; 
}
