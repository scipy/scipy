#include <algorithm>

#include "HConfig.h"
#include "Highs.h"
#include "HighsRandom.h"
#include "catch.hpp"

bool GetBasisSolvesSolutionNzOk(int numRow, double* pass_solution_vector,
                                int* solution_num_nz, int* solution_indices) {
  if (solution_num_nz == NULL) return true;
  double* solution_vector = (double*)malloc(sizeof(double) * numRow);
  bool solution_nz_ok = true;
  for (int row = 0; row < numRow; row++)
    solution_vector[row] = pass_solution_vector[row];
  // Check that the indexed entries are nonzero
  for (int ix = 0; ix < *solution_num_nz; ix++) {
    int row = solution_indices[ix];
    if (!solution_vector[row]) {
      printf("SolutionNzOk: Indexed entry solution_vector[%2d] = %11.4g\n", row,
             solution_vector[row]);
      solution_nz_ok = false;
    } else {
      solution_vector[row] = 0;
    }
  }
  // Solution should now be zero
  for (int row = 0; row < numRow; row++) {
    if (solution_vector[row]) {
      printf("SolutionNzOk: Non-indexed entry solution_vector[%2d] = %11.4g\n",
             row, solution_vector[row]);
      solution_nz_ok = false;
    }
  }
  free(solution_vector);
  return solution_nz_ok;
}
double GetBasisSolvesCheckSolution(HighsLp& lp, int* basic_variables,
                                   double* rhs, double* solution,
                                   const bool transpose = false) {
  const double residual_tolerance = 1e-8;
  double residual_norm = 0;
  //  for (int k=0; k<lp.numRow_; k++) printf("solution[%2d]=%11.4g\n", k,
  //  solution[k]);
  if (transpose) {
    for (int k = 0; k < lp.numRow_; k++) {
      double residual = 0;
      int var = basic_variables[k];
      if (var < 0) {
        int row = -(1 + var);
        residual = fabs(rhs[k] - solution[row]);
        if (residual > residual_tolerance)
          printf("Row |[B^Tx-b]_{%2d}| = %11.4g\n", k, residual);
      } else {
        int col = var;
        for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
          int row = lp.Aindex_[el];
          residual += lp.Avalue_[el] * solution[row];
          //	  printf("k=%1d; col=%1d; el=%1d; row=%1d;
          // lp.Avalue_[col]=%11.4g; solution[row]=%11.4g; residual=%1.4g\n", k,
          // col, el, row, lp.Avalue_[col], solution[row], residual);
        }
        residual = fabs(rhs[k] - residual);
        if (residual > residual_tolerance)
          printf("Col |[B^Tx-b]_{%2d}| = %11.4g\n", k, residual);
      }
      residual_norm += residual;
    }
  } else {
    vector<double> basis_matrix_times_solution;
    basis_matrix_times_solution.assign(lp.numRow_, 0);
    for (int k = 0; k < lp.numRow_; k++) {
      int var = basic_variables[k];
      if (var < 0) {
        int row = -(1 + var);
        basis_matrix_times_solution[row] += solution[k];
      } else {
        int col = var;
        for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
          int row = lp.Aindex_[el];
          basis_matrix_times_solution[row] += lp.Avalue_[el] * solution[k];
        }
      }
    }
    for (int k = 0; k < lp.numRow_; k++) {
      double residual = fabs(rhs[k] - basis_matrix_times_solution[k]);
      if (residual > residual_tolerance)
        printf("|[B^Tx-b]_{%2d}| = %11.4g\n", k, residual);
      residual_norm += residual;
    }
  }
  return residual_norm;
}

void GetBasisSolvesFormRHS(HighsLp& lp, int* basic_variables, double* solution,
                           double* rhs, const bool transpose = false) {
  if (transpose) {
    for (int k = 0; k < lp.numRow_; k++) {
      rhs[k] = 0;
      int var = basic_variables[k];
      if (var < 0) {
        int row = -(1 + var);
        rhs[k] = solution[row];
      } else {
        int col = var;
        for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
          int row = lp.Aindex_[el];
          rhs[k] += lp.Avalue_[el] * solution[row];
        }
      }
    }
  } else {
    for (int k = 0; k < lp.numRow_; k++) rhs[k] = 0;
    for (int k = 0; k < lp.numRow_; k++) {
      int var = basic_variables[k];
      if (var < 0) {
        int row = -(1 + var);
        rhs[row] += solution[k];
      } else {
        int col = var;
        for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
          int row = lp.Aindex_[el];
          rhs[row] += lp.Avalue_[el] * solution[k];
        }
      }
    }
  }
}

// No commas in test case name.
TEST_CASE("Basis-solves", "[highs_basis_solves]") {
  std::cout << std::string(HIGHS_DIR) << std::endl;

  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/chip.mps";
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  //  filename = std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";

  // printf("CMAKE %s\n", HIGHS_DIR);

  Highs highs;

  int* basic_variables = nullptr;
  double* rhs = nullptr;
  double* known_solution;
  double* solution_vector = nullptr;
  int solution_num_nz;
  int* solution_indices = (int*)malloc(sizeof(int) * 1);

  HighsStatus highs_status;

  highs_status = highs.getBasicVariables(basic_variables);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getBasisInverseRow(0, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getBasisInverseCol(0, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getBasisSolve(rhs, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getBasisTransposeSolve(rhs, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getReducedRow(0, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getReducedColumn(0, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.readModel(filename);
  REQUIRE(highs_status == HighsStatus::OK);

  HighsLp lp = highs.getLp();
  REQUIRE(highs_status == HighsStatus::OK);

  highs_status = highs.writeModel("");
  REQUIRE(highs_status == HighsStatus::OK);

  int numRow = lp.numRow_;
  int numCol = lp.numCol_;
  int check_row = 0;
  int check_col = 0;

  basic_variables = (int*)malloc(sizeof(int) * numRow);
  known_solution = (double*)malloc(sizeof(double) * numRow);
  solution_vector = (double*)malloc(sizeof(double) * numRow);
  solution_indices = (int*)malloc(sizeof(int) * numRow);
  rhs = (double*)malloc(sizeof(double) * numRow);

  highs_status = highs.getBasicVariables(basic_variables);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getBasisInverseRow(check_row, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getBasisInverseCol(0, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getBasisSolve(rhs, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getBasisTransposeSolve(rhs, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getReducedRow(0, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.getReducedColumn(0, solution_vector);
  REQUIRE(highs_status == HighsStatus::Error);

  highs_status = highs.run();
  REQUIRE(highs_status == HighsStatus::OK);

  highs_status = highs.getBasicVariables(basic_variables);
  REQUIRE(highs_status == HighsStatus::OK);

  /*
  for (int row=0; row < numRow; row++) {
    printf("Basic variable %3d is ", row);
    int var = basic_variables[row];
    if (var<0) {
      printf("row %d\n", -(1+var));
    } else {
      printf("col %d\n", var);
    }
  }
  */

  double residual_norm;
  const double residual_norm_tolerance = 1e-8;
  const double solution_error_tolerance = 1e-8;
  HighsRandom random;

  int basic_col;

  for (int ix = 0; ix < numRow; ix++) known_solution[ix] = 0;
  bool transpose = true;
  int num_ix = 3;
  int col;
  col = 6;
  basic_col = basic_variables[col];
  known_solution[col] = random.fraction();
  //  printf("Known solution col %2d is basic_col %2d\n", col, basic_col);

  if (num_ix > 1) {
    col = 15;
    basic_col = basic_variables[col];
    known_solution[col] = random.fraction();
    //    printf("Known solution col %2d is basic_col %2d\n", col, basic_col);
  }

  if (num_ix > 2) {
    col = 12;
    basic_col = basic_variables[col];
    known_solution[col] = random.fraction();
    //    printf("Known solution col %2d is basic_col %2d\n", col, basic_col);
  }

  GetBasisSolvesFormRHS(lp, basic_variables, known_solution, rhs, transpose);
  if (transpose) {
    highs_status = highs.getBasisTransposeSolve(rhs, solution_vector);
  } else {
    highs_status = highs.getBasisSolve(rhs, solution_vector);
  }
  REQUIRE(highs_status == HighsStatus::OK);
  residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                              solution_vector, transpose);
  REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
  double solution_error_norm = 0;
  for (int ix = 0; ix < numRow; ix++) {
    double solution_error = fabs(known_solution[ix] - solution_vector[ix]);
    if (solution_error > solution_error_tolerance)
      printf("Row %2d: |x-x^|_i = %11.4g\n", ix, solution_error);
    solution_error_norm += solution_error;
  }
  printf(
      "Test 0:     residual_norm = %11.4g\n      solution_error_norm = %11.4g "
      "(Known solution)\n",
      residual_norm, solution_error_norm);

  double max_residual_norm;
  int max_k = min(numRow, 9);
  int k;

  k = 0;
  max_residual_norm = 0;
  for (int row = 0; row < numRow; row++) {
    int var = basic_variables[row];
    if (var >= 0) {
      basic_col = var;
      // int rhs_nnz = lp.Astart_[basic_col+1]-lp.Astart_[basic_col];
      //      printf("Row %2d; Var %3d; RHS nnz = %d\n", row, basic_col,
      //      rhs_nnz);
      for (int ix = 0; ix < numRow; ix++) rhs[ix] = 0;
      for (int el = lp.Astart_[basic_col]; el < lp.Astart_[basic_col + 1]; el++)
        rhs[lp.Aindex_[el]] = lp.Avalue_[el];

      highs_status = highs.getBasisSolve(rhs, solution_vector, &solution_num_nz,
                                         solution_indices);
      REQUIRE(highs_status == HighsStatus::OK);
      bool solution_nz_ok = GetBasisSolvesSolutionNzOk(
          numRow, solution_vector, &solution_num_nz, solution_indices);
      REQUIRE(solution_nz_ok == true);
      residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                  solution_vector, false);
      max_residual_norm = std::max(residual_norm, max_residual_norm);
      if (residual_norm > residual_norm_tolerance)
        printf("getBasisSolve(%d): residual_norm = %g\n", k, residual_norm);
      REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
      if (k < max_k)
        k++;
      else
        k *= 2;
    }
    if (k >= numRow) break;
  }
  printf("Test 1: max_residual_norm = %11.4g (Basic column)\n",
         max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  for (;;) {
    check_row = k;
    // Determine row check_row of B^{-1}
    highs_status = highs.getBasisInverseRow(check_row, solution_vector,
                                            &solution_num_nz, solution_indices);
    REQUIRE(highs_status == HighsStatus::OK);
    bool solution_nz_ok = GetBasisSolvesSolutionNzOk(
        numRow, solution_vector, &solution_num_nz, solution_indices);
    REQUIRE(solution_nz_ok == true);
    // Check solution
    // Set up RHS as e_{check_row}
    for (int row = 0; row < numRow; row++) rhs[row] = 0;
    rhs[check_row] = 1;
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_vector, true);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance)
      printf("getBasisInverseRow(%d): residual_norm = %g\n", k, residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  printf("Test 2: max_residual_norm = %11.4g (getBasisInverseRow)\n",
         max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  for (;;) {
    check_col = k;
    // Determine col check_col of B^{-1}
    highs_status = highs.getBasisInverseCol(check_col, solution_vector,
                                            &solution_num_nz, solution_indices);
    REQUIRE(highs_status == HighsStatus::OK);
    bool solution_nz_ok = GetBasisSolvesSolutionNzOk(
        numRow, solution_vector, &solution_num_nz, solution_indices);
    REQUIRE(solution_nz_ok == true);
    // Check solution
    // Set up RHS as e_{check_col}
    for (int row = 0; row < numRow; row++) rhs[row] = 0;
    rhs[check_col] = 1;
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_vector, false);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance)
      printf("getBasisInverseCol(%d): residual_norm = %g\n", k, residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  printf("Test 3: max_residual_norm = %11.4g (getBasisInverseCol)\n",
         max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  for (;;) {
    for (int row = 0; row < numRow; row++) rhs[row] = random.fraction();
    highs_status = highs.getBasisSolve(rhs, solution_vector);
    REQUIRE(highs_status == HighsStatus::OK);
    // Check solution
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_vector, false);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance)
      printf("getBasisSolve(%d): residual_norm = %g\n", k, residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  printf("Test 4: max_residual_norm = %11.4g (getBasisSolve)\n",
         max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  for (;;) {
    for (int row = 0; row < numRow; row++) rhs[row] = random.fraction();
    highs_status = highs.getBasisTransposeSolve(rhs, solution_vector);
    REQUIRE(highs_status == HighsStatus::OK);
    // Check solution
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_vector, true);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance)
      printf("getBasisTransposeSolve(%d): residual_norm = %g\n", k,
             residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  printf("Test 5: max_residual_norm = %11.4g (getBasisTransposeSolve)\n",
         max_residual_norm);

  if (numCol > numRow) {
    solution_vector = (double*)malloc(sizeof(double) * numCol);
    solution_indices = (int*)malloc(sizeof(int) * numCol);
  }

  k = 0;
  max_residual_norm = 0;
  max_k = min(numRow, 9);
  for (;;) {
    check_row = k;
    highs_status = highs.getReducedRow(check_row, solution_vector,
                                       &solution_num_nz, solution_indices);
    REQUIRE(highs_status == HighsStatus::OK);
    bool solution_nz_ok = GetBasisSolvesSolutionNzOk(
        numCol, solution_vector, &solution_num_nz, solution_indices);
    REQUIRE(solution_nz_ok == true);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  printf("Test 6: max_residual_norm = %11.4g (getReducedRow)\n",
         max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  max_k = min(numCol, 9);
  for (;;) {
    check_col = k;
    highs_status = highs.getReducedColumn(check_col, solution_vector,
                                          &solution_num_nz, solution_indices);
    REQUIRE(highs_status == HighsStatus::OK);
    // Check solution
    for (int row = 0; row < numRow; row++) rhs[row] = 0;
    for (int el = lp.Astart_[check_col]; el < lp.Astart_[check_col + 1]; el++)
      rhs[lp.Aindex_[el]] = lp.Avalue_[el];
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_vector, false);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance)
      printf("getBasisTransposeSolve(%d): residual_norm = %g\n", k,
             residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numCol) break;
  }
  printf("Test 7: max_residual_norm = %11.4g (getReducedColumn)\n",
         max_residual_norm);
}
