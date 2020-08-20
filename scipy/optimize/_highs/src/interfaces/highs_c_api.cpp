/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "highs_c_api.h"

#include "Highs.h"

int Highs_call(int numcol, int numrow, int numnz, double* colcost,
               double* collower, double* colupper, double* rowlower,
               double* rowupper, int* astart, int* aindex, double* avalue,
               double* colvalue, double* coldual, double* rowvalue,
               double* rowdual, int* colbasisstatus, int* rowbasisstatus,
               int* modelstatus) {
  Highs highs;

  int status =
      Highs_passLp(&highs, numcol, numrow, numnz, colcost, collower, colupper,
                   rowlower, rowupper, astart, aindex, avalue);
  if (status != 0) {
    return status;
  }

  status = (int)highs.run();

  if (status == 0) {
    HighsSolution solution;
    HighsBasis basis;
    solution = highs.getSolution();
    basis = highs.getBasis();
    *modelstatus = (int)highs.getModelStatus();

    for (int i = 0; i < numcol; i++) {
      colvalue[i] = solution.col_value[i];
      coldual[i] = solution.col_dual[i];

      colbasisstatus[i] = (int)basis.col_status[i];
    }

    for (int i = 0; i < numrow; i++) {
      rowvalue[i] = solution.row_value[i];
      rowdual[i] = solution.row_dual[i];

      rowbasisstatus[i] = (int)basis.row_status[i];
    }
  }

  return status;
}

void* Highs_create() { return new Highs(); }

void Highs_destroy(void* highs) { delete (Highs*)highs; }

int Highs_run(void* highs) { return (int)((Highs*)highs)->run(); }

int Highs_readModel(void* highs, const char* filename) {
  return (int)((Highs*)highs)->readModel(std::string(filename));
}

int Highs_writeModel(void* highs, const char* filename) {
  return (int)((Highs*)highs)->writeModel(std::string(filename));
}

int Highs_writeSolution(void* highs, const char* filename) {
  return (int)((Highs*)highs)->writeSolution(std::string(filename));
}

int Highs_passLp(void* highs, int numcol, int numrow, int numnz,
                 double* colcost, double* collower, double* colupper,
                 double* rowlower, double* rowupper, int* astart, int* aindex,
                 double* avalue) {
  HighsLp lp;
  lp.numCol_ = numcol;
  lp.numRow_ = numrow;

  lp.colCost_.assign(colcost, colcost + numcol);
  lp.colLower_.assign(collower, collower + numcol);
  lp.colUpper_.assign(colupper, colupper + numcol);

  lp.rowLower_.assign(rowlower, rowlower + numrow);
  lp.rowUpper_.assign(rowupper, rowupper + numrow);

  lp.Astart_.assign(astart, astart + numcol);
  lp.Astart_.resize(numcol + 1);
  lp.Astart_[numcol] = numnz;
  lp.Aindex_.assign(aindex, aindex + numnz);
  lp.Avalue_.assign(avalue, avalue + numnz);

  return (int)((Highs*)highs)->passModel(lp);
}

int Highs_clearModel(void* highs) { return (int)((Highs*)highs)->clearModel(); }

int Highs_runQuiet(void* highs) {
  int return_status = (int)((Highs*)highs)->setHighsLogfile(NULL);
  if (return_status) return return_status;
  return (int)((Highs*)highs)->setHighsOutput(NULL);
}

int Highs_setHighsBoolOptionValue(void* highs, const char* option,
                                  const int value) {
  return (int)((Highs*)highs)
      ->setHighsOptionValue(std::string(option), (bool)value);
}

int Highs_setHighsIntOptionValue(void* highs, const char* option,
                                 const int value) {
  return (int)((Highs*)highs)->setHighsOptionValue(std::string(option), value);
}

int Highs_setHighsDoubleOptionValue(void* highs, const char* option,
                                    const double value) {
  return (int)((Highs*)highs)->setHighsOptionValue(std::string(option), value);
}

int Highs_setHighsStringOptionValue(void* highs, const char* option,
                                    const char* value) {
  return (int)((Highs*)highs)
      ->setHighsOptionValue(std::string(option), std::string(value));
}

int Highs_setHighsOptionValue(void* highs, const char* option,
                              const char* value) {
  return (int)((Highs*)highs)
      ->setHighsOptionValue(std::string(option), std::string(value));
}

int Highs_getHighsBoolOptionValue(void* highs, const char* option, int* value) {
  bool v;
  int retcode =
      (int)((Highs*)highs)->getHighsOptionValue(std::string(option), v);
  *value = (int)v;
  return retcode;
}

int Highs_getHighsIntOptionValue(void* highs, const char* option, int* value) {
  return (int)((Highs*)highs)->getHighsOptionValue(std::string(option), *value);
}

int Highs_getHighsDoubleOptionValue(void* highs, const char* option,
                                    double* value) {
  return (int)((Highs*)highs)->getHighsOptionValue(std::string(option), *value);
}

int Highs_getHighsStringOptionValue(void* highs, const char* option,
                                    char* value) {
  std::string v;
  int retcode =
      (int)((Highs*)highs)->getHighsOptionValue(std::string(option), v);
  strcpy(value, v.c_str());
  return retcode;
}

int Highs_resetHighsOptions(void* highs) {
  return (int)((Highs*)highs)->resetHighsOptions();
}

int Highs_getHighsIntInfoValue(void* highs, const char* info, int* value) {
  return (int)((Highs*)highs)->getHighsInfoValue(info, *value);
}

int Highs_getHighsDoubleInfoValue(void* highs, const char* info,
                                  double* value) {
  return (int)((Highs*)highs)->getHighsInfoValue(info, *value);
}

void Highs_getSolution(void* highs, double* colvalue, double* coldual,
                       double* rowvalue, double* rowdual) {
  HighsSolution solution = ((Highs*)highs)->getSolution();

  for (int i = 0; i < (int)solution.col_value.size(); i++) {
    colvalue[i] = solution.col_value[i];
  }

  for (int i = 0; i < (int)solution.col_dual.size(); i++) {
    coldual[i] = solution.col_dual[i];
  }

  for (int i = 0; i < (int)solution.row_value.size(); i++) {
    rowvalue[i] = solution.row_value[i];
  }

  for (int i = 0; i < (int)solution.row_dual.size(); i++) {
    rowdual[i] = solution.row_dual[i];
  }
}

void Highs_getBasis(void* highs, int* colstatus, int* rowstatus) {
  HighsBasis basis = ((Highs*)highs)->getBasis();
  for (int i = 0; i < (int)basis.col_status.size(); i++) {
    colstatus[i] = (int)basis.col_status[i];
  }

  for (int i = 0; i < (int)basis.row_status.size(); i++) {
    rowstatus[i] = (int)basis.row_status[i];
  }
}

int Highs_getModelStatus(void* highs, const int scaled_model) {
  return (int)((Highs*)highs)->getModelStatus(scaled_model);
}

double Highs_getObjectiveValue(void* highs) {
  return ((Highs*)highs)->getObjectiveValue();
}

int Highs_getIterationCount(void* highs) {
  return Highs_getSimplexIterationCount(highs);
}

int Highs_getSimplexIterationCount(void* highs) {
  return (int)((Highs*)highs)->getSimplexIterationCount();
}

int Highs_getBasicVariables(void* highs, int* basic_variables) {
  return (int)((Highs*)highs)->getBasicVariables(basic_variables);
}

int Highs_getBasisInverseRow(void* highs, const int row, double* row_vector,
                             int* row_num_nz, int* row_indices) {
  return (int)((Highs*)highs)
      ->getBasisInverseRow(row, row_vector, row_num_nz, row_indices);
}

int Highs_getBasisInverseCol(void* highs, const int col, double* col_vector,
                             int* col_num_nz, int* col_indices) {
  return (int)((Highs*)highs)
      ->getBasisInverseCol(col, col_vector, col_num_nz, col_indices);
}

int Highs_getBasisSolve(void* highs, const double* rhs, double* solution_vector,
                        int* solution_num_nz, int* solution_indices) {
  return (int)((Highs*)highs)
      ->getBasisSolve(rhs, solution_vector, solution_num_nz, solution_indices);
}

int Highs_getBasisTransposeSolve(void* highs, const double* rhs,
                                 double* solution_vector, int* solution_nz,
                                 int* solution_indices) {
  return (int)((Highs*)highs)
      ->getBasisTransposeSolve(rhs, solution_vector, solution_nz,
                               solution_indices);
}

int Highs_getReducedRow(void* highs, const int row, double* row_vector,
                        int* row_num_nz, int* row_indices) {
  return (int)((Highs*)highs)
      ->getReducedRow(row, row_vector, row_num_nz, row_indices);
}

int Highs_getReducedColumn(void* highs, const int col, double* col_vector,
                           int* col_num_nz, int* col_indices) {
  return (int)((Highs*)highs)
      ->getReducedColumn(col, col_vector, col_num_nz, col_indices);
}

int Highs_setBasis(void* highs, const int* colstatus, const int* rowstatus) {
  HighsBasis basis;
  const int num_col = Highs_getNumCols(highs);
  basis.col_status.resize(num_col);
  for (int i = 0; i < num_col; i++) {
    if (colstatus[i] == (int)HighsBasisStatus::LOWER) {
      basis.col_status[i] = HighsBasisStatus::LOWER;
    } else if (colstatus[i] == (int)HighsBasisStatus::BASIC) {
      basis.col_status[i] = HighsBasisStatus::BASIC;
    } else if (colstatus[i] == (int)HighsBasisStatus::UPPER) {
      basis.col_status[i] = HighsBasisStatus::UPPER;
    } else if (colstatus[i] == (int)HighsBasisStatus::ZERO) {
      basis.col_status[i] = HighsBasisStatus::ZERO;
    } else if (colstatus[i] == (int)HighsBasisStatus::NONBASIC) {
      basis.col_status[i] = HighsBasisStatus::NONBASIC;
    } else if (colstatus[i] == (int)HighsBasisStatus::SUPER) {
      basis.col_status[i] = HighsBasisStatus::SUPER;
    } else {
      return (int)HighsStatus::Error;
    }
  }
  const int num_row = Highs_getNumRows(highs);
  basis.row_status.resize(num_row);
  for (int i = 0; i < num_row; i++) {
    if (rowstatus[i] == (int)HighsBasisStatus::LOWER) {
      basis.row_status[i] = HighsBasisStatus::LOWER;
    } else if (rowstatus[i] == (int)HighsBasisStatus::BASIC) {
      basis.row_status[i] = HighsBasisStatus::BASIC;
    } else if (rowstatus[i] == (int)HighsBasisStatus::UPPER) {
      basis.row_status[i] = HighsBasisStatus::UPPER;
    } else if (rowstatus[i] == (int)HighsBasisStatus::ZERO) {
      basis.row_status[i] = HighsBasisStatus::ZERO;
    } else if (rowstatus[i] == (int)HighsBasisStatus::NONBASIC) {
      basis.row_status[i] = HighsBasisStatus::NONBASIC;
    } else if (rowstatus[i] == (int)HighsBasisStatus::SUPER) {
      basis.row_status[i] = HighsBasisStatus::SUPER;
    } else {
      return (int)HighsStatus::Error;
    }
  }
  return (int)((Highs*)highs)->setBasis(basis);
}

int Highs_setLogicalBasis(void* highs) {
  return (int)((Highs*)highs)->setBasis();
}

double Highs_getHighsRunTime(void* highs) {
  return (double)((Highs*)highs)->getHighsRunTime();
}

int Highs_addRow(void* highs, const double lower, const double upper,
                 const int num_new_nz, const int* indices,
                 const double* values) {
  return ((Highs*)highs)->addRow(lower, upper, num_new_nz, indices, values);
}

int Highs_addRows(void* highs, const int num_new_row, const double* lower,
                  const double* upper, const int num_new_nz, const int* starts,
                  const int* indices, const double* values) {
  return ((Highs*)highs)
      ->addRows(num_new_row, lower, upper, num_new_nz, starts, indices, values);
}

int Highs_addCol(void* highs, const double cost, const double lower,
                 const double upper, const int num_new_nz, const int* indices,
                 const double* values) {
  return ((Highs*)highs)
      ->addCol(cost, lower, upper, num_new_nz, indices, values);
}

int Highs_addCols(void* highs, const int num_new_col, const double* costs,
                  const double* lower, const double* upper,
                  const int num_new_nz, const int* starts, const int* indices,
                  const double* values) {
  return ((Highs*)highs)
      ->addCols(num_new_col, costs, lower, upper, num_new_nz, starts, indices,
                values);
}

int Highs_changeObjectiveSense(void* highs, const int sense) {
  ObjSense pass_sense = ObjSense::MINIMIZE;
  if (sense == (int)ObjSense::MAXIMIZE) pass_sense = ObjSense::MAXIMIZE;
  return ((Highs*)highs)->changeObjectiveSense(pass_sense);
}

int Highs_changeColCost(void* highs, const int col, const double cost) {
  return ((Highs*)highs)->changeColCost(col, cost);
}

int Highs_changeColsCostBySet(void* highs, const int num_set_entries,
                              const int* set, const double* cost) {
  return ((Highs*)highs)->changeColsCost(num_set_entries, set, cost);
}

int Highs_changeColsCostByMask(void* highs, const int* mask,
                               const double* cost) {
  return ((Highs*)highs)->changeColsCost(mask, cost);
}

int Highs_changeColBounds(void* highs, const int col, const double lower,
                          const double upper) {
  return ((Highs*)highs)->changeColBounds(col, lower, upper);
}

int Highs_changeColsBoundsByRange(void* highs, const int from_col,
                                  const int to_col, const double* lower,
                                  const double* upper) {
  return ((Highs*)highs)->changeColsBounds(from_col, to_col, lower, upper);
}

int Highs_changeColsBoundsBySet(void* highs, const int num_set_entries,
                                const int* set, const double* lower,
                                const double* upper) {
  return ((Highs*)highs)->changeColsBounds(num_set_entries, set, lower, upper);
}

int Highs_changeColsBoundsByMask(void* highs, const int* mask,
                                 const double* lower, const double* upper) {
  return ((Highs*)highs)->changeColsBounds(mask, lower, upper);
}

int Highs_changeRowBounds(void* highs, const int row, const double lower,
                          const double upper) {
  return ((Highs*)highs)->changeRowBounds(row, lower, upper);
}

int Highs_changeRowsBoundsBySet(void* highs, const int num_set_entries,
                                const int* set, const double* lower,
                                const double* upper) {
  return ((Highs*)highs)->changeRowsBounds(num_set_entries, set, lower, upper);
}

int Highs_changeRowsBoundsByMask(void* highs, const int* mask,
                                 const double* lower, const double* upper) {
  return ((Highs*)highs)->changeRowsBounds(mask, lower, upper);
}

int Highs_getObjectiveSense(void* highs, int* sense) {
  ObjSense get_sense;
  int status = ((Highs*)highs)->getObjectiveSense(get_sense);
  *sense = (int)get_sense;
  return status;
}

int Highs_getColsByRange(void* highs, const int from_col, const int to_col,
                         int* num_col, double* costs, double* lower,
                         double* upper, int* num_nz, int* matrix_start,
                         int* matrix_index, double* matrix_value) {
  int numcol, numnz;
  int status = ((Highs*)highs)
                   ->getCols(from_col, to_col, numcol, costs, lower, upper,
                             numnz, matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return status;
}

int Highs_getColsBySet(void* highs, const int num_set_entries, const int* set,
                       int* num_col, double* costs, double* lower,
                       double* upper, int* num_nz, int* matrix_start,
                       int* matrix_index, double* matrix_value) {
  int numcol, numnz;
  int status = ((Highs*)highs)
                   ->getCols(num_set_entries, set, numcol, costs, lower, upper,
                             numnz, matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return status;
}

int Highs_getColsByMask(void* highs, const int* mask, int* num_col,
                        double* costs, double* lower, double* upper,
                        int* num_nz, int* matrix_start, int* matrix_index,
                        double* matrix_value) {
  int numcol, numnz;
  int status = ((Highs*)highs)
                   ->getCols(mask, numcol, costs, lower, upper, numnz,
                             matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return status;
}

int Highs_getRowsByRange(void* highs, const int from_row, const int to_row,
                         int* num_row, double* lower, double* upper,
                         int* num_nz, int* matrix_start, int* matrix_index,
                         double* matrix_value) {
  int numrow, numnz;
  int status = ((Highs*)highs)
                   ->getRows(from_row, to_row, numrow, lower, upper, numnz,
                             matrix_start, matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return status;
}

int Highs_getRowsBySet(void* highs, const int num_set_entries, const int* set,
                       int* num_row, double* lower, double* upper, int* num_nz,
                       int* matrix_start, int* matrix_index,
                       double* matrix_value) {
  int numrow, numnz;
  int status = ((Highs*)highs)
                   ->getRows(num_set_entries, set, numrow, lower, upper, numnz,
                             matrix_start, matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return status;
}

int Highs_getRowsByMask(void* highs, const int* mask, int* num_row,
                        double* lower, double* upper, int* num_nz,
                        int* matrix_start, int* matrix_index,
                        double* matrix_value) {
  int numrow, numnz;
  int status = ((Highs*)highs)
                   ->getRows(mask, numrow, lower, upper, numnz, matrix_start,
                             matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return status;
}

int Highs_deleteColsByRange(void* highs, const int from_col, const int to_col) {
  return ((Highs*)highs)->deleteCols(from_col, to_col);
}

int Highs_deleteColsBySet(void* highs, const int num_set_entries,
                          const int* set) {
  return ((Highs*)highs)->deleteCols(num_set_entries, set);
}

int Highs_deleteColsByMask(void* highs, int* mask) {
  return ((Highs*)highs)->deleteCols(mask);
}

int Highs_deleteRowsByRange(void* highs, const int from_row, const int to_row) {
  return ((Highs*)highs)->deleteRows(from_row, to_row);
}

int Highs_deleteRowsBySet(void* highs, const int num_set_entries,
                          const int* set) {
  return ((Highs*)highs)->deleteRows(num_set_entries, set);
}

int Highs_deleteRowsByMask(void* highs, int* mask) {
  return ((Highs*)highs)->deleteRows(mask);
}

double Highs_getHighsInfinity(void* highs) {
  return ((Highs*)highs)->getHighsInfinity();
}

int Highs_getNumCols(void* highs) { return ((Highs*)highs)->getLp().numCol_; }

int Highs_getNumRows(void* highs) { return ((Highs*)highs)->getLp().numRow_; }

int Highs_getNumNz(void* highs) {
  int numCol = Highs_getNumCols(highs);
  if (numCol <= 0) return 0;
  return ((Highs*)highs)->getLp().Astart_[numCol];
}

const char* Highs_highsModelStatusToChar(void* highs,
                                         int int_highs_model_status) {
  const char* illegal_highs_model_status = "Model status out of range";
  if (int_highs_model_status < (int)HighsModelStatus::HIGHS_MODEL_STATUS_MIN ||
      int_highs_model_status > (int)HighsModelStatus::HIGHS_MODEL_STATUS_MAX)
    return illegal_highs_model_status;
  return ((Highs*)highs)
      ->highsModelStatusToString(
          static_cast<HighsModelStatus>(int_highs_model_status))
      .c_str();
}

const char* Highs_primalDualStatusToChar(void* highs,
                                         int int_primal_dual_status) {
  const char* illegal_primal_dual_status = "Primal/Dual status out of range";
  if (int_primal_dual_status < PrimalDualStatus::STATUS_MIN ||
      int_primal_dual_status > PrimalDualStatus::STATUS_MAX)
    return illegal_primal_dual_status;
  return ((Highs*)highs)
      ->primalDualStatusToString(int_primal_dual_status)
      .c_str();
}
