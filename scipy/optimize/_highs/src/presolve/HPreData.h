/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPreData.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_HPREDATA_H_
#define PRESOLVE_HPREDATA_H_

#include <cstring>
#include <list>
#include <stack>
#include <utility>
#include <vector>

#include "lp_data/HConst.h"
#include "test/KktCh2.h"

using std::pair;
using std::stack;
using std::string;
using std::vector;

namespace presolve {
struct change {
  int type;
  int row;
  int col;
};

class HPreData {
 public:
  HPreData();
  virtual ~HPreData() {}

  // Model data
  int numCol;
  int numRow;
  int numRowOriginal;
  int numColOriginal;
  int numTot;

  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<double> rowLower;
  vector<double> rowUpper;

  // during postsolve hold the reduced solution, then at the end of postsolve
  // they hold the recovered. passed to dev kkt checker.
  vector<double> colValue;
  vector<double> colDual;
  vector<double> rowValue;
  vector<double> rowDual;

  // Row wise copy of matrix.
  vector<int> ARstart;
  vector<int> ARindex;
  vector<double> ARvalue;

  vector<int> Aend;

  // Solution
  // The first numColOriginal elements are the primal variables, slacks after
  vector<double> valuePrimal;
  vector<double> valueColDual;
  vector<double> valueRowDual;

  vector<int> nzCol;  // nonzeros in columns and rows
  vector<int> nzRow;
  vector<int> flagCol;
  vector<int> flagRow;

  const bool use_simplex_basis_logic = false;  // true;//
  vector<int> nonbasicFlag;

  // Record of whether a column or row is basic or nonbasic
  vector<HighsBasisStatus> col_status;
  vector<HighsBasisStatus> row_status;

  vector<double> colCostAtEl;
  vector<double> rowLowerAtEl;
  vector<double> rowUpperAtEl;

  void makeARCopy();
  void makeACopy();
  double getaij(int i, int j);
  bool isZeroA(int i, int j);
  double getRowValue(int i);

  stack<double> postValue;

  // to match reduced solution to original
  vector<int> rIndex;
  vector<int> cIndex;

  dev_kkt_check::KktChStep chk2;

  stack<change> chng;
  stack<pair<int, vector<double>>> oldBounds;  //(j, l, u)
};

struct MainLoop {
  int rows;
  int cols;
  int nnz;
};

struct DevStats {
  int n_loops = 0;
  std::vector<MainLoop> loops;
};

struct PresolveStats {
  DevStats dev;

  int n_rows_removed = 0;
  int n_cols_removed = 0;
  int n_nnz_removed = 0;
};

void initPresolve(PresolveStats& stats);

}  // namespace presolve

#endif /* PRESOLVE_HPREDATA_H_ */
