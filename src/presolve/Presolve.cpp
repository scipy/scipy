/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/Presolve.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "presolve/Presolve.h"

#include "io/HighsIO.h"
#include "lp_data/HConst.h"

//#include "simplex/HFactor.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <queue>
#include <sstream>

#include "test/KktChStep.h"

using std::cout;
using std::endl;
using std::flush;
using std::get;
using std::ios;
using std::list;
using std::make_pair;
using std::max;
using std::min;
using std::ofstream;
using std::setprecision;
using std::setw;
using std::stringstream;

void Presolve::load(const HighsLp& lp) {
  timer.recordStart(MATRIX_COPY);
  numCol = lp.numCol_;
  numRow = lp.numRow_;
  numTot = numTot;
  Astart = lp.Astart_;
  Aindex = lp.Aindex_;
  Avalue = lp.Avalue_;

  colCost = lp.colCost_;
  if (lp.sense_ == ObjSense::MAXIMIZE) {
    for (unsigned int col = 0; col < lp.colCost_.size(); col++)
      colCost[col] = -colCost[col];
  }

  colLower = lp.colLower_;
  colUpper = lp.colUpper_;
  rowLower = lp.rowLower_;
  rowUpper = lp.rowUpper_;

  modelName = lp.model_name_;
  timer.recordFinish(MATRIX_COPY);
}

void PresolveInfo::negateColDuals(bool reduced) {
  if (reduced)
    for (unsigned int col = 0; col < reduced_solution_.col_dual.size(); col++)
      reduced_solution_.col_dual[col] = reduced_solution_.col_dual[col];
  else
    for (unsigned int col = 0; col < recovered_solution_.col_dual.size(); col++)
      recovered_solution_.col_dual[col] = recovered_solution_.col_dual[col];
  return;
}

void PresolveInfo::negateReducedCosts() {
  for (unsigned int col = 0; col < reduced_lp_.colCost_.size(); col++)
    reduced_lp_.colCost_[col] = -reduced_lp_.colCost_[col];
}

HighsLp& PresolveInfo::getReducedProblem() {
  if (presolve_.size() == 0) {
    std::cout << "Error during presolve. No presolve initialized." << std::endl;
  } else if (presolve_[0].status != presolve_[0].stat::Reduced) {
    std::cout << "Error during presolve. No reduced LP. status: "
              << presolve_[0].status << std::endl;
  } else {
    if (presolve_[0].numRow == 0 && presolve_[0].numCol == 0) {
      // Reduced problem has been already moved to this.reduced_lp_;
      return reduced_lp_;
    } else {
      // Move vectors so no copying happens. presolve does not need that lp
      // any more.
      reduced_lp_.numCol_ = presolve_[0].numCol;
      reduced_lp_.numRow_ = presolve_[0].numRow;
      reduced_lp_.Astart_ = std::move(presolve_[0].Astart);
      reduced_lp_.Aindex_ = std::move(presolve_[0].Aindex);
      reduced_lp_.Avalue_ = std::move(presolve_[0].Avalue);
      reduced_lp_.colCost_ = std::move(presolve_[0].colCost);
      reduced_lp_.colLower_ = std::move(presolve_[0].colLower);
      reduced_lp_.colUpper_ = std::move(presolve_[0].colUpper);
      reduced_lp_.rowLower_ = std::move(presolve_[0].rowLower);
      reduced_lp_.rowUpper_ = std::move(presolve_[0].rowUpper);

      reduced_lp_.sense_ = ObjSense::MINIMIZE;
      reduced_lp_.offset_ = 0;
      reduced_lp_.model_name_ =
          std::move(presolve_[0].modelName);  //"Presolved model";

      reduced_lp_.nnz_ = reduced_lp_.Avalue_.size();
    }
  }
  return reduced_lp_;
}

void Presolve::setBasisInfo(
    const std::vector<HighsBasisStatus>& pass_col_status,
    const std::vector<HighsBasisStatus>& pass_row_status) {
  col_status = pass_col_status;
  row_status = pass_row_status;
}

int Presolve::presolve(int print) {
  iPrint = print;

  // iPrint = 1;
  // iKKTcheck = 1;
  // chk.print = 1;

  if (iPrint > 0) {
    cout << "Presolve started ..." << endl;
    cout << "Original problem ... N=" << numCol << "  M=" << numRow << endl;
  }

  initializeVectors();
  if (status) return status;

  int iter = 1;
  // print(0);

  timer.recordStart(FIXED_COL);
  for (int j = 0; j < numCol; ++j)
    if (flagCol.at(j)) {
      removeIfFixed(j);
      if (status) return status;
    }
  timer.recordFinish(FIXED_COL);

  while (hasChange == 1) {
    hasChange = false;
    if (iPrint > 0) cout << "PR: main loop " << iter << ":" << endl;
    //***************** main loop ******************
    checkBoundsAreConsistent();
    if (status) return status;

    removeRowSingletons();
    if (status) return status;
    removeForcingConstraints(iter);
    if (status) return status;

    removeRowSingletons();
    if (status) return status;
    removeDoubletonEquations();
    if (status) return status;

    removeRowSingletons();
    if (status) return status;
    removeColumnSingletons();
    if (status) return status;

    removeDominatedColumns();
    if (status) return status;

    //***************** main loop ******************
    iter++;
  }

  timer.recordStart(RESIZE_MATRIX);
  checkForChanges(iter);
  timer.recordFinish(RESIZE_MATRIX);

  timer.updateInfo();

  return status;
}

HighsPresolveStatus Presolve::presolve() {
  timer.recordStart(TOTAL_PRESOLVE_TIME);
  HighsPresolveStatus presolve_status = HighsPresolveStatus::NotReduced;
  int result = presolve(0);
  switch (result) {
    case stat::Unbounded:
      presolve_status = HighsPresolveStatus::Unbounded;
      break;
    case stat::Infeasible:
      presolve_status = HighsPresolveStatus::Infeasible;
      break;
    case stat::Reduced:
      if (numCol > 0 || numRow > 0)
        presolve_status = HighsPresolveStatus::Reduced;
      else
        presolve_status = HighsPresolveStatus::ReducedToEmpty;
      break;
    case stat::Empty:
      presolve_status = HighsPresolveStatus::Empty;
      break;
    case stat::Optimal:
      // reduced problem solution indicated as optimal by
      // the solver.
      break;
  }
  timer.recordFinish(TOTAL_PRESOLVE_TIME);

  return presolve_status;
}

void Presolve::checkBoundsAreConsistent() {
  for (int col = 0; col < numCol; col++) {
    if (flagCol[col]) {
      if (colUpper[col] - colLower[col] < -tol) {
        status = Infeasible;
        return;
      }
    }
  }

  for (int row = 0; row < numRow; row++) {
    if (flagRow[row]) {
      if (rowUpper[row] - rowLower[row] < -tol) {
        status = Infeasible;
        return;
      }
    }
  }
}

/**
 * returns <x, y>
 * 		   <x, -1> if we need to skip row
 *
 * 		   row is of form akx_x + aky_y = b,
 */
pair<int, int> Presolve::getXYDoubletonEquations(const int row) {
  pair<int, int> colIndex;
  // row is of form akx_x + aky_y = b, where k=row and y is present in fewer
  // constraints

  int col1 = -1;
  int col2 = -1;
  int kk = ARstart.at(row);
  while (kk < ARstart.at(row + 1)) {
    if (flagCol.at(ARindex.at(kk))) {
      if (col1 == -1)
        col1 = ARindex.at(kk);
      else if (col2 == -1)
        col2 = ARindex.at(kk);
      else {
        cout << "ERROR: doubleton eq row" << row
             << " has more than two variables. \n";
        col2 = -2;
        break;
      }
      ++kk;
    } else
      ++kk;
  }
  if (col2 == -1)
    cout << "ERROR: doubleton eq row" << row
         << " has less than two variables. \n";
  if (col2 < 0) {
    colIndex.second = -1;
    return colIndex;
  }

  int x, y;
  if (nzCol.at(col1) <= nzCol.at(col2)) {
    y = col1;
    x = col2;
  } else {
    x = col1;
    y = col2;
  }

  //	if (nzCol.at(y) == 1 && nzCol.at(x) == 1) { //two singletons case
  // handled elsewhere 		colIndex.second = -1; 		return colIndex;
  //	}

  colIndex.first = x;
  colIndex.second = y;
  return colIndex;
}

void Presolve::processRowDoubletonEquation(const int row, const int x,
                                           const int y, const double akx,
                                           const double aky, const double b) {
  postValue.push(akx);
  postValue.push(aky);
  postValue.push(b);

  // modify bounds on variable x (j), variable y (col,k) is substituted out
  // double aik = Avalue.at(k);
  // double aij = Avalue.at(kk);
  pair<double, double> p = getNewBoundsDoubletonConstraint(row, y, x, aky, akx);
  double low = p.first;
  double upp = p.second;

  // add old bounds of x to checker and for postsolve
  if (iKKTcheck == 1) {
    vector<pair<int, double>> bndsL, bndsU, costS;
    bndsL.push_back(make_pair(x, colLower.at(x)));
    bndsU.push_back(make_pair(x, colUpper.at(x)));
    costS.push_back(make_pair(x, colCost.at(x)));
    chk.cLowers.push(bndsL);
    chk.cUppers.push(bndsU);
    chk.costs.push(costS);
  }

  vector<double> bnds({colLower.at(y), colUpper.at(y), colCost.at(y)});
  vector<double> bnds2({colLower.at(x), colUpper.at(x), colCost.at(x)});
  oldBounds.push(make_pair(y, bnds));
  oldBounds.push(make_pair(x, bnds2));

  if (low > colLower.at(x)) colLower.at(x) = low;
  if (upp < colUpper.at(x)) colUpper.at(x) = upp;

  // modify cost of xj
  colCost.at(x) = colCost.at(x) - colCost.at(y) * akx / aky;

  // for postsolve: need the new bounds too
  vector<double> bnds3({colLower.at(x), colUpper.at(x), colCost.at(x)});
  oldBounds.push(make_pair(x, bnds3));

  addChange(DOUBLETON_EQUATION, row, y);

  // remove y (col) and the row
  if (iPrint > 0)
    cout << "PR: Doubleton equation removed. Row " << row << ", column " << y
         << ", column left is " << x << "    nzy=" << nzCol.at(y) << endl;

  flagRow.at(row) = 0;
  nzCol.at(x)--;

  countRemovedRows(DOUBLETON_EQUATION);
  countRemovedCols(DOUBLETON_EQUATION);

  //----------------------------
  flagCol.at(y) = 0;
  if (!hasChange) hasChange = true;
}

void Presolve::removeDoubletonEquations() {
  timer.recordStart(DOUBLETON_EQUATION);
  // flagCol should have one more element at end which is zero
  // needed for AR matrix manipulation
  if ((int)flagCol.size() == numCol) flagCol.push_back(0);

  double b, akx, aky;
  int x, y;
  int iter = 0;

  for (int row = 0; row < numRow; row++)
    if (flagRow.at(row))
      if (nzRow.at(row) == 2 && rowLower[row] > -HIGHS_CONST_INF &&
          rowUpper[row] < HIGHS_CONST_INF &&
          fabs(rowLower[row] - rowUpper[row]) < tol) {
        // row is of form akx_x + aky_y = b, where k=row and y is present in
        // fewer constraints
        b = rowLower.at(row);
        pair<int, int> colIndex = getXYDoubletonEquations(row);
        x = colIndex.first;
        y = colIndex.second;

        // two singletons case handled elsewhere
        if (y < 0 || ((nzCol.at(y) == 1 && nzCol.at(x) == 1))) continue;

        akx = getaij(row, x);
        aky = getaij(row, y);
        processRowDoubletonEquation(row, x, y, akx, aky, b);

        for (int k = Astart.at(y); k < Aend.at(y); ++k)
          if (flagRow.at(Aindex.at(k)) && Aindex.at(k) != row) {
            int i = Aindex.at(k);
            double aiy = Avalue.at(k);

            // update row bounds
            if (iKKTcheck == 1) {
              vector<pair<int, double>> bndsL, bndsU;
              bndsL.push_back(make_pair(i, rowLower.at(i)));
              bndsU.push_back(make_pair(i, rowUpper.at(i)));
              chk.rLowers.push(bndsL);
              chk.rUppers.push(bndsU);
              addChange(DOUBLETON_EQUATION_ROW_BOUNDS_UPDATE, i, y);
            }

            if (rowLower.at(i) > -HIGHS_CONST_INF)
              rowLower.at(i) -= b * aiy / aky;
            if (rowUpper.at(i) < HIGHS_CONST_INF)
              rowUpper.at(i) -= b * aiy / aky;

            if (implRowValueLower.at(i) > -HIGHS_CONST_INF)
              implRowValueLower.at(i) -= b * aiy / aky;
            if (implRowValueUpper.at(i) < HIGHS_CONST_INF)
              implRowValueUpper.at(i) -= b * aiy / aky;

            // update matrix coefficients
            if (isZeroA(i, x))
              UpdateMatrixCoeffDoubletonEquationXzero(i, x, y, aiy, akx, aky);
            else
              UpdateMatrixCoeffDoubletonEquationXnonZero(i, x, y, aiy, akx,
                                                         aky);
          }
        if (Avalue.size() > 40000000) {
          trimA();
        }

        iter++;
      }
  timer.recordFinish(DOUBLETON_EQUATION);
}

void Presolve::UpdateMatrixCoeffDoubletonEquationXzero(const int i, const int x,
                                                       const int y,
                                                       const double aiy,
                                                       const double akx,
                                                       const double aky) {
  // case x is zero initially
  // row nonzero count doesn't change here
  // cout<<"case: x not present "<<i<<" "<<endl;

  // update AR
  int ind;
  for (ind = ARstart.at(i); ind < ARstart.at(i + 1); ++ind)
    if (ARindex.at(ind) == y) {
      break;
    }

  postValue.push(ARvalue.at(ind));
  postValue.push(y);
  addChange(DOUBLETON_EQUATION_X_ZERO_INITIALLY, i, x);

  ARindex.at(ind) = x;
  ARvalue.at(ind) = -aiy * akx / aky;

  // just row rep in checker
  if (iKKTcheck == 1) {
    chk.ARvalue.at(ind) = ARvalue.at(ind);
    chk.ARindex.at(ind) = ARindex.at(ind);
  }

  // update A: append X column to end of array
  int st = Avalue.size();
  for (int ind = Astart.at(x); ind < Aend.at(x); ++ind) {
    Avalue.push_back(Avalue.at(ind));
    Aindex.push_back(Aindex.at(ind));
  }
  Avalue.push_back(-aiy * akx / aky);
  Aindex.push_back(i);
  Astart.at(x) = st;
  Aend.at(x) = Avalue.size();

  nzCol.at(x)++;
  // nzRow does not change here.
  if (nzCol.at(x) == 2) singCol.remove(x);
}

void Presolve::UpdateMatrixCoeffDoubletonEquationXnonZero(
    const int i, const int x, const int y, const double aiy, const double akx,
    const double aky) {
  int ind;

  // update nonzeros: for removal of
  nzRow.at(i)--;
  if (nzRow.at(i) == 1) singRow.push_back(i);

  if (nzRow.at(i) == 0) {
    singRow.remove(i);
    removeEmptyRow(i);
    countRemovedRows(DOUBLETON_EQUATION);
  }

  double xNew;
  for (ind = ARstart.at(i); ind < ARstart.at(i + 1); ++ind)
    if (ARindex.at(ind) == x) break;

  xNew = ARvalue.at(ind) - (aiy * akx) / aky;
  if (fabs(xNew) > tol) {
    // case new x != 0
    // cout<<"case: x still there row "<<i<<" "<<endl;

    postValue.push(ARvalue.at(ind));
    addChange(DOUBLETON_EQUATION_NEW_X_NONZERO, i, x);
    ARvalue.at(ind) = xNew;

    if (iKKTcheck == 1) chk.ARvalue.at(ind) = xNew;

    // update A:
    for (ind = Astart.at(x); ind < Aend.at(x); ++ind)
      if (Aindex.at(ind) == i) {
        break;
      }
    Avalue.at(ind) = xNew;
  } else if (xNew < tol) {
    // case new x == 0
    // cout<<"case: x also disappears from row "<<i<<" "<<endl;
    // update nz row
    nzRow.at(i)--;
    // update singleton row list
    if (nzRow.at(i) == 1) singRow.push_back(i);

    if (nzRow.at(i) == 0) {
      singRow.remove(i);
      removeEmptyRow(i);
      countRemovedRows(DOUBLETON_EQUATION);
    }

    if (nzRow.at(i) > 0) {
      // AR update
      // set ARindex of element for x to numCol
      // flagCol[numCol] = false
      // mind when resizing: should be OK
      postValue.push(ARvalue.at(ind));

      ARindex.at(ind) = numCol;
      if (iKKTcheck == 1) {
        chk.ARindex.at(ind) = ARindex.at(ind);
        chk.ARvalue.at(ind) = ARvalue.at(ind);
      }

      addChange(DOUBLETON_EQUATION_NEW_X_ZERO_AR_UPDATE, i, x);
    }

    if (nzCol.at(x) > 0) {
      // A update for case when x is zero: move x entry to end and set
      // Aend to be Aend - 1;
      int indi;
      for (indi = Astart.at(x); indi < Aend.at(x); ++indi)
        if (Aindex.at(indi) == i) break;

      postValue.push(Avalue.at(indi));

      // if indi is not Aend-1 swap elements indi and Aend-1
      if (indi != Aend.at(x) - 1) {
        double tmp = Avalue.at(Aend.at(x) - 1);
        int tmpi = Aindex.at(Aend.at(x) - 1);
        Avalue.at(Aend.at(x) - 1) = Avalue.at(indi);
        Aindex.at(Aend.at(x) - 1) = Aindex.at(indi);
        Avalue.at(indi) = tmp;
        Aindex.at(indi) = tmpi;
      }
      Aend.at(x)--;
      addChange(DOUBLETON_EQUATION_NEW_X_ZERO_A_UPDATE, i, x);
    }

    // update nz col
    nzCol.at(x)--;
    // update singleton col list
    if (nzCol.at(x) == 1) singCol.push_back(x);
    if (nzCol.at(x) == 0) {
      removeEmptyColumn(x);
    }
  }
  if (y) {
  }  // surpress warning.
}

void Presolve::trimA() {
  int cntEl = 0;
  for (int j = 0; j < numCol; ++j)
    if (flagCol.at(j)) cntEl += nzCol.at(j);

  vector<pair<int, size_t>> vp;
  vp.reserve(numCol);

  for (int i = 0; i != numCol; ++i) {
    vp.push_back(make_pair(Astart.at(i), i));
  }

  // Sorting will put lower values ahead of larger ones,
  // resolving ties using the original index
  sort(vp.begin(), vp.end());

  vector<int> Aendtmp;
  Aendtmp = Aend;

  int iPut = 0;
  for (size_t i = 0; i != vp.size(); ++i) {
    int col = vp.at(i).second;
    if (flagCol.at(col)) {
      int k = vp.at(i).first;
      Astart.at(col) = iPut;
      while (k < Aendtmp.at(col)) {
        if (flagRow.at(Aindex.at(k))) {
          Avalue[iPut] = Avalue.at(k);
          Aindex[iPut] = Aindex.at(k);
          iPut++;
        }
        k++;
      }
      Aend.at(col) = iPut;
    }
  }
  Avalue.resize(iPut);
  Aindex.resize(iPut);
}

void Presolve::resizeProblem() {
  int i, j, k;

  int nz = 0;
  int nR = 0;
  int nC = 0;

  // arrays to keep track of indices
  rIndex.assign(numRow, -1);
  cIndex.assign(numCol, -1);

  for (i = 0; i < numRow; ++i)
    if (flagRow.at(i)) {
      nz += nzRow.at(i);
      rIndex.at(i) = nR;
      nR++;
    }

  for (i = 0; i < numCol; ++i)
    if (flagCol.at(i)) {
      cIndex.at(i) = nC;
      nC++;
    }

  // counts
  numRowOriginal = numRow;
  numColOriginal = numCol;
  numRow = nR;
  numCol = nC;
  numTot = nR + nC;

  if (nR + nC == 0) {
    status = Empty;
    return;
  }

  // matrix
  vector<int> iwork(numCol, 0);
  Astart.assign(numCol + 1, 0);
  Aend.assign(numCol + 1, 0);
  Aindex.resize(nz);
  Avalue.resize(nz);

  for (i = 0; i < numRowOriginal; ++i)
    if (flagRow.at(i))
      for (int k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
        j = ARindex.at(k);
        if (flagCol.at(j)) iwork.at(cIndex.at(j))++;
      }
  for (i = 1; i <= numCol; ++i)
    Astart.at(i) = Astart.at(i - 1) + iwork.at(i - 1);
  for (i = 0; i < numCol; ++i) iwork.at(i) = Aend.at(i) = Astart.at(i);
  for (i = 0; i < numRowOriginal; ++i) {
    if (flagRow.at(i)) {
      int iRow = rIndex.at(i);
      for (k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
        j = ARindex.at(k);
        if (flagCol.at(j)) {
          int iCol = cIndex.at(j);
          int iPut = iwork.at(iCol)++;
          Aindex.at(iPut) = iRow;
          Avalue.at(iPut) = ARvalue.at(k);
        }
      }
    }
  }
  // For KKT checker: pass vectors before you trim them
  if (iKKTcheck == 1) {
    chk.setFlags(flagRow, flagCol);
    chk.setBoundsCostRHS(colUpper, colLower, colCost, rowLower, rowUpper);
  }

  // also call before trimming
  resizeImpliedBounds();

  // cost, bounds
  colCostAtEl = colCost;
  vector<double> tempCost = colCost;
  vector<double> temp = colLower;
  vector<double> teup = colUpper;

  colCost.resize(numCol);
  colLower.resize(numCol);
  colUpper.resize(numCol);

  k = 0;
  for (i = 0; i < numColOriginal; ++i)
    if (flagCol.at(i)) {
      colCost.at(k) = tempCost.at(i);
      colLower.at(k) = temp.at(i);
      colUpper.at(k) = teup.at(i);
      k++;
    }

  // RHS and bounds
  rowLowerAtEl = rowLower;
  rowUpperAtEl = rowUpper;
  temp = rowLower;
  teup = rowUpper;
  rowLower.resize(numRow);
  rowUpper.resize(numRow);
  k = 0;
  for (i = 0; i < numRowOriginal; ++i)
    if (flagRow.at(i)) {
      rowLower.at(k) = temp.at(i);
      rowUpper.at(k) = teup.at(i);
      k++;
    }

  if (chk.print == 3) {
    ofstream myfile;
    myfile.open("../experiments/out", ios::app);
    myfile << " eliminated rows " << (numRowOriginal - numRow) << " cols "
           << (numColOriginal - numCol);
    myfile.close();

    myfile.open("../experiments/t3", ios::app);
    myfile << (numRowOriginal) << "  &  " << (numColOriginal) << "  & ";
    myfile << (numRowOriginal - numRow) << "  &  " << (numColOriginal - numCol)
           << "  & " << endl;

    myfile.close();
  }
}

void Presolve::initializeVectors() {
  // copy original bounds
  colCostOriginal = colCost;
  rowUpperOriginal = rowUpper;
  rowLowerOriginal = rowLower;
  colUpperOriginal = colUpper;
  colLowerOriginal = colLower;

  makeARCopy();

  valueRowDual.resize(numRow);
  valuePrimal.resize(numCol);
  valueColDual.resize(numCol);

  flagCol.assign(numCol, 1);
  flagRow.assign(numRow, 1);

  if (iKKTcheck) setKKTcheckerData();

  nzCol.assign(numCol, 0);
  nzRow.assign(numRow, 0);

  for (int i = 0; i < numRow; ++i) {
    nzRow.at(i) = ARstart.at(i + 1) - ARstart.at(i);
    if (nzRow.at(i) == 1) singRow.push_back(i);
    if (nzRow.at(i) == 0) {
      timer.recordStart(EMPTY_ROW);
      removeEmptyRow(i);
      countRemovedRows(EMPTY_ROW);
      timer.recordFinish(EMPTY_ROW);
    }
  }

  Aend.resize(numCol + 1);
  for (int i = 0; i < numCol; ++i) {
    Aend.at(i) = Astart.at(i + 1);
    nzCol.at(i) = Aend.at(i) - Astart.at(i);
    if (nzCol.at(i) == 1) singCol.push_back(i);
  }
  objShift = 0;

  implColUpper = colUpper;  // working copies of primal variable bounds
  implColLower = colLower;
  implColLowerRowIndex.assign(numCol, -1);
  implColUpperRowIndex.assign(numCol, -1);

  implRowDualLowerSingColRowIndex.assign(numRow, -1);
  implRowDualUpperSingColRowIndex.assign(numRow, -1);
  implRowDualLower.assign(numRow, -HIGHS_CONST_INF);
  implRowDualUpper.assign(numRow, HIGHS_CONST_INF);

  implColDualLower.assign(numCol, -HIGHS_CONST_INF);
  implColDualUpper.assign(numCol, HIGHS_CONST_INF);
  implRowValueLower = rowLower;
  implRowValueUpper = rowUpper;

  for (int i = 0; i < numRow; ++i) {
    if (rowLower.at(i) == -HIGHS_CONST_INF) implRowDualUpper.at(i) = 0;
    if (rowUpper.at(i) == HIGHS_CONST_INF) implRowDualLower.at(i) = 0;
  }

  for (int i = 0; i < numCol; ++i) {
    if (colLower.at(i) == -HIGHS_CONST_INF) implColDualUpper.at(i) = 0;
    if (colUpper.at(i) == HIGHS_CONST_INF) implColDualLower.at(i) = 0;
  }

  colCostAtEl = colCost;
  rowLowerAtEl = rowLower;
  rowUpperAtEl = rowUpper;
}

void Presolve::removeIfFixed(int j) {
  if (colLower.at(j) == colUpper.at(j)) {
    setPrimalValue(j, colUpper.at(j));
    addChange(FIXED_COL, 0, j);
    if (iPrint > 0)
      cout << "PR: Fixed variable " << j << " = " << colUpper.at(j)
           << ". Column eliminated." << endl;

    countRemovedCols(FIXED_COL);

    for (int k = Astart.at(j); k < Aend.at(j); ++k) {
      if (flagRow.at(Aindex.at(k))) {
        int i = Aindex.at(k);

        if (nzRow.at(i) == 0) {
          removeEmptyRow(i);
          countRemovedRows(FIXED_COL);
        }
      }
    }
  }
}

void Presolve::removeEmptyRow(int i) {
  if (rowLower.at(i) <= tol && rowUpper.at(i) >= -tol) {
    if (iPrint > 0) cout << "PR: Empty row " << i << " removed. " << endl;
    flagRow.at(i) = 0;
    valueRowDual.at(i) = 0;
    addChange(EMPTY_ROW, i, 0);
  } else {
    if (iPrint > 0) cout << "PR: Problem infeasible." << endl;
    status = Infeasible;
    return;
  }
}

void Presolve::removeEmptyColumn(int j) {
  flagCol.at(j) = 0;
  singCol.remove(j);
  double value;
  if ((colCost.at(j) < 0 && colUpper.at(j) == HIGHS_CONST_INF) ||
      (colCost.at(j) > 0 && colLower.at(j) == -HIGHS_CONST_INF)) {
    if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
    status = Unbounded;
    return;
  }

  if (colCost.at(j) > 0)
    value = colLower.at(j);
  else if (colCost.at(j) < 0)
    value = colUpper.at(j);
  else if (colUpper.at(j) >= 0 && colLower.at(j) <= 0)
    value = 0;
  else if (colUpper.at(j) < 0)
    value = colUpper.at(j);
  else
    value = colLower.at(j);

  setPrimalValue(j, value);
  valueColDual.at(j) = colCost.at(j);

  addChange(EMPTY_COL, 0, j);

  if (iPrint > 0)
    cout << "PR: Column: " << j
         << " eliminated: all nonzero rows have been removed. Cost = "
         << colCost.at(j) << ", value = " << value << endl;

  countRemovedCols(EMPTY_COL);
}

void Presolve::rowDualBoundsDominatedColumns() {
  int col, i, k;

  // for each row calc yihat and yibar and store in implRowDualLower and
  // implRowDualUpper
  for (list<int>::iterator it = singCol.begin(); it != singCol.end(); ++it)
    if (flagCol.at(*it)) {
      col = *it;
      k = getSingColElementIndexInA(col);
      i = Aindex.at(k);

      if (!flagRow.at(i)) {
        cout << "ERROR: column singleton " << col << " is in row " << i
             << " which is already mapped off\n";
        exit(-1);
      }

      if (colLower.at(col) == -HIGHS_CONST_INF ||
          colUpper.at(col) == HIGHS_CONST_INF) {
        if (colLower.at(col) > -HIGHS_CONST_INF &&
            colUpper.at(col) == HIGHS_CONST_INF) {
          if (Avalue.at(k) > 0)
            if ((colCost.at(col) / Avalue.at(k)) < implRowDualUpper.at(i))
              implRowDualUpper.at(i) = colCost.at(col) / Avalue.at(k);
          if (Avalue.at(k) < 0)
            if ((colCost.at(col) / Avalue.at(k)) > implRowDualLower.at(i))
              implRowDualLower.at(i) = colCost.at(col) / Avalue.at(k);
        } else if (colLower.at(col) == -HIGHS_CONST_INF &&
                   colUpper.at(col) < HIGHS_CONST_INF) {
          if (Avalue.at(k) > 0)
            if ((colCost.at(col) / Avalue.at(k)) > implRowDualLower.at(i))
              implRowDualUpper.at(i) = -colCost.at(col) / Avalue.at(k);
          if (Avalue.at(k) < 0)
            if ((colCost.at(col) / Avalue.at(k)) < implRowDualUpper.at(i))
              implRowDualUpper.at(i) = colCost.at(col) / Avalue.at(k);
        } else if (colLower.at(col) == -HIGHS_CONST_INF &&
                   colUpper.at(col) == HIGHS_CONST_INF) {
          // all should be removed earlier but use them
          if ((colCost.at(col) / Avalue.at(k)) > implRowDualLower.at(i))
            implRowDualLower.at(i) = colCost.at(col) / Avalue.at(k);
          if ((colCost.at(col) / Avalue.at(k)) < implRowDualUpper.at(i))
            implRowDualUpper.at(i) = colCost.at(col) / Avalue.at(k);
        }

        if (implRowDualLower.at(i) > implRowDualUpper.at(i)) {
          cout << "Error: inconstistent bounds for Lagrange multiplier for row "
               << i << " detected after column singleton " << col
               << ". In presolve::dominatedColumns" << endl;
          exit(0);
        }
      }
    }
}

pair<double, double> Presolve::getImpliedColumnBounds(int j) {
  pair<double, double> out;
  double e = 0;
  double d = 0;

  int i;
  for (int k = Astart.at(j); k < Aend.at(j); ++k) {
    i = Aindex.at(k);
    if (flagRow.at(i)) {
      if (Avalue.at(k) < 0) {
        if (implRowDualUpper.at(i) < HIGHS_CONST_INF)
          e += Avalue.at(k) * implRowDualUpper.at(i);
        else {
          e = -HIGHS_CONST_INF;
          break;
        }
      } else {
        if (implRowDualLower.at(i) > -HIGHS_CONST_INF)
          e += Avalue.at(k) * implRowDualLower.at(i);
        else {
          e = -HIGHS_CONST_INF;
          break;
        }
      }
    }
  }

  for (int k = Astart.at(j); k < Aend.at(j); ++k) {
    i = Aindex.at(k);
    if (flagRow.at(i)) {
      if (Avalue.at(k) < 0) {
        if (implRowDualLower.at(i) > -HIGHS_CONST_INF)
          d += Avalue.at(k) * implRowDualLower.at(i);
        else {
          d = HIGHS_CONST_INF;
          break;
        }
      } else {
        if (implRowDualUpper.at(i) < HIGHS_CONST_INF)
          d += Avalue.at(k) * implRowDualUpper.at(i);
        else {
          d = HIGHS_CONST_INF;
          break;
        }
      }
    }
  }

  if (e > d) {
    cout << "Error: inconstistent bounds for Lagrange multipliers for column "
         << j << ": e>d. In presolve::dominatedColumns" << endl;
    exit(-1);
  }
  out.first = d;
  out.second = e;
  return out;
}

void Presolve::removeDominatedColumns() {
  // for each column j calculate e and d and check:
  double e, d;
  pair<double, double> p;
  for (int j = 0; j < numCol; ++j)
    if (flagCol.at(j)) {
      timer.recordStart(DOMINATED_COLS);

      p = getImpliedColumnBounds(j);
      d = p.first;
      e = p.second;

      // check if it is dominated
      if (colCost.at(j) - d > tol) {
        if (colLower.at(j) == -HIGHS_CONST_INF) {
          if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
          status = Unbounded;
          return;
        }
        setPrimalValue(j, colLower.at(j));
        addChange(DOMINATED_COLS, 0, j);
        if (iPrint > 0)
          cout << "PR: Dominated column " << j
               << " removed. Value := " << valuePrimal.at(j) << endl;
        countRemovedCols(DOMINATED_COLS);
      } else if (colCost.at(j) - e < -tol) {
        if (colUpper.at(j) == HIGHS_CONST_INF) {
          if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
          status = Unbounded;
          return;
        }
        setPrimalValue(j, colUpper.at(j));
        addChange(DOMINATED_COLS, 0, j);
        if (iPrint > 0)
          cout << "PR: Dominated column " << j
               << " removed. Value := " << valuePrimal.at(j) << endl;
        countRemovedCols(DOMINATED_COLS);
      } else {
        // update implied bounds
        if (implColDualLower.at(j) < (colCost.at(j) - d))
          implColDualLower.at(j) = colCost.at(j) - d;
        if (implColDualUpper.at(j) > (colCost.at(j) - e))
          implColDualUpper.at(j) = colCost.at(j) - e;
        if (implColDualLower.at(j) > implColDualUpper.at(j))
          cout << "INCONSISTENT\n";

        timer.recordFinish(DOMINATED_COLS);

        removeIfWeaklyDominated(j, d, e);
        continue;
      }
      timer.recordFinish(DOMINATED_COLS);
    }
}

void Presolve::removeIfWeaklyDominated(const int j, const double d,
                                       const double e) {
  timer.recordStart(WEAKLY_DOMINATED_COLS);

  int i;
  // check if it is weakly dominated: Excluding singletons!
  if (nzCol.at(j) > 1) {
    if (d < HIGHS_CONST_INF && fabs(colCost.at(j) - d) < tol &&
        colLower.at(j) > -HIGHS_CONST_INF) {
      setPrimalValue(j, colLower.at(j));
      addChange(WEAKLY_DOMINATED_COLS, 0, j);
      if (iPrint > 0)
        cout << "PR: Weakly Dominated column " << j
             << " removed. Value := " << valuePrimal.at(j) << endl;

      countRemovedCols(WEAKLY_DOMINATED_COLS);
    } else if (e > -HIGHS_CONST_INF && fabs(colCost.at(j) - e) < tol &&
               colUpper.at(j) < HIGHS_CONST_INF) {
      setPrimalValue(j, colUpper.at(j));
      addChange(WEAKLY_DOMINATED_COLS, 0, j);
      if (iPrint > 0)
        cout << "PR: Weakly Dominated column " << j
             << " removed. Value := " << valuePrimal.at(j) << endl;

      countRemovedCols(WEAKLY_DOMINATED_COLS);
    } else {
      double bnd;

      // calculate new bounds
      if (colLower.at(j) > -HIGHS_CONST_INF ||
          colUpper.at(j) == HIGHS_CONST_INF)
        for (int kk = Astart.at(j); kk < Aend.at(j); ++kk)
          if (flagRow.at(Aindex.at(kk)) && d < HIGHS_CONST_INF) {
            i = Aindex.at(kk);
            if (Avalue.at(kk) > 0 &&
                implRowDualLower.at(i) > -HIGHS_CONST_INF) {
              bnd =
                  -(colCost.at(j) + d) / Avalue.at(kk) + implRowDualLower.at(i);
              if (bnd < implRowDualUpper.at(i) &&
                  !(bnd < implRowDualLower.at(i)))
                implRowDualUpper.at(i) = bnd;
            } else if (Avalue.at(kk) < 0 &&
                       implRowDualUpper.at(i) < HIGHS_CONST_INF) {
              bnd =
                  -(colCost.at(j) + d) / Avalue.at(kk) + implRowDualUpper.at(i);
              if (bnd > implRowDualLower.at(i) &&
                  !(bnd > implRowDualUpper.at(i)))
                implRowDualLower.at(i) = bnd;
            }
          }

      if (colLower.at(j) == -HIGHS_CONST_INF ||
          colUpper.at(j) < HIGHS_CONST_INF)
        for (int kk = Astart.at(j); kk < Aend.at(j); ++kk)
          if (flagRow.at(Aindex.at(kk)) && e > -HIGHS_CONST_INF) {
            i = Aindex.at(kk);
            if (Avalue.at(kk) > 0 && implRowDualUpper.at(i) < HIGHS_CONST_INF) {
              bnd =
                  -(colCost.at(j) + e) / Avalue.at(kk) + implRowDualUpper.at(i);
              if (bnd > implRowDualLower.at(i) &&
                  !(bnd > implRowDualUpper.at(i)))
                implRowDualLower.at(i) = bnd;
            } else if (Avalue.at(kk) < 0 &&
                       implRowDualLower.at(i) > -HIGHS_CONST_INF) {
              bnd =
                  -(colCost.at(j) + e) / Avalue.at(kk) + implRowDualLower.at(i);
              if (bnd < implRowDualUpper.at(i) &&
                  !(bnd < implRowDualLower.at(i)))
                implRowDualUpper.at(i) = bnd;
            }
          }
    }
  }
  timer.recordFinish(WEAKLY_DOMINATED_COLS);
}

void Presolve::setProblemStatus(const int s) {
  if (s == Infeasible)
    cout << "NOT-OPT status = 1, returned from solver after presolve: Problem "
            "infeasible.\n";
  else if (s == Unbounded)
    cout << "NOT-OPT status = 2, returned from solver after presolve: Problem "
            "unbounded.\n";
  else if (s == 0) {
    status = Optimal;
    return;
  } else
    cout << "unknown problem status returned from solver after presolve: " << s
         << endl;
  status = s;
}

void Presolve::setKKTcheckerData() {
  // after initializing equations.
  chk.setMatrixAR(numCol, numRow, ARstart, ARindex, ARvalue);
  chk.setFlags(flagRow, flagCol);
  chk.setBoundsCostRHS(colUpper, colLower, colCost, rowLower, rowUpper);
}

pair<double, double> Presolve::getNewBoundsDoubletonConstraint(int row, int col,
                                                               int j,
                                                               double aik,
                                                               double aij) {
  int i = row;

  double upp = HIGHS_CONST_INF;
  double low = -HIGHS_CONST_INF;

  if (aij > 0 && aik > 0) {
    if (colLower.at(col) > -HIGHS_CONST_INF)
      upp = (rowUpper.at(i) - aik * colLower.at(col)) / aij;
    if (colUpper.at(col) < HIGHS_CONST_INF)
      low = (rowLower.at(i) - aik * colUpper.at(col)) / aij;
  } else if (aij > 0 && aik < 0) {
    if (colLower.at(col) > -HIGHS_CONST_INF)
      low = (rowLower.at(i) - aik * colLower.at(col)) / aij;
    if (colUpper.at(col) < HIGHS_CONST_INF)
      upp = (rowUpper.at(i) - aik * colUpper.at(col)) / aij;
  } else if (aij < 0 && aik > 0) {
    if (colLower.at(col) > -HIGHS_CONST_INF)
      low = (rowUpper.at(i) - aik * colLower.at(col)) / aij;
    if (colUpper.at(col) < HIGHS_CONST_INF)
      upp = (rowLower.at(i) - aik * colUpper.at(col)) / aij;
  } else {
    if (colLower.at(col) > -HIGHS_CONST_INF)
      upp = (rowLower.at(i) - aik * colLower.at(col)) / aij;
    if (colUpper.at(col) < HIGHS_CONST_INF)
      low = (rowUpper.at(i) - aik * colUpper.at(col)) / aij;
  }

  if (j) {
  }  // surpress warning.

  return make_pair(low, upp);
}

void Presolve::removeFreeColumnSingleton(const int col, const int row,
                                         const int k) {
  if (iPrint > 0)
    cout << "PR: Free column singleton " << col << " removed. Row " << row
         << " removed." << endl;

  // modify costs
  vector<pair<int, double>> newCosts;
  int j;
  for (int kk = ARstart.at(row); kk < ARstart.at(row + 1); ++kk) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col) {
      newCosts.push_back(make_pair(j, colCost.at(j)));
      colCost.at(j) =
          colCost.at(j) - colCost.at(col) * ARvalue.at(kk) / Avalue.at(k);
    }
  }
  if (iKKTcheck == 1) chk.costs.push(newCosts);

  flagCol.at(col) = 0;
  postValue.push(colCost.at(col));
  fillStackRowBounds(row);

  valueColDual.at(col) = 0;
  valueRowDual.at(row) = -colCost.at(col) / Avalue.at(k);

  addChange(FREE_SING_COL, row, col);
  removeRow(row);

  countRemovedCols(FREE_SING_COL);
  countRemovedRows(FREE_SING_COL);
}

bool Presolve::removeColumnSingletonInDoubletonInequality(const int col,
                                                          const int i,
                                                          const int k) {
  // second column index j
  // second column row array index kk
  int j = -1;

  // count
  int kk = ARstart.at(i);
  while (kk < ARstart.at(i + 1)) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col)
      break;
    else
      ++kk;
  }
  if (kk == ARstart.at(i + 1))
    cout << "ERROR: nzRow[" << i << "]=2, but no second variable in row. \n";

  // only inequality case and case two singletons here,
  // others handled in doubleton equation
  if ((fabs(rowLower.at(i) - rowUpper.at(i)) < tol) && (nzCol.at(j) > 1))
    return false;

  // additional check if it is indeed implied free
  // needed since we handle inequalities and it may not be true
  // low and upp to be tighter than original bounds for variable col
  // so it is indeed implied free and we can remove it
  pair<double, double> p =
      getNewBoundsDoubletonConstraint(i, j, col, ARvalue.at(kk), Avalue.at(k));
  if (!(colLower.at(col) <= p.first && colUpper.at(col) >= p.second)) {
    return false;
  }

  postValue.push(ARvalue.at(kk));
  postValue.push(Avalue.at(k));

  // modify bounds on variable j, variable col (k) is substituted out
  // double aik = Avalue.at(k);
  // double aij = Avalue.at(kk);
  p = getNewBoundsDoubletonConstraint(i, col, j, Avalue.at(k), ARvalue.at(kk));
  double low = p.first;
  double upp = p.second;

  // add old bounds of xj to checker and for postsolve
  if (iKKTcheck == 1) {
    vector<pair<int, double>> bndsL, bndsU, costS;
    bndsL.push_back(make_pair(j, colLower.at(j)));
    bndsU.push_back(make_pair(j, colUpper.at(j)));
    costS.push_back(make_pair(j, colCost.at(j)));
    chk.cLowers.push(bndsL);
    chk.cUppers.push(bndsU);
    chk.costs.push(costS);
  }

  vector<double> bndsCol({colLower.at(col), colUpper.at(col), colCost.at(col)});
  vector<double> bndsJ({colLower.at(j), colUpper.at(j), colCost.at(j)});
  oldBounds.push(make_pair(col, bndsCol));
  oldBounds.push(make_pair(j, bndsJ));

  // modify bounds of xj
  if (low > colLower.at(j)) colLower.at(j) = low;
  if (upp < colUpper.at(j)) colUpper.at(j) = upp;

  // modify cost of xj
  colCost.at(j) =
      colCost.at(j) - colCost.at(col) * ARvalue.at(kk) / Avalue.at(k);

  // for postsolve: need the new bounds too
  // oldBounds.push_back(colLower.at(j)); oldBounds.push_back(colUpper.at(j));
  bndsJ.at(0) = (colLower.at(j));
  bndsJ.at(1) = (colUpper.at(j));
  bndsJ.at(2) = (colCost.at(j));
  oldBounds.push(make_pair(j, bndsJ));

  // remove col as free column singleton
  if (iPrint > 0)
    cout << "PR: Column singleton " << col
         << " in a doubleton inequality constraint removed. Row " << i
         << " removed. variable left is " << j << endl;

  flagCol.at(col) = 0;
  fillStackRowBounds(i);
  countRemovedCols(SING_COL_DOUBLETON_INEQ);
  countRemovedRows(SING_COL_DOUBLETON_INEQ);

  valueColDual.at(col) = 0;
  valueRowDual.at(i) =
      -colCost.at(col) /
      Avalue.at(k);  // may be changed later, depending on bounds.
  addChange(SING_COL_DOUBLETON_INEQ, i, col);

  // if not special case two column singletons
  if (nzCol.at(j) > 1)
    removeRow(i);
  else if (nzCol.at(j) == 1)
    removeSecondColumnSingletonInDoubletonRow(j, i);

  return true;
}

void Presolve::removeSecondColumnSingletonInDoubletonRow(const int j,
                                                         const int i) {
  // case two singleton columns
  // when we get here bounds on xj are updated so we can choose low/upper one
  // depending on the cost of xj
  flagRow.at(i) = 0;
  double value;
  if (colCost.at(j) > 0) {
    if (colLower.at(j) == -HIGHS_CONST_INF) {
      if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
      status = Unbounded;
      return;
    }
    value = colLower.at(j);
  } else if (colCost.at(j) < 0) {
    if (colUpper.at(j) == HIGHS_CONST_INF) {
      if (iPrint > 0) cout << "PR: Problem unbounded." << endl;
      status = Unbounded;
      return;
    }
    value = colUpper.at(j);
  } else {  //(colCost.at(j) == 0)
    if (colUpper.at(j) >= 0 && colLower.at(j) <= 0)
      value = 0;
    else if (fabs(colUpper.at(j)) < fabs(colLower.at(j)))
      value = colUpper.at(j);
    else
      value = colLower.at(j);
  }
  setPrimalValue(j, value);
  addChange(SING_COL_DOUBLETON_INEQ_SECOND_SING_COL, 0, j);
  if (iPrint > 0)
    cout << "PR: Second singleton column " << j << " in doubleton row " << i
         << " removed.\n";
  countRemovedCols(SING_COL_DOUBLETON_INEQ);
  singCol.remove(j);
}

void Presolve::removeColumnSingletons() {
  int i, k, col;
  list<int>::iterator it = singCol.begin();

  while (it != singCol.end()) {
    if (flagCol[*it]) {
      col = *it;
      k = getSingColElementIndexInA(col);
      i = Aindex.at(k);

      // free
      if (colLower.at(col) == -HIGHS_CONST_INF &&
          colUpper.at(col) == HIGHS_CONST_INF) {
        timer.recordStart(FREE_SING_COL);
        removeFreeColumnSingleton(col, i, k);
        it = singCol.erase(it);
        timer.recordFinish(FREE_SING_COL);
        continue;
      }
      // singleton column in a doubleton inequality
      // case two column singletons
      else if (nzRow.at(i) == 2) {
        timer.recordStart(SING_COL_DOUBLETON_INEQ);
        bool result = removeColumnSingletonInDoubletonInequality(col, i, k);
        timer.recordFinish(SING_COL_DOUBLETON_INEQ);
        if (result) {
          it = singCol.erase(it);
          continue;
        }
      }
      // implied free
      else {
        timer.recordStart(IMPLIED_FREE_SING_COL);
        bool result = removeIfImpliedFree(col, i, k);
        timer.recordFinish(IMPLIED_FREE_SING_COL);
        if (result) {
          it = singCol.erase(it);
          continue;
        }
      }
      it++;
    } else
      it = singCol.erase(it);
  }
}

pair<double, double> Presolve::getBoundsImpliedFree(double lowInit,
                                                    double uppInit,
                                                    const int col, const int i,
                                                    const int k) {
  double low = lowInit;
  double upp = uppInit;

  // use implied bounds with original bounds
  int j;
  double l, u;
  // if at any stage low becomes  or upp becomes inf break loop
  // can't use bounds for variables generated by the same row.
  // low
  for (int kk = ARstart.at(i); kk < ARstart.at(i + 1); ++kk) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col) {
      // check if new bounds are precisely implied bounds from same row
      if (i != implColLowerRowIndex.at(j))
        l = max(colLower.at(j), implColLower.at(j));
      else
        l = colLower.at(j);
      if (i != implColUpperRowIndex.at(j))
        u = min(colUpper.at(j), implColUpper.at(j));
      else
        u = colUpper.at(j);

      if ((Avalue.at(k) < 0 && ARvalue.at(kk) > 0) ||
          (Avalue.at(k) > 0 && ARvalue.at(kk) < 0))
        if (l == -HIGHS_CONST_INF) {
          low = -HIGHS_CONST_INF;
          break;
        } else
          low -= ARvalue.at(kk) * l;
      else if (u == HIGHS_CONST_INF) {
        low = -HIGHS_CONST_INF;
        break;
      } else
        low -= ARvalue.at(kk) * u;
    }
  }
  // upp
  for (int kk = ARstart.at(i); kk < ARstart.at(i + 1); ++kk) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col) {
      // check if new bounds are precisely implied bounds from same row
      if (i != implColLowerRowIndex.at(j))
        l = max(colLower.at(j), implColLower.at(j));
      else
        l = colLower.at(j);
      if (i != implColUpperRowIndex.at(j))
        u = min(colUpper.at(j), implColUpper.at(j));
      else
        u = colUpper.at(j);
      // if at any stage low becomes  or upp becomes inf it's not implied free
      // low::
      if ((Avalue.at(k) < 0 && ARvalue.at(kk) > 0) ||
          (Avalue.at(k) > 0 && ARvalue.at(kk) < 0))
        if (u == HIGHS_CONST_INF) {
          upp = HIGHS_CONST_INF;
          break;
        } else
          upp -= ARvalue.at(kk) * u;
      else if (l == -HIGHS_CONST_INF) {
        upp = HIGHS_CONST_INF;
        break;
      } else
        upp -= ARvalue.at(kk) * l;
    }
  }
  return make_pair(low, upp);
}

void Presolve::removeImpliedFreeColumn(const int col, const int i,
                                       const int k) {
  if (iPrint > 0)
    cout << "PR: Implied free column singleton " << col << " removed.  Row "
         << i << " removed." << endl;

  countRemovedCols(IMPLIED_FREE_SING_COL);
  countRemovedRows(IMPLIED_FREE_SING_COL);

  // modify costs
  int j;
  vector<pair<int, double>> newCosts;
  for (int kk = ARstart.at(i); kk < ARstart.at(i + 1); ++kk) {
    j = ARindex.at(kk);
    if (flagCol.at(j) && j != col) {
      newCosts.push_back(make_pair(j, colCost.at(j)));
      colCost.at(j) =
          colCost.at(j) - colCost.at(col) * ARvalue.at(kk) / Avalue.at(k);
    }
  }
  if (iKKTcheck == 1) chk.costs.push(newCosts);

  flagCol.at(col) = 0;
  postValue.push(colCost.at(col));
  fillStackRowBounds(i);

  valueColDual.at(col) = 0;
  valueRowDual.at(i) = -colCost.at(col) / Avalue.at(k);
  addChange(IMPLIED_FREE_SING_COL, i, col);
  removeRow(i);
}

bool Presolve::removeIfImpliedFree(int col, int i, int k) {
  // first find which bound is active for row i
  // A'y + c = z so yi = -ci/aij
  double aij = getaij(i, col);
  if (aij != Avalue.at(k)) cout << "ERROR during implied free";
  double yi = -colCost.at(col) / aij;
  double low, upp;

  if (yi > 0) {
    if (rowUpper.at(i) == HIGHS_CONST_INF) return false;
    low = rowUpper.at(i);
    upp = rowUpper.at(i);
  } else if (yi < 0) {
    if (rowLower.at(i) == -HIGHS_CONST_INF) return false;
    low = rowLower.at(i);
    upp = rowLower.at(i);
  } else {
    low = rowLower.at(i);
    upp = rowUpper.at(i);
  }

  pair<double, double> p = getBoundsImpliedFree(low, upp, col, i, k);
  low = p.first;
  upp = p.second;

  if (low > -HIGHS_CONST_INF) low = low / Avalue.at(k);
  if (upp < HIGHS_CONST_INF) upp = upp / Avalue.at(k);

  // if implied free
  if (colLower.at(col) <= low && low <= upp && upp <= colUpper.at(col)) {
    removeImpliedFreeColumn(col, i, k);
    return true;
  }
  // else calculate implied bounds
  else if (colLower.at(col) <= low && low <= upp) {
    if (implColLower.at(col) < low) {
      implColLower.at(col) = low;
      implColUpperRowIndex.at(col) = i;
    }
    // JAJH(190419): Segfault here since i is a row index and
    // segfaults (on dcp1) due to i=4899 exceeding column dimension of
    // 3007. Hence correction. Also pattern-matches the next case :-)
    //    implColDualUpper.at(i) = 0;
    implColDualUpper.at(col) = 0;
  } else if (low <= upp && upp <= colUpper.at(col)) {
    if (implColUpper.at(col) > upp) {
      implColUpper.at(col) = upp;
      implColUpperRowIndex.at(col) = i;
    }
    implColDualLower.at(col) = 0;
  }

  return false;
}

// used to remove column too, now possible to just modify bounds
void Presolve::removeRow(int i) {
  hasChange = true;
  flagRow.at(i) = 0;
  for (int k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
    int j = ARindex.at(k);
    if (flagCol.at(j)) {
      nzCol.at(j)--;
      // if now singleton add to list
      if (nzCol.at(j) == 1) {
        int index = getSingColElementIndexInA(j);
        if (index >= 0)
          singCol.push_back(j);
        else
          cout << "Warning: Column " << j
               << " with 1 nz but not in singCol or? Row removing of " << i
               << ". Ignored.\n";
      }
      // if it was a singleton column remove from list and problem
      if (nzCol.at(j) == 0) removeEmptyColumn(j);
    }
  }
}

void Presolve::fillStackRowBounds(int row) {
  postValue.push(rowUpper.at(row));
  postValue.push(rowLower.at(row));
}

pair<double, double> Presolve::getImpliedRowBounds(int row) {
  double g = 0;
  double h = 0;

  int col;
  for (int k = ARstart.at(row); k < ARstart.at(row + 1); ++k) {
    col = ARindex.at(k);
    if (flagCol.at(col)) {
      if (ARvalue.at(k) < 0) {
        if (colUpper.at(col) < HIGHS_CONST_INF)
          g += ARvalue.at(k) * colUpper.at(col);
        else {
          g = -HIGHS_CONST_INF;
          break;
        }
      } else {
        if (colLower.at(col) > -HIGHS_CONST_INF)
          g += ARvalue.at(k) * colLower.at(col);
        else {
          g = -HIGHS_CONST_INF;
          break;
        }
      }
    }
  }

  for (int k = ARstart.at(row); k < ARstart.at(row + 1); ++k) {
    col = ARindex.at(k);
    if (flagCol.at(col)) {
      if (ARvalue.at(k) < 0) {
        if (colLower.at(col) > -HIGHS_CONST_INF)
          h += ARvalue.at(k) * colLower.at(col);
        else {
          h = HIGHS_CONST_INF;
          break;
        }
      } else {
        if (colUpper.at(col) < HIGHS_CONST_INF)
          h += ARvalue.at(k) * colUpper.at(col);
        else {
          h = HIGHS_CONST_INF;
          break;
        }
      }
    }
  }
  return make_pair(g, h);
}

void Presolve::setVariablesToBoundForForcingRow(const int row,
                                                const bool isLower) {
  int k, col;
  if (iPrint > 0)
    cout << "PR: Forcing row " << row
         << " removed. Following variables too:   nzRow=" << nzRow.at(row)
         << endl;

  flagRow.at(row) = 0;
  addChange(FORCING_ROW, row, 0);
  k = ARstart.at(row);
  while (k < ARstart.at(row + 1)) {
    col = ARindex.at(k);
    if (flagCol.at(col)) {
      double value;
      if ((ARvalue.at(k) < 0 && isLower) || (ARvalue.at(k) > 0 && !isLower))
        value = colUpper.at(col);
      else
        value = colLower.at(col);

      setPrimalValue(col, value);
      valueColDual.at(col) = colCost.at(col);
      vector<double> bnds({colLower.at(col), colUpper.at(col)});
      oldBounds.push(make_pair(col, bnds));
      addChange(FORCING_ROW_VARIABLE, 0, col);

      if (iPrint > 0)
        cout << "PR:      Variable  " << col << " := " << value << endl;
      countRemovedCols(FORCING_ROW);
    }
    ++k;
  }

  if (nzRow.at(row) == 1) singRow.remove(row);

  countRemovedRows(FORCING_ROW);
}

void Presolve::dominatedConstraintProcedure(const int i, const double g,
                                            const double h) {
  int j;
  double val;
  if (h < HIGHS_CONST_INF) {
    // fill in implied bounds arrays
    if (h < implRowValueUpper.at(i)) {
      implRowValueUpper.at(i) = h;
    }
    if (h <= rowUpper.at(i)) implRowDualLower.at(i) = 0;

    // calculate implied bounds for discovering free column singletons
    for (int k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
      j = ARindex.at(k);
      if (flagCol.at(j)) {
        if (ARvalue.at(k) < 0 && colLower.at(j) > -HIGHS_CONST_INF) {
          val = (rowLower.at(i) - h) / ARvalue.at(k) + colLower.at(j);
          if (val < implColUpper.at(j)) {
            implColUpper.at(j) = val;
            implColUpperRowIndex.at(j) = i;
          }
        } else if (ARvalue.at(k) > 0 && colUpper.at(j) < HIGHS_CONST_INF) {
          val = (rowLower.at(i) - h) / ARvalue.at(k) + colUpper.at(j);
          if (val > implColLower.at(j)) {
            implColLower.at(j) = val;
            implColLowerRowIndex.at(j) = i;
          }
        }
      }
    }
  }
  if (g > -HIGHS_CONST_INF) {
    // fill in implied bounds arrays
    if (g > implRowValueLower.at(i)) {
      implRowValueLower.at(i) = g;
    }
    if (g >= rowLower.at(i)) implRowDualUpper.at(i) = 0;

    // calculate implied bounds for discovering free column singletons
    for (int k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
      int j = ARindex.at(k);
      if (flagCol.at(j)) {
        if (ARvalue.at(k) < 0 && colUpper.at(j) < HIGHS_CONST_INF) {
          val = (rowUpper.at(i) - g) / ARvalue.at(k) + colUpper.at(j);
          if (val > implColLower.at(j)) {
            implColLower.at(j) = val;
            implColLowerRowIndex.at(j) = i;
          }
        } else if (ARvalue.at(k) > 0 && colLower.at(j) > -HIGHS_CONST_INF) {
          val = (rowUpper.at(i) - g) / ARvalue.at(k) + colLower.at(j);
          if (val < implColUpper.at(j)) {
            implColUpper.at(j) = val;
            implColUpperRowIndex.at(j) = i;
          }
        }
      }
    }
  }
}

void Presolve::removeForcingConstraints(int mainIter) {
  double g, h;
  pair<double, double> implBounds;

  for (int i = 0; i < numRow; ++i)
    if (flagRow.at(i)) {
      if (nzRow.at(i) == 0) {
        removeEmptyRow(i);
        countRemovedRows(EMPTY_ROW);
        continue;
      }

      // removeRowSingletons will handle just after removeForcingConstraints
      if (nzRow.at(i) == 1) continue;

      timer.recordStart(FORCING_ROW);
      implBounds = getImpliedRowBounds(i);

      g = implBounds.first;
      h = implBounds.second;

      // Infeasible row
      if (g > rowUpper.at(i) || h < rowLower.at(i)) {
        if (iPrint > 0) cout << "PR: Problem infeasible." << endl;
        status = Infeasible;
        timer.recordFinish(FORCING_ROW);
        return;
      }
      // Forcing row
      else if (g == rowUpper.at(i)) {
        setVariablesToBoundForForcingRow(i, true);
      } else if (h == rowLower.at(i)) {
        setVariablesToBoundForForcingRow(i, false);
      }
      // Redundant row
      else if (g >= rowLower.at(i) && h <= rowUpper.at(i)) {
        removeRow(i);
        addChange(REDUNDANT_ROW, i, 0);
        if (iPrint > 0)
          cout << "PR: Redundant row " << i << " removed." << endl;
        countRemovedRows(REDUNDANT_ROW);
      }
      // Dominated constraints
      else {
        timer.recordFinish(FORCING_ROW);
        timer.recordStart(DOMINATED_ROW_BOUNDS);
        dominatedConstraintProcedure(i, g, h);
        timer.recordFinish(DOMINATED_ROW_BOUNDS);
        continue;
      }
      timer.recordFinish(FORCING_ROW);
    }
  if (mainIter) {
  }  // surpress warning.
}

void Presolve::removeRowSingletons() {
  timer.recordStart(SING_ROW);
  int i;
  int singRowZ = singRow.size();
  /*
  if (singRowZ == 36) {
    printf("JAJH: singRow.size() = %d\n", singRowZ);fflush(stdout);
  }
  */
  while (!(singRow.empty())) {
    i = singRow.front();
    singRow.pop_front();

    assert(flagRow[i]);

    int k = getSingRowElementIndexInAR(i);
    // JAJH(190419): This throws a segfault with greenbea and greenbeb since
    // k=-1
    if (k < 0) {
      printf("In removeRowSingletons: %d = k < 0\n", k);
      printf("   Occurs for case when initial singRow.size() = %d\n", singRowZ);
      fflush(stdout);
    }
    int j = ARindex.at(k);

    // add old bounds OF X to checker and for postsolve
    if (iKKTcheck == 1) {
      vector<pair<int, double>> bndsL, bndsU, costS;
      bndsL.push_back(make_pair(j, colLower.at(j)));
      bndsU.push_back(make_pair(j, colUpper.at(j)));
      chk.cLowers.push(bndsL);
      chk.cUppers.push(bndsU);
    }

    vector<double> bnds(
        {colLower.at(j), colUpper.at(j), rowLower.at(i), rowUpper.at(i)});
    oldBounds.push(make_pair(j, bnds));

    double aij = ARvalue.at(k);
    /*		//before update bounds of x take it out of rows with implied row
    bounds for (int r = Astart.at(j); r<Aend.at(j); r++) { if
    (flagRow[Aindex[r]]) { int rr = Aindex[r]; if (implRowValueLower[rr] >
    -HIGHS_CONST_INF) { if (aij > 0) implRowValueLower[rr] =
    implRowValueLower[rr] - aij*colLower.at(j); else implRowValueLower[rr] =
    implRowValueLower[rr] - aij*colUpper.at(j);
                    }
                    if (implRowValueUpper[rr] < HIGHS_CONST_INF) {
                            if (aij > 0)
                                    implRowValueUpper[rr] =
    implRowValueUpper[rr] - aij*colUpper.at(j); else implRowValueUpper[rr] =
    implRowValueUpper[rr] - aij*colLower.at(j);
                    }
            }
    }*/

    // update bounds of X
    if (aij > 0) {
      if (rowLower.at(i) != -HIGHS_CONST_INF)
        colLower.at(j) =
            max(max(rowLower.at(i) / aij, -HIGHS_CONST_INF), colLower.at(j));
      if (rowUpper.at(i) != HIGHS_CONST_INF)
        colUpper.at(j) =
            min(min(rowUpper.at(i) / aij, HIGHS_CONST_INF), colUpper.at(j));
    } else if (aij < 0) {
      if (rowLower.at(i) != -HIGHS_CONST_INF)
        colUpper.at(j) =
            min(min(rowLower.at(i) / aij, HIGHS_CONST_INF), colUpper.at(j));
      if (rowUpper.at(i) != HIGHS_CONST_INF)
        colLower.at(j) =
            max(max(rowUpper.at(i) / aij, -HIGHS_CONST_INF), colLower.at(j));
    }

    /*		//after update bounds of x add to rows with implied row bounds
    for (int r = Astart.at(j); r<Aend.at(j); r++) {
            if (flagRow[r]) {
                    int rr = Aindex[r];
                    if (implRowValueLower[rr] > -HIGHS_CONST_INF) {
                            if (aij > 0)
                                    implRowValueLower[rr] =
    implRowValueLower[rr] + aij*colLower.at(j); else implRowValueLower[rr] =
    implRowValueLower[rr] + aij*colUpper.at(j);
                    }
                    if (implRowValueUpper[rr] < HIGHS_CONST_INF) {
                            if (aij > 0)
                                    implRowValueUpper[rr] =
    implRowValueUpper[rr] + aij*colUpper.at(j); else implRowValueUpper[rr] =
    implRowValueUpper[rr] + aij*colLower.at(j);
                    }
            }
    }*/

    // check for feasibility
    if (colLower.at(j) > colUpper.at(j) + tol) {
      status = Infeasible;
      timer.recordFinish(SING_ROW);
      return;
    }

    if (iPrint > 0)
      cout << "PR: Singleton row " << i << " removed. Bounds of variable  " << j
           << " modified: l= " << colLower.at(j) << " u=" << colUpper.at(j)
           << ", aij = " << aij << endl;

    addChange(SING_ROW, i, j);
    postValue.push(colCost.at(j));
    removeRow(i);

    if (flagCol.at(j) && colLower.at(j) == colUpper.at(j)) removeIfFixed(j);

    countRemovedRows(SING_ROW);
  }
  timer.recordFinish(SING_ROW);
}

void Presolve::addChange(PresolveRule type, int row, int col) {
  change ch;
  ch.type = type;
  ch.row = row;
  ch.col = col;
  chng.push(ch);

  if (type < PRESOLVE_RULES_COUNT) timer.addChange(type);
}

// when setting a value to a primal variable and eliminating row update b,
// singleton Rows linked list, number of nonzeros in rows
void Presolve::setPrimalValue(int j, double value) {
  flagCol.at(j) = 0;
  if (!hasChange) hasChange = true;
  valuePrimal.at(j) = value;

  // update nonzeros
  for (int k = Astart.at(j); k < Aend.at(j); ++k) {
    int row = Aindex.at(k);
    if (flagRow.at(row)) {
      nzRow.at(row)--;

      // update singleton row list
      if (nzRow.at(row) == 1)
        singRow.push_back(row);
      else if (nzRow.at(row) == 0)
        singRow.remove(row);
    }
  }

  // update values if necessary
  if (fabs(value) > 0) {
    // RHS
    vector<pair<int, double>> bndsL, bndsU;

    for (int k = Astart.at(j); k < Aend.at(j); ++k)
      if (flagRow.at(Aindex.at(k))) {
        if (iKKTcheck == 1) {
          bndsL.push_back(make_pair(Aindex.at(k), rowLower.at(Aindex.at(k))));
          bndsU.push_back(make_pair(Aindex.at(k), rowUpper.at(Aindex.at(k))));
        }
        if (rowLower.at(Aindex.at(k)) > -HIGHS_CONST_INF)
          rowLower.at(Aindex.at(k)) -= Avalue.at(k) * value;
        if (rowUpper.at(Aindex.at(k)) < HIGHS_CONST_INF)
          rowUpper.at(Aindex.at(k)) -= Avalue.at(k) * value;

        if (implRowValueLower.at(Aindex.at(k)) > -HIGHS_CONST_INF)
          implRowValueLower.at(Aindex.at(k)) -= Avalue.at(k) * value;
        if (implRowValueUpper.at(Aindex.at(k)) < HIGHS_CONST_INF)
          implRowValueUpper.at(Aindex.at(k)) -= Avalue.at(k) * value;
      }

    if (iKKTcheck == 1) {
      chk.rLowers.push(bndsL);
      chk.rUppers.push(bndsU);
    }

    // shift objective
    if (colCost.at(j) != 0) objShift += colCost.at(j) * value;
  }
}

void Presolve::checkForChanges(int iteration) {
  if (iteration <= 2) {
    // flagCol has one more element at end which is zero
    // from removeDoubletonEquatoins, needed for AR matrix manipulation
    if (none_of(flagCol.begin(), flagCol.begin() + numCol,
                [](int i) { return i == 0; }) &&
        none_of(flagRow.begin(), flagRow.begin() + numRow,
                [](int i) { return i == 0; })) {
      if (iPrint > 0)
        cout << "PR: No variables were eliminated at presolve." << endl;
      noPostSolve = true;
      return;
    }
  }
  resizeProblem();
  status = stat::Reduced;
}

// void Presolve::reportTimes() {
//   int reportList[] = {EMPTY_ROW,
//                       FIXED_COL,
//                       SING_ROW,
//                       DOUBLETON_EQUATION,
//                       FORCING_ROW,
//                       REDUNDANT_ROW,
//                       FREE_SING_COL,
//                       SING_COL_DOUBLETON_INEQ,
//                       IMPLIED_FREE_SING_COL,
//                       DOMINATED_COLS,
//                       WEAKLY_DOMINATED_COLS};
//   int reportCount = sizeof(reportList) / sizeof(int);

//   printf("Presolve rules ");
//   for (int i = 0; i < reportCount; ++i) {
//     printf(" %s", timer.itemNames[reportList[i]].c_str());
//     cout << flush;
//   }

//   printf("\n");
//   cout << "Time spent     " << flush;
//   for (int i = 0; i < reportCount; ++i) {
//     float f = (float)timer.itemTicks[reportList[i]];
//     if (f < 0.01)
//       cout << setw(4) << " <.01 ";
//     else
//       printf(" %3.2f ", f);
//   }
//   printf("\n");
// }

// void Presolve::recordCounts(const string fileName) {
//   ofstream myfile;
//   myfile.open(fileName.c_str(), ios::app);
//   int reportList[] = {EMPTY_ROW,
//                       FIXED_COL,
//                       SING_ROW,
//                       DOUBLETON_EQUATION,
//                       FORCING_ROW,
//                       REDUNDANT_ROW,
//                       FREE_SING_COL,
//                       SING_COL_DOUBLETON_INEQ,
//                       IMPLIED_FREE_SING_COL,
//                       DOMINATED_COLS,
//                       WEAKLY_DOMINATED_COLS,
//                       EMPTY_COL};
//   int reportCount = sizeof(reportList) / sizeof(int);

//   myfile << "Problem " << modelName << ":\n";
//   myfile << "Rule   , removed rows , removed cols , time  \n";

//   int cRows = 0, cCols = 0;
//   for (int i = 0; i < reportCount; ++i) {
//     float f = (float)timer.itemTicks[reportList[i]];

//     myfile << setw(7) << timer.itemNames[reportList[i]].c_str() << ", "
//            << setw(7) << countRemovedRows[reportList[i]] << ", " << setw(7)
//            << countRemovedCols[reportList[i]] << ", ";
//     if (f < 0.001)
//       myfile << setw(7) << " <.001 ";
//     else
//       myfile << setw(7) << setprecision(3) << f;
//     myfile << endl;

//     cRows += countRemovedRows[reportList[i]];
//     cCols += countRemovedCols[reportList[i]];
//   }

//   if (!noPostSolve) {
//     if (cRows != numRowOriginal - numRow) cout << "Wrong row reduction
//     count\n"; if (cCols != numColOriginal - numCol) cout << "Wrong col
//     reduction count\n";

//     myfile << setw(7) << "Total "
//            << ", " << setw(7) << numRowOriginal - numRow << ", " << setw(7)
//            << numColOriginal - numCol;
//   } else {
//     myfile << setw(7) << "Total "
//            << ", " << setw(7) << 0 << ", " << setw(7) << 0;
//   }
//   myfile << endl << " \\\\ " << endl;
//   myfile.close();
// }

void Presolve::resizeImpliedBounds() {
  // implied bounds for crashes
  // row duals
  vector<double> temp = implRowDualLower;
  vector<double> teup = implRowDualUpper;
  implRowDualLower.resize(numRow);
  implRowDualUpper.resize(numRow);

  int k = 0;
  for (int i = 0; i < numRowOriginal; ++i)
    if (flagRow.at(i)) {
      implRowDualLower.at(k) = temp.at(i);
      implRowDualUpper.at(k) = teup.at(i);
      k++;
    }

  // row value
  temp = implRowValueLower;
  teup = implRowValueUpper;
  implRowValueLower.resize(numRow);
  implRowValueUpper.resize(numRow);
  k = 0;
  for (int i = 0; i < numRowOriginal; ++i)
    if (flagRow.at(i)) {
      if (temp.at(i) < rowLower.at(i)) temp.at(i) = rowLower.at(i);
      implRowValueLower.at(k) = temp.at(i);
      if (teup.at(i) > rowUpper.at(i)) teup.at(i) = rowUpper.at(i);
      implRowValueUpper.at(k) = teup.at(i);
      k++;
    }

  // column dual
  temp = implColDualLower;
  teup = implColDualUpper;
  implColDualLower.resize(numCol);
  implColDualUpper.resize(numCol);

  k = 0;
  for (int i = 0; i < numColOriginal; ++i)
    if (flagCol.at(i)) {
      implColDualLower.at(k) = temp.at(i);
      implColDualUpper.at(k) = teup.at(i);
      k++;
    }

  // column value
  temp = implColLower;
  teup = implColUpper;
  implColLower.resize(numCol);
  implColUpper.resize(numCol);

  k = 0;
  for (int i = 0; i < numColOriginal; ++i)
    if (flagCol.at(i)) {
      if (temp.at(i) < colLower.at(i)) temp.at(i) = colLower.at(i);
      implColLower.at(k) = temp.at(i);
      if (teup.at(i) > colUpper.at(i)) teup.at(i) = colUpper.at(i);
      implColUpper.at(k) = teup.at(i);
      k++;
    }
}

int Presolve::getSingRowElementIndexInAR(int i) {
  int k = ARstart.at(i);
  while (!flagCol.at(ARindex.at(k))) ++k;
  if (k >= ARstart.at(i + 1)) {
    cout << "Error during presolve: no variable found in singleton row " << i
         << endl;
    return -1;
  }
  int rest = k + 1;
  while (rest < ARstart.at(i + 1) && !flagCol.at(ARindex.at(rest))) ++rest;
  if (rest < ARstart.at(i + 1)) {
    cout << "Error during presolve: more variables found in singleton row " << i
         << endl;
    return -1;
  }
  return k;
}

int Presolve::getSingColElementIndexInA(int j) {
  int k = Astart.at(j);
  while (!flagRow.at(Aindex.at(k))) ++k;
  if (k >= Aend.at(j)) {
    cout << "Error during presolve: no variable found in singleton col " << j
         << ".";
    return -1;
  }
  int rest = k + 1;
  while (rest < Aend.at(j) && !flagRow.at(Aindex.at(rest))) ++rest;
  if (rest < Aend.at(j)) {
    cout << "Error during presolve: more variables found in singleton col " << j
         << ".";
    return -1;
  }
  return k;
}

void Presolve::testAnAR(int post) {
  int rows = numRow;
  int cols = numCol;
  int i, j, k;

  double valueA = 0;
  double valueAR = 0;
  bool hasValueA, hasValueAR;

  if (post) {
    rows = numRowOriginal;
    cols = numColOriginal;
  }

  // check that A = AR
  for (i = 0; i < rows; ++i) {
    for (j = 0; j < cols; ++j) {
      if (post == 0)
        if (!flagRow.at(i) || !flagCol.at(j)) continue;
      hasValueA = false;
      for (k = Astart.at(j); k < Aend.at(j); ++k)
        if (Aindex.at(k) == i) {
          hasValueA = true;
          valueA = Avalue.at(k);
        }

      hasValueAR = false;
      for (k = ARstart.at(i); k < ARstart.at(i + 1); ++k)
        if (ARindex.at(k) == j) {
          hasValueAR = true;
          valueAR = ARvalue.at(k);
        }

      if (hasValueA != hasValueAR)
        cout << "    MATRIX is0   DIFF row=" << i << " col=" << j
             << "           ------------A: " << hasValueA
             << "  AR: " << hasValueAR << endl;
      else if (hasValueA && valueA != valueAR)
        cout << "    MATRIX VAL  DIFF row=" << i << " col=" << j
             << "           ------------A: " << valueA << "  AR: " << valueAR
             << endl;
    }
  }

  if (post == 0) {
    // check nz
    int nz = 0;
    for (i = 0; i < rows; ++i) {
      if (!flagRow.at(i)) continue;
      nz = 0;
      for (k = ARstart.at(i); k < ARstart.at(i + 1); ++k)
        if (flagCol.at(ARindex.at(k))) nz++;
      if (nz != nzRow.at(i))
        cout << "    NZ ROW      DIFF row=" << i << " nzRow=" << nzRow.at(i)
             << " actually " << nz << "------------" << endl;
    }

    for (j = 0; j < cols; ++j) {
      if (!flagCol.at(j)) continue;
      nz = 0;
      for (k = Astart.at(j); k < Aend.at(j); ++k)
        if (flagRow.at(Aindex.at(k))) nz++;
      if (nz != nzCol.at(j))
        cout << "    NZ COL      DIFF col=" << j << " nzCol=" << nzCol.at(j)
             << " actually " << nz << "------------" << endl;
    }
  }
}

// todo: error reporting.
HighsPostsolveStatus Presolve::postsolve(const HighsSolution& reduced_solution,
                                         HighsSolution& recovered_solution) {
  colValue = reduced_solution.col_value;
  colDual = reduced_solution.col_dual;
  rowDual = reduced_solution.row_dual;

  // todo: add nonbasic flag to Solution.
  // todo: change to new basis info structure later or keep.
  // basis info and solution should be somehow connected to each other.

  // here noPostSolve is always false. If the problem has not been reduced
  // Presolve::postsolve(..) is never called. todo: delete block below. For now
  // left just as legacy.
  // if (noPostSolve) {
  //   // set valuePrimal
  //   for (int i = 0; i < numCol; ++i) {
  //     valuePrimal.at(i) = colValue.at(i);
  //     valueColDual.at(i) = colDual.at(i);
  //   }
  //   for (int i = 0; i < numRow; ++i) valueRowDual.at(i) = rowDual.at(i);
  //   // For KKT check: first check solverz` results before we do any postsolve
  //   if (iKKTcheck == 1) {
  //     chk.passSolution(colValue, colDual, rowDual);
  //     chk.passBasis(col_status, row_status);
  //     chk.makeKKTCheck();
  //   }
  //   // testBasisMatrixSingularity();
  //   return HighsPostsolveStatus::NoPostsolve;
  // }

  // For KKT check: first check solver results before we do any postsolve
  if (iKKTcheck == 1) {
    cout << "----KKT check on HiGHS solution-----\n";

    chk.passSolution(colValue, colDual, rowDual);
    chk.passBasis(col_status, row_status);
    chk.makeKKTCheck();
  }
  // So there have been changes definitely ->
  makeACopy();  // so we can efficiently calculate primal and dual values

  //	iKKTcheck = false;
  // set corresponding parts of solution vectors:
  int j = 0;
  vector<int> eqIndexOfReduced(numCol, -1);
  vector<int> eqIndexOfReduROW(numRow, -1);
  for (int i = 0; i < numColOriginal; ++i)
    if (cIndex.at(i) > -1) {
      eqIndexOfReduced.at(j) = i;
      ++j;
    }
  j = 0;
  for (int i = 0; i < numRowOriginal; ++i)
    if (rIndex.at(i) > -1) {
      eqIndexOfReduROW.at(j) = i;
      ++j;
    }

  vector<HighsBasisStatus> temp_col_status = col_status;
  vector<HighsBasisStatus> temp_row_status = row_status;

  nonbasicFlag.assign(numColOriginal + numRowOriginal, 1);
  col_status.assign(numColOriginal, HighsBasisStatus::NONBASIC);  // Was LOWER
  row_status.assign(numRowOriginal, HighsBasisStatus::NONBASIC);  // Was LOWER

  for (int i = 0; i < numCol; ++i) {
    int iCol = eqIndexOfReduced.at(i);
    valuePrimal[iCol] = colValue.at(i);
    valueColDual[iCol] = colDual.at(i);
    col_status.at(iCol) = temp_col_status.at(i);
  }

  for (int i = 0; i < numRow; ++i) {
    int iRow = eqIndexOfReduROW.at(i);
    valueRowDual[iRow] = rowDual.at(i);
    row_status.at(iRow) = temp_row_status.at(i);
  }

  // cmpNBF(-1, -1);

  double z;
  vector<int> fRjs;
  while (!chng.empty()) {
    change c = chng.top();
    chng.pop();
    // cout<<"chng.pop:       "<<c.col<<"       "<<c.row << endl;

    setBasisElement(c);
    if (iKKTcheck == 1) chk.replaceBasis(col_status, row_status);

    switch (c.type) {
      case DOUBLETON_EQUATION: {  // Doubleton equation row
        getDualsDoubletonEquation(c.row, c.col);

        if (iKKTcheck == 1) {
          if (chk.print == 1)
            cout
                << "----KKT check after doubleton equation re-introduced. Row: "
                << c.row << ", column " << c.col << " -----\n";
          chk.addChange(17, c.row, c.col, valuePrimal[c.col],
                        valueColDual[c.col], valueRowDual[c.row]);
          chk.replaceBasis(col_status, row_status);
          chk.makeKKTCheck();
        }
        // exit(2);
        break;
      }
      case DOUBLETON_EQUATION_ROW_BOUNDS_UPDATE: {
        // new bounds from doubleton equation, retrieve old ones
        // just for KKT check, not called otherwise
        chk.addChange(171, c.row, c.col, 0, 0, 0);
        break;
      }
      case DOUBLETON_EQUATION_NEW_X_NONZERO: {
        // matrix transformation from doubleton equation, case x still there
        // case new x is not 0
        // just change value of entry in row for x

        int indi;
        for (indi = ARstart[c.row]; indi < ARstart[c.row + 1]; ++indi)
          if (ARindex.at(indi) == c.col) break;
        ARvalue.at(indi) = postValue.top();
        for (indi = Astart[c.col]; indi < Aend[c.col]; ++indi)
          if (Aindex.at(indi) == c.row) break;
        Avalue.at(indi) = postValue.top();

        if (iKKTcheck == 1)
          chk.addChange(172, c.row, c.col, postValue.top(), 0, 0);
        postValue.pop();

        break;
      }
      case DOUBLETON_EQUATION_X_ZERO_INITIALLY: {
        // matrix transformation from doubleton equation, retrieve old value
        // case when row does not have x initially: entries for row i swap x and
        // y cols

        int indi, yindex;
        yindex = (int)postValue.top();
        postValue.pop();

        // reverse AR for case when x is zero and y entry has moved
        for (indi = ARstart[c.row]; indi < ARstart[c.row + 1]; ++indi)
          if (ARindex.at(indi) == c.col) break;
        ARvalue.at(indi) = postValue.top();
        ARindex.at(indi) = yindex;

        // reverse A for case when x is zero and y entry has moved
        for (indi = Astart[c.col]; indi < Aend[c.col]; ++indi)
          if (Aindex.at(indi) == c.row) break;

        // recover x: column decreases by 1
        // if indi is not Aend-1 swap elements indi and Aend-1
        if (indi != Aend[c.col] - 1) {
          double tmp = Avalue[Aend[c.col] - 1];
          int tmpi = Aindex[Aend[c.col] - 1];
          Avalue[Aend[c.col] - 1] = Avalue.at(indi);
          Aindex[Aend[c.col] - 1] = Aindex.at(indi);
          Avalue.at(indi) = tmp;
          Aindex.at(indi) = tmpi;
        }
        Aend[c.col]--;

        // recover y: column increases by 1
        // update A: append X column to end of array
        int st = Avalue.size();
        for (int ind = Astart[yindex]; ind < Aend[yindex]; ++ind) {
          Avalue.push_back(Avalue.at(ind));
          Aindex.push_back(Aindex.at(ind));
        }
        Avalue.push_back(postValue.top());
        Aindex.push_back(c.row);
        Astart[yindex] = st;
        Aend[yindex] = Avalue.size();

        if (iKKTcheck == 1)
          chk.addChange(173, c.row, c.col, postValue.top(), (double)yindex, 0);
        postValue.pop();

        break;
      }
      case DOUBLETON_EQUATION_NEW_X_ZERO_AR_UPDATE: {
        // sp case x disappears row representation change
        int indi;
        for (indi = ARstart[c.row]; indi < ARstart[c.row + 1]; ++indi)
          if (ARindex.at(indi) == numColOriginal) break;
        ARindex.at(indi) = c.col;
        ARvalue.at(indi) = postValue.top();

        if (iKKTcheck == 1) {
          chk.ARindex.at(indi) = c.col;
          chk.ARvalue.at(indi) = postValue.top();
        }

        postValue.pop();

        break;
      }
      case DOUBLETON_EQUATION_NEW_X_ZERO_A_UPDATE: {
        // sp case x disappears column representation change
        // here A is copied from AR array at end of presolve so need to expand x
        // column  Aend[c.col]++; wouldn't do because old value is overriden
        double oldXvalue = postValue.top();
        postValue.pop();
        int x = c.col;

        // update A: append X column to end of array
        int st = Avalue.size();
        for (int ind = Astart.at(x); ind < Aend.at(x); ++ind) {
          Avalue.push_back(Avalue.at(ind));
          Aindex.push_back(Aindex.at(ind));
        }
        Avalue.push_back(oldXvalue);
        Aindex.push_back(c.row);
        Astart.at(x) = st;
        Aend.at(x) = Avalue.size();

        break;
      }
      case EMPTY_ROW: {
        valueRowDual[c.row] = 0;
        flagRow[c.row] = 1;
        if (iKKTcheck == 1) {
          if (chk.print == 1)
            cout << "----KKT check after empty row " << c.row
                 << " re-introduced-----\n";
          chk.addChange(0, c.row, 0, 0, 0, 0);
          chk.makeKKTCheck();
        }
        break;
      }
      case SING_ROW: {
        // valuePrimal is already set for this one, colDual also, we need
        // rowDual. AR copy keeps full matrix.  col dual maybe infeasible, we
        // need to check.  recover old bounds and see
        getDualsSingletonRow(c.row, c.col);

        if (iKKTcheck == 1) {
          if (chk.print == 1)
            cout << "----KKT check after singleton row " << c.row
                 << " re-introduced. Variable: " << c.col << " -----\n";
          chk.addChange(1, c.row, c.col, valuePrimal[c.col],
                        valueColDual[c.col], valueRowDual[c.row]);
          chk.replaceBasis(col_status, row_status);
          chk.makeKKTCheck();
        }
        break;
      }
      case FORCING_ROW_VARIABLE:
        fRjs.push_back(c.col);
        flagCol[c.col] = 1;
        if (iKKTcheck == 1 && valuePrimal[c.col] != 0)
          chk.addChange(22, c.row, c.col, 0, 0, 0);
        break;
      case FORCING_ROW: {
        string str = getDualsForcingRow(c.row, fRjs);

        if (iKKTcheck == 1) {
          if (chk.print == 1)
            cout << "----KKT check after forcing row " << c.row
                 << " re-introduced. Variable(s): " << str << " -----\n";
          chk.replaceBasis(col_status, row_status);
          chk.addChange(3, c.row, 0, 0, 0, valueRowDual[c.row]);
          chk.makeKKTCheck();
        }
        fRjs.clear();
        break;
      }
      case REDUNDANT_ROW: {
        // this is not zero if the row bounds got relaxed and transferred to a
        // column which then had a nonzero dual.
        valueRowDual[c.row] = 0;

        flagRow[c.row] = 1;

        if (iKKTcheck == 1) {
          if (chk.print == 1)
            cout << "----KKT check after redundant row " << c.row
                 << " re-introduced.----------------\n";
          chk.addChange(0, c.row, 0, 0, 0, 0);
          chk.makeKKTCheck();
        }
        break;
      }
      case FREE_SING_COL:
      case IMPLIED_FREE_SING_COL: {
        // colDual rowDual already set.
        // calculate row value without xj
        double aij = getaij(c.row, c.col);
        double sum = 0;
        for (int k = ARstart[c.row]; k < ARstart[c.row + 1]; ++k)
          if (flagCol.at(ARindex.at(k)))
            sum += valuePrimal.at(ARindex.at(k)) * ARvalue.at(k);

        double rowlb = postValue.top();
        postValue.pop();
        double rowub = postValue.top();
        postValue.pop();

        // calculate xj
        if (valueRowDual[c.row] < 0) {
          // row is at lower bound
          valuePrimal[c.col] = (rowlb - sum) / aij;
        } else if (valueRowDual[c.row] > 0) {
          // row is at upper bound
          valuePrimal[c.col] = (rowub - sum) / aij;
        } else if (rowlb == rowub)
          valuePrimal[c.col] = (rowlb - sum) / aij;
        else if (colCostAtEl[c.col] > 0) {
          // we are interested in the lowest possible value of x:
          // max { l_j, bound implied by row i }
          double bndL;
          if (aij > 0)
            bndL = (rowlb - sum) / aij;
          else
            bndL = (rowub - sum) / aij;
          valuePrimal[c.col] = max(colLowerOriginal[c.col], bndL);
        } else if (colCostAtEl[c.col] < 0) {
          // we are interested in the highest possible value of x:
          // min { u_j, bound implied by row i }
          double bndU;
          if (aij < 0)
            bndU = (rowlb - sum) / aij;
          else
            bndU = (rowub - sum) / aij;
          valuePrimal[c.col] = min(colUpperOriginal[c.col], bndU);
        } else {  // cost is zero
          double bndL, bndU;
          if (aij > 0) {
            bndL = (rowlb - sum) / aij;
            bndU = (rowub - sum) / aij;
          } else {
            bndL = (rowub - sum) / aij;
            bndU = (rowlb - sum) / aij;
          }
          double valuePrimalUB = min(colUpperOriginal[c.col], bndU);
          double valuePrimalLB = max(colLowerOriginal[c.col], bndL);
          if (valuePrimalUB < valuePrimalLB - tol) {
            cout << "Postsolve error: inconsistent bounds for implied free "
                    "column singleton "
                 << c.col << endl;
          }

          if (fabs(valuePrimalLB) < fabs(valuePrimalUB))
            valuePrimal[c.col] = valuePrimalLB;
          else
            valuePrimal[c.col] = valuePrimalUB;
        }
        sum = sum + valuePrimal[c.col] * aij;

        double costAtTimeOfElimination = postValue.top();
        postValue.pop();
        objShift += (costAtTimeOfElimination * sum) / aij;

        flagRow[c.row] = 1;
        flagCol[c.col] = 1;
        // valueRowDual[c.row] = 0;

        if (iKKTcheck == 1) {
          chk.addCost(c.col, costAtTimeOfElimination);
          if (c.type == FREE_SING_COL && chk.print == 1)
            cout << "----KKT check after free col singleton " << c.col
                 << " re-introduced. Row: " << c.row << " -----\n";
          else if (c.type == IMPLIED_FREE_SING_COL && chk.print == 1)
            cout << "----KKT check after implied free col singleton " << c.col
                 << " re-introduced. Row: " << c.row << " -----\n";
          chk.addChange(4, c.row, c.col, valuePrimal[c.col],
                        valueColDual[c.col], valueRowDual[c.row]);
          chk.makeKKTCheck();
        }
        break;
      }
      case SING_COL_DOUBLETON_INEQ: {
        // column singleton in a doubleton equation.
        // colDual already set. need valuePrimal from stack. maybe change
        // rowDual depending on bounds. old bounds kept in oldBounds. variables
        // j,k : we eliminated j and are left with changed bounds on k and no
        // row. c.col is column COL (K) - eliminated, j is with new bounds
        pair<int, vector<double>> p = oldBounds.top();
        oldBounds.pop();
        vector<double> v = get<1>(p);
        int j = get<0>(p);
        // double ubNew = v[1];
        // double lbNew = v[0];
        double cjNew = v[2];
        p = oldBounds.top();
        oldBounds.pop();
        v = get<1>(p);
        double ubOld = v[1];
        double lbOld = v[0];
        double cjOld = v[2];
        p = oldBounds.top();
        oldBounds.pop();
        v = get<1>(p);
        double ubCOL = v[1];
        double lbCOL = v[0];
        double ck = v[2];

        double rowlb = postValue.top();
        postValue.pop();
        double rowub = postValue.top();
        postValue.pop();
        double aik = postValue.top();
        postValue.pop();
        double aij = postValue.top();
        postValue.pop();
        double xj = valuePrimal.at(j);

        // calculate xk, depending on signs of coeff and cost
        double upp = HIGHS_CONST_INF;
        double low = -HIGHS_CONST_INF;

        if ((aij > 0 && aik > 0) || (aij < 0 && aik < 0)) {
          upp = (rowub - aij * xj) / aik;
          low = (rowlb - aij * xj) / aik;
        } else {
          upp = (rowub - aij * xj) / aik;
          low = (rowlb - aij * xj) / aik;
        }

        double xkValue = 0;
        if (ck == 0) {
          if (low < 0 && upp > 0)
            xkValue = 0;
          else if (fabs(low) < fabs(upp))
            xkValue = low;
          else
            xkValue = upp;
        }

        else if ((ck > 0 && aik > 0) || (ck < 0 && aik < 0)) {
          if (low == -HIGHS_CONST_INF)
            cout << "ERROR UNBOUNDED? unnecessary check";
          xkValue = low;
        } else if ((ck > 0 && aik < 0) || (ck < 0 && aik > 0)) {
          if (upp == HIGHS_CONST_INF)
            cout << "ERROR UNBOUNDED? unnecessary check";
          xkValue = upp;
        }

        // primal value and objective shift
        valuePrimal[c.col] = xkValue;
        objShift += -cjNew * xj + cjOld * xj + ck * xkValue;

        // fix duals

        double rowVal = aij * xj + aik * xkValue;
        if (rowub - rowVal > tol && rowVal - rowlb > tol) {
          row_status.at(c.row) = HighsBasisStatus::BASIC;
          col_status.at(c.col) = HighsBasisStatus::NONBASIC;
          valueRowDual[c.row] = 0;
          flagRow[c.row] = 1;
          valueColDual[c.col] = getColumnDualPost(c.col);
        } else {
          double lo, up;
          if (fabs(rowlb - rowub) < tol) {
            lo = -HIGHS_CONST_INF;
            up = HIGHS_CONST_INF;
          } else if (fabs(rowub - rowVal) <= tol) {
            lo = 0;
            up = HIGHS_CONST_INF;
          } else if (fabs(rowlb - rowVal) <= tol) {
            lo = -HIGHS_CONST_INF;
            ;
            up = 0;
          }

          colCostAtEl.at(j) = cjOld;  // revert cost before calculating duals
          getBoundOnLByZj(c.row, j, &lo, &up, lbOld, ubOld);
          getBoundOnLByZj(c.row, c.col, &lo, &up, lbCOL, ubCOL);

          // calculate yi
          if (lo - up > tol)
            cout << "PR: Error in postsolving doubleton inequality " << c.row
                 << " : inconsistent bounds for its dual value.\n";

          // WARNING: bound_row_dual not used. commented out to surpress warning
          // but maybe this causes trouble. Look into when you do dual postsolve
          // again (todo)
          //
          //
          // double bound_row_dual = 0;
          // if (lo > 0) {
          //   bound_row_dual = lo;
          // } else if (up < 0) {
          //   bound_row_dual = up;
          // }

          if (lo > 0 || up < 0) {
            // row is nonbasic, since dual value zero for it is infeasible.
            row_status.at(c.row) = HighsBasisStatus::NONBASIC;
            col_status.at(c.col) = HighsBasisStatus::BASIC;
            valueColDual[c.col] = 0;
            flagRow[c.row] = 1;
            valueRowDual[c.row] = getRowDualPost(c.row, c.col);
            valueColDual[j] = getColumnDualPost(j);
          } else {
            // zero row dual is feasible, set row to basic and column to
            // nonbasic.
            row_status.at(c.row) = HighsBasisStatus::BASIC;
            col_status.at(c.col) = HighsBasisStatus::NONBASIC;
            valueRowDual[c.row] = 0;
            flagRow[c.row] = 1;
            valueColDual[c.col] = getColumnDualPost(c.col);
          }

          if (iKKTcheck == 1) chk.colDual.at(j) = valueColDual.at(j);
        }

        flagCol[c.col] = 1;

        if (iKKTcheck == 1) {
          if (chk.print == 1)
            cout << "----KKT check after col singleton " << c.col
                 << " in doubleton eq re-introduced. Row: " << c.row
                 << " -----\n";

          chk.addChange(5, c.row, c.col, valuePrimal[c.col],
                        valueColDual[c.col], valueRowDual[c.row]);
          chk.replaceBasis(col_status, row_status);
          chk.makeKKTCheck();
        }
        // exit(2);
        break;
      }
      case EMPTY_COL:
      case DOMINATED_COLS:
      case WEAKLY_DOMINATED_COLS: {
        // got valuePrimal, need colDual
        if (c.type != EMPTY_COL) {
          z = colCostAtEl[c.col];
          for (int k = Astart[c.col]; k < Astart[c.col + 1]; ++k)
            if (flagRow.at(Aindex.at(k)))
              z = z + valueRowDual.at(Aindex.at(k)) * Avalue.at(k);
          valueColDual[c.col] = z;
        }

        flagCol[c.col] = 1;
        if (iKKTcheck == 1) {
          if (c.type == EMPTY_COL && chk.print == 1)
            cout << "----KKT check after empty column " << c.col
                 << " re-introduced.-----------\n";
          else if (c.type == DOMINATED_COLS && chk.print == 1)
            cout << "----KKT check after dominated column " << c.col
                 << " re-introduced.-----------\n";
          else if (c.type == WEAKLY_DOMINATED_COLS && chk.print == 1)
            cout << "----KKT check after weakly dominated column " << c.col
                 << " re-introduced.-----------\n";

          chk.addChange(6, 0, c.col, valuePrimal[c.col], valueColDual[c.col],
                        0);
          chk.makeKKTCheck();
        }
        break;
      }

      case FIXED_COL: {
        // got valuePrimal, need colDual
        valueColDual[c.col] = getColumnDualPost(c.col);

        flagCol[c.col] = 1;
        if (iKKTcheck == 1) {
          if (chk.print == 1)
            cout << "----KKT check after fixed variable " << c.col
                 << " re-introduced.-----------\n";
          chk.addChange(7, 0, c.col, valuePrimal[c.col], valueColDual[c.col],
                        0);
          chk.makeKKTCheck();
        }
        break;
      }
    }
    // cmpNBF(c.row, c.col);
  }

  // cmpNBF();

  // Check number of basic variables
  int num_basic_var = 0;
  for (int iCol = 0; iCol < numColOriginal; iCol++) {
    if (col_status[iCol] == HighsBasisStatus::BASIC) {
      assert(num_basic_var < numRowOriginal);
      if (num_basic_var == numRowOriginal) {
        printf("Error in postsolve: more basic variables than rows\n");
        break;
      }
      num_basic_var++;
    }
  }
  for (int iRow = 0; iRow < numRowOriginal; iRow++) {
    // int iVar = numColOriginal + iRow;
    if (row_status[iRow] == HighsBasisStatus::BASIC) {
      assert(num_basic_var < numRowOriginal);
      if (num_basic_var == numRowOriginal) {
        printf("Error from postsolve: more basic variables than rows\n");
        break;
      }
      num_basic_var++;
    }
  }
  // Return error if the number of basic variables does not equal the
  // number of rows in the original LP
  assert(num_basic_var == numRowOriginal);
  if (num_basic_var != numRowOriginal) {
    printf(
        "Error from postsolve: number of basic variables = %d != %d = number "
        "of rows\n",
        num_basic_var, numRowOriginal);
    return HighsPostsolveStatus::BasisError;
  }

  // cout<<"Singularity check at end of postsolve: ";
  // testBasisMatrixSingularity();

  if (iKKTcheck == 2) {
    if (chk.print == 3) chk.print = 2;
    chk.passSolution(valuePrimal, valueColDual, valueRowDual);
    chk.makeKKTCheck();
  }

  // now recover original model data to pass back to HiGHS
  // A is already recovered!
  // however, A is expressed in terms of Astart, Aend and columns are in
  // different order so
  makeACopy();

  numRow = numRowOriginal;
  numCol = numColOriginal;
  numTot = numRow + numCol;

  rowUpper = rowUpperOriginal;
  rowLower = rowLowerOriginal;

  colUpper = colUpperOriginal;
  colLower = colLowerOriginal;

  colCost = colCostOriginal;

  /*
  nonbasicMove.resize(numTot, 0);
  for (int i = 0; i < numColOriginal; ++i) {
    if (colLower.at(i) != colUpper.at(i) && colLower.at(i) != -HIGHS_CONST_INF)
      nonbasicMove.at(i) = 1;
    else if (colUpper.at(i) != HIGHS_CONST_INF)
      nonbasicMove.at(i) = -1;
    else
      nonbasicMove.at(i) = 0;
  }
  */
  colValue = valuePrimal;
  colDual = valueColDual;
  rowDual = valueRowDual;
  rowValue.assign(numRow, 0);
  for (int i = 0; i < numRowOriginal; ++i) {
    for (int k = ARstart.at(i); k < ARstart.at(i + 1); ++k)
      rowValue.at(i) += valuePrimal.at(ARindex.at(k)) * ARvalue.at(k);
  }
  // JAJH(120519) Added following four lines so that recovered solution is
  // returned
  recovered_solution.col_value = colValue;
  recovered_solution.col_dual = colDual;
  recovered_solution.row_value = rowValue;
  recovered_solution.row_dual = rowDual;
  return HighsPostsolveStatus::SolutionRecovered;
}

void Presolve::setBasisElement(change c) {
  // col_status starts off as [numCol] and has already been increased to
  // [numColOriginal] and row_status starts off as [numRow] and has already been
  // increased to [numRowOriginal] so fill fill in gaps in both

  switch (c.type) {
    case EMPTY_ROW: {
      if (report_postsolve) {
        printf("2.1 : Recover row %3d as %3d (basic): empty row\n", c.row,
               numColOriginal + c.row);
      }
      row_status.at(c.row) = HighsBasisStatus::BASIC;
      break;
    }
    case SING_ROW:
    case FORCING_ROW_VARIABLE:
    case FORCING_ROW:
    case SING_COL_DOUBLETON_INEQ:
      break;
    case REDUNDANT_ROW: {
      if (report_postsolve) {
        printf("2.3 : Recover row %3d as %3d (basic): redundant\n", c.row,
               numColOriginal + c.row);
      }
      row_status.at(c.row) = HighsBasisStatus::BASIC;
      break;
    }
    case FREE_SING_COL:
    case IMPLIED_FREE_SING_COL: {
      if (report_postsolve) {
        printf(
            "2.4a: Recover col %3d as %3d (basic): implied free singleton "
            "column\n",
            c.col, numColOriginal + c.row);
      }
      col_status.at(c.col) = HighsBasisStatus::BASIC;

      if (report_postsolve) {
        printf(
            "2.5b: Recover row %3d as %3d (nonbasic): implied free singleton "
            "column\n",
            c.row, numColOriginal + c.row);
      }
      row_status.at(c.row) = HighsBasisStatus::NONBASIC;  // Was LOWER
      break;
    }
    case EMPTY_COL:
    case DOMINATED_COLS:
    case WEAKLY_DOMINATED_COLS: {
      if (report_postsolve) {
        printf("2.7 : Recover column %3d (nonbasic): weakly dominated column\n",
               c.col);
      }
      col_status.at(c.col) = HighsBasisStatus::NONBASIC;  // Was LOWER
      break;
    }
    case FIXED_COL: {  // fixed variable:
      // check if it was NOT after singRow
      if (chng.size() > 0)
        if (chng.top().type != SING_ROW) {
          if (report_postsolve) {
            printf(
                "2.8 : Recover column %3d (nonbasic): weakly dominated "
                "column\n",
                c.col);
          }
          col_status.at(c.col) = HighsBasisStatus::NONBASIC;  // Was LOWER
        }
      break;
    }
  }
}

/* testing and dev
int Presolve::testBasisMatrixSingularity() {

        HFactor factor;

        //resize matrix in M so we can pass to factor
        int i, j, k;
        int nz = 0;
        int nR = 0;
        int nC = 0;

        numRowOriginal = rowLowerOriginal.size();
        numColOriginal = colLowerOriginal.size();
        //arrays to keep track of indices
        vector<int> rIndex_(numRowOriginal, -1);
        vector<int> cIndex_(numColOriginal, -1);

        for (i=0;i<numRowOriginal;++i)
                if (flagRow.at(i)) {
                        for (j = ARstart.at(i); j<ARstart.at(i+1); ++j)
                                if (flagCol[ARindex.at(j)])
                                        nz ++;
                        rIndex_.at(i) = nR;
                        nR++;
                        }

        for (i=0;i<numColOriginal;++i)
                if (flagCol.at(i)) {
                        cIndex_.at(i) = nC;
                        nC++;
                }


        //matrix
        vector<int>    Mstart(nC + 1, 0);
        vector<int>    Mindex(nz);
        vector<double> Mvalue(nz);

    vector<int> iwork(nC, 0);

    for (i = 0;i<numRowOriginal; ++i)
        if (flagRow.at(i))
            for (int k = ARstart.at(i); k < ARstart.at(i+1);++k ) {
                j = ARindex.at(k);
                if (flagCol.at(j))
                                iwork[cIndex_.at(j)]++;
                        }
    for (i = 1; i <= nC; ++i)
        Mstart.at(i) = Mstart[i - 1] + iwork[i - 1];
   for (i = 0; i < numColOriginal; ++i)
        iwork.at(i) = Mstart.at(i);

   for (i = 0; i < numRowOriginal; ++i) {
        if (flagRow.at(i)) {
                        int iRow = rIndex_.at(i);
                    for (k = ARstart.at(i); k < ARstart[i + 1];++k ) {
                        j = ARindex.at(k);
                        if (flagCol.at(j)) {
                                int iCol = cIndex_.at(j);
                                    int iPut = iwork[iCol]++;
                                    Mindex[iPut] = iRow;
                                    Mvalue[iPut] = ARvalue.at(k);
                                }
                    }
                }
    }

    vector<int>  bindex(nR);
    int countBasic=0;

    printf("To recover this test need to use col/row_status\n");
     for (int i=0; i< nonbasicFlag.size();++i) {
         if (nonbasicFlag.at(i) == 0)
                         countBasic++;
     }

     if (countBasic != nR)
         cout<<" Wrong count of basic variables: != numRow"<<endl;

     int c=0;
     for (int i=0; i< nonbasicFlag.size();++i) {
         if (nonbasicFlag.at(i) == 0) {
                        if (i < numColOriginal)
                                bindex[c] = cIndex_.at(i);
                        else
                                bindex[c] = nC + rIndex_[i - numColOriginal];
                        c++;
         }
    }

        factor.setup(nC, nR, &Mstart[0], &Mindex[0], &Mvalue[0],  &bindex[0]);
/ *	if (1) // for this check both A and M are the full matrix again
        {
                if (nC - numColOriginal != 0)
                        cout<<"columns\n";
                if (nR - numRowOriginal != 0)
                        cout<<"rows\n";
                for (int i=0; i< Mstart.size();++i)
                        if (Mstart.at(i) - Astart.at(i) != 0)
                                cout<<"Mstart "<<i<<"\n";
                for (int i=0; i< Mindex.size();++i)
                        if (Mindex.at(i) - Aindex.at(i) != 0)
                                cout<<"Mindex "<<i<<"\n";
                for (int i=0; i< Mvalue.size();++i)
                        if (Mvalue.at(i) - Avalue.at(i) != 0)
                                cout<<"Mvalue "<<i<<"\n";
                for (int i=0; i< bindex.size();++i)
                        if (nonbasicFlag.at(i) - nbffull.at(i) != 0)
                                cout<<"nbf "<<i<<"\n";
        } * /

        try {
        factor.build();
    } catch (runtime_error& error) {
        cout << error.what() << endl;
        cout << "Postsolve: could not factorize basis matrix." << endl;
        return 0;
    }
    cout << "Postsolve: basis matrix successfully factorized." << endl;

    return 1;
}*/

/***
 * lo and up refer to the place storing the current bounds on y_row
 *
 */
void Presolve::getBoundOnLByZj(int row, int j, double* lo, double* up,
                               double colLow, double colUpp) {
  double cost = colCostAtEl.at(j);  // valueColDual.at(j);
  double x = -cost;

  double sum = 0;
  for (int kk = Astart.at(j); kk < Aend.at(j); ++kk)
    if (flagRow.at(Aindex.at(kk))) {
      sum = sum + Avalue.at(kk) * valueRowDual.at(Aindex.at(kk));
    }
  x = x - sum;

  double aij = getaij(row, j);
  x = x / aij;

  if (fabs(colLow - colUpp) < tol)
    return;  // here there is no restriction on zj so no bound on y

  if ((valuePrimal.at(j) - colLow) > tol &&
      (colUpp - valuePrimal.at(j)) > tol) {
    // set both bounds
    if (x < *up) *up = x;
    if (x > *lo) *lo = x;
  }

  else if ((valuePrimal.at(j) == colLow && aij < 0) ||
           (valuePrimal.at(j) == colUpp && aij > 0)) {
    if (x < *up) *up = x;
  } else if ((valuePrimal.at(j) == colLow && aij > 0) ||
             (valuePrimal.at(j) == colUpp && aij < 0)) {
    if (x > *lo) *lo = x;
  }
}

/**
 * returns z_col
 * z = A'y + c
 */
double Presolve::getColumnDualPost(int col) {
  int row;
  double z;
  double sum = 0;
  for (int cnt = Astart.at(col); cnt < Aend.at(col); cnt++)
    if (flagRow.at(Aindex.at(cnt))) {
      row = Aindex.at(cnt);
      sum = sum + valueRowDual.at(row) * Avalue.at(cnt);
    }
  z = sum + colCostAtEl.at(col);
  return z;
}

/***
 * A'y + c = z
 *
 * returns y_row = -(A'y      +   c   - z )/a_rowcol
 *               (except row)  (at el)
 */
double Presolve::getRowDualPost(int row, int col) {
  double x = 0;

  for (int kk = Astart.at(col); kk < Aend.at(col); ++kk)
    if (flagRow.at(Aindex.at(kk)) && Aindex.at(kk) != row)
      x = x + Avalue.at(kk) * valueRowDual.at(Aindex.at(kk));

  x = x + colCostAtEl.at(col) - valueColDual.at(col);

  double y = getaij(row, col);
  return -x / y;
}

string Presolve::getDualsForcingRow(int row, vector<int>& fRjs) {
  double z;
  stringstream ss;
  int j;

  double lo = -HIGHS_CONST_INF;
  double up = HIGHS_CONST_INF;
  int lo_col = -1;
  int up_col = -1;

  double cost, sum;

  for (size_t jj = 0; jj < fRjs.size(); ++jj) {
    j = fRjs[jj];

    pair<int, vector<double>> p = oldBounds.top();
    vector<double> v = get<1>(p);
    oldBounds.pop();
    double colLow = v[0];
    double colUpp = v[1];

    // calculate bound x imposed by zj
    double save_lo = lo;
    double save_up = up;
    getBoundOnLByZj(row, j, &lo, &up, colLow, colUpp);
    if (lo > save_lo) lo_col = j;
    if (up < save_up) up_col = j;
  }

  // calculate yi
  if (lo > up)
    cout << "PR: Error in postsolving forcing row " << row
         << " : inconsistent bounds for its dual value.\n";

  if (lo <= 0 && up >= 0) {
    valueRowDual.at(row) = 0;
    row_status[row] = HighsBasisStatus::BASIC;
  } else if (lo > 0) {
    // row is set to basic and column to non-basic but that should change
    row_status[row] = HighsBasisStatus::NONBASIC;
    col_status.at(lo_col) = HighsBasisStatus::BASIC;
    valueRowDual.at(row) = lo;
    valueColDual.at(lo_col) = 0;
    // valueColDual[lo_col] should be zero since it imposed the lower bound.
  } else if (up < 0) {
    // row is set to basic and column to non-basic but that should change
    row_status[row] = HighsBasisStatus::NONBASIC;
    col_status.at(up_col) = HighsBasisStatus::BASIC;
    valueRowDual.at(row) = up;
    valueColDual.at(up_col) = 0;
  }

  flagRow.at(row) = 1;

  for (size_t jj = 0; jj < fRjs.size(); ++jj) {
    j = fRjs[jj];
    if (lo > 0 && j == lo_col) continue;
    if (up < 0 && j == up_col) continue;

    col_status[j] = HighsBasisStatus::NONBASIC;

    cost = valueColDual.at(j);
    sum = 0;
    for (int k = Astart.at(j); k < Aend.at(j); ++k)
      if (flagRow.at(Aindex.at(k))) {
        sum = sum + valueRowDual.at(Aindex.at(k)) * Avalue.at(k);
        // cout<<" row "<<Aindex.at(k)<<" dual
        // "<<valueRowDual.at(Aindex.at(k))<<" a_"<<Aindex.at(k)<<"_"<<j<<"\n";
      }
    z = cost + sum;

    valueColDual.at(j) = z;

    if (iKKTcheck == 1) {
      ss << j;
      ss << " ";
      chk.addChange(2, 0, j, valuePrimal.at(j), valueColDual.at(j), cost);
    }
  }

  return ss.str();
}

void Presolve::getDualsSingletonRow(int row, int col) {
  pair<int, vector<double>> bnd = oldBounds.top();
  oldBounds.pop();

  valueRowDual.at(row) = 0;
  //   double cost = postValue.top();
  postValue.pop();
  double aij = getaij(row, col);
  double l = (get<1>(bnd))[0];
  double u = (get<1>(bnd))[1];
  double lrow = (get<1>(bnd))[2];
  double urow = (get<1>(bnd))[3];

  flagRow.at(row) = 1;

  HighsBasisStatus local_status;
  local_status = col_status.at(col);
  if (local_status != HighsBasisStatus::BASIC) {
    // x was not basic but is now
    // if x is strictly between original bounds or a_ij*x_j is at a bound.
    if (fabs(valuePrimal.at(col) - l) > tol &&
        fabs(valuePrimal.at(col) - u) > tol) {
      if (report_postsolve) {
        printf("3.1 : Make column %3d basic and row %3d nonbasic\n", col, row);
      }
      col_status.at(col) = HighsBasisStatus::BASIC;
      row_status.at(row) = HighsBasisStatus::NONBASIC;  // Was LOWER
      valueColDual[col] = 0;
      valueRowDual[row] = getRowDualPost(row, col);
    } else {
      // column is at bound
      bool isRowAtLB = fabs(aij * valuePrimal[col] - lrow) < tol;
      bool isRowAtUB = fabs(aij * valuePrimal[col] - urow) < tol;

      double save_dual = valueColDual[col];
      valueColDual[col] = 0;
      double row_dual = getRowDualPost(row, col);

      if ((isRowAtLB && !isRowAtUB && row_dual > 0) ||
          (!isRowAtLB && isRowAtUB && row_dual < 0) ||
          (!isRowAtLB && !isRowAtUB)) {
        // make row basic
        row_status.at(row) = HighsBasisStatus::BASIC;
        valueRowDual[row] = 0;
        valueColDual[col] = save_dual;
      } else {
        // column is basic
        col_status.at(col) = HighsBasisStatus::BASIC;
        row_status.at(row) = HighsBasisStatus::NONBASIC;
        valueColDual[col] = 0;
        valueRowDual[row] = getRowDualPost(row, col);
      }
    }
  } else {
    // x is basic
    if (report_postsolve) {
      printf("3.3 : Make row %3d basic\n", row);
    }
    row_status.at(row) = HighsBasisStatus::BASIC;
    valueRowDual[row] = 0;
    // if the row dual is zero it does not contribute to the column dual.
  }

  if (iKKTcheck == 1) {
    chk.colDual.at(col) = valueColDual.at(col);
    chk.rowDual.at(row) = valueRowDual.at(row);
  }
}

void Presolve::getDualsDoubletonEquation(int row, int col) {
  // colDual already set. need valuePrimal from stack. maybe change rowDual
  // depending on bounds. old bounds kept in oldBounds. variables j,k : we
  // eliminated col(k)(c.col) and are left with changed bounds on j and no row.
  //                               y x

  pair<int, vector<double>> p = oldBounds.top();
  oldBounds.pop();
  vector<double> v = get<1>(p);
  int x = get<0>(p);
  double ubxNew = v[1];
  double lbxNew = v[0];
  double cxNew = v[2];
  p = oldBounds.top();
  oldBounds.pop();
  v = get<1>(p);
  double ubxOld = v[1];
  double lbxOld = v[0];
  double cxOld = v[2];
  p = oldBounds.top();
  oldBounds.pop();
  v = get<1>(p);
  double uby = v[1];
  double lby = v[0];
  double cy = v[2];

  int y = col;

  double b = postValue.top();
  postValue.pop();
  double aky = postValue.top();
  postValue.pop();
  double akx = postValue.top();
  postValue.pop();
  double valueX = valuePrimal.at(x);

  // primal value and objective shift
  valuePrimal.at(y) = (b - akx * valueX) / aky;
  objShift += -cxNew * valueX + cxOld * valueX + cy * valuePrimal.at(y);

  // column cost of x
  colCostAtEl.at(x) = cxOld;

  flagRow.at(row) = 1;

  HighsBasisStatus local_status;
  if (x < numColOriginal) {
    local_status = col_status.at(x);
  } else {
    local_status = row_status.at(x - numColOriginal);
  }
  if ((local_status != HighsBasisStatus::BASIC && valueX == ubxNew &&
       ubxNew < ubxOld) ||
      (local_status != HighsBasisStatus::BASIC && valueX == lbxNew &&
       lbxNew > lbxOld)) {
    if (x < numColOriginal) {
      col_status.at(x) = HighsBasisStatus::BASIC;
      // transfer dual of x to dual of row
      valueColDual.at(x) = 0;
      valueRowDual.at(row) = getRowDualPost(row, x);
      valueColDual.at(y) = getColumnDualPost(y);

      if (report_postsolve) printf("4.1 : Make column %3d basic\n", x);
    } else {
      row_status.at(x - numColOriginal) = HighsBasisStatus::BASIC;
      if (report_postsolve)
        printf("4.1 : Make row    %3d basic\n", x - numColOriginal);

      valueRowDual.at(row) = 0;
      valueColDual.at(x) = getColumnDualPost(x);
      valueColDual.at(y) = getColumnDualPost(y);
    }
  } else {
    // row becomes basic unless y is between bounds, in which case y is basic
    if (valuePrimal.at(y) - lby > tol && uby - valuePrimal.at(y) > tol) {
      if (y < numColOriginal) {
        col_status.at(y) = HighsBasisStatus::BASIC;
        if (report_postsolve) printf("4.2 : Make column %3d basic\n", y);

        valueColDual.at(y) = 0;
        valueRowDual.at(row) = getRowDualPost(row, y);
      } else {
        row_status.at(y - numColOriginal) = HighsBasisStatus::BASIC;
        if (report_postsolve)
          printf("4.2 : Make row    %3d basic\n", y - numColOriginal);

        valueRowDual.at(row) = 0;
        valueColDual.at(x) = getColumnDualPost(x);
        valueColDual.at(y) = getColumnDualPost(y);
      }
    } else if (fabs(valueX - ubxNew) < tol || fabs(valueX - lbxNew) < tol) {
      if (y < numColOriginal) {
        col_status.at(y) = HighsBasisStatus::BASIC;
        if (report_postsolve) printf("4.3 : Make column %3d basic\n", y);

        valueColDual.at(y) = 0;
        valueRowDual.at(row) = getRowDualPost(row, y);
      } else {
        row_status.at(y - numColOriginal) = HighsBasisStatus::BASIC;
        if (report_postsolve)
          printf("4.3 : Make row    %3d basic\n", y - numColOriginal);

        valueRowDual.at(row) = 0;
        valueColDual.at(x) = getColumnDualPost(x);
        valueColDual.at(y) = getColumnDualPost(y);
      }
    } else {
      if (report_postsolve) {
        printf("4.4 : Make row    %3d basic\n", row);
      }
      row_status.at(row) = HighsBasisStatus::BASIC;
      valueRowDual.at(row) = 0;
      valueColDual.at(x) = getColumnDualPost(x);
      valueColDual.at(y) = getColumnDualPost(y);
    }
  }
  if (iKKTcheck == 1) {
    chk.colDual.at(x) = valueColDual.at(x);
    chk.colDual.at(y) = valueColDual.at(y);
    chk.rowDual.at(row) = valueRowDual.at(row);
  }

  flagCol.at(y) = 1;
}

void Presolve::countRemovedRows(PresolveRule rule) {
  timer.increaseCount(true, rule);
}

void Presolve::countRemovedCols(PresolveRule rule) {
  timer.increaseCount(false, rule);
}
