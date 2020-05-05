/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file test/KktChStep.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "test/KktChStep.h"

#include <utility>

#include "test/KktCheck.h"

using std::cout;
using std::endl;
using std::get;
using std::pair;
using std::setw;
using std::vector;

void KktChStep::replaceBasis(const vector<HighsBasisStatus>& columns,
                             const vector<HighsBasisStatus>& rows) {
  col_status = columns;
  row_status = rows;
}

void KktChStep::passBasis(const vector<HighsBasisStatus>& columns,
                          const vector<HighsBasisStatus>& rows) {
  col_status.resize(RnumCol);
  row_status.resize(RnumRow);

  // This can be done much faster but it's just for checking and
  // re-used the code from passSolution.

  numRow = rows.size();
  numCol = columns.size();

  // arrays to keep track of indices
  vector<int> rIndex(RnumRow, -1);
  vector<int> cIndex(RnumCol, -1);
  int nR = 0;
  int nC = 0;

  for (int i = 0; i < RnumRow; i++)
    if (flagRow[i]) {
      rIndex[i] = nR;
      nR++;
    }

  for (int i = 0; i < RnumCol; i++)
    if (flagCol[i]) {
      cIndex[i] = nC;
      nC++;
    }

  // set corresponding parts of solution vectors:
  int j = 0;
  vector<int> eqIndexOfReduced(numCol, -1);
  vector<int> eqIndexOfReduROW(numRow, -1);
  for (int i = 0; i < RnumCol; i++)
    if (cIndex[i] > -1) {
      eqIndexOfReduced[j] = i;
      j++;
    }
  j = 0;
  for (int i = 0; i < RnumRow; i++)
    if (rIndex[i] > -1) {
      eqIndexOfReduROW[j] = i;
      j++;
    }

  for (int i = 0; i < numCol; i++) col_status[eqIndexOfReduced[i]] = columns[i];

  for (int i = 0; i < numRow; i++) row_status[eqIndexOfReduROW[i]] = rows[i];
}

void KktChStep::passSolution(const vector<double>& colVal,
                             const vector<double>& colDu,
                             const vector<double>& rDu) {
  colValue.resize(RnumCol);
  colDual.resize(RnumCol);
  rowDual.resize(RnumRow);

  numRow = rDu.size();
  numCol = colDu.size();
  // arrays to keep track of indices
  vector<int> rIndex(RnumRow, -1);
  vector<int> cIndex(RnumCol, -1);
  int nR = 0;
  int nC = 0;

  for (int i = 0; i < RnumRow; i++)
    if (flagRow[i]) {
      rIndex[i] = nR;
      nR++;
    }

  for (int i = 0; i < RnumCol; i++)
    if (flagCol[i]) {
      cIndex[i] = nC;
      nC++;
    }

  // set corresponding parts of solution vectors:
  int j = 0;
  vector<int> eqIndexOfReduced(numCol, -1);
  vector<int> eqIndexOfReduROW(numRow, -1);
  for (int i = 0; i < RnumCol; i++)
    if (cIndex[i] > -1) {
      eqIndexOfReduced[j] = i;
      j++;
    }
  j = 0;
  for (int i = 0; i < RnumRow; i++)
    if (rIndex[i] > -1) {
      eqIndexOfReduROW[j] = i;
      j++;
    }

  for (int i = 0; i < numCol; i++) {
    colValue[eqIndexOfReduced[i]] = colVal[i];
    colDual[eqIndexOfReduced[i]] = colDu[i];
  }

  for (int i = 0; i < numRow; i++) rowDual[eqIndexOfReduROW[i]] = rDu[i];
}

// full matrix
void KktChStep::setMatrixAR(int nCol, int nRow, const vector<int>& ARstart_,
                            const vector<int>& ARindex_,
                            const vector<double>& ARvalue_) {
  RnumCol = nCol;
  RnumRow = nRow;
  ARstart = ARstart_;
  ARindex = ARindex_;
  ARvalue = ARvalue_;
}

void KktChStep::setBoundsCostRHS(const vector<double>& colUpper_,
                                 const vector<double>& colLower_,
                                 const vector<double>& cost,
                                 const vector<double>& rowLower_,
                                 const vector<double>& rowUpper_) {
  RcolLower = colLower_;
  RcolUpper = colUpper_;
  RrowLower = rowLower_;
  RrowUpper = rowUpper_;
  RcolCost = cost;
}

void KktChStep::addCost(int col, double value) { RcolCost[col] = value; }

/*
1  SING_ROW
17 DOUBLETON_EQUATION

171 DOUBLETON_EQUATION_ROW_BOUNDS_UPDATE
173 DOUBLETON_EQUATION_X_ZERO_INITIALLY
172 DOUBLETON_EQUATION_NEW_X_NONZERO
174 DOUBLETON_EQUATION_NEW_X_ZERO_AR_UPDATE
175 DOUBLETON_EQUATION_NEW_X_ZERO_A_UPDATE

7 FIXED_COL
0 EMPTY_ROW
6 EMPTY_COL
9 DOMINATED_COLS
10 WEAKLY_DOMINATED_COLS
4 FREE_SING_COL
5 SING_COL_DOUBLETON_INEQ
19 SING_COL_DOUBLETON_INEQ_SECOND_SING_COL
8 IMPLIED_FREE_SING_COL
3 FORCING_ROW
2 FORCING_ROW_VARIABLE
16 REDUNDANT_ROW

*/
void KktChStep::addChange(int type, int row, int col, double valC, double dualC,
                          double dualR) {
  // when updating fill new values for b, c, bounds in Rb RcolCost RcolUpper,
  // RcolLower
  vector<pair<int, double>> upd;

  switch (type) {
    case 171:  // new bounds from doubleton equation, retrieve old ones
      upd = rLowers.top();
      rLowers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RrowLower[ind] = get<1>(upd[i]);
      }
      upd = rUppers.top();
      rUppers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RrowUpper[ind] = get<1>(upd[i]);
      }
      break;
    case 172: {  // matrix transformation from doubleton equation, retrieve old
                 // value
      // here col is x, the value of which we are changing back.
      // valC is the old value of the x entry in the matrix
      int i = ARstart[row];
      while (i < ARstart[row + 1] && ARindex[i] != col) i++;
      ARvalue[i] = valC;
      break;
    }
    case 173: {  // matrix transformation from doubleton equation, retrieve old
                 // value
      // here col is x, the value of which we are changing back.
      // valC is the old value of the x entry in the matrix
      int i = ARstart[row];
      while (i < ARstart[row + 1] && ARindex[i] != col) i++;
      ARvalue[i] = valC;
      ARindex[i] = (int)dualC;
      break;
    }
    case 0:  // empty row
      flagRow[row] = 1;
      break;
    case 1:  // row singleton
      flagRow[row] = 1;
      colValue[col] = valC;
      colDual[col] = dualC;
      rowDual[row] = dualR;
      if (valC != 0) {
        upd = cLowers.top();
        cLowers.pop();
        for (size_t i = 0; i < upd.size(); i++) {
          int ind = get<0>(upd[i]);
          RcolLower[ind] = get<1>(upd[i]);
        }
        upd = cUppers.top();
        cUppers.pop();
        for (size_t i = 0; i < upd.size(); i++) {
          int ind = get<0>(upd[i]);
          RcolUpper[ind] = get<1>(upd[i]);
        }
      }
      break;
    case 2:  // each variable at forcing row: rowDual is cost here
      colValue[col] = valC;
      colDual[col] = dualC;
      flagCol[col] = 1;
      RcolCost[col] = dualR;
      break;
    case 21:  //
      rowDual[row] = dualR;
      break;
    case 22:  //
      upd = rLowers.top();
      rLowers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RrowLower[ind] = get<1>(upd[i]);
      }
      upd = rUppers.top();
      rUppers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RrowUpper[ind] = get<1>(upd[i]);
      }
      break;
    case 3:  // the row that is forcing
      rowDual[row] = dualR;
      flagRow[row] = 1;
      if (valC != 0) {
        upd = rLowers.top();
        rLowers.pop();
        for (size_t i = 0; i < upd.size(); i++) {
          int ind = get<0>(upd[i]);
          RrowLower[ind] = get<1>(upd[i]);
        }
        upd = rUppers.top();
        rUppers.pop();
        for (size_t i = 0; i < upd.size(); i++) {
          int ind = get<0>(upd[i]);
          RrowUpper[ind] = get<1>(upd[i]);
        }
      }
      break;
    case 4:  // implied free column singleton (also from duplicate row)
      flagRow[row] = 1;
      flagCol[col] = 1;
      colValue[col] = valC;
      colDual[col] = dualC;
      rowDual[row] = dualR;
      upd = costs.top();
      costs.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolCost[ind] = get<1>(upd[i]);
      }
      break;
    case 5:  // doubleton eq with singleton col
      flagRow[row] = 1;
      flagCol[col] = 1;
      colValue[col] = valC;
      colDual[col] = dualC;
      rowDual[row] = dualR;
      upd = cLowers.top();
      cLowers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolLower[ind] = get<1>(upd[i]);
      }
      upd = cUppers.top();
      cUppers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolUpper[ind] = get<1>(upd[i]);
      }
      upd = costs.top();
      costs.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolCost[ind] = get<1>(upd[i]);
      }
      break;
    case 17: {  // doubleton equation
      flagRow[row] = 1;
      flagCol[col] = 1;
      colValue[col] = valC;
      colDual[col] = dualC;
      rowDual[row] = dualR;
      upd = cLowers.top();
      cLowers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolLower[ind] = get<1>(upd[i]);
      }
      upd = cUppers.top();
      cUppers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolUpper[ind] = get<1>(upd[i]);
      }
      upd = costs.top();
      costs.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolCost[ind] = get<1>(upd[i]);
      }
      break;
    }
    case 6:  // empty column, dominated column or weakly dominated
      colValue[col] = valC;
      colDual[col] = dualC;
      flagCol[col] = 1;
      if (valC != 0) {
        upd = rLowers.top();
        rLowers.pop();
        for (size_t i = 0; i < upd.size(); i++) {
          int ind = get<0>(upd[i]);
          RrowLower[ind] = get<1>(upd[i]);
        }
        upd = rUppers.top();
        rUppers.pop();
        for (size_t i = 0; i < upd.size(); i++) {
          int ind = get<0>(upd[i]);
          RrowUpper[ind] = get<1>(upd[i]);
        }
      }
      break;
    case 7:  // fixed variable
      colValue[col] = valC;
      colDual[col] = dualC;
      flagCol[col] = 1;
      if (valC != 0) {
        upd = rLowers.top();
        rLowers.pop();
        for (size_t i = 0; i < upd.size(); i++) {
          int ind = get<0>(upd[i]);
          RrowLower[ind] = get<1>(upd[i]);
        }
        upd = rUppers.top();
        rUppers.pop();
        for (size_t i = 0; i < upd.size(); i++) {
          int ind = get<0>(upd[i]);
          RrowUpper[ind] = get<1>(upd[i]);
        }
      }
      break;
    case 11:  // empty row from duplucate rows
      upd = rLowers.top();
      rLowers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RrowLower[ind] = get<1>(upd[i]);
      }
      upd = rUppers.top();
      rUppers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RrowUpper[ind] = get<1>(upd[i]);
      }
      flagRow[row] = 1;
      rowDual[row] = dualR;
      break;
    case 12:  // doubleton eq from dupliocate rows;
      flagRow[row] = 1;
      flagCol[col] = 1;
      colValue[col] = valC;
      colDual[col] = dualC;
      rowDual[row] = dualR;
      upd = cLowers.top();
      cLowers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolLower[ind] = get<1>(upd[i]);
      }
      upd = cUppers.top();
      cUppers.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolUpper[ind] = get<1>(upd[i]);
      }
      upd = costs.top();
      costs.pop();
      for (size_t i = 0; i < upd.size(); i++) {
        int ind = get<0>(upd[i]);
        RcolCost[ind] = get<1>(upd[i]);
      }
      break;
    case 121:  //
      colDual[col] = dualC;
      break; /*
case 14: //two duplicate columns by one
       colValue[col] = valC;
       colDual[col] = dualC;
       RcolLower = cLowers.top(); cLowers.pop();
       RcolUpper = cUppers.top(); cUppers.pop();
       break;
case 15: //sing variable on initEq
       flagRow[row] = true;
       rowDual[row] = dualR;
       break;*/
  }
}

void KktChStep::setFlags(vector<int>& r, vector<int>& c) {
  flagRow = r;
  flagCol = c;
}

void KktChStep::resizeProblemMatrix(KktCheck& checker) {
  // full matrix data in AR copy, A copy made to be passed to checker
  int i, j, k;
  int nz = 0;
  int nR = 0;
  int nC = 0;

  // nzRow is the current # on nonzeros, not the general one!
  vector<int> nzRow(RnumRow, 0);
  for (i = 0; i < RnumRow; i++)
    if (flagRow[i])
      for (k = ARstart[i]; k < ARstart[i + 1]; k++)
        if (flagCol[ARindex[k]]) nzRow[i]++;

  numRow = RnumRow;
  numCol = RnumCol;
  // arrays to keep track of indices
  vector<int> rIndex(numRow, -1);
  vector<int> cIndex(numCol, -1);

  for (i = 0; i < numRow; i++)
    if (flagRow[i]) {
      nz += nzRow[i];
      rIndex[i] = nR;
      nR++;
    }

  for (i = 0; i < numCol; i++)
    if (flagCol[i]) {
      cIndex[i] = nC;
      nC++;
    }

  // set indices for checker
  checker.setIndexVectors(rIndex, cIndex);

  // counts
  int numRowEq = numRow;
  int numColEq = numCol;
  numRow = nR;
  numCol = nC;

  // matrix
  vector<int> iwork(numCol, 0);
  Astart.resize(numCol + 1, 0);
  Aindex.resize(nz);
  Avalue.resize(nz);

  for (i = 0; i < numRowEq; i++)
    if (flagRow[i])
      for (int k = ARstart[i]; k < ARstart[i + 1]; k++) {
        j = ARindex[k];
        if (flagCol[j]) iwork[cIndex[j]]++;
      }
  for (i = 1; i <= numCol; i++) Astart[i] = Astart[i - 1] + iwork[i - 1];
  for (i = 0; i < numCol; i++) iwork[i] = Astart[i];
  for (i = 0; i < numRowEq; i++) {
    if (flagRow[i]) {
      int iRow = rIndex[i];
      for (k = ARstart[i]; k < ARstart[i + 1]; k++) {
        j = ARindex[k];
        if (flagCol[j]) {
          int iCol = cIndex[j];
          int iPut = iwork[iCol]++;
          Aindex[iPut] = iRow;
          Avalue[iPut] = ARvalue[k];
        }
      }
    }
  }

  // resize cost, RHS, bounds
  colCost.resize(numCol);
  colLower.resize(numCol);
  colUpper.resize(numCol);

  k = 0;
  for (i = 0; i < numColEq; i++)
    if (flagCol[i]) {
      colCost[k] = RcolCost[i];
      colLower[k] = RcolLower[i];
      colUpper[k] = RcolUpper[i];
      k++;
    }

  rowLower.resize(numRow);
  rowUpper.resize(numRow);

  k = 0;
  for (i = 0; i < numRowEq; i++)
    if (flagRow[i]) {
      rowUpper[k] = RrowUpper[i];
      rowLower[k] = RrowLower[i];

      // b[k] = Rb[i];
      k++;
    }
}

void KktChStep::makeKKTCheck() {
  KktCheck checker;
  resizeProblemMatrix(checker);
  // KktCheck checker;
  checker.setMatrix(Astart, Aindex, Avalue);
  checker.setBounds(colUpper, colLower);
  checker.setNumbersCostRHS(numCol, numRow, rowLower, rowUpper, colCost);
  if (print) checker.print = print;

  // resize and pass solutions
  vector<double> cV;
  vector<double> cD;
  vector<double> rD;

  std::vector<HighsBasisStatus> col_stat;
  std::vector<HighsBasisStatus> row_stat;

  cV.resize(numCol);
  cD.resize(numCol);
  col_stat.resize(numCol);

  int k = 0;
  for (int i = 0; i < RnumCol; i++)
    if (flagCol[i]) {
      cV[k] = colValue[i];
      cD[k] = colDual[i];
      col_stat[k] = col_status[i];
      k++;
    }

  rD.resize(numRow);
  row_stat.resize(numRow);

  k = 0;
  for (int i = 0; i < RnumRow; i++) {
    if (flagRow[i]) {
      rD[k] = rowDual[i];
      row_stat[k] = row_status[i];
      k++;
    }
  }

  checker.passSolution(cV, cD, rD);
  checker.col_status = col_stat;
  checker.row_status = row_stat;
  checker.checkKKT();

  int count = 0;
  for (int i = 0; i < (int)col_status.size(); i++)
    if (col_status[i] == HighsBasisStatus::BASIC) count++;
  for (int i = 0; i < (int)row_status.size(); i++)
    if (row_status[i] == HighsBasisStatus::BASIC) count++;

  if (count != numRow)
    std::cout << "Wrong number of basic variables: numRow = " << numRow
              << ", count = " << count << std::endl;
}

void KktChStep::printA() {
  char buff[7];
  cout << "\n-----cost-----\n";

  for (int i = 0; i < numCol; i++) {
    cout << colCost[i] << " ";
  }
  cout << endl;
  cout << "------A-|-b-----\n";
  for (int i = 0; i < numRow; i++) {
    for (int j = 0; j < numCol; j++) {
      int ind = Astart[j];
      while (Aindex[ind] != i && ind < Astart[j + 1]) ind++;
      // if a_ij is nonzero print
      if (Aindex[ind] == i && ind < Astart[j + 1]) {
        cout << Avalue[ind] << " ";
      } else
        cout << " ";
    }
    cout << "  |   " << RrowLower[i] << " < < " << RrowUpper[i] << endl;
  }
  cout << "------l------\n";
  for (int i = 0; i < numCol; i++) {
    if (colLower[i] > -HIGHS_CONST_INF)
      cout << colLower[i] << " ";
    else
      cout << "-inf ";
    cout << setw(9) << buff;
  }
  cout << endl;
  cout << "------u------\n";
  for (int i = 0; i < numCol; i++) {
    if (colUpper[i] < HIGHS_CONST_INF)
      cout << colUpper[i] << " ";
    else
      cout << "inf ";
  }
  cout << endl;
}

void KktChStep::printAR() {
  cout << "\n-----cost-----\n";

  for (int i = 0; i < numCol; i++) {
    cout << colCost[i] << " ";
  }
  cout << endl;
  cout << "------AR-|-b-----\n";
  for (int i = 0; i < RnumRow; i++) {
    for (int j = 0; j < RnumCol; j++) {
      int ind = ARstart[i];
      while (ARindex[ind] != j && ind < ARstart[i + 1]) ind++;
      // if a_ij is nonzero print
      if (ARindex[ind] == j && ind < ARstart[i + 1]) {
        cout << ARvalue[ind] << " ";
      } else
        cout << " ";
    }
    cout << "  |   " << RrowLower[i] << " < < " << RrowUpper[i] << endl;
  }

  cout << endl;
}
