/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file test/KktCheck.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "test/KktCheck.h"

#include <cassert>
#include <vector>

void KktCheck::printAR() {
  std::cout << "N=" << numCol << ",  M=" << numRow
            << ",  NZ= " << ARstart[numRow] << '\n';

  std::cout << "\n-----cost-----\n";
  for (size_t i = 0; i < colCost.size(); i++) {
    std::cout << colCost[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "------AR | b----KktCheck-\n";
  for (i = 0; i < numRow; i++) {
    for (j = 0; j < numCol; j++) {
      int ind = ARstart[i];
      while (ARindex[ind] != j && ind < ARstart[i + 1]) ind++;
      // if a_ij is nonzero print
      if (ARindex[ind] == j && ind < ARstart[i + 1]) {
        std::cout << ARvalue[ind] << " ";
      } else
        std::cout << "   ";
    }
    std::cout << "  |   " << rowLower[i] << " < < " << rowUpper[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "------l------\n";
  for (int i = 0; i < numCol; i++) {
    if (colLower[i] > -HIGHS_CONST_INF)
      std::cout << colLower[i] << " ";
    else
      std::cout << "-inf ";
  }
  std::cout << std::endl;
  std::cout << "------u------\n";
  for (int i = 0; i < numCol; i++) {
    if (colUpper[i] < HIGHS_CONST_INF)
      std::cout << colUpper[i] << " ";
    else
      std::cout << "inf ";
  }
  std::cout << std::endl;
}

void KktCheck::makeARCopy() {
  tol = 0.00001;
  // Make a AR copy
  std::vector<int> iwork(numRow, 0);
  ARstart.resize(numRow + 1, 0);
  int AcountX = Aindex.size();
  ARindex.resize(AcountX);
  ARvalue.resize(AcountX);
  for (k = 0; k < AcountX; k++) iwork[Aindex[k]]++;
  for (i = 1; i <= numRow; i++) ARstart[i] = ARstart[i - 1] + iwork[i - 1];
  for (i = 0; i < numRow; i++) iwork[i] = ARstart[i];
  for (int iCol = 0; iCol < numCol; iCol++) {
    for (k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      int iRow = Aindex[k];
      int iPut = iwork[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = Avalue[k];
    }
  }
}

void KktCheck::chPrimalBounds() {
  for (i = 0; i < numCol; i++) {
    if ((colLower[i] - colValue[i] > tol) ||
        (colValue[i] - colUpper[i] > tol)) {
      if (print == 1)
        std::cout << "Variable " << cIndexRev[i]
                  << " infeasible: lb=" << colLower[i]
                  << ", vaule=" << colValue[i] << ",  ub=" << colUpper[i]
                  << std::endl;
      // std::cout<<"Variable "<<i<<" infeasible: lb="<<colLower[i]<<",
      // vaule="<<colValue[i]<<",  ub="<<colUpper[i]<<std::endl;
      istrueGlb = true;
    }
  }
}

void KktCheck::chPrimalFeas() {
  bool istrue = true;
  double rowV;
  // Ax = b
  for (i = 0; i < numRow; i++) {
    rowV = 0;
    for (k = ARstart[i]; k < ARstart[i + 1]; k++)
      rowV = rowV + colValue[ARindex[k]] * ARvalue[k];

    if (((rowV - rowLower[i]) < 0) && (fabs(rowV - rowLower[i]) > tol)) {
      if (print == 1)
        std::cout << "Row " << rIndexRev[i] << " infeasible: Row value=" << rowV
                  << "  L=" << rowLower[i] << "  U=" << rowUpper[i]
                  << std::endl;
      // std::cout<<"Row "<<i<<" infeasible: Row value="<<rowV<<"
      // L="<<rowLower[i]<<"  U="<<rowUpper[i]<<std::endl;
      istrue = false;
    }

    if (((rowV - rowUpper[i]) > 0) && (fabs(rowV - rowUpper[i]) > tol)) {
      if (print == 1)
        std::cout << "Row " << rIndexRev[i] << " infeasible: Row value=" << rowV
                  << "  L=" << rowLower[i] << "  U=" << rowUpper[i]
                  << std::endl;
      // std::cout<<"Row "<<i<<" infeasible: Row value="<<rowV<<"
      // L="<<rowLower[i]<<"  U="<<rowUpper[i]<<std::endl;
      istrue = false;
    }
  }

  if (istrue) {
    if (print == 1) std::cout << "Primal feasible.\n";
  } else {
    if (print == 1) std::cout << "KKT check error: Primal infeasible.\n";
    istrueGlb = true;
  }
}

void KktCheck::chDualFeas() {
  bool istrue = true;

  // check values of z_j are dual feasible
  for (i = 0; i < numCol; i++) {
    // j not in L or U
    if (colLower[i] == -HIGHS_CONST_INF && colUpper[i] == HIGHS_CONST_INF) {
      if (fabs(colDual[i]) > tol) {
        if (print == 1)
          std::cout << "Dual feasibility fail: l=-inf, x[" << cIndexRev[i]
                    << "]=" << colValue[i] << ", u=inf, z[" << i
                    << "]=" << colDual[i] << std::endl;
        // std::cout<<"Dual feasibility fail: l=-inf,
        // x["<<i<<"]="<<colValue[i]<<", u=inf,
        // z["<<i<<"]="<<colDual[i]<<std::endl;
        istrue = false;
      }
    }
    // j in L: x=l and l<u
    else if (colValue[i] == colLower[i] && colLower[i] < colUpper[i]) {
      if (colDual[i] < 0 && fabs(colDual[i]) > tol) {
        if (print == 1)
          std::cout << "Dual feasibility fail: l[" << cIndexRev[i]
                    << "]=" << colLower[i] << " = x[" << cIndexRev[i]
                    << "]=" << colValue[i] << ", z[" << cIndexRev[i]
                    << "]=" << colDual[i] << std::endl;
        // std::cout<<"Dual feasibility fail: l["<<i<<"]="<<colLower[i]<<" =
        // x["<<i<<"]="<<colValue[i]<<", z["<<i<<"]="<<colDual[i]<<std::endl;
        istrue = false;
      }
    }
    // j in U: x=u and l<u
    else if (colValue[i] == colUpper[i] && colLower[i] < colUpper[i]) {
      if (colDual[i] > tol) {
        if (print == 1)
          std::cout << "Dual feasibility fail: x[" << cIndexRev[i]
                    << "]=" << colValue[i] << "=u[" << cIndexRev[i] << "], z["
                    << cIndexRev[i] << "]=" << colDual[i] << std::endl;
        // std::cout<<"Dual feasibility fail:
        // x["<<i<<"]="<<colValue[i]<<"=u["<<i<<"],
        // z["<<i<<"]="<<colDual[i]<<std::endl;
        istrue = false;
      }
    }
  }

  // check values of y_i are dual feasible
  for (i = 0; i < numRow; i++) {
    double rowV = 0;
    for (k = ARstart[i]; k < ARstart[i + 1]; k++)
      rowV = rowV + colValue[ARindex[k]] * ARvalue[k];

    // L = Ax = U can be any sign
    if (fabs(rowLower[i] - rowV) < tol && fabs(rowUpper[i] - rowV) < tol) {
    }
    // L = Ax < U
    else if (fabs(rowLower[i] - rowV) < tol && rowV < rowUpper[i]) {
      if (rowDual[i] > tol) {
        if (print == 1)
          std::cout << "Dual feasibility fail for row " << rIndexRev[i]
                    << ": L= " << rowLower[i] << ", Ax=" << rowV
                    << ", U=" << rowUpper[i] << ", y=" << rowDual[i]
                    << std::endl;
        // std::cout<<"Dual feasibility fail for row "<<i<<": L= "<<rowLower[i]
        // <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<std::endl;
        istrue = false;
      }
    }
    // L < Ax = U
    else if (rowLower[i] < rowV && fabs(rowV - rowUpper[i]) < tol) {
      // std::cout<<"Dual feasibility fail for row "<<i<<": L= "<<rowLower[i]
      // <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<std::endl;
      if (rowDual[i] < -tol) {
        if (print == 1)
          std::cout << "Dual feasibility fail for row " << rIndexRev[i]
                    << ": L= " << rowLower[i] << ", Ax=" << rowV
                    << ", U=" << rowUpper[i] << ", y=" << rowDual[i]
                    << std::endl;
        // std::cout<<"Dual feasibility fail for row "<<i<<": L= "<<rowLower[i]
        // <<", Ax="<<rowV<<", U="<<rowUpper[i]<<", y="<<rowDual[i]<<std::endl;
        istrue = false;
      }
    }
    // L < Ax < U
    else if ((rowLower[i] < (rowV + tol)) && (rowV < (rowUpper[i] + tol))) {
      if (fabs(rowDual[i]) > tol) {
        if (print == 1)
          std::cout << "Dual feasibility fail for row " << rIndexRev[i]
                    << ": L= " << rowLower[i] << ", Ax=" << rowV
                    << ", U=" << rowUpper[i] << ", y=" << rowDual[i]
                    << std::endl;
        // std::cout<<"Dual feasibility fail for row "<<i<<": L= "<<rowLower[i]
        // <<", Ax="<<rowV<<", U="<<rowUpper[i]<<",
        // y="<<rowDual[i]<<std::endl;istrue = false;
        istrue = false;
      }
    }
  }

  if (istrue) {
    if (print == 1) std::cout << "Dual feasible.\n";
  } else {
    if (print == 1) std::cout << "KKT check error: Dual feasibility fail.\n";
    istrueGlb = true;
  }
}

void KktCheck::chComplementarySlackness() {
  bool istrue = true;

  for (i = 0; i < numCol; i++) {
    if (colLower[i] > -HIGHS_CONST_INF)
      if (fabs((colValue[i] - colLower[i]) * (colDual[i])) > tol &&
          colValue[i] != colUpper[i] && fabs(colDual[i]) > tol) {
        if (print == 1)
          std::cout << "Comp. slackness fail: "
                    << "l[" << cIndexRev[i] << "]=" << colLower[i] << ", x["
                    << i << "]=" << colValue[i] << ", z[" << i
                    << "]=" << colDual[i] << std::endl;
        // std::cout<<"Comp. slackness fail: "<<"l["<<i<<"]="<<colLower[i]<<",
        // x["<<i<<"]="<<colValue[i]<<", z["<<i<<"]="<<colDual[i]<<std::endl;
        istrue = false;
      }
    if (colUpper[i] < HIGHS_CONST_INF)
      if (fabs((colUpper[i] - colValue[i]) * (colDual[i])) > tol &&
          colValue[i] != colLower[i] && fabs(colDual[i]) > tol) {
        if (print == 1)
          std::cout << "Comp. slackness fail: x[" << cIndexRev[i]
                    << "]=" << colValue[i] << ", u[" << i << "]=" << colUpper[i]
                    << ", z[" << i << "]=" << colDual[i] << std::endl;
        // std::cout<<"Comp. slackness fail: x["<<i<<"]="<<colValue[i]<<",
        // u["<<i<<"]="<<colUpper[i]<<", z["<<i<<"]="<<colDual[i]<<std::endl;
        istrue = false;
      }
  }

  if (istrue) {
    if (print == 1) std::cout << "Complementary Slackness.\n";
  } else {
    if (print == 1) std::cout << "KKT check error: Comp slackness fail.\n";
    istrueGlb = true;
  }
}

void KktCheck::printSol() {
  char buff[10];
  std::cout << std::endl << "Col value: ";
  for (size_t i = 0; i < colValue.size(); i++) {
    sprintf(buff, "%2.2f ", colValue[i]);
    std::cout << std::setw(5) << buff;
  }
  std::cout << std::endl << "Col dual:  ";
  for (size_t i = 0; i < colDual.size(); i++) {
    sprintf(buff, "%2.2f ", colDual[i]);
    std::cout << std::setw(5) << buff;
  }
  /*	cout<<std::endl<<"Row value: ";
          for (i=0;i<numRow;i++) {
                  sprintf(buff, "%2.2f ", rowValue[i]);
                  std::cout<<setw(5)<<buff;
                  }*/
  std::cout << std::endl << "Row dual:  ";
  for (size_t i = 0; i < rowDual.size(); i++) {
    sprintf(buff, "%2.2f ", rowDual[i]);
    std::cout << std::setw(5) << buff;
  }
  std::cout << std::endl << std::endl;
}

void KktCheck::chStOfLagrangian() {
  bool istrue = true;
  double lagrV;
  // A'y + c - z = 0
  for (j = 0; j < numCol; j++) {
    lagrV = colCost[j] - colDual[j];
    for (k = Astart[j]; k < Astart[j + 1]; k++)
      lagrV = lagrV + rowDual[Aindex[k]] * Avalue[k];

    if (fabs(lagrV) > tol) {
      if (print == 1)
        std::cout << "Column " << cIndexRev[j]
                  << " fails stationary of Lagrangian: dL/dx" << j << " = "
                  << lagrV << ", rather than zero." << std::endl;
      // std::cout<<"Column "<<j<<" fails stationary of Lagrangian: dL/dx"<<j<<"
      // =
      // "<<lagrV<<", rather than zero."<<std::endl;
      istrue = false;
    }
  }

  if (istrue) {
    if (print == 1) std::cout << "Stationarity of Lagrangian.\n";
  } else {
    if (print == 1)
      std::cout << "KKT check error: Lagrangian is not stationary.\n";
    istrueGlb = true;
  }
}

void KktCheck::checkBFS() {
  // Go over cols and check that the duals of basic values are zero.
  assert((int)col_status.size() == numCol);
  assert((int)colDual.size() == numCol);
  for (int j = 0; j < numCol; j++) {
    if (col_status[j] == HighsBasisStatus::BASIC && colDual[j] != 0) {
      if (print == 1)
        std::cout << "Col " << cIndexRev[j] << " is basic but has nonzero dual."
                  << std::endl;
    }
  }

  // Go over rows and check that the duals of basic values are zero.
  assert((int)row_status.size() == numRow);
  assert((int)rowDual.size() == numRow);
  for (int i = 0; i < numRow; i++) {
    if (row_status[i] == HighsBasisStatus::BASIC && rowDual[i] != 0) {
      if (print == 1)
        std::cout << "Row " << rIndexRev[i] << " is basic but has nonzero dual."
                  << std::endl;
    }
  }
}

void KktCheck::checkKKT() {
  if (numCol == 0) return;

  istrueGlb = false;

  makeARCopy();
  // printAR();printSol();
  chPrimalBounds();
  chPrimalFeas();
  chDualFeas();
  chComplementarySlackness();
  chStOfLagrangian();

  checkBFS();

  // if (print == 2) {
  //   std::ofstream myfile;
  //   myfile.open("../experiments/out", std::ios::app);
  //   if (istrueGlb)
  //     myfile << "           KKT fail      ";
  //   else
  //     myfile << "           KKT pass      ";
  //   myfile.close();
  // }
}

void KktCheck::passSolution(const std::vector<double>& colVal,
                            const std::vector<double>& colDu,
                            const std::vector<double>& rDu) {
  colValue = colVal;
  colDual = colDu;
  rowDual = rDu;
}
// get DATA
void KktCheck::setMatrix(const std::vector<int>& Astart_,
                         const std::vector<int>& Aindex_,
                         const std::vector<double>& Avalue_) {
  Astart = Astart_;
  Aindex = Aindex_;
  Avalue = Avalue_;
}

void KktCheck::setBounds(const std::vector<double>& colUpper_,
                         const std::vector<double>& colLower_) {
  colLower = colLower_;
  colUpper = colUpper_;
}

void KktCheck::setNumbersCostRHS(int nCol, int nRow,
                                 const std::vector<double>& rowLower_,
                                 const std::vector<double>& rowUpper_,
                                 const std::vector<double>& cost) {
  numCol = nCol;
  numRow = nRow;
  colCost = cost;
  rowLower = rowLower_;
  rowUpper = rowUpper_;
}

void KktCheck::setIndexVectors(std::vector<int>& rIndex,
                               std::vector<int>& cIndex) {
  rIndexRev.clear();
  cIndexRev.clear();

  for (size_t i = 0; i < rIndex.size(); i++) {
    if (rIndex[i] != -1) {
      rIndexRev.push_back(i);
    }
  }
  for (size_t i = 0; i < cIndex.size(); i++) {
    if (cIndex[i] != -1) {
      cIndexRev.push_back(i);
    }
  }
}
