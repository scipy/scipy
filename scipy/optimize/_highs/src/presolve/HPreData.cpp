/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPreData.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "presolve/HPreData.h"

using std::cout;
using std::endl;
using std::setw;

namespace presolve {

HPreData::HPreData() {}

double HPreData::getRowValue(int i) {
  double sum = 0;
  for (int k = ARstart[i]; k < ARstart[i + 1]; k++)
    if (flagRow[ARindex[k]]) sum += ARvalue[k] * valuePrimal[ARindex[k]];
  return sum;
}

double HPreData::getaij(int i, int j) {
  int k = ARstart[i];
  while (j != ARindex[k] && k <= ARstart[i + 1]) k++;
  if (k == ARstart[i + 1]) {
    // cout<<"Error: No such element in A: row "<<i<<", column "<<j<<endl;
    // exit(0);
  }
  return ARvalue[k];
}

bool HPreData::isZeroA(int i, int j) {
  int k = ARstart[i];
  while (k < ARstart[i + 1] && j != ARindex[k]) k++;
  if (k == ARstart[i + 1]) {
    return true;
  }
  return false;
}

void HPreData::makeARCopy() {
  // Make a AR copy
  vector<int> iwork(numRow, 0);
  ARstart.resize(numRow + 1, 0);
  int AcountX = Aindex.size();
  ARindex.resize(AcountX);
  ARvalue.resize(AcountX);
  for (int k = 0; k < AcountX; k++) iwork.at(Aindex.at(k))++;
  for (int i = 1; i <= numRow; i++)
    ARstart.at(i) = ARstart.at(i - 1) + iwork.at(i - 1);
  for (int i = 0; i < numRow; i++) iwork.at(i) = ARstart.at(i);
  for (int iCol = 0; iCol < numCol; iCol++) {
    for (int k = Astart.at(iCol); k < Astart.at(iCol + 1); k++) {
      int iRow = Aindex.at(k);
      int iPut = iwork.at(iRow)++;
      ARindex.at(iPut) = iCol;
      ARvalue.at(iPut) = Avalue[k];
    }
  }
}

void HPreData::makeACopy() {
  // Make a A copy

  vector<int> iwork(numColOriginal, 0);
  Astart.assign(numColOriginal + 1, 0);
  const int AcountX = ARindex.size();
  Aindex.resize(AcountX);
  Avalue.resize(AcountX);
  for (int k = 0; k < AcountX; k++)
    if (ARindex[k] < numColOriginal) iwork[ARindex[k]]++;
  for (int i = 1; i <= numColOriginal; i++)
    Astart[i] = Astart[i - 1] + iwork[i - 1];
  for (int i = 0; i < numColOriginal; i++) iwork[i] = Astart[i];
  for (int iRow = 0; iRow < numRowOriginal; iRow++) {
    for (int k = ARstart[iRow]; k < ARstart[iRow + 1]; k++) {
      const int iColumn = ARindex[k];
      if (iColumn != numColOriginal) {
        int iPut = iwork[iColumn]++;
        Aindex[iPut] = iRow;
        Avalue[iPut] = ARvalue[k];
      }
    }
  }

  Aend.resize(numColOriginal + 1, 0);
  for (int i = 0; i < numColOriginal; i++) Aend[i] = Astart[i + 1];
}

void initPresolve(PresolveStats& stats) {
  std::cout << "Init Presolve form HiGHS" << std::endl;
}

}  // namespace presolve