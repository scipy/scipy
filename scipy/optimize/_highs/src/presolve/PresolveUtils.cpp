#include "presolve/PresolveUtils.h"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include "lp_data/HConst.h"

namespace presolve {

using std::setw;

void printRowOneLine(
    const int row, const int numRow, const int numCol,
    const std::vector<int>& flagRow, const std::vector<int>& flagCol,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<double>& values, const std::vector<int>& ARstart,
    const std::vector<int>& ARindex, const std::vector<double>& ARvalue) {
  assert(row >= 0 && row < numRow);

  // go over row and sum
  // col
  // if flagCol[]
  // a_ij * value_j

  double sum = 0.0;
  for (int k = ARstart[row]; k < ARstart[row + 1]; k++) {
    const int col = ARindex[k];
    assert(col >= 0 && col <= numCol);
    sum += ARvalue[k] * values[col];
  }

  std::cout << "row " << row << ": " << flagRow[row] << "   " << rowLower[row]
            << " <= " << sum << " <= " << rowUpper[row] << std::endl;
}

void printRow(const int row, const int numRow, const int numCol,
              const std::vector<int>& flagRow, const std::vector<int>& flagCol,
              const std::vector<double>& rowLower,
              const std::vector<double>& rowUpper,
              const std::vector<double>& values,
              const std::vector<int>& ARstart, const std::vector<int>& ARindex,
              const std::vector<double>& ARvalue) {
  assert(row >= 0 && row < numRow);

  std::cout << "row " << row << ": " << flagRow[row] << "   " << rowLower[row]
            << " <= ... <= " << rowUpper[row] << std::endl
            << "..." << std::endl;
  // go over row and print
  // col
  // flagCol[] ..... next col in row
  // a_ij
  // x_j
  for (int k = ARstart[row]; k < ARstart[row + 1]; k++) {
    const int col = ARindex[k];
    assert(col >= 0 && col <= numCol);
    (void)col;
  }

  for (int k = ARstart[row]; k < ARstart[row + 1]; k++)
    std::cout << setw(3) << ARindex[k] << " ";

  std::cout << std::endl;

  for (int k = ARstart[row]; k < ARstart[row + 1]; k++)
    std::cout << setw(3) << flagCol[ARindex[k]] << " ";

  std::cout << std::endl;
  for (int k = ARstart[row]; k < ARstart[row + 1]; k++)
    std::cout << setw(3) << ARvalue[k] << " ";

  std::cout << std::endl;
  for (int k = ARstart[row]; k < ARstart[row + 1]; k++)
    std::cout << setw(3) << values[ARindex[k]] << " ";

  std::cout << std::endl;
}

void printCol(const int col, const int numRow, const int numCol,
              const std::vector<int>& flagRow, const std::vector<int>& flagCol,
              const std::vector<double>& colLower,
              const std::vector<double>& colUpper,
              const std::vector<double>& row_values,
              const std::vector<int>& Astart, const std::vector<int>& Aend,
              const std::vector<int>& Aindex,
              const std::vector<double>& Avalue) {
  assert(col >= 0 && col < numCol);

  std::cout << "col" << col << ": " << flagCol[col] << "   " << colLower[col]
            << " <= ... <= " << colUpper[col] << std::endl
            << "..." << std::endl;

  // go over col and print
  // row flagRow[] a_ij x_j
  // ...
  // next row in column

  for (int k = Astart[col]; k < Aend[col]; k++) {
    const int row = Aindex[k];
    assert(row >= 0 && row <= numRow);
    std::cout << setw(3) << row << " ";
    std::cout << setw(3) << flagRow[row] << " ";
    std::cout << setw(3) << Avalue[k] << " ";
    std::cout << setw(3) << row_values[row] << " ";
    std::cout << std::endl;
  }

  std::cout << std::endl;
}

void printRowWise(
    const int numRow, const int numCol, const std::vector<double>& colCost,
    const std::vector<double>& colLower, const std::vector<double>& colUpper,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<int>& ARstart, const std::vector<int>& ARindex,
    const std::vector<double>& ARvalue) {
  const int rows = numRow;
  const int cols = numCol;

  std::cout << "\n-----cost-----\n";

  for (unsigned int i = 0; i < colCost.size(); i++) {
    std::cout << colCost[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "------AR-|-L-U-----\n";
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      int ind = ARstart[i];
      while (ARindex[ind] != j && ind < ARstart[i + 1]) ind++;
      // if a_ij is nonzero print
      if (ARindex[ind] == j && ind < ARstart[i + 1])
        std::cout << ARvalue[ind];
      else
        std::cout << "   ";
    }
    std::cout << "  |   " << rowLower[i] << " < < " << rowUpper[i] << std::endl;
  }
  std::cout << "------l------\n";
  for (int i = 0; i < cols; i++) {
    if (colLower[i] > -HIGHS_CONST_INF)
      std::cout << colLower[i] << " ";
    else
      std::cout << "-inf";
  }
  std::cout << std::endl;
  std::cout << "------u------\n";
  for (int i = 0; i < cols; i++) {
    if (colUpper[i] < HIGHS_CONST_INF)
      std::cout << colUpper[i] << " ";
    else
      std::cout << "inf ";
  }
  std::cout << std::endl;
}

void printA(const int numRow, const int numCol,
            const std::vector<double>& colCost,
            const std::vector<double>& rowLower,
            const std::vector<double>& rowUpper,
            const std::vector<double>& colLower,
            const std::vector<double>& colUpper, const std::vector<int>& Astart,
            const std::vector<int>& Aindex, std::vector<double>& Avalue) {
  char buff[7];
  std::cout << "\n-----cost-----\n";

  for (int i = 0; i < numCol; i++) {
    std::cout << colCost[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "------A-|-b-----\n";
  for (int i = 0; i < numRow; i++) {
    for (int j = 0; j < numCol; j++) {
      int ind = Astart[j];
      while (Aindex[ind] != i && ind < Astart[j + 1]) ind++;
      // if a_ij is nonzero print
      if (Aindex[ind] == i && ind < Astart[j + 1]) {
        std::cout << Avalue[ind] << " ";
      } else
        std::cout << " ";
    }
    std::cout << "  |   " << rowLower[i] << " < < " << rowUpper[i] << std::endl;
  }
  std::cout << "------l------\n";
  for (int i = 0; i < numCol; i++) {
    if (colLower[i] > -HIGHS_CONST_INF)
      std::cout << colLower[i] << " ";
    else
      std::cout << "-inf ";
    std::cout << setw(9) << buff;
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

void printAR(const int numRow, const int numCol,
             const std::vector<double>& colCost,
             const std::vector<double>& rowLower,
             const std::vector<double>& rowUpper,
             const std::vector<int>& ARstart, const std::vector<int>& ARindex,
             std::vector<double>& ARvalue) {
  std::cout << "\n-----cost-----\n";

  for (int i = 0; i < numCol; i++) {
    std::cout << colCost[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "------AR-|-b-----\n";
  for (int i = 0; i < numRow; i++) {
    for (int j = 0; j < numCol; j++) {
      int ind = ARstart[i];
      while (ARindex[ind] != j && ind < ARstart[i + 1]) ind++;
      // if a_ij is nonzero print
      if (ARindex[ind] == j && ind < ARstart[i + 1]) {
        std::cout << ARvalue[ind] << " ";
      } else
        std::cout << " ";
    }
    std::cout << "  |   " << rowLower[i] << " < < " << rowUpper[i] << std::endl;
  }

  std::cout << std::endl;
}

}  // namespace presolve