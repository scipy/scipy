#ifndef PRESOLVE_PRESOLVE_UTILS_H_
#define PRESOLVE_PRESOLVE_UTILS_H_

#include <vector>

namespace presolve {

void printRowWise(
    const int numRow, const int numCol, const std::vector<double>& colCost,
    const std::vector<double>& colLower, const std::vector<double>& colUpper,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<int>& ARstart, const std::vector<int>& ARindex,
    const std::vector<double>& ARvalue);

void printRow(const int row, const int numRow, const int numCol,
              const std::vector<int>& flagRow, const std::vector<int>& flagCol,
              const std::vector<double>& rowLower,
              const std::vector<double>& rowUpper,
              const std::vector<double>& values,
              const std::vector<int>& ARstart, const std::vector<int>& ARindex,
              const std::vector<double>& ARvalue);

void printCol(const int col, const int numRow, const int numCol,
              const std::vector<int>& flagRow, const std::vector<int>& flagCol,
              const std::vector<double>& colLower,
              const std::vector<double>& colUpper,
              const std::vector<double>& values, const std::vector<int>& Astart,
              const std::vector<int>& Aend, const std::vector<int>& Aindex,
              const std::vector<double>& Avalue);

void printCol(const int col, const int numRow, const int numCol,
              const std::vector<int>& flagRow, const std::vector<int>& flagCol,
              const std::vector<double>& colLower,
              const std::vector<double>& colUpper,
              const std::vector<double>& values, const std::vector<int>& Astart,
              const std::vector<int>& Aindex,
              const std::vector<double>& Avalue);

void printRowOneLine(
    const int row, const int numRow, const int numCol,
    const std::vector<int>& flagRow, const std::vector<int>& flagCol,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<double>& values, const std::vector<int>& ARstart,
    const std::vector<int>& ARindex, const std::vector<double>& ARvalue);

void printAR(const int numRow, const int numCol,
             const std::vector<double>& colCost,
             const std::vector<double>& rowLower,
             const std::vector<double>& rowUpper,
             const std::vector<int>& ARstart, const std::vector<int>& ARindex,
             std::vector<double>& ARvalue);

void printA(const int numRow, const int numCol,
            const std::vector<double>& colCost,
            const std::vector<double>& rowLower,
            const std::vector<double>& rowUpper,
            const std::vector<double>& colLower,
            const std::vector<double>& colUpper, const std::vector<int>& Astart,
            const std::vector<int>& Aindex, std::vector<double>& Avalue);

}  // namespace presolve

#endif