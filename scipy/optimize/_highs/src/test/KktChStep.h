/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file test/KktChStep.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef TEST_KKTCHSTEP_H_
#define TEST_KKTCHSTEP_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stack>
#include <string>
#include <vector>

#include "lp_data/HConst.h"

class KktCheck;

class KktChStep {
 public:
  KktChStep() { print = 0; }

  // model: full matrix in AR (row-wise) and working copy(column-wise)
  std::vector<int> ARstart;
  std::vector<int> ARindex;
  std::vector<double> ARvalue;

 private:
  int RnumCol;
  int RnumRow;

  // the 4 std::vectors below always of full length
  std::vector<double> RcolCost;
  std::vector<double> RcolLower;
  std::vector<double> RcolUpper;
  // std::vector<double> Rb;
  std::vector<double> RrowLower;
  std::vector<double> RrowUpper;

  std::vector<int> flagCol;
  std::vector<int> flagRow;

  // testing
  void printA();
  void printAR();

 public:
  // data for actual check
  int numCol;
  int numRow;
  std::vector<int> Astart;
  std::vector<int> Aindex;
  std::vector<double> Avalue;
  std::vector<double> colCost;
  std::vector<double> colLower;
  std::vector<double> colUpper;
  std::vector<double> rowLower;
  std::vector<double> rowUpper;
  int print;

  // basis
  std::vector<HighsBasisStatus> col_status;
  std::vector<HighsBasisStatus> row_status;

  // solution
  std::vector<double> colValue;
  std::vector<double> colDual;
  std::vector<double> rowDual;

  // std::stack<std::vector<double> > bs;
  std::stack<std::vector<std::pair<int, double> > > rLowers;
  std::stack<std::vector<std::pair<int, double> > > rUppers;
  std::stack<std::vector<std::pair<int, double> > > cLowers;
  std::stack<std::vector<std::pair<int, double> > > cUppers;
  std::stack<std::vector<std::pair<int, double> > > costs;
  // std::stack<double> M;

  void passBasis(const std::vector<HighsBasisStatus>& columns,
                 const std::vector<HighsBasisStatus>& rows);

  void replaceBasis(const std::vector<HighsBasisStatus>& columns,
                    const std::vector<HighsBasisStatus>& rows);

  void passSolution(const std::vector<double>& colVal,
                    const std::vector<double>& colDu,
                    const std::vector<double>& rDu);
  // full matrix
  void setMatrixAR(int nCol, int nRow, const std::vector<int>& ARstart_,
                   const std::vector<int>& ARindex_,
                   const std::vector<double>& ARvalue_);
  void setBoundsCostRHS(const std::vector<double>& colUpper_,
                        const std::vector<double>& colLower_,
                        const std::vector<double>& cost,
                        const std::vector<double>& rowLower_,
                        const std::vector<double>& rowUpper_);
  void addChange(int type, int row, int col, double valC, double dualC,
                 double dualR);
  void setFlags(std::vector<int>& r, std::vector<int>& c);
  void makeKKTCheck();
  void resizeProblemMatrix(KktCheck& checker);
  void addCost(int col, double value);
};
#endif /* TEST_KKTCHSTEP_H_ */
