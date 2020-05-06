/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file test/KktCheck.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef TEST_KKTCHECK_H_
#define TEST_KKTCHECK_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "lp_data/HConst.h"

class KktCheck {
  // model
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

  // Row wise
  std::vector<int> ARstart;
  std::vector<int> ARindex;
  std::vector<double> ARvalue;
  int i, j, k;
  double tol;

  bool istrueGlb;

  // index std::vectors
  std::vector<int> rIndexRev;
  std::vector<int> cIndexRev;

 public:
  int print;
  // solution
  std::vector<double> colValue;
  std::vector<double> colDual;  // lambda
  std::vector<double> rowDual;  // mu

  // basis
  std::vector<HighsBasisStatus> col_status;
  std::vector<HighsBasisStatus> row_status;

  void printAR();
  void makeARCopy();

  void chPrimalBounds();
  void chPrimalFeas();
  void chDualFeas();
  void chComplementarySlackness();
  void chStOfLagrangian();
  void checkBFS();

  void checkKKT();
  void printSol();
  void setIndexVectors(std::vector<int>& rows, std::vector<int>& cols);

  void passSolution(const std::vector<double>& colVal,

                    const std::vector<double>& colDu,
                    const std::vector<double>& rDu);
  void setMatrix(const std::vector<int>& Astart_,
                 const std::vector<int>& Aindex_,
                 const std::vector<double>& Avalue_);
  void setBounds(const std::vector<double>& colUpper_,
                 const std::vector<double>& colLower_);
  void setNumbersCostRHS(int nCol, int nRow,
                         const std::vector<double>& rowLower_,
                         const std::vector<double>& rowUpper_,
                         const std::vector<double>& cost);
};
#endif /* TEST_KKTCHECK_H_ */
