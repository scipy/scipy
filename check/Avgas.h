/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/Avgas.h
 * @brief Utilities for tests with AVGAS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_AVGAS_H_
#define SIMPLEX_AVGAS_H_

#include <vector>

/**
 * @brief Utilities for tests with AVGAS
 */
class Avgas {
 public:
  void row(int row, int& num_row, int& num_row_nz,
           std::vector<double>& rowLower, std::vector<double>& rowUpper,
           std::vector<int>& ARstart, std::vector<int>& ARindex,
           std::vector<double>& ARvalue);

  void col(int col, int& num_col, int& num_col_nz, std::vector<double>& colCost,
           std::vector<double>& colLower, std::vector<double>& colUpper,
           std::vector<int>& Astart, std::vector<int>& Aindex,
           std::vector<double>& Avalue);
};
#endif /* SIMPLEX_AVGAS_H_ */
