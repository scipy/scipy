/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsLp.h"

bool isSolutionConsistent(const HighsLp& lp, const HighsSolution& solution) {
  if (solution.col_dual.size() == (size_t)lp.numCol_ ||
      solution.col_value.size() == (size_t)lp.numCol_ ||
      solution.row_dual.size() == (size_t)lp.numRow_ ||
      solution.row_value.size() == (size_t)lp.numRow_)
    return true;
  return false;
}
