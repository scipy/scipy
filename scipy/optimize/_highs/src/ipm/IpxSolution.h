/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/IpxSolution.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IPM_IPX_SOLUTION_H_
#define IPM_IPX_SOLUTION_H_

#include <stdint.h>

#include <vector>
typedef int64_t ipxint;

struct IpxSolution {
  ipxint num_col;
  ipxint num_row;
  std::vector<double> ipx_col_value;
  std::vector<double> ipx_row_value;
  std::vector<double> ipx_col_dual;
  std::vector<double> ipx_row_dual;
  std::vector<ipxint> ipx_col_status;
  std::vector<ipxint> ipx_row_status;
};

#endif
