/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HToyIO.h
 * @brief Read an LP or MIP model in toy format
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IO_HTOYIO_H_
#define IO_HTOYIO_H_

//#include <cmath>
//#include <cstring>
//#include <cstdio>
//#include <fstream>
//#include <iostream>

int readToy_LP_cpp(const char* filename, int* m_p, int* n_p, int* maxmin,
                   double* offset, double** A, double** b, double** c,
                   double** lb, double** ub);

int readToy_MIP_cpp(const char* filename, int* numRow_p, int* numCol_p,
                    int* objSense_p, double* objOffset_p, double** A_cw_p,
                    double** rhs_p, double** cost_p, double** lb_p,
                    double** ub_p, int** integerColumn);

#endif /* IO_HTOYIO_H_ */
