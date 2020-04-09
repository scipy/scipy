/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HToyIO_C.h
 * @brief C wrapper to read an LP model in toy format
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IO_HTOYIO_C_H_
#define IO_HTOYIO_C_H_

int readToy_LP_c(const char* filename, int* m_p, int* n_p, int* maxmin,
                 double* offset, double** A, double** b, double** c,
                 double** lb, double** ub);

#endif /* IO_HTOYIO_C_H_ */
