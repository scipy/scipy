/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/LoadOptions.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#ifndef IO_LOAD_OPTIONS_H_
#define IO_LOAD_OPTIONS_H_

#include "lp_data/HighsOptions.h"

// For extended options to be parsed from a file. Assuming options file is
// specified.
bool loadOptionsFromFile(HighsOptions& options);

#endif
