/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderEms.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#ifndef IO_FILEREADER_EMS_H_
#define IO_FILEREADER_EMS_H_

#include <list>

#include "io/Filereader.h"
#include "io/HighsIO.h"  // For messages.

class FilereaderEms : public Filereader {
 public:
  FilereaderRetcode readModelFromFile(const HighsOptions& options,
                                      HighsLp& model);
  HighsStatus writeModelToFile(const HighsOptions& options,
                               const std::string filename, HighsLp& model);
};

#endif
