/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/Filereader.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IO_FILEREADER_H_
#define IO_FILEREADER_H_

#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"

enum class FilereaderRetcode {
  OK = 0,
  FILENOTFOUND = 1,
  PARSERERROR = 2,
  NOT_IMPLEMENTED = 3,
  TIMEOUT
};

void interpretFilereaderRetcode(FILE* logfile, const std::string filename,
                                const FilereaderRetcode code);
std::string extractModelName(const std::string filename);

class Filereader {
 public:
  virtual FilereaderRetcode readModelFromFile(const HighsOptions& options,
                                              HighsLp& model) = 0;
  virtual HighsStatus writeModelToFile(const HighsOptions& options,
                                       const std::string filename,
                                       HighsLp& model) = 0;
  static Filereader* getFilereader(const std::string filename);

  virtual ~Filereader(){};
};
#endif
