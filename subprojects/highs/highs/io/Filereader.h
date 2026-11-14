/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/Filereader.h
 * @brief
 */
#ifndef IO_FILEREADER_H_
#define IO_FILEREADER_H_

#include "io/HighsIO.h"
#include "lp_data/HighsOptions.h"
#include "model/HighsModel.h"

enum class FilereaderRetcode {
  kOk = 0,
  kWarning = 1,
  kFileNotFound = 2,
  kParserError = 3,
  kNotImplemented = 4,
  kTimeout
};

void interpretFilereaderRetcode(const HighsLogOptions& log_options,
                                const std::string& filename,
                                const FilereaderRetcode code);
std::string extractModelName(const std::string& filename);

class Filereader {
 public:
  virtual FilereaderRetcode readModelFromFile(const HighsOptions& options,
                                              const std::string filename,
                                              HighsModel& model) = 0;
  virtual HighsStatus writeModelToFile(const HighsOptions& options,
                                       const std::string filename,
                                       const HighsModel& model) = 0;
  static Filereader* getFilereader(const HighsLogOptions& log_options,
                                   const std::string& filename);

  virtual ~Filereader(){};
};
#endif
