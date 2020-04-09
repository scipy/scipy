/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/LoadProblem.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "io/LoadProblem.h"

#include <sys/stat.h>

#include "io/Filereader.h"
#include "lp_data/HighsLpUtils.h"

HighsStatus loadLpFromFile(const HighsOptions& options, HighsLp& lp) {
  if (options.model_file.size() == 0) return HighsStatus::Error;

  // Make sure it is not a folder.

  struct stat info;
  const char* pathname = options.model_file.c_str();
  printf("loadLpFromFile: %s\n", pathname);
  if (stat(pathname, &info) != 0) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Cannot access %s", pathname);
    return HighsStatus::Error;
  } else if (info.st_mode & S_IFDIR) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "%s is a directory: please specify a file", pathname);
    return HighsStatus::Error;
  }

  Filereader* reader = Filereader::getFilereader(options.model_file.c_str());
  FilereaderRetcode success = reader->readModelFromFile(options, lp);
  delete reader;

  switch (success) {
    case FilereaderRetcode::FILENOTFOUND:
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "File not found");
      return HighsStatus::Error;
    case FilereaderRetcode::PARSERERROR:
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Error when parsing file");
      return HighsStatus::Error;
    default:
      break;
  }

  lp.nnz_ = lp.Avalue_.size();

  // Extract model name.
  std::string name = options.model_file;
  std::size_t found = name.find_last_of("/\\");
  if (found < name.size()) name = name.substr(found + 1);
  found = name.find_last_of(".");
  if (found < name.size()) name.erase(found, name.size() - found);
  lp.model_name_ = name;

  lp.numInt_ = 0;
  for (unsigned int i = 0; i < lp.integrality_.size(); i++)
    if (lp.integrality_[i]) lp.numInt_++;

  // Don't check validity of the LP here: do it when calling highs.initializeLp
  //  return assessLp(lp, options);
  return HighsStatus::OK;
}
