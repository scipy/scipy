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

#include "lp_data/HighsLp.h"
#include "lp_data/HighsModelBuilder.h"
#include "lp_data/HighsOptions.h"

enum class HighsInputStatus {
  OK,
  FileNotFound,
  ErrorMatrixDimensions,
  ErrorMatrixIndices,
  ErrorMatrixStart,
  ErrorMatrixValue,
  ErrorColBounds,
  ErrorRowBounds,
  ErrorObjective
};

enum class FilereaderRetcode {
  OK = 0,
  FILENOTFOUND = 1,
  PARSERERROR = 2,
  NOT_IMPLEMENTED = 3
};

class Filereader {
 public:
  virtual FilereaderRetcode readModelFromFile(const HighsOptions& options,
                                              HighsLp& model) = 0;
  virtual FilereaderRetcode readModelFromFile(const char* filename,
                                              HighsModelBuilder& model) = 0;
  virtual HighsStatus writeModelToFile(const HighsOptions& options,
                                       const char* filename,
                                       HighsLp& model) = 0;
  static Filereader* getFilereader(const char* filename);

  virtual ~Filereader(){};
};

// Return a string representation of ParseStatus.
std::string HighsInputStatusToString(HighsInputStatus status);

#endif
