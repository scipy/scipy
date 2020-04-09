/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "io/Filereader.h"

#include <cstring>  // For strrchr
#include <stdexcept>

#include "io/FilereaderEms.h"
#include "io/FilereaderLp.h"
#include "io/FilereaderMps.h"

static const char* getFilenameExt(const char* filename) {
  const char* dot = strrchr(filename, '.');
  if (!dot || dot == filename) return "";
  return dot + 1;
}

Filereader* Filereader::getFilereader(const char* filename) {
  Filereader* reader;
  const char* extension = getFilenameExt(filename);
  if (strcmp(extension, "mps") == 0) {
    reader = new FilereaderMps();
  } else if (strcmp(extension, "lp") == 0) {
    reader = new FilereaderLp();
  } else if (strcmp(extension, "ems") == 0) {
    reader = new FilereaderEms();
  } else {
    // use .mps filereader by default
    reader = new FilereaderMps();
  }
  return reader;
}

// Return a string representation of ParseStatus.}
std::string HighsInputStatusToString(HighsInputStatus status) {
  switch (status) {
    case HighsInputStatus::OK:
      return "OK";
      break;
    case HighsInputStatus::FileNotFound:
      return "Error: File not found";
      break;
    case HighsInputStatus::ErrorMatrixDimensions:
      return "Error Matrix Dimensions";
      break;
    case HighsInputStatus::ErrorMatrixIndices:
      return "Error Matrix Indices";
      break;
    case HighsInputStatus::ErrorMatrixStart:
      return "Error Matrix Start";
      break;
    case HighsInputStatus::ErrorMatrixValue:
      return "Error Matrix Value";
      break;
    case HighsInputStatus::ErrorColBounds:
      return "Error Col Bound";
      break;
    case HighsInputStatus::ErrorRowBounds:
      return "Error Row Bounds";
      break;
    case HighsInputStatus::ErrorObjective:
      return "Error Objective";
      break;
  }
  return "";
}
