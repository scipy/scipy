#ifndef HIPO_LOGGER_H
#define HIPO_LOGGER_H

#include <sstream>

#include "io/HighsIO.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

// Class for Highs logging
// Use setOptions to pass the log_options.

namespace hipo {

void logger_printf(const char* format, ...);

class Logger {
  const HighsLogOptions* log_options_ = nullptr;
  const bool use_printf_ = false;

 public:
  Logger(bool use_printf = false);
  void setOptions(const HighsLogOptions* log_options);
  bool debug(Int level) const;

  template <typename... Args>
  void print(const char* format, Args... args) const {
    if (log_options_)
      highsLogUser(*log_options_, HighsLogType::kInfo, format, args...);
    else if (use_printf_)
      logger_printf(format, args...);
  }

  template <typename... Args>
  void printw(const char* format, Args... args) const {
    if (log_options_)
      highsLogUser(*log_options_, HighsLogType::kWarning, format, args...);
    else if (use_printf_)
      logger_printf(format, args...);
  }

  template <typename... Args>
  void printe(const char* format, Args... args) const {
    if (log_options_)
      highsLogUser(*log_options_, HighsLogType::kError, format, args...);
    else if (use_printf_)
      logger_printf(format, args...);
  }

  template <typename... Args>
  void printInfo(const char* format, Args... args) const {
    if (log_options_)
      highsLogDev(*log_options_, HighsLogType::kInfo, format, args...);
    else if (use_printf_)
      logger_printf(format, args...);
  }

  template <typename... Args>
  void printDetailed(const char* format, Args... args) const {
    if (log_options_)
      highsLogDev(*log_options_, HighsLogType::kDetailed, format, args...);
    else if (use_printf_)
      logger_printf(format, args...);
  }

  template <typename... Args>
  void printVerbose(const char* format, Args... args) const {
    if (log_options_)
      highsLogDev(*log_options_, HighsLogType::kVerbose, format, args...);
    else if (use_printf_)
      logger_printf(format, args...);
  }
};

// Functions to print using streams, taken from IPX.
std::string format(double d, Int width, Int prec,
                   std::ios_base::fmtflags floatfield);
std::string integer(Int i, Int width = 0);
inline std::string sci(double d, Int width, Int prec) {
  return format(d, width, prec, std::ios_base::scientific);
}
inline std::string fix(double d, Int width, Int prec) {
  return format(d, width, prec, std::ios_base::fixed);
}
template <typename T>
std::string textline(const T& text) {
  std::ostringstream s;
  s.setf(std::ios_base::left);
  s.width(32);
  s << text;
  return s.str();
}

}  // namespace hipo
#endif