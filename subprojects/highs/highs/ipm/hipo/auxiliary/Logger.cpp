#include "Logger.h"

#include <stdarg.h>

namespace hipo {

Logger::Logger(bool use_printf) : use_printf_{use_printf} {}

void Logger::setOptions(const HighsLogOptions* log_options) {
  log_options_ = log_options;
}

bool Logger::debug(Int level) const {
  // Return true if level agrees with log_dev_level

  if (log_options_ && log_options_->log_dev_level) {
    if (*log_options_->log_dev_level == kHighsLogDevLevelInfo)
      return level == 1;

    if (*log_options_->log_dev_level == kHighsLogDevLevelDetailed)
      return level == 1 || level == 2;

    if (*log_options_->log_dev_level == kHighsLogDevLevelVerbose)
      return level == 1 || level == 2 || level == 3;
  }
  return false;
}

std::string format(double d, Int width, Int prec,
                   std::ios_base::fmtflags floatfield) {
  std::ostringstream s;
  s.precision(prec);
  s.width(width);
  s.setf(floatfield, std::ios_base::floatfield);
  s << d;
  return s.str();
}

std::string integer(Int i, Int width) {
  std::ostringstream s;
  s.width(width);
  s << i;
  return s.str();
}

void logger_printf(const char* format, ...) {
  va_list argptr;
  va_start(argptr, format);
  vprintf(format, argptr);
  va_end(argptr);
}

}  // namespace hipo
