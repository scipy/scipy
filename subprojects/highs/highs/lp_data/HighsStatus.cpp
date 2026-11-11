/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "lp_data/HighsStatus.h"

#include <cassert>

std::string highsStatusToString(HighsStatus status) {
  switch (status) {
    case HighsStatus::kOk:
      return "OK";
    case HighsStatus::kWarning:
      return "Warning";
    case HighsStatus::kError:
      return "Error";
    default:
      assert(1 == 0);
      return "Unrecognised HiGHS status";
  }
}

HighsStatus interpretCallStatus(const HighsLogOptions& log_options,
                                const HighsStatus call_status,
                                const HighsStatus from_return_status,
                                const std::string& message) {
  HighsStatus to_return_status;
  to_return_status = worseStatus(call_status, from_return_status);
  if (call_status != HighsStatus::kOk)
    highsLogDev(log_options, HighsLogType::kWarning,
                "%s return of HighsStatus::%s\n", message.c_str(),
                highsStatusToString(call_status).c_str());
  return to_return_status;
}

HighsStatus worseStatus(const HighsStatus status0, const HighsStatus status1) {
  HighsStatus return_status = HighsStatus::kError;
  if (status0 == HighsStatus::kError || status1 == HighsStatus::kError)
    return_status = HighsStatus::kError;
  else if (status0 == HighsStatus::kWarning || status1 == HighsStatus::kWarning)
    return_status = HighsStatus::kWarning;
  else
    return_status = HighsStatus::kOk;
  return return_status;
}
