/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef LP_DATA_HIGHS_STATUS_H_
#define LP_DATA_HIGHS_STATUS_H_

#include <string>

// HiGHS status
enum class HighsStatus { OK = 0, Warning, Error };

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status);

// Return the maximum of two HighsStatus and possibly report on
// call_status not being HighsStatus::OK
HighsStatus interpretCallStatus(const HighsStatus call_status,
                                const HighsStatus from_return_status,
                                const std::string& message = "");

// Return the maximum of two HighsStatus
HighsStatus worseStatus(HighsStatus status0, HighsStatus status1);
#endif
