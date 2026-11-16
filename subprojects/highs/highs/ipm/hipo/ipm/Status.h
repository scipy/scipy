#ifndef HIPO_STATUS_H
#define HIPO_STATUS_H

#include <map>
#include <string>

#include "ipm/ipx/ipx_status.h"

namespace hipo {

enum Error {
  kOk = 0,
  kErrorModel,
  kErrorOverflow,
  kErrorAnalyse,
  kErrorFactorise,
  kErrorSolve,
  kErrorIpx,
  kErrorNan,
  kErrorNegativeComponent,
  kErrorInvalidPointer,
  kErrorFailedAllocation,
};

enum Status {
  kStatusTypeStopped = 100,
  kStatusNotSet,
  kStatusNoProgress,
  kStatusMaxIter,
  kStatusTimeLimit,
  kStatusUserInterrupt,

  kStatusTypeFailed,
  kStatusError,
  kStatusUnknown,

  kStatusTypeSolved,
  kStatusImprecise,
  kStatusPrimalInfeasible,
  kStatusDualInfeasible,
  kStatusOptimal,
};

inline Status IpxToHipoStatus(Int ipx_status) {
  static const std::map<Int, Status> status_map{
      {IPX_STATUS_not_run, kStatusNotSet},
      {IPX_STATUS_no_progress, kStatusNoProgress},
      {IPX_STATUS_iter_limit, kStatusMaxIter},
      {IPX_STATUS_time_limit, kStatusTimeLimit},
      {IPX_STATUS_user_interrupt, kStatusUserInterrupt},
      {IPX_STATUS_failed, kStatusError},
      {IPX_STATUS_imprecise, kStatusImprecise},
      {IPX_STATUS_primal_infeas, kStatusPrimalInfeasible},
      {IPX_STATUS_dual_infeas, kStatusDualInfeasible},
      {IPX_STATUS_optimal, kStatusOptimal},
  };

  auto found = status_map.find(ipx_status);
  if (found != status_map.end()) return found->second;
  return kStatusUnknown;
}

inline std::string errorString(Error error) {
  static const std::map<Error, std::string> status_map{
      {kErrorModel, "Error with model"},
      {kErrorOverflow, "Integer overflow"},
      {kErrorAnalyse, "Error in analyse phase"},
      {kErrorFactorise, "Error in factorise phase"},
      {kErrorSolve, "Error in solve phase"},
      {kErrorIpx, "Error in IPX"},
      {kErrorNan, "NaN detected"},
      {kErrorNegativeComponent, "Negative component detected"},
      {kErrorInvalidPointer, "Invalid pointer"},
      {kErrorFailedAllocation, "Failed allocation"},
  };

  auto found = status_map.find(error);
  if (found != status_map.end()) return found->second;
  return "Unknown error";
}

inline std::string statusString(Status status) {
  static const std::map<Status, std::string> status_map{
      {kStatusNotSet, "Not run"},
      {kStatusNoProgress, "No progress"},
      {kStatusMaxIter, "Reached maximum iterations"},
      {kStatusTimeLimit, "Time limit"},
      {kStatusUserInterrupt, "User interrupt"},
      {kStatusError, "Internal error"},
      {kStatusImprecise, "Imprecise solution"},
      {kStatusPrimalInfeasible, "Primal infeasible"},
      {kStatusDualInfeasible, "Dual infeasible"},
      {kStatusOptimal, "Optimal"}};

  auto found = status_map.find(status);
  if (found != status_map.end()) return found->second;
  return "Unknown status";
}

}  // namespace hipo

#endif